"""Exact local repair passes for already-valid Subway Challenge routes.

The LNS route is a legal time-expanded path, but some local slices may still be
non-shortest between their exact endpoints because the route was optimizing
station coverage rather than endpoint-to-endpoint travel. This module searches
for faster replacement subpaths, splices them back into the whole route, and
keeps only candidates that still validate against all 472 stations.
"""
from __future__ import annotations

import argparse
import collections
import glob
import heapq
import json
import pickle
import random
import time
from pathlib import Path

from .run_layer import RunLayer
from .search import (_node_tables, _terminal_run_layer, realize_from,
                     regret_insert, static_station_metric)
from .solver import GRAPH_PKL, Problem, hms, result_line, validate
from .stations import StationIndex

INF = float("inf")


def _load_solution(path_json):
    data = json.loads(Path(path_json).read_text())
    return data.get("path", data), data.get("meta", {})


def _nodes_from_path(prob, path):
    nodes = []
    for stop, t, *_rest in path:
        nid = prob.node_of(stop, int(t))
        if nid is None:
            raise ValueError(f"node {stop}@{int(t)} not found in graph")
        nodes.append(nid)
    return nodes


def _path_from_nodes(G, nodes):
    return [[G.nodes[n]["stop"], int(G.nodes[n]["t"])] for n in nodes]


def _path_cost_modes(prob, nodes):
    total = 0
    modes = collections.Counter()
    for u, v in zip(nodes, nodes[1:]):
        w, mode = prob.transition(u, v)
        if w is None:
            du, dv = prob.G.nodes[u], prob.G.nodes[v]
            raise ValueError(
                f"illegal transition {du['stop']}@{int(du['t'])} -> "
                f"{dv['stop']}@{int(dv['t'])}"
            )
        total += int(w)
        modes[mode] += 1
    return total, modes


def _covered_stations(prob, nodes):
    return {prob.stations.resolve(prob.G.nodes[n]["stop"]) for n in nodes} & prob.canonical


def _reconstruct(prev, src, target):
    path = [target]
    x = target
    while x != src:
        x = prev[x]
        path.append(x)
    path.reverse()
    return path


def _dijkstra_to_node(adj, src, target, max_cost, cap, runs=None):
    if src == target:
        return [src], 0, 0
    dist = {src: 0}
    prev = {}
    pq = [(0, src)]
    popped = 0
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, INF):
            continue
        if d > max_cost:
            break
        popped += 1
        if popped > cap:
            return None, INF, popped
        if u == target:
            return _reconstruct(prev, src, target), int(d), popped
        for v, ed in adj[u].items():
            nd = d + int(ed["weight"])
            if nd <= max_cost and nd < dist.get(v, INF):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
        if runs is not None:
            for v, w, _info in runs.run_successors(u):
                nd = d + int(w)
                if nd <= max_cost and nd < dist.get(v, INF):
                    dist[v] = nd
                    prev[v] = u
                    heapq.heappush(pq, (nd, v))
    return None, INF, popped


def _parse_ints(value):
    return [int(x) for x in str(value).split(",") if x.strip()]


def _run_metric_edges(runs, si):
    from .walk_transfers import load_complexes

    cid2off = {cid: {si.resolve(p) for p in c.parents}
               for cid, c in load_complexes().items()}
    return [(a, b, secs)
            for cid, nb in runs.adjacency.items() for ocid, secs, _m in nb
            for a in cid2off.get(cid, ()) for b in cid2off.get(ocid, ())]


def _run_layer_and_edges(G, si, node_off, mode, radius_m):
    if mode == "none":
        return None, []
    if mode == "terminal":
        return _terminal_run_layer(G, si, node_off, radius_m)
    if mode == "all":
        runs = RunLayer.from_graph(G, radius_m=radius_m)
        return runs, _run_metric_edges(runs, si)
    raise ValueError(f"unknown run mode: {mode!r}")


def _run_layer(G, si, node_off, mode, radius_m):
    runs, _extra_edges = _run_layer_and_edges(G, si, node_off, mode, radius_m)
    return runs


def _window_stats(prob, nodes, a, b):
    segment = nodes[a:b + 1]
    cost, modes = _path_cost_modes(prob, segment)
    covered = _covered_stations(prob, segment)
    score = cost + 180 * modes.get("run", 0) + 90 * modes.get("wait", 0)
    return {
        "a": a,
        "b": b,
        "edges": b - a,
        "cost": cost,
        "modes": dict(modes),
        "covered": len(covered),
        "score": score,
    }


def _candidate_windows(prob, nodes, widths, stride, max_windows):
    windows = []
    n = len(nodes)
    for width in widths:
        if width < 1 or width >= n:
            continue
        for a in range(0, n - width, stride):
            b = a + width
            windows.append(_window_stats(prob, nodes, a, b))
    windows.sort(key=lambda r: (r["score"], r["cost"], r["edges"]), reverse=True)
    if max_windows:
        windows = windows[:max_windows]
    return windows


def _anchor_sequence_with_positions(prob, nodes, mode):
    if mode not in {"first", "revisit"}:
        raise ValueError(f"unknown anchor mode: {mode!r}")
    anchors = []
    positions = []
    seen = set()
    last = None
    for i, n in enumerate(nodes):
        st = prob.stations.resolve(prob.G.nodes[n]["stop"])
        if st == last:
            continue
        last = st
        if st not in prob.canonical:
            continue
        if st not in seen:
            seen.add(st)
            anchors.append(st)
            positions.append(i)
        elif mode == "revisit":
            anchors.append(st)
            positions.append(i)
    return anchors, positions


def _try_window(prob, adj, nodes, window, runs, args):
    a, b = window["a"], window["b"]
    src, target = nodes[a], nodes[b]
    original_cost = window["cost"]
    max_cost = original_cost - args.min_saving
    if max_cost <= 0:
        return None

    alt, alt_cost, popped = _dijkstra_to_node(
        adj, src, target, max_cost=max_cost, cap=args.cap, runs=runs)
    if alt is None:
        return {"status": "no_path", "popped": popped}
    if alt == nodes[a:b + 1]:
        return {"status": "same_path", "popped": popped}

    candidate_nodes = nodes[:a] + alt + nodes[b + 1:]
    candidate_path = _path_from_nodes(prob.G, candidate_nodes)
    r = validate(prob, candidate_path)
    saving = original_cost - alt_cost
    return {
        "status": "valid" if r["valid"] else "invalid",
        "result": r,
        "saving": saving,
        "alt_cost": alt_cost,
        "alt_nodes": len(alt),
        "candidate_nodes": candidate_nodes if r["valid"] else None,
        "reason": r.get("reason"),
        "missing": r.get("total_stations", 0) - r.get("covered", 0),
        "popped": popped,
    }


def cmd_exact_window(args):
    with open(args.graph, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    node_off, _node_stop, _node_t = _node_tables(G, si)

    seed_path, seed_meta = _load_solution(args.seed)
    validation_radius = args.validation_radius or seed_meta.get("radius_m", 5000)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=validation_radius), si)
    current_nodes = _nodes_from_path(prob, seed_path)
    current_path = _path_from_nodes(G, current_nodes)
    current_result = validate(prob, current_path)
    if not current_result["valid"]:
        print(result_line(current_result))
        return 1

    runs = _run_layer(G, si, node_off, args.run_mode, args.run_radius)
    widths = _parse_ints(args.widths)
    best_nodes = current_nodes
    best_result = current_result
    started = time.time()
    attempts = 0
    accepted = 0
    invalid = 0

    print(f"exact-window: seed {args.seed} {hms(best_result['elapsed_s'])} "
          f"nodes={len(best_nodes)} run_mode={args.run_mode} radius={args.run_radius}",
          flush=True)

    for pass_no in range(1, args.passes + 1):
        windows = _candidate_windows(prob, best_nodes, widths, args.stride, args.max_windows)
        pass_best = None
        print(f"pass {pass_no}: testing {len(windows)} windows "
              f"widths={widths} stride={args.stride}", flush=True)
        for i, window in enumerate(windows, start=1):
            attempts += 1
            out = _try_window(prob, G._adj, best_nodes, window, runs, args)
            if out is None:
                continue
            if out["status"] == "invalid":
                invalid += 1
            if out["status"] == "valid":
                r = out["result"]
                delta = best_result["elapsed_s"] - r["elapsed_s"]
                if delta >= args.min_saving:
                    accepted += 1
                    candidate = (r["elapsed_s"], delta, window, out)
                    if pass_best is None or candidate[0] < pass_best[0]:
                        pass_best = candidate
                    print(f"  improvement candidate {hms(r['elapsed_s'])} "
                          f"delta={delta}s window={window['a']}:{window['b']} "
                          f"orig={hms(window['cost'])} alt={hms(out['alt_cost'])} "
                          f"nodes={out['alt_nodes']} popped={out['popped']}",
                          flush=True)
            if args.progress_every and i % args.progress_every == 0:
                elapsed = time.time() - started
                print(f"  checked {i}/{len(windows)} in pass {pass_no}; "
                      f"best={hms(best_result['elapsed_s'])} "
                      f"accepted={accepted} invalid={invalid} wall={elapsed:.1f}s",
                      flush=True)
            if args.time_budget and time.time() - started >= args.time_budget:
                print("exact-window: time budget reached", flush=True)
                break

        if pass_best is None:
            print(f"pass {pass_no}: no accepted improvement", flush=True)
            break
        elapsed_s, delta, window, out = pass_best
        best_nodes = out["candidate_nodes"]
        best_result = out["result"]
        print(f"pass {pass_no}: accepted {hms(elapsed_s)} "
              f"delta={delta}s window={window['a']}:{window['b']}", flush=True)
        if args.time_budget and time.time() - started >= args.time_budget:
            break

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    improved = best_result["elapsed_s"] < current_result["elapsed_s"]
    meta = {
        "radius_m": validation_radius,
        "elapsed_s": best_result["elapsed_s"],
        "improved": improved,
        "source": str(args.seed),
        "notes": (
            "exact-window repair "
            f"run_mode={args.run_mode} run_radius={args.run_radius} "
            f"widths={args.widths} stride={args.stride} attempts={attempts}"
        ),
    }
    out_path.write_text(json.dumps({"meta": meta, "path": _path_from_nodes(G, best_nodes)},
                                   indent=2))
    print(result_line(best_result))
    print(f"exact-window: wrote {out_path} improved={improved} "
          f"attempts={attempts} accepted={accepted} invalid={invalid} "
          f"wall={time.time() - started:.1f}s", flush=True)
    return 0


def _anchor_window_stats(prob, nodes, positions, i, width):
    start_pos = 0 if i == 0 else positions[i - 1]
    end_pos = positions[min(len(positions) - 1, i + width - 1)]
    segment = nodes[start_pos:end_pos + 1]
    cost, modes = _path_cost_modes(prob, segment)
    return {
        "i": i,
        "width": width,
        "start_pos": start_pos,
        "end_pos": end_pos,
        "cost": cost,
        "modes": dict(modes),
        "score": cost + 180 * modes.get("run", 0) + 90 * modes.get("wait", 0),
    }


def _anchor_windows(prob, nodes, positions, widths, stride, max_windows):
    out = []
    n = len(positions)
    for width in widths:
        if width < 1 or width >= n:
            continue
        for i in range(0, n - width, stride):
            out.append(_anchor_window_stats(prob, nodes, positions, i, width))
    out.sort(key=lambda r: (r["score"], r["cost"], r["width"]), reverse=True)
    if max_windows:
        out = out[:max_windows]
    return out


def _block_insert_tail(kept, removed, rng):
    tail = list(kept)
    block = list(removed)
    rng.shuffle(block)
    pos = rng.randrange(len(tail) + 1) if tail else 0
    return tail[:pos] + block + tail[pos:]


def cmd_anchor_reinsert(args):
    with open(args.graph, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    seed_path, seed_meta = _load_solution(args.seed)
    validation_radius = args.validation_radius or seed_meta.get("radius_m", 5000)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=validation_radius), si)
    seed_nodes = _nodes_from_path(prob, seed_path)
    seed_result = validate(prob, _path_from_nodes(G, seed_nodes))
    if not seed_result["valid"]:
        print(result_line(seed_result))
        return 1

    runs, extra_edges = _run_layer_and_edges(G, si, node_off, args.run_mode, args.run_radius)
    D, _metric_nodes = static_station_metric(G, node_off, extra_edges)
    anchors, positions = _anchor_sequence_with_positions(prob, seed_nodes, args.anchor_mode)
    widths = _parse_ints(args.widths)
    windows = _anchor_windows(prob, seed_nodes, positions, widths, args.stride, args.max_windows)
    rng = random.Random(args.seed_value)
    best_nodes = seed_nodes
    best_result = seed_result
    attempts = 0
    feasible = 0
    improved = 0
    started = time.time()

    print(f"anchor-reinsert: seed {args.seed} {hms(seed_result['elapsed_s'])} "
          f"anchors={len(anchors)} mode={args.anchor_mode} windows={len(windows)} "
          f"run_mode={args.run_mode} radius={args.run_radius}",
          flush=True)

    for wi, window in enumerate(windows, start=1):
        i, width = window["i"], window["width"]
        cut_pos = 0 if i == 0 else positions[i - 1]
        prefix_nodes = seed_nodes[:cut_pos + 1]
        current = prefix_nodes[-1]
        visited0 = _covered_stations(prob, prefix_nodes)
        removed = anchors[i:i + width]
        kept_tail = anchors[i + width:]
        trial_tails = [list(removed) + list(kept_tail)]
        if kept_tail:
            trial_tails.append(regret_insert(list(kept_tail), list(removed), D, rng, 0.0))
        while len(trial_tails) < args.attempts_per_window:
            if kept_tail and rng.random() < args.regret_probability:
                trial_tails.append(regret_insert(
                    list(kept_tail), list(removed), D, rng, args.jitter))
            else:
                trial_tails.append(_block_insert_tail(kept_tail, removed, rng))

        for tail in trial_tails[:args.attempts_per_window]:
            attempts += 1
            suffix, end, visited = realize_from(
                G._adj, tables, current, visited0, tail, cap=args.cap, runs=runs,
                skip_visited=args.skip_visited_anchors)
            if suffix is None or len(visited & si.canonical_stations) < len(si.canonical_stations):
                continue
            feasible += 1
            cand_nodes = prefix_nodes + suffix
            r = validate(prob, _path_from_nodes(G, cand_nodes))
            if not r["valid"]:
                continue
            delta = best_result["elapsed_s"] - r["elapsed_s"]
            if delta >= args.min_saving:
                improved += 1
                best_nodes = cand_nodes
                best_result = r
                print(f"  improvement {hms(r['elapsed_s'])} delta={delta}s "
                      f"anchor_window={i}:{i + width} cut_node={cut_pos} "
                      f"tail={len(tail)} end={node_stop[end]}@{node_t[end]}",
                      flush=True)
        if args.progress_every and wi % args.progress_every == 0:
            print(f"  checked {wi}/{len(windows)} windows; "
                  f"best={hms(best_result['elapsed_s'])} attempts={attempts} "
                  f"feasible={feasible} improved={improved} "
                  f"wall={time.time() - started:.1f}s",
                  flush=True)
        if args.time_budget and time.time() - started >= args.time_budget:
            print("anchor-reinsert: time budget reached", flush=True)
            break

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    meta = {
        "radius_m": validation_radius,
        "elapsed_s": best_result["elapsed_s"],
        "improved": best_result["elapsed_s"] < seed_result["elapsed_s"],
        "source": str(args.seed),
        "notes": (
            "anchor-reinsert "
            f"anchor_mode={args.anchor_mode} run_mode={args.run_mode} "
            f"run_radius={args.run_radius} widths={args.widths} "
            f"attempts={attempts}"
        ),
    }
    out_path.write_text(json.dumps({"meta": meta, "path": _path_from_nodes(G, best_nodes)},
                                   indent=2))
    print(result_line(best_result))
    print(f"anchor-reinsert: wrote {out_path} attempts={attempts} feasible={feasible} "
          f"improved={improved} wall={time.time() - started:.1f}s", flush=True)
    return 0


def _expand_patterns(patterns):
    files = []
    seen = set()
    for raw in str(patterns).split(","):
        pat = raw.strip()
        if not pat:
            continue
        matches = sorted(glob.glob(pat))
        if not matches and Path(pat).exists():
            matches = [pat]
        for path in matches:
            if path not in seen:
                seen.add(path)
                files.append(path)
    return files


def _route_record(prob, path_json):
    path, meta = _load_solution(path_json)
    nodes = _nodes_from_path(prob, path)
    result = validate(prob, _path_from_nodes(prob.G, nodes))
    if not result["valid"]:
        return None

    prefix_cost = [0]
    for u, v in zip(nodes, nodes[1:]):
        w, _mode = prob.transition(u, v)
        prefix_cost.append(prefix_cost[-1] + int(w))

    prefix_cov = []
    seen = set()
    for n in nodes:
        st = prob.stations.resolve(prob.G.nodes[n]["stop"])
        if st in prob.canonical:
            seen.add(st)
        prefix_cov.append(set(seen))

    suffix_cov = [set() for _ in nodes]
    seen = set()
    for i in range(len(nodes) - 1, -1, -1):
        st = prob.stations.resolve(prob.G.nodes[nodes[i]]["stop"])
        if st in prob.canonical:
            seen.add(st)
        suffix_cov[i] = set(seen)

    positions = collections.defaultdict(list)
    for i, n in enumerate(nodes):
        positions[n].append(i)

    return {
        "path_json": str(path_json),
        "meta": meta,
        "nodes": nodes,
        "result": result,
        "prefix_cost": prefix_cost,
        "prefix_cov": prefix_cov,
        "suffix_cov": suffix_cov,
        "positions": positions,
    }


def cmd_splice(args):
    with open(args.graph, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)
    route_files = _expand_patterns(args.routes)
    if not route_files:
        raise ValueError(f"no route files matched {args.routes!r}")

    routes = []
    for path_json in route_files:
        rec = _route_record(prob, path_json)
        if rec is None:
            print(f"skip invalid {path_json}", flush=True)
            continue
        routes.append(rec)
        print(f"route {len(routes):02d}: {hms(rec['result']['elapsed_s'])} "
              f"nodes={len(rec['nodes'])} {path_json}", flush=True)
    if len(routes) < 2:
        print("splice: need at least two valid routes")
        return 1

    baseline = min(r["result"]["elapsed_s"] for r in routes)
    candidates = []
    checked = 0
    coverage_hits = 0
    validated = 0
    started = time.time()

    for ai, a in enumerate(routes):
        for bi, b in enumerate(routes):
            if ai == bi and not args.self_splice:
                continue
            shared = set(a["positions"]) & set(b["positions"])
            for node in shared:
                for ix in a["positions"][node]:
                    if ix == len(a["nodes"]) - 1:
                        continue
                    for jx in b["positions"][node]:
                        if jx == 0:
                            continue
                        checked += 1
                        cov = a["prefix_cov"][ix] | b["suffix_cov"][jx]
                        if len(cov) < len(prob.canonical):
                            continue
                        coverage_hits += 1
                        estimate = (a["prefix_cost"][ix]
                                    + b["result"]["elapsed_s"] - b["prefix_cost"][jx])
                        if estimate >= baseline - args.min_saving:
                            continue
                        cand_nodes = a["nodes"][:ix + 1] + b["nodes"][jx + 1:]
                        result = validate(prob, _path_from_nodes(G, cand_nodes))
                        validated += 1
                        if result["valid"]:
                            candidates.append({
                                "elapsed_s": result["elapsed_s"],
                                "nodes": cand_nodes,
                                "result": result,
                                "prefix_route": a["path_json"],
                                "suffix_route": b["path_json"],
                                "prefix_index": ix,
                                "suffix_index": jx,
                                "shared_node": node,
                            })
                            print(f"  splice candidate {hms(result['elapsed_s'])} "
                                  f"prefix={Path(a['path_json']).name}:{ix} "
                                  f"suffix={Path(b['path_json']).name}:{jx}",
                                  flush=True)
                        if args.time_budget and time.time() - started >= args.time_budget:
                            break
                    if args.time_budget and time.time() - started >= args.time_budget:
                        break
                if args.time_budget and time.time() - started >= args.time_budget:
                    break
            if args.time_budget and time.time() - started >= args.time_budget:
                break
        if args.time_budget and time.time() - started >= args.time_budget:
            break

    candidates.sort(key=lambda c: c["elapsed_s"])
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for rank, cand in enumerate(candidates[:args.top_k], start=1):
        out = out_dir / f"splice_top{rank:02d}_{cand['elapsed_s']}.json"
        meta = {
            "radius_m": args.validation_radius,
            "elapsed_s": cand["elapsed_s"],
            "source": "route-splice",
            "prefix_route": cand["prefix_route"],
            "suffix_route": cand["suffix_route"],
            "prefix_index": cand["prefix_index"],
            "suffix_index": cand["suffix_index"],
            "notes": "exact shared-event route splice",
        }
        out.write_text(json.dumps({"meta": meta, "path": _path_from_nodes(G, cand["nodes"])},
                                  indent=2))
        print(f"  wrote {hms(cand['elapsed_s'])} -> {out}", flush=True)

    best = candidates[0]["elapsed_s"] if candidates else baseline
    print(f"splice: checked={checked} coverage_hits={coverage_hits} "
          f"validated={validated} candidates={len(candidates)} "
          f"baseline={hms(baseline)} best={hms(best)} "
          f"wall={time.time() - started:.1f}s", flush=True)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description="Exact local route repair tools.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    ew = sub.add_parser("exact-window",
                        help="Replace exact route windows with faster exact-endpoint paths.")
    ew.add_argument("--seed", default="solutions/best.json")
    ew.add_argument("--out", default="reports/optimization_runs/exact_window_repair.json")
    ew.add_argument("--graph", default=GRAPH_PKL)
    ew.add_argument("--run-mode", choices=("none", "terminal", "all"), default="all")
    ew.add_argument("--run-radius", type=float, default=5000)
    ew.add_argument("--validation-radius", type=float, default=5000)
    ew.add_argument("--widths", default="6,10,16,24,32,48")
    ew.add_argument("--stride", type=int, default=2)
    ew.add_argument("--max-windows", type=int, default=300)
    ew.add_argument("--min-saving", type=int, default=30)
    ew.add_argument("--cap", type=int, default=120000,
                    help="Maximum popped nodes per window Dijkstra.")
    ew.add_argument("--passes", type=int, default=3)
    ew.add_argument("--time-budget", type=float, default=0,
                    help="Optional wall-clock budget in seconds.")
    ew.add_argument("--progress-every", type=int, default=25)
    ew.set_defaults(func=cmd_exact_window)

    ar = sub.add_parser("anchor-reinsert",
                        help="Freeze an exact prefix, reinsert anchor windows, and re-realize suffix.")
    ar.add_argument("--seed", default="solutions/best.json")
    ar.add_argument("--out", default="reports/optimization_runs/anchor_reinsert.json")
    ar.add_argument("--graph", default=GRAPH_PKL)
    ar.add_argument("--run-mode", choices=("none", "terminal", "all"), default="all")
    ar.add_argument("--run-radius", type=float, default=5000)
    ar.add_argument("--validation-radius", type=float, default=5000)
    ar.add_argument("--anchor-mode", choices=("first", "revisit"), default="first")
    ar.add_argument("--skip-visited-anchors", dest="skip_visited_anchors",
                    action="store_true", default=True)
    ar.add_argument("--no-skip-visited-anchors", dest="skip_visited_anchors",
                    action="store_false")
    ar.add_argument("--widths", default="4,6,8,12,16,24")
    ar.add_argument("--stride", type=int, default=4)
    ar.add_argument("--max-windows", type=int, default=120)
    ar.add_argument("--attempts-per-window", type=int, default=4)
    ar.add_argument("--regret-probability", type=float, default=0.75)
    ar.add_argument("--jitter", type=float, default=0.2)
    ar.add_argument("--seed-value", type=int, default=0)
    ar.add_argument("--min-saving", type=int, default=30)
    ar.add_argument("--cap", type=int, default=300000)
    ar.add_argument("--time-budget", type=float, default=0)
    ar.add_argument("--progress-every", type=int, default=10)
    ar.set_defaults(func=cmd_anchor_reinsert)

    sp = sub.add_parser("splice",
                        help="Combine prefixes/suffixes of valid routes at shared exact events.")
    sp.add_argument("--routes", default="solutions/*.json",
                    help="Comma-separated route globs/files.")
    sp.add_argument("--out-dir", default="reports/optimization_runs/splices")
    sp.add_argument("--graph", default=GRAPH_PKL)
    sp.add_argument("--validation-radius", type=float, default=5000)
    sp.add_argument("--min-saving", type=int, default=30)
    sp.add_argument("--top-k", type=int, default=10)
    sp.add_argument("--self-splice", action="store_true")
    sp.add_argument("--time-budget", type=float, default=0)
    sp.set_defaults(func=cmd_splice)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
