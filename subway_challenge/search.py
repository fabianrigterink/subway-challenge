"""Route optimizer for the Subway Challenge: time-dependent Large-Neighborhood
Search (LNS) with out-of-system "runs" between dead-end terminals.

The route is represented as an ordered list of target stations (anchors).
``realize_from`` prices an anchor order on the real timetable (earliest-arrival
per leg, optionally using runs). LNS then repeatedly **ruins** a window of the
order and **recreates** it with regret-2 insertion, accepting via simulated
annealing. The recreate is guided by a static station-to-station metric that is
made *run-aware* (terminal runs folded in). Iterate with ``--seed-from``.

This is the method that produced ``solutions/best.json`` (24:24:30). Run:

    python -m subway_challenge.search lns --seed-from solutions/best.json --terminal-runs

For broader exploration before exploitation, ``start-grid`` sweeps start
platforms, service days, and times of day, then writes checkpointed results
under ``reports/``.
"""
from __future__ import annotations

import argparse
import bisect
import collections
import csv
import heapq
import json
import math
import pickle
import random
import sys
import time
from pathlib import Path

from .build_graph import WEEK
from .solver import GRAPH_PKL, hms
from .stations import StationIndex

INF = float("inf")
DAY_NAMES = {
    "mon": 0,
    "monday": 0,
    "tue": 1,
    "tues": 1,
    "tuesday": 1,
    "wed": 2,
    "wednesday": 2,
    "thu": 3,
    "thur": 3,
    "thurs": 3,
    "thursday": 3,
    "fri": 4,
    "friday": 4,
    "sat": 5,
    "saturday": 5,
    "sun": 6,
    "sunday": 6,
}


def _node_tables(G, si):
    """Per-node-id arrays (ids are contiguous): official station id, stop, t."""
    n = G.number_of_nodes()
    off, stop, t = [None] * n, [None] * n, [0] * n
    p2s = si.parent_to_station
    for nid, d in G.nodes(data=True):
        stop[nid] = d["stop"]
        t[nid] = int(d["t"])
        off[nid] = p2s.get(d["station"], d["station"])
    return off, stop, t


def _elapsed(start_t, end_t):
    """Elapsed seconds in the cyclic GTFS week."""
    return int((int(end_t) - int(start_t)) % WEEK)


def _parse_time_of_day(value):
    """Parse HH:MM[:SS] or raw seconds since midnight."""
    value = str(value).strip()
    if ":" not in value:
        return int(value)
    parts = [int(x) for x in value.split(":")]
    if len(parts) == 2:
        h, m = parts
        s = 0
    elif len(parts) == 3:
        h, m, s = parts
    else:
        raise ValueError(f"bad time-of-day: {value!r}")
    return h * 3600 + m * 60 + s


def _parse_days(value):
    days = []
    for raw in str(value).split(","):
        item = raw.strip().lower()
        if not item:
            continue
        if "-" in item and all(x.strip().isdigit() for x in item.split("-", 1)):
            a, b = [int(x.strip()) for x in item.split("-", 1)]
            days.extend(range(a, b + 1))
        elif item in DAY_NAMES:
            days.append(DAY_NAMES[item])
        else:
            days.append(int(item))
    bad = [d for d in days if d < 0 or d > 6]
    if bad:
        raise ValueError(f"days must be in 0..6, got {bad}")
    return days


def _parse_times(value):
    return [_parse_time_of_day(x) for x in str(value).split(",") if x.strip()]


def _stop_index(node_stop, node_t):
    by_stop = collections.defaultdict(list)
    for nid, stop in enumerate(node_stop):
        if stop is not None:
            by_stop[stop].append((node_t[nid], nid))
    return {stop: (tuple(t for t, _ in rows), tuple(n for _, n in rows))
            for stop, rows in ((s, sorted(v)) for s, v in by_stop.items())}


def _first_event_at_or_after(stop_events, stop, target_t):
    """First event for stop at/after target_t, wrapping within the cyclic week."""
    times, nodes = stop_events.get(stop, ((), ()))
    if not times:
        return None
    target_t %= WEEK
    i = bisect.bisect_left(times, target_t)
    if i == len(times):
        i = 0
    return nodes[i]


def _seed_anchors_from_path(spath, si, mode="first"):
    """Extract station anchors from a route path.

    ``first`` is the historical behavior: keep only first visits. ``all`` keeps
    every station transition, including revisits. ``revisit`` keeps first visits
    plus later revisits, which preserves branch/turn structure without anchoring
    every intermediate stop.
    """
    if mode not in {"first", "all", "revisit"}:
        raise ValueError(f"unknown anchor mode: {mode!r}")
    seen, anchors = set(), []
    last = None
    for stop, _t in spath:
        st = si.resolve(stop)
        if st == last:
            continue
        last = st
        if mode == "all":
            anchors.append(st)
        elif st not in seen:
            seen.add(st)
            anchors.append(st)
        elif mode == "revisit":
            anchors.append(st)
    return anchors


def _load_seed_anchors(path_json, si, mode="first"):
    spath = json.loads(Path(path_json).read_text())["path"]
    anchors = _seed_anchors_from_path(spath, si, mode)
    return spath, anchors


def _rotate_anchors(anchors, start_station):
    if start_station not in anchors:
        return list(anchors)
    i = anchors.index(start_station)
    return list(anchors[i:] + anchors[:i])


def _terminal_stops(gtfs_dir="data/gtfs"):
    """Platform stop ids that are the first or last stop of some trip (line termini)."""
    import pandas as pd
    st = pd.read_csv(Path(gtfs_dir) / "stop_times.txt", dtype=str)
    st["stop_sequence"] = st["stop_sequence"].astype(int)
    g = st.sort_values("stop_sequence").groupby("trip_id")["stop_id"]
    return set(g.first()) | set(g.last())


# -- earliest-arrival primitives ---------------------------------------------

def _dijkstra_until_unvisited(adj, node_off, visited, src):
    """From src, return (node, prev) at the first popped node at an unvisited
    station (used by the greedy seed construction)."""
    dist = {src: 0.0}
    prev = {}
    pq = [(0.0, src)]
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, INF):
            continue
        if u != src and node_off[u] not in visited:
            return u, prev
        for v, ed in adj[u].items():
            nd = d + ed["weight"]
            if nd < dist.get(v, INF):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
    return None, prev


def _dijkstra_to_station(adj, node_off, src, target, cap, runs=None):
    """Earliest-arrival from src to the nearest node of ``target`` station (capped).
    With ``runs``, may use out-of-system run shortcuts."""
    dist = {src: 0.0}
    prev = {}
    pq = [(0.0, src)]
    popped = 0
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, INF):
            continue
        popped += 1
        if popped > cap:
            return None, prev
        if u != src and node_off[u] == target:
            return u, prev
        for v, ed in adj[u].items():
            nd = d + ed["weight"]
            if nd < dist.get(v, INF):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
        if runs is not None:
            for v, w, _info in runs.run_successors(u):
                nd = d + w
                if nd < dist.get(v, INF):
                    dist[v] = nd
                    prev[v] = u
                    heapq.heappush(pq, (nd, v))
    return None, prev


def greedy_anchors(adj, tables, canonical, start_node):
    """Greedy nearest-unvisited construction -> (path, anchors, elapsed). Anchors
    are the ordered target stations; the LNS optimizes this order."""
    node_off, _stop, node_t = tables
    visited = {node_off[start_node]}
    path = [start_node]
    anchors = []
    current = start_node
    while len(visited & canonical) < len(canonical):
        tgt, prev = _dijkstra_until_unvisited(adj, node_off, visited, current)
        if tgt is None:
            return None, None, INF
        seg = []
        x = tgt
        while x != current:
            seg.append(x)
            x = prev[x]
        for nid in reversed(seg):
            path.append(nid)
            visited.add(node_off[nid])
        anchors.append(node_off[tgt])
        current = tgt
    return path, anchors, _elapsed(node_t[start_node], node_t[current])


def realize_from(adj, tables, current, visited0, anchors, cap=300000, runs=None,
                 skip_visited=True):
    """Realize an anchor order on the schedule from node ``current`` (with
    ``visited0`` already covered). Anchors covered incidentally are skipped.
    Returns (suffix_path_nodes, end_node, visited) or (None, None, None)."""
    node_off = tables[0]
    visited = set(visited0)
    suffix = []
    for a in anchors:
        if skip_visited and a in visited:
            continue
        tgt, prev = _dijkstra_to_station(adj, node_off, current, a, cap, runs)
        if tgt is None:
            return None, None, None
        seg = []
        x = tgt
        while x != current:
            seg.append(x)
            x = prev[x]
        for nid in reversed(seg):
            suffix.append(nid)
            visited.add(node_off[nid])
        current = tgt
    return suffix, current, visited


# -- static metric & regret recreate -----------------------------------------

def static_station_metric(G, node_off, extra_edges=None):
    """Station-to-station best-case time matrix D (dict of dicts) from min edge
    weights. ``extra_edges`` = [(a, b, w)] (terminal runs) folded in -> run-aware."""
    import networkx as nx
    edge_min = {}
    for u, nbrs in G._adj.items():
        a = node_off[u]
        for v, ed in nbrs.items():
            b = node_off[v]
            if a == b:
                continue
            w = ed["weight"]
            if w < edge_min.get((a, b), INF):
                edge_min[(a, b)] = w
    for a, b, w in (extra_edges or []):
        if w < edge_min.get((a, b), INF):
            edge_min[(a, b)] = w
    SG = nx.DiGraph()
    for (a, b), w in edge_min.items():
        SG.add_edge(a, b, weight=w)
    return dict(nx.all_pairs_dijkstra_path_length(SG)), list(SG.nodes())


def regret_insert(order, removed, D, rng=None, jitter=0.0):
    """Regret-2 insertion: insert each removed station at its best position by D,
    highest-regret first. ``jitter`` occasionally picks the 2nd-best position to
    diversify (lets SA escape local optima)."""
    order = list(order)
    remaining = list(removed)
    while remaining:
        pick = None
        for s in remaining:
            costs = []
            for pos in range(len(order) + 1):
                prev = order[pos - 1] if pos > 0 else None
                nxt = order[pos] if pos < len(order) else None
                c = ((D[prev].get(s, INF) if prev else 0)
                     + (D[s].get(nxt, INF) if nxt else 0)
                     - (D[prev].get(nxt, INF) if (prev and nxt) else 0))
                costs.append((c, pos))
            costs.sort()
            best_c, best_pos = costs[0]
            if rng is not None and jitter and len(costs) > 1 and rng.random() < jitter:
                best_pos = costs[1][1]
            regret = (costs[1][0] - best_c) if len(costs) > 1 else 0
            if pick is None or regret > pick[0]:
                pick = (regret, s, best_pos)
        _, s, pos = pick
        order.insert(pos, s)
        remaining.remove(s)
    return order


# -- LNS ----------------------------------------------------------------------

def lns_run(adj, tables, canonical, D, snode, split, wmin, wmax, budget, t0_temp,
            tend_temp, seed, jitter, anchors, runs, skip_visited=True):
    """One LNS run: freeze a prefix of the order, then ruin+regret-recreate+SA on
    the tail. Seeds from ``anchors`` if given, else greedy. Returns (path, elapsed)."""
    node_off, _stop, node_t = tables
    start_t = node_t[snode]
    if anchors is None:
        _, anchors, _ = greedy_anchors(adj, tables, canonical, snode)
    if not anchors:
        return None, INF
    p = int(len(anchors) * split)
    prefix_anchors, tail = anchors[:p], anchors[p:]
    pre_suffix, pre_end, pre_visited = realize_from(
        adj, tables, snode, {node_off[snode]}, prefix_anchors,
        runs=runs, skip_visited=skip_visited)
    if pre_suffix is None:
        return None, INF
    prefix_path = [snode] + pre_suffix
    base_path, end, _ = realize_from(
        adj, tables, pre_end, pre_visited, tail, runs=runs, skip_visited=skip_visited)
    if base_path is None:
        return None, INF
    cur_tail, cur_e = tail[:], _elapsed(start_t, node_t[end])
    best_tail, best_e, best_path = tail[:], cur_e, base_path
    rng = random.Random(seed)
    t0 = time.time()
    while time.time() - t0 < budget and len(tail) > 4:
        T = t0_temp * (tend_temp / t0_temp) ** ((time.time() - t0) / budget)
        W = rng.randint(wmin, wmax)
        i = rng.randrange(max(1, len(cur_tail) - W))
        removed = cur_tail[i:i + W]
        kept = cur_tail[:i] + cur_tail[i + W:]
        if not removed or not kept:
            continue
        cand = regret_insert(kept, removed, D, rng, jitter)
        tp, tend, tvis = realize_from(
            adj, tables, pre_end, pre_visited, cand, runs=runs, skip_visited=skip_visited)
        if tp is None or len(tvis & canonical) < len(canonical):
            continue
        e = _elapsed(start_t, node_t[tend])
        if e - cur_e < 0 or rng.random() < math.exp(-(e - cur_e) / T):
            cur_tail, cur_e = cand, e
            if e < best_e:
                best_tail, best_e, best_path = cand, e, tp
    return prefix_path + best_path, best_e


def _terminal_run_layer(G, si, node_off, radius_m):
    """Run layer restricted to runs *from* line termini (the high-value
    dead-end-to-line shortcuts), plus the run-aware metric edges."""
    from .run_layer import RunLayer
    from .walk_transfers import load_complexes
    runs = RunLayer.from_graph(G, radius_m=radius_m)
    cid2off = {cid: {si.resolve(p) for p in c.parents} for cid, c in load_complexes().items()}
    deg = collections.defaultdict(set)
    for u, nbrs in G._adj.items():
        a = node_off[u]
        for v, ed in nbrs.items():
            if ed["mode"] == "train" and node_off[v] != a:
                deg[a].add(node_off[v])
    terminals = {s for s, nb in deg.items() if len(nb) == 1} | {si.resolve(s) for s in _terminal_stops()}
    off2cid = {s: cid for cid, ss in cid2off.items() for s in ss}
    tcids = {off2cid[s] for s in terminals if s in off2cid}
    runs.adjacency = {cid: nb for cid, nb in runs.adjacency.items() if cid in tcids}
    extra_edges = [(a, b, secs)
                   for cid, nb in runs.adjacency.items() for ocid, secs, _m in nb
                   for a in cid2off.get(cid, ()) for b in cid2off.get(ocid, ())]
    print(f"terminal runs: {len(terminals)} termini -> {len(runs.adjacency)} run-source complexes")
    return runs, extra_edges


def _run_metric_edges(runs, si):
    from .walk_transfers import load_complexes

    cid2off = {cid: {si.resolve(p) for p in c.parents}
               for cid, c in load_complexes().items()}
    return [(a, b, secs)
            for cid, nb in runs.adjacency.items() for ocid, secs, _m in nb
            for a in cid2off.get(cid, ()) for b in cid2off.get(ocid, ())]


def _run_layer_and_metric_edges(G, si, node_off, mode, radius_m):
    if mode == "none":
        return None, []
    if mode == "terminal":
        return _terminal_run_layer(G, si, node_off, radius_m)
    if mode == "all":
        from .run_layer import RunLayer
        runs = RunLayer.from_graph(G, radius_m=radius_m)
        print(f"all runs: {len(runs.adjacency)} run-source complexes")
        return runs, _run_metric_edges(runs, si)
    raise ValueError(f"unknown run mode: {mode!r}")


def _effective_run_mode(args):
    return "terminal" if getattr(args, "terminal_runs", False) else args.run_mode


def _start_stops(spec, node_off, node_stop, seed_stop=None):
    """Resolve a start-grid spec into platform stop ids.

    Supported presets are ``seed``, ``terminals``, and ``all``. Explicit tokens
    may be platform stop ids (``A02S``) or official station ids, in which case
    all platforms observed for that station are included.
    """
    all_stops = sorted({s for s in node_stop if s})
    stops_by_station = collections.defaultdict(set)
    for off, stop in zip(node_off, node_stop):
        if off and stop:
            stops_by_station[off].add(stop)

    out = []
    for raw in str(spec).split(","):
        item = raw.strip()
        if not item:
            continue
        key = item.lower()
        if key == "seed":
            if seed_stop:
                out.append(seed_stop)
        elif key in {"terminal", "terminals"}:
            terminal_set = _terminal_stops()
            out.extend(s for s in all_stops if s in terminal_set)
        elif key == "all":
            out.extend(all_stops)
        elif item in all_stops:
            out.append(item)
        elif item in stops_by_station:
            out.extend(sorted(stops_by_station[item]))
        else:
            raise ValueError(f"unknown start stop/station: {item!r}")

    seen, unique = set(), []
    for stop in out:
        if stop not in seen:
            seen.add(stop)
            unique.append(stop)
    return unique


def cmd_start_grid(args) -> int:
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    stop_events = _stop_index(node_stop, node_t)

    seed_stop = None
    seed_anchors = None
    if args.seed_from:
        spath, seed_anchors = _load_seed_anchors(args.seed_from, si, args.anchor_mode)
        seed_stop = spath[0][0]
        print(f"start-grid: seeded from {args.seed_from} "
              f"({len(seed_anchors)} {args.anchor_mode} anchors)")

    starts = _start_stops(args.starts, node_off, node_stop, seed_stop)
    if args.max_starts:
        starts = starts[:args.max_starts]
    days = _parse_days(args.days)
    times = _parse_times(args.times)
    splits = [float(x) for x in args.splits.split(",") if x.strip()]

    run_mode = _effective_run_mode(args)
    runs, extra_edges = _run_layer_and_metric_edges(G, si, node_off, run_mode, args.run_radius)
    D, _ = static_station_metric(G, node_off, extra_edges)

    if args.config_csv:
        seen_configs = set()
        base_configs = []
        with Path(args.config_csv).open() as f:
            source_rows = [r for r in csv.DictReader(f) if r.get("elapsed_s")]
        source_rows.sort(key=lambda r: int(r["elapsed_s"]))
        for row in source_rows:
            key = (row["start_stop"], int(row["day"]), int(row["time_of_day_s"]))
            if key in seen_configs:
                continue
            seen_configs.add(key)
            base_configs.append(key)
            if args.config_top_n and len(base_configs) >= args.config_top_n:
                break
    else:
        base_configs = [(st, day, tod) for st in starts for day in days for tod in times]

    all_configs = [(st, day, tod, sp, sd)
                   for st, day, tod in base_configs
                   for sp in splits
                   for sd in range(args.seeds)]
    indexed_configs = list(enumerate(all_configs, start=1))
    if args.skip_configs:
        indexed_configs = indexed_configs[args.skip_configs:]
    if args.max_configs:
        indexed_configs = indexed_configs[:args.max_configs]
    per = args.per_config if args.per_config is not None else args.time_budget / max(1, len(indexed_configs))
    if args.config_csv:
        print(f"start-grid: {len(base_configs)} configs from {args.config_csv} "
              f"x {len(splits)} splits x {args.seeds} seeds = {len(all_configs)} configs",
              flush=True)
    else:
        print(f"start-grid: {len(starts)} starts x {len(days)} days x {len(times)} times "
              f"x {len(splits)} splits x {args.seeds} seeds = {len(all_configs)} configs",
              flush=True)
    if args.skip_configs or args.max_configs:
        print(f"start-grid: running {len(indexed_configs)} configs "
              f"(skip={args.skip_configs}, max={args.max_configs})", flush=True)
    print(f"start-grid: {per:.1f}s per config", flush=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report_csv = out_dir / "start_grid_results.csv"
    csv_fields = [
        "rank_seen", "start_stop", "start_station", "day", "time_of_day_s",
        "target_t", "actual_t", "split", "seed", "elapsed_s", "elapsed", "steps"]
    rows = []
    kept = []
    t0 = time.time()
    interrupted = False
    with report_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields)
        writer.writeheader()
        try:
            for idx, (st, day, tod, sp, sd) in indexed_configs:
                target_t = day * 86400 + tod
                snode = _first_event_at_or_after(stop_events, st, target_t)
                if snode is None:
                    continue
                start_station = node_off[snode]
                anchors = None
                if seed_anchors is not None:
                    anchors = (_rotate_anchors(seed_anchors, start_station)
                               if args.rotate_seed else list(seed_anchors))
                full, elapsed = lns_run(
                    G._adj, tables, si.canonical_stations, D, snode, sp, args.wmin,
                    args.wmax, per, args.t0, args.tend, sd, args.jitter, anchors, runs,
                    skip_visited=args.skip_visited_anchors)
                row = {
                    "rank_seen": idx,
                    "start_stop": st,
                    "start_station": start_station,
                    "day": day,
                    "time_of_day_s": tod,
                    "target_t": target_t,
                    "actual_t": node_t[snode],
                    "split": sp,
                    "seed": sd,
                    "elapsed_s": None if elapsed == INF else int(elapsed),
                    "elapsed": None if elapsed == INF else hms(elapsed),
                    "steps": 0 if not full else len(full),
                }
                rows.append(row)
                writer.writerow(row)
                f.flush()
                if full and elapsed < INF:
                    kept.append((elapsed, row, full))
                    kept.sort(key=lambda x: x[0])
                    del kept[args.top_k:]
                    if elapsed <= kept[0][0]:
                        print(f"  {idx:4d}/{len(all_configs)} best {hms(elapsed)} "
                              f"start={st}@day{day} {hms(tod)} split={sp} seed={sd}",
                              flush=True)
                elif idx % max(1, args.progress_every) == 0:
                    print(f"  {idx:4d}/{len(all_configs)} no route for "
                          f"start={st}@day{day} {hms(tod)}", flush=True)
        except KeyboardInterrupt:
            interrupted = True
            print("start-grid interrupted; writing partial results", flush=True)

    report_json = out_dir / "start_grid_results.json"
    report_json.write_text(json.dumps(rows, indent=2))

    for rank, (elapsed, row, full) in enumerate(kept, start=1):
        out = out_dir / f"{args.prefix}_top{rank:02d}_{int(elapsed)}.json"
        meta = {
            "radius_m": args.run_radius if run_mode != "none" else 5000,
            "run_mode": run_mode,
            "elapsed_s": int(elapsed),
            "notes": ("start-grid "
                      f"start={row['start_stop']} day={row['day']} "
                      f"time={hms(row['time_of_day_s'])} split={row['split']} seed={row['seed']}"),
        }
        out.write_text(json.dumps({"meta": meta,
                                   "path": [[node_stop[n], node_t[n]] for n in full]}))
        print(f"  top {rank}: {hms(elapsed)} -> {out}", flush=True)

    status = "interrupted" if interrupted else "done"
    print(f"start-grid {status} in {time.time() - t0:.1f}s; wrote {report_json} and {report_csv}",
          flush=True)
    if kept:
        print(f"start-grid best: {hms(kept[0][0])} {kept[0][1]}", flush=True)
    else:
        print("start-grid best: no valid realized routes", flush=True)
    return 0


def cmd_lns(args) -> int:
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables

    run_mode = _effective_run_mode(args)
    runs, extra_edges = _run_layer_and_metric_edges(G, si, node_off, run_mode, args.run_radius)
    D, _ = static_station_metric(G, node_off, extra_edges)

    splits = [float(x) for x in args.splits.split(",")]
    seed_anchors = None
    if args.seed_from:                                # iterated LNS: seed from a solution
        spath = json.loads(Path(args.seed_from).read_text())["path"]
        starts = [spath[0][0]]
        seed_anchors = _seed_anchors_from_path(spath, si, args.anchor_mode)
        print(f"iterated LNS: seeded from {args.seed_from} "
              f"({len(seed_anchors)} {args.anchor_mode} anchors)", flush=True)
    else:
        starts = args.start.split(",")

    configs = [(st, sp, args.seed_offset + sd)
               for st in starts for sp in splits for sd in range(args.seeds)]
    per = args.time_budget / len(configs)
    print(f"LNS sweep: {len(configs)} configs x {per:.0f}s each "
          f"(seed_offset={args.seed_offset})", flush=True)

    best = (INF, None, None)
    for i, (st, sp, sd) in enumerate(configs, start=1):
        cands = [n for n, d in G.nodes(data=True) if d["stop"] == st and d["t"] >= args.after]
        if not cands:
            continue
        snode = min(cands, key=lambda n: node_t[n])
        t0_cfg = time.time()
        full, e = lns_run(G._adj, tables, si.canonical_stations, D, snode, sp, args.wmin,
                          args.wmax, per, args.t0, args.tend, sd, args.jitter,
                          seed_anchors, runs,
                          skip_visited=args.skip_visited_anchors)
        status = hms(e) if full and e < INF else "failed"
        print(f"  config {i:3d}/{len(configs)} elapsed={status} "
              f"start={st} split={sp} seed={sd} wall={time.time() - t0_cfg:.1f}s",
              flush=True)
        if full and e < best[0]:
            best = (e, full, (st, sp, sd))
            print(f"  new best {hms(e)} from start={st} split={sp} seed={sd}", flush=True)

    e, full, key = best
    print(f"LNS sweep best: {hms(e)} {key}", flush=True)
    if full is None:
        return 1
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": args.run_radius if run_mode != "none" else 5000,
                                        "run_mode": run_mode,
                                        "elapsed_s": int(e),
                                        "notes": f"LNS {key}"},
                               "path": [[node_stop[n], node_t[n]] for n in full]}))
    print(f"wrote {len(full)} nodes -> {out}", flush=True)
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Subway Challenge optimizer (time-dependent LNS).")
    sub = p.add_subparsers(dest="cmd", required=True)
    ln = sub.add_parser("lns", help="Time-dependent LNS: ruin + regret recreate + SA.")
    ln.add_argument("--start", default="A02S", help="Comma start stops (greedy seed).")
    ln.add_argument("--seed-from", default=None, help="Seed anchor order from a solution JSON (iterate).")
    ln.add_argument("--anchor-mode", choices=("first", "all", "revisit"), default="first",
                    help="How to extract anchors from --seed-from.")
    ln.add_argument("--skip-visited-anchors", dest="skip_visited_anchors",
                    action="store_true", default=True,
                    help="Skip target anchors already covered incidentally.")
    ln.add_argument("--no-skip-visited-anchors", dest="skip_visited_anchors",
                    action="store_false",
                    help="Force realization to obey repeated/revisited anchors.")
    ln.add_argument("--splits", default="0.1,0.15,0.2,0.3", help="Comma freeze fractions.")
    ln.add_argument("--seeds", type=int, default=4, help="SA seeds per (start, split).")
    ln.add_argument("--seed-offset", type=int, default=0, help="First SA seed value to try.")
    ln.add_argument("--after", type=int, default=21600, help="Earliest start t (sec in week).")
    ln.add_argument("--run-mode", choices=("none", "terminal", "all"), default="none",
                    help="Run-transfer policy used in realization and static metric.")
    ln.add_argument("--terminal-runs", action="store_true",
                    help="Alias for --run-mode terminal.")
    ln.add_argument("--run-radius", type=float, default=2500, help="Run layer radius (m).")
    ln.add_argument("--wmin", type=int, default=3)
    ln.add_argument("--wmax", type=int, default=22)
    ln.add_argument("--jitter", type=float, default=0.15)
    ln.add_argument("--time-budget", type=float, default=600)
    ln.add_argument("--t0", type=float, default=600, help="SA start temperature (sec).")
    ln.add_argument("--tend", type=float, default=10, help="SA end temperature (sec).")
    ln.add_argument("--out", default="solutions/lns.json")
    ln.set_defaults(func=cmd_lns)

    sg = sub.add_parser("start-grid", help="Grid-search start stops, service days, and times.")
    sg.add_argument("--seed-from", default="solutions/best.json",
                    help="Anchor order seed to rotate/reuse; use empty string to disable.")
    sg.add_argument("--anchor-mode", choices=("first", "all", "revisit"), default="first",
                    help="How to extract anchors from --seed-from.")
    sg.add_argument("--skip-visited-anchors", dest="skip_visited_anchors",
                    action="store_true", default=True,
                    help="Skip target anchors already covered incidentally.")
    sg.add_argument("--no-skip-visited-anchors", dest="skip_visited_anchors",
                    action="store_false",
                    help="Force realization to obey repeated/revisited anchors.")
    sg.add_argument("--starts", default="terminals",
                    help="Comma list of stop ids/station ids, or presets: seed, terminals, all.")
    sg.add_argument("--max-starts", type=int, default=0,
                    help="Optional cap after resolving --starts, useful for smoke runs.")
    sg.add_argument("--days", default="0-4", help="Comma/range days, Mon=0 ... Sun=6.")
    sg.add_argument("--times", default="05:00,06:00,07:00,08:00",
                    help="Comma times of day as HH:MM[:SS] or seconds since midnight.")
    sg.add_argument("--config-csv", default=None,
                    help="Optional previous start-grid CSV; refines exact start/day/time rows.")
    sg.add_argument("--config-top-n", type=int, default=0,
                    help="With --config-csv, use only the N best unique start/day/time rows.")
    sg.add_argument("--splits", default="0.0,0.1,0.2",
                    help="Comma freeze fractions for the rotated/greedy anchor order.")
    sg.add_argument("--seeds", type=int, default=1, help="SA seeds per grid point.")
    sg.add_argument("--per-config", type=float, default=None,
                    help="Seconds per config; defaults to --time-budget / configs.")
    sg.add_argument("--time-budget", type=float, default=600,
                    help="Total budget if --per-config is not set.")
    sg.add_argument("--max-configs", type=int, default=0,
                    help="Optional cap after building the grid, useful for smoke runs.")
    sg.add_argument("--skip-configs", type=int, default=0,
                    help="Skip this many grid configs before running; useful for resuming chunks.")
    sg.add_argument("--run-mode", choices=("none", "terminal", "all"), default="none",
                    help="Run-transfer policy used in realization and static metric.")
    sg.add_argument("--terminal-runs", action="store_true",
                    help="Alias for --run-mode terminal.")
    sg.add_argument("--run-radius", type=float, default=2500, help="Run layer radius (m).")
    sg.add_argument("--wmin", type=int, default=3)
    sg.add_argument("--wmax", type=int, default=18)
    sg.add_argument("--jitter", type=float, default=0.2)
    sg.add_argument("--t0", type=float, default=600, help="SA start temperature (sec).")
    sg.add_argument("--tend", type=float, default=10, help="SA end temperature (sec).")
    sg.add_argument("--rotate-seed", dest="rotate_seed", action="store_true", default=True,
                    help="Rotate seeded anchor order so the selected start station is first.")
    sg.add_argument("--no-rotate-seed", dest="rotate_seed", action="store_false")
    sg.add_argument("--top-k", type=int, default=5, help="Number of top route JSONs to write.")
    sg.add_argument("--progress-every", type=int, default=25)
    sg.add_argument("--out-dir", default="reports/start_grid")
    sg.add_argument("--prefix", default="start_grid")
    sg.set_defaults(func=cmd_start_grid)
    args = p.parse_args(argv)
    if getattr(args, "seed_from", None) == "":
        args.seed_from = None
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
