"""Route construction/search for the Subway Challenge.

v1: greedy nearest-unvisited construction on the time-expanded graph (graph edges
only -- train/wait/transfer -- which already connect all 475 stations). From the
current node, an earliest-arrival Dijkstra runs until it reaches the nearest node
at an unvisited station; we commit that segment, mark every station passed as
visited, and repeat until all 472 official stations are covered.

Runs (the on-demand layer) are intentionally NOT used for construction -- they
are an optimization for the improvement phase. Output is a solutions JSON that
``solver.py`` validates and scores.

    python -m subway_challenge.search greedy --start 101S --after 18000 --out solutions/greedy.json
"""
from __future__ import annotations

import argparse
import heapq
import json
import pickle
import sys
import time
from pathlib import Path

from .solver import GRAPH_PKL, hms
from .stations import StationIndex

INF = float("inf")


def _node_tables(G, si):
    """Per-node-id arrays: official station id, stop, t (ids are contiguous)."""
    n = G.number_of_nodes()
    off, stop, t = [None] * n, [None] * n, [0] * n
    p2s = si.parent_to_station
    for nid, d in G.nodes(data=True):
        stop[nid] = d["stop"]
        t[nid] = int(d["t"])
        off[nid] = p2s.get(d["station"], d["station"])
    return off, stop, t


def _pick_start(G, stop_pref: str, after: int) -> int:
    """Earliest node at ``stop_pref`` with t >= after that has a train out-edge."""
    adj = G._adj
    best, best_t = None, INF
    for nid, d in G.nodes(data=True):
        if d["stop"] == stop_pref and after <= d["t"] < best_t:
            if any(e["mode"] == "train" for e in adj[nid].values()):
                best, best_t = nid, d["t"]
    if best is None:
        raise SystemExit(f"no start node at {stop_pref} after {after}")
    return best


def _dijkstra_until_unvisited(adj, node_off, visited, src):
    """Earliest-arrival search from src; return (node, prev) at the first popped
    node whose station is unvisited, or (None, prev)."""
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


def _dijkstra_until_unvisited_runs(adj, node_off, visited, src, runs):
    """Earliest-arrival search using graph edges AND on-demand run successors."""
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
        for v, w, _info in runs.run_successors(u):
            nd = d + w
            if nd < dist.get(v, INF):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
    return None, prev


def greedy_from_runs(adj, tables, canonical, start_node, runs):
    node_off, _stop, node_t = tables
    visited = {node_off[start_node]}
    path = [start_node]
    current = start_node
    while len(visited & canonical) < len(canonical):
        tgt, prev = _dijkstra_until_unvisited_runs(adj, node_off, visited, current, runs)
        if tgt is None:
            return None, INF
        seg = []
        x = tgt
        while x != current:
            seg.append(x)
            x = prev[x]
        for nid in reversed(seg):
            path.append(nid)
            visited.add(node_off[nid])
        current = tgt
    return path, node_t[current] - node_t[start_node]


def greedy_from(adj, tables, canonical, start_node):
    """Greedy nearest-unvisited from start_node. Returns (path_node_ids, elapsed_s)
    or (None, INF) if it gets stuck. Reuses precomputed node tables."""
    node_off, node_stop, node_t = tables
    visited = {node_off[start_node]}
    path = [start_node]
    current = start_node
    while len(visited & canonical) < len(canonical):
        tgt, prev = _dijkstra_until_unvisited(adj, node_off, visited, current)
        if tgt is None:
            return None, INF
        seg = []
        x = tgt
        while x != current:
            seg.append(x)
            x = prev[x]
        for nid in reversed(seg):
            path.append(nid)
            visited.add(node_off[nid])
        current = tgt
    return path, node_t[current] - node_t[start_node]


def greedy_anchors(adj, tables, canonical, start_node):
    """Greedy construction returning (path, anchors, elapsed). Anchors = the
    ordered list of target stations the greedy explicitly headed to."""
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
    return path, anchors, node_t[current] - node_t[start_node]


def _dijkstra_to_station(adj, node_off, src, target_station):
    """Earliest-arrival from src to the nearest node of target_station."""
    dist = {src: 0.0}
    prev = {}
    pq = [(0.0, src)]
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, INF):
            continue
        if u != src and node_off[u] == target_station:
            return u, prev
        for v, ed in adj[u].items():
            nd = d + ed["weight"]
            if nd < dist.get(v, INF):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
    return None, prev


def realize(adj, tables, canonical, start_node, anchors):
    """Realize an anchor order on the schedule. Anchors already covered
    incidentally are skipped. Returns (path, elapsed) or (None, INF) if the order
    fails to cover all stations."""
    node_off, _stop, node_t = tables
    visited = {node_off[start_node]}
    path = [start_node]
    current = start_node
    for a in anchors:
        if a in visited:
            continue
        tgt, prev = _dijkstra_to_station(adj, node_off, current, a)
        if tgt is None:
            return None, INF
        seg = []
        x = tgt
        while x != current:
            seg.append(x)
            x = prev[x]
        for nid in reversed(seg):
            path.append(nid)
            visited.add(node_off[nid])
        current = tgt
    if len(visited & canonical) < len(canonical):
        return None, INF
    return path, node_t[current] - node_t[start_node]


def _dijkstra_to_station_capped(adj, node_off, src, target_station, cap, runs=None):
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
        if u != src and node_off[u] == target_station:
            return u, prev
        for v, ed in adj[u].items():
            nd = d + ed["weight"]
            if nd < dist.get(v, INF):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
        if runs is not None:                         # out-of-system run shortcuts
            for v, w, _info in runs.run_successors(u):
                nd = d + w
                if nd < dist.get(v, INF):
                    dist[v] = nd
                    prev[v] = u
                    heapq.heappush(pq, (nd, v))
    return None, prev


def realize_from(adj, tables, current, visited0, anchors, cap=300000, runs=None):
    """Realize anchors starting from node ``current`` with ``visited0`` already
    covered. Returns (suffix_path_nodes, end_node, visited) or (None, None, None).
    With ``runs``, the realizer may use out-of-system run/walk shortcuts."""
    node_off = tables[0]
    visited = set(visited0)
    suffix = []
    for a in anchors:
        if a in visited:
            continue
        tgt, prev = _dijkstra_to_station_capped(adj, node_off, current, a, cap, runs)
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


def static_station_metric(G, node_off, extra_edges=None):
    """Static station-to-station best-case time matrix D (dict of dicts) and the
    ordered station list, from min edge weights in the time-expanded graph.
    ``extra_edges`` = optional [(a, b, w)] (e.g. terminal runs) folded in so the
    metric is run-aware."""
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
    D = dict(nx.all_pairs_dijkstra_path_length(SG))
    return D, list(SG.nodes())


def _nn_order(D, stations, start):
    """Nearest-neighbor tour over stations starting at `start`, using D."""
    remaining = set(stations)
    remaining.discard(start)
    order = [start]
    cur = start
    while remaining:
        nxt = min(remaining, key=lambda s: D[cur].get(s, INF))
        order.append(nxt)
        remaining.discard(nxt)
        cur = nxt
    return order


def _oropt_improve(D, order, passes=4):
    """Relocate single stations to cheaper positions (asymmetric Or-opt)."""
    def cost_at(a, b):
        return D[a].get(b, INF)
    n = len(order)
    for _ in range(passes):
        improved = False
        for i in range(1, n):
            s = order[i]
            prev, nxt = order[i - 1], order[i + 1] if i + 1 < n else None
            removed = cost_at(prev, s) + (cost_at(s, nxt) if nxt else 0) \
                - (cost_at(prev, nxt) if nxt else 0)
            best_delta, best_j = 0, None
            for j in range(n - 1):
                if j in (i - 1, i):
                    continue
                a, b = order[j], order[j + 1]
                added = cost_at(a, s) + cost_at(s, b) - cost_at(a, b)
                if added - removed < best_delta - 1e-9:
                    best_delta, best_j = added - removed, j
            if best_j is not None:
                order.pop(i)
                order.insert(best_j + 1 if best_j < i else best_j, s)
                improved = True
        if not improved:
            break
    return order


def _dijkstra_to_any(adj, node_off, src, targets):
    """Earliest-arrival from src to the nearest node whose station is in `targets`."""
    dist = {src: 0.0}
    prev = {}
    pq = [(0.0, src)]
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, INF):
            continue
        if u != src and node_off[u] in targets:
            return u, prev
        for v, ed in adj[u].items():
            nd = d + ed["weight"]
            if nd < dist.get(v, INF):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
    return None, prev


def realize_windowed(adj, tables, canonical, start_node, tsp_order, window):
    """Sweep tsp_order as a frontier: at each step go to the schedule-nearest of
    the next `window` unvisited stations in order. Coherent + schedule-aware."""
    node_off, _stop, node_t = tables
    visited = {node_off[start_node]}
    path = [start_node]
    current = start_node
    order = [s for s in tsp_order if s != node_off[start_node]]
    while len(visited & canonical) < len(canonical):
        win = []
        for s in order:
            if s not in visited:
                win.append(s)
                if len(win) >= window:
                    break
        if not win:
            break
        tgt, prev = _dijkstra_to_any(adj, node_off, current, set(win))
        if tgt is None:
            return None, INF
        seg = []
        x = tgt
        while x != current:
            seg.append(x)
            x = prev[x]
        for nid in reversed(seg):
            path.append(nid)
            visited.add(node_off[nid])
        current = tgt
    if len(visited & canonical) < len(canonical):
        return None, INF
    return path, node_t[current] - node_t[start_node]


def cmd_sweep(args) -> int:
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    adj = G._adj
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    canonical = si.canonical_stations

    D, _ = static_station_metric(G, node_off)
    windows = [int(x) for x in args.windows.split(",")]
    starts = args.starts.split(",")

    best = (INF, None, None)
    for st in starts:
        snode = min((n for n, d in G.nodes(data=True)
                     if d["stop"] == st and d["t"] >= args.after), key=lambda n: node_t[n])
        order = _oropt_improve(D, _nn_order(D, canonical, node_off[snode]), passes=args.passes)
        for w in windows:
            path, elapsed = realize_windowed(adj, tables, canonical, snode, order, w)
            if path:
                print(f"  start {st} W={w}: {hms(elapsed)}")
                if elapsed < best[0]:
                    best = (elapsed, path, (st, w))
    elapsed, path, key = best
    print(f"sweep best: {hms(elapsed)} from {key}")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": 5000, "notes": f"windowed sweep {key}"},
                               "path": [[node_stop[n], node_t[n]] for n in path]}))
    print(f"wrote {len(path)} nodes -> {out}")
    return 0


def cmd_tsp(args) -> int:
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    adj = G._adj
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    canonical = si.canonical_stations

    t0 = time.time()
    D, stations = static_station_metric(G, node_off)
    print(f"static metric: {len(stations)} stations, {time.time()-t0:.0f}s")

    snode = min((n for n, d in G.nodes(data=True)
                 if d["stop"] == args.start and d["t"] >= args.after),
                key=lambda n: node_t[n])
    start_station = node_off[snode]

    order = _nn_order(D, canonical, start_station)
    static_cost = sum(D[order[i]].get(order[i + 1], INF) for i in range(len(order) - 1))
    print(f"NN order static cost {hms(static_cost)}")
    order = _oropt_improve(D, order, passes=args.passes)
    static_cost = sum(D[order[i]].get(order[i + 1], INF) for i in range(len(order) - 1))
    print(f"Or-opt order static cost {hms(static_cost)}")

    path, elapsed = realize(adj, tables, canonical, snode, order)
    if path is None:
        print("realize failed (incomplete coverage)")
        return 1
    print(f"realized makespan: {hms(elapsed)} from {args.start}, {time.time()-t0:.0f}s total")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": 5000, "notes": f"tsp order from {args.start}"},
                               "path": [[node_stop[n], node_t[n]] for n in path]}))
    print(f"wrote {len(path)} nodes -> {out}")
    return 0


def tail_sa(adj, tables, canonical, snode, split, maxseg, budget, t0_temp, tend_temp, seed):
    """Greedy from snode, then simulated-annealing or-opt on the tail anchors.
    Returns (full_path_nodes, elapsed) or (None, INF)."""
    import math
    import random
    node_off, _stop, node_t = tables
    path, anchors, elapsed = greedy_anchors(adj, tables, canonical, snode)
    if path is None:
        return None, INF
    start_t = node_t[snode]
    p = int(len(anchors) * split)
    prefix_anchors, tail = anchors[:p], anchors[p:]
    pre_suffix, pre_end, pre_visited = realize_from(adj, tables, snode, {node_off[snode]}, prefix_anchors)
    prefix_path = [snode] + pre_suffix
    base_tail, end, _ = realize_from(adj, tables, pre_end, pre_visited, tail)
    best_order, best_elapsed, best_tail_path = tail[:], node_t[end] - start_t, base_tail

    rng = random.Random(seed)
    mt = len(tail)
    cur_order, cur_elapsed = best_order[:], best_elapsed
    t0 = time.time()
    while time.time() - t0 < budget and mt > 2:
        frac = (time.time() - t0) / budget
        T = t0_temp * (tend_temp / t0_temp) ** frac
        cand = cur_order[:]
        if rng.random() < 0.5:
            i = rng.randrange(mt)
            cand.insert(rng.randrange(len(cand) + 1), cand.pop(i))
        else:
            i = rng.randrange(mt - 1)
            j = min(i + 2 + rng.randrange(maxseg), mt)
            cand[i:j] = cand[i:j][::-1]
        tp, tend, tvis = realize_from(adj, tables, pre_end, pre_visited, cand)
        if tp is None or len(tvis & canonical) < len(canonical):
            continue
        e = node_t[tend] - start_t
        delta = e - cur_elapsed
        if delta < 0 or rng.random() < math.exp(-delta / T):
            cur_order, cur_elapsed = cand, e
            if e < best_elapsed:
                best_order, best_elapsed, best_tail_path = cand, e, tp
    return prefix_path + best_tail_path, best_elapsed


def cmd_optimize_tail(args) -> int:
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    snode = min((n for n, d in G.nodes(data=True)
                 if d["stop"] == args.start and d["t"] >= args.after), key=lambda n: node_t[n])
    full, elapsed = tail_sa(G._adj, tables, si.canonical_stations, snode,
                            args.split, args.maxseg, args.time_budget, args.t0, args.tend, args.seed)
    print(f"tail-SA {args.start}: {hms(elapsed)}")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": 5000, "notes": "tail or-opt local search"},
                               "path": [[node_stop[n], node_t[n]] for n in full]}))
    print(f"wrote {len(full)} nodes -> {out}")
    return 0


def realize_sequence(adj, tables, canonical, start_node, seq):
    """Realize a station visit sequence on the schedule: at each step go to the
    next station in seq via earliest arrival. Returns (path, elapsed)."""
    node_off, _stop, node_t = tables
    current = start_node
    path = [start_node]
    visited = {node_off[start_node]}
    for s in seq:
        if node_off[current] == s:
            continue
        tgt, prev = _dijkstra_to_station(adj, node_off, current, s)
        if tgt is None:
            return None, INF
        seg = []
        x = tgt
        while x != current:
            seg.append(x)
            x = prev[x]
        for nid in reversed(seg):
            path.append(nid)
            visited.add(node_off[nid])
        current = tgt
    if len(visited & canonical) < len(canonical):
        return None, INF
    return path, node_t[current] - node_t[start_node]


def regret_insert(order, removed, D, rng=None, jitter=0.0):
    """Insert each station in `removed` into `order` at its best position by the
    static metric D (highest-regret first). With `rng` and `jitter>0`, sometimes
    insert at the 2nd-best position to diversify (lets SA escape local optima)."""
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


def lns_run(adj, tables, canonical, D, snode, split, wmin, wmax, budget, t0_temp,
            tend_temp, seed, jitter=0.15, anchors=None, runs=None, coords=None,
            cluster_frac=0.4):
    """One LNS run (ruin + randomized-regret recreate + SA) from snode. If
    `anchors` is given, seed from them (iterated LNS); else from greedy. With
    `runs`, the realizer may use out-of-system run/walk shortcuts. Returns
    (full_path, elapsed)."""
    import math
    import random
    from .walk_transfers import _haversine_m
    node_off, _stop, node_t = tables
    start_t = node_t[snode]
    if anchors is None:
        _, anchors, _ = greedy_anchors(adj, tables, canonical, snode)
    p = int(len(anchors) * split)
    prefix_anchors, tail = anchors[:p], anchors[p:]
    pre_suffix, pre_end, pre_visited = realize_from(adj, tables, snode, {node_off[snode]}, prefix_anchors, runs=runs)
    prefix_path = [snode] + pre_suffix
    base_path, end, _ = realize_from(adj, tables, pre_end, pre_visited, tail, runs=runs)
    cur_tail, cur_e = tail[:], node_t[end] - start_t
    best_tail, best_e, best_path = tail[:], cur_e, base_path
    rng = random.Random(seed)
    t0 = time.time()
    while time.time() - t0 < budget and len(tail) > 4:
        frac = (time.time() - t0) / budget
        T = t0_temp * (tend_temp / t0_temp) ** frac
        if coords is not None and rng.random() < cluster_frac:
            # geographic-cluster ruin: tear out all tail anchors near a random center
            center = cur_tail[rng.randrange(len(cur_tail))]
            if center not in coords:
                continue
            clat, clon = coords[center]
            rad = rng.uniform(1000, 3000)
            removed = [s for s in cur_tail if s in coords
                       and _haversine_m(coords[s][0], coords[s][1], clat, clon) <= rad]
            if len(removed) > 3 * wmax:
                removed = rng.sample(removed, 3 * wmax)
            rset = set(removed)
            kept = [s for s in cur_tail if s not in rset]
        else:
            W = rng.randint(wmin, wmax)
            i = rng.randrange(max(1, len(cur_tail) - W))
            removed = cur_tail[i:i + W]
            kept = cur_tail[:i] + cur_tail[i + W:]
        if not removed or not kept:
            continue
        cand = regret_insert(kept, removed, D, rng, jitter)
        tp, tend, tvis = realize_from(adj, tables, pre_end, pre_visited, cand, runs=runs)
        if tp is None or len(tvis & canonical) < len(canonical):
            continue
        e = node_t[tend] - start_t
        if e - cur_e < 0 or rng.random() < math.exp(-(e - cur_e) / T):
            cur_tail, cur_e = cand, e
            if e < best_e:
                best_tail, best_e, best_path = cand, e, tp
    return prefix_path + best_path, best_e


def cmd_lns(args) -> int:
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables

    runs = None
    extra_edges = []
    if args.use_runs or args.terminal_runs:
        from .run_layer import RunLayer
        from .walk_transfers import load_complexes
        runs = RunLayer.from_graph(G, radius_m=args.run_radius)
        cid2off = {cid: {si.resolve(p) for p in c.parents}
                   for cid, c in load_complexes().items()}
        if args.terminal_runs:
            # restrict runs to originate only at dead-end terminals (deg-1 in the
            # track graph) -- ride a line to its end, run to a nearby line.
            import collections
            deg = collections.defaultdict(set)
            for u, nbrs in G._adj.items():
                a = node_off[u]
                for v, ed in nbrs.items():
                    if ed["mode"] == "train" and node_off[v] != a:
                        deg[a].add(node_off[v])
            terminals = {s for s, nb in deg.items() if len(nb) == 1}
            terminals |= {si.resolve(s) for s in _terminal_stops()}   # service line termini
            off2cid = {s: cid for cid, ss in cid2off.items() for s in ss}
            tcids = {off2cid[s] for s in terminals if s in off2cid}
            runs.adjacency = {cid: nb for cid, nb in runs.adjacency.items() if cid in tcids}
            print(f"terminal runs: {len(terminals)} terminals -> "
                  f"{len(runs.adjacency)} run-source complexes (radius {args.run_radius}m)")
        else:
            print(f"runs enabled (radius {args.run_radius}m)")
        # make the static metric run-aware so recreate proposes run-friendly orders
        for cid, nbrs in runs.adjacency.items():
            for ocid, secs, _m in nbrs:
                for a in cid2off.get(cid, ()):
                    for b in cid2off.get(ocid, ()):
                        extra_edges.append((a, b, secs))

    D, _ = static_station_metric(G, node_off, extra_edges)

    coords = None                                    # official station -> (lat, lon)
    if args.cluster_ruin:
        import pandas as pd
        sdf = pd.read_csv(STATIONS_CSV := "data/official/mta_subway_stations.csv", dtype=str)
        coords = {r["Station ID"]: (float(r["GTFS Latitude"]), float(r["GTFS Longitude"]))
                  for _, r in sdf.iterrows()}
        print("geographic-cluster ruin enabled")

    splits = [float(x) for x in args.splits.split(",")]
    seed_anchors = None
    if args.seed_from:                                # iterated LNS: seed from a solution
        data = json.loads(Path(args.seed_from).read_text())
        spath = data["path"]
        s0, t0 = spath[0][0], int(spath[0][1])
        starts = [s0]
        seen, seed_anchors = set(), []
        for stop, _t in spath:
            st = si.resolve(stop)
            if st not in seen:
                seen.add(st)
                seed_anchors.append(st)
        print(f"iterated LNS: seeded from {args.seed_from} ({len(seed_anchors)} anchors)")
    else:
        starts = args.start.split(",")
    configs = [(st, sp, sd) for st in starts for sp in splits for sd in range(args.seeds)]
    per = args.time_budget / len(configs)
    print(f"LNS sweep: {len(configs)} configs x {per:.0f}s each")

    best = (INF, None, None)
    for st, sp, sd in configs:
        cands = [n for n, d in G.nodes(data=True) if d["stop"] == st and d["t"] >= args.after]
        if not cands:
            continue
        snode = min(cands, key=lambda n: node_t[n])
        full, e = lns_run(G._adj, tables, si.canonical_stations, D, snode, sp,
                          args.wmin, args.wmax, per, args.t0, args.tend, sd,
                          anchors=seed_anchors, runs=runs, coords=coords)
        if full and e < best[0]:
            best = (e, full, (st, sp, sd))
            print(f"  new best {hms(e)} from start={st} split={sp} seed={sd}")
    e, full, key = best
    print(f"LNS sweep best: {hms(e)} {key}")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": 5000, "notes": f"LNS sweep {key}"},
                               "path": [[node_stop[n], node_t[n]] for n in full]}))
    print(f"wrote {len(full)} nodes -> {out}")
    return 0


def cmd_postman(args) -> int:
    import networkx as nx
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    adj = G._adj
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    canonical = si.canonical_stations

    # undirected required-edge (track) graph with best-case hop times
    Greq = nx.Graph()
    for u, nbrs in G._adj.items():
        a = node_off[u]
        for v, ed in nbrs.items():
            if ed["mode"] != "train" or node_off[v] == a:
                continue
            b, w = node_off[v], ed["weight"]
            if not Greq.has_edge(a, b) or w < Greq[a][b]["weight"]:
                Greq.add_edge(a, b, weight=w)

    comps = sorted(nx.connected_components(Greq), key=len, reverse=True)
    print(f"track components: {[len(c) for c in comps]}")
    seqs = []
    for comp in comps:
        Gc = Greq.subgraph(comp).copy()
        Hc = nx.eulerize(Gc)                         # Chinese-Postman within component
        start = min(comp, key=lambda s: (Gc.degree(s) != 1, s))  # prefer a terminus
        circuit = list(nx.eulerian_circuit(Hc, source=start))
        seqs.append([circuit[0][0]] + [v for _, v in circuit])
    full_seq = [s for seq in seqs for s in seq]      # components concatenated (realize bridges)

    snode = min((n for n, d in G.nodes(data=True)
                 if node_off[n] == full_seq[0] and d["t"] >= args.after),
                key=lambda n: node_t[n])
    path, elapsed = realize_sequence(adj, tables, canonical, snode, full_seq)
    if path is None:
        print("realize failed")
        return 1
    print(f"postman realized: {hms(elapsed)} from {node_stop[snode]}")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": 5000, "notes": "postman (eulerized) route"},
                               "path": [[node_stop[n], node_t[n]] for n in path]}))
    print(f"wrote {len(path)} nodes -> {out}")
    return 0


def cmd_portfolio(args) -> int:
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    adj = G._adj
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    canonical = si.canonical_stations
    starts = args.starts.split(",")
    per = args.time_budget / len(starts)
    best = (INF, None, None)
    for st in starts:
        snode = min((n for n, d in G.nodes(data=True)
                     if d["stop"] == st and d["t"] >= args.after), key=lambda n: node_t[n])
        full, elapsed = tail_sa(adj, tables, canonical, snode, args.split, args.maxseg,
                                per, args.t0, args.tend, args.seed)
        print(f"  {st}: {hms(elapsed)}")
        if full and elapsed < best[0]:
            best = (elapsed, full, st)
    elapsed, full, st = best
    print(f"portfolio best: {hms(elapsed)} from {st}")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": 5000, "notes": f"portfolio tail-SA {st}"},
                               "path": [[node_stop[n], node_t[n]] for n in full]}))
    print(f"wrote {len(full)} nodes -> {out}")
    return 0


def cmd_optimize(args) -> int:
    import random
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    adj = G._adj
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    canonical = si.canonical_stations

    # initial route + anchors from the best terminal start
    snode = min((n for n, d in G.nodes(data=True)
                 if d["stop"] == args.start and d["t"] >= args.after),
                key=lambda n: node_t[n])
    path, anchors, elapsed = greedy_anchors(adj, tables, canonical, snode)
    print(f"initial: {hms(elapsed)} from {args.start}, {len(anchors)} anchors")

    rng = random.Random(args.seed)
    best_anchors, best_elapsed, best_path = anchors[:], elapsed, path
    m = len(anchors)
    t0 = time.time()
    tries = accepts = 0
    while time.time() - t0 < args.time_budget:
        tries += 1
        cand = best_anchors[:]
        if rng.random() < 0.5:                       # relocate one anchor (or-opt)
            i = rng.randrange(m)
            x = cand.pop(i)
            cand.insert(rng.randrange(len(cand) + 1), x)
        else:                                        # reverse a short segment (2-opt)
            i = rng.randrange(m - 1)
            j = min(i + 1 + rng.randrange(args.maxseg), m)
            cand[i:j] = cand[i:j][::-1]
        p, e = realize(adj, tables, canonical, snode, cand)
        if p is not None and e < best_elapsed:
            best_anchors, best_elapsed, best_path = cand, e, p
            accepts += 1
            print(f"  {time.time()-t0:5.0f}s try {tries}: improved -> {hms(e)}")
    print(f"optimize done: {hms(best_elapsed)} ({accepts} accepts / {tries} tries)")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": 5000, "notes": "or-opt local search"},
                               "path": [[node_stop[n], node_t[n]] for n in best_path]}))
    print(f"wrote {len(best_path)} nodes -> {out}")
    return 0


def greedy_construct(G, si, start_node, log_every=50):
    tables = _node_tables(G, si)
    t0 = time.time()
    path, elapsed = greedy_from(G._adj, tables, si.canonical_stations, start_node)
    print(f"done: route time {hms(elapsed)}, {time.time()-t0:.0f}s wall")
    node_stop, node_t = tables[1], tables[2]
    return [[node_stop[n], node_t[n]] for n in path]


def _dijkstra_k_unvisited(adj, node_off, visited, src, k):
    """Earliest-arrival search returning up to k nearest nodes at *distinct*
    unvisited stations, plus the prev map for reconstruction."""
    dist = {src: 0.0}
    prev = {}
    pq = [(0.0, src)]
    found = []
    found_stations = set()
    while pq and len(found) < k:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, INF):
            continue
        s = node_off[u]
        if u != src and s not in visited and s not in found_stations:
            found.append(u)
            found_stations.add(s)
        for v, ed in adj[u].items():
            nd = d + ed["weight"]
            if nd < dist.get(v, INF):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
    return found, prev


def greedy_rand(adj, tables, canonical, start_node, k, rng):
    """Randomized greedy: at each step pick among the k nearest unvisited
    stations, biased toward the nearest. Returns (path, elapsed_s)."""
    node_off, _stop, node_t = tables
    visited = {node_off[start_node]}
    path = [start_node]
    current = start_node
    while len(visited & canonical) < len(canonical):
        found, prev = _dijkstra_k_unvisited(adj, node_off, visited, current, k)
        if not found:
            return None, INF
        # weight choices toward the nearest (rank 0 heaviest)
        n = len(found)
        weights = list(range(n, 0, -1))
        tgt = rng.choices(found, weights=weights, k=1)[0]
        seg = []
        x = tgt
        while x != current:
            seg.append(x)
            x = prev[x]
        for nid in reversed(seg):
            path.append(nid)
            visited.add(node_off[nid])
        current = tgt
    return path, node_t[current] - node_t[start_node]


def cmd_grasp(args) -> int:
    import random
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    adj = G._adj
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    canonical = si.canonical_stations

    # seed starts: best from a multi pass (top terminals at one time)
    times = [21600]
    terminals = _terminal_stops() & set(node_stop)
    starts = {}
    for nid, d in G.nodes(data=True):
        if d["stop"] in terminals and d["t"] >= times[0] \
                and any(e["mode"] == "train" for e in adj[nid].values()):
            key = d["stop"]
            if key not in starts or d["t"] < node_t[starts[key]]:
                starts[key] = nid
    start_nodes = list(starts.values())

    rng = random.Random(args.seed)
    best_path, best_elapsed, best_src = None, INF, None
    # warm start from current best.json if present
    t0 = time.time()
    for it in range(args.iters):
        snode = rng.choice(start_nodes)
        path, elapsed = greedy_rand(adj, tables, canonical, snode, args.k, rng)
        if path and elapsed < best_elapsed:
            best_path, best_elapsed, best_src = path, elapsed, node_stop[snode]
            print(f"  it {it}: new best {hms(elapsed)} from {best_src}")
        if time.time() - t0 > args.time_budget:
            print(f"  time budget reached at it {it}")
            break
    print(f"GRASP best: {hms(best_elapsed)} from {best_src}, {time.time()-t0:.0f}s")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": 5000,
                                        "notes": f"grasp k={args.k} from {best_src}"},
                               "path": [[node_stop[n], node_t[n]] for n in best_path]}))
    print(f"wrote {len(best_path)} nodes -> {out}")
    return 0


def _terminal_stops(gtfs_dir="data/gtfs"):
    """Platform stop ids that are the first or last stop of some trip (terminals
    -- natural route start points that minimize backtracking)."""
    import pandas as pd
    st = pd.read_csv(Path(gtfs_dir) / "stop_times.txt", dtype=str)
    st["stop_sequence"] = st["stop_sequence"].astype(int)
    g = st.sort_values("stop_sequence").groupby("trip_id")["stop_id"]
    return set(g.first()) | set(g.last())


def cmd_multi(args) -> int:
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    adj = G._adj
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    canonical = si.canonical_stations

    # candidate start nodes: each terminal platform, at each candidate start time
    times = [int(x) for x in args.times.split(",")]
    terminals = _terminal_stops() & set(node_stop)
    # earliest node at each (terminal, time) with a train out-edge
    starts = {}
    for nid, d in G.nodes(data=True):
        stop = d["stop"]
        if stop not in terminals:
            continue
        for ct in times:
            if d["t"] >= ct and any(e["mode"] == "train" for e in adj[nid].values()):
                key = (stop, ct)
                if key not in starts or d["t"] < node_t[starts[key]]:
                    starts[key] = nid
    print(f"trying {len(starts)} starts ({len(terminals)} terminals x {len(times)} times)")

    t0 = time.time()
    best_path, best_elapsed, best_key = None, INF, None
    for i, (key, snode) in enumerate(sorted(starts.items())):
        path, elapsed = greedy_from(adj, tables, canonical, snode)
        if path and elapsed < best_elapsed:
            best_path, best_elapsed, best_key = path, elapsed, key
            print(f"  [{i+1}/{len(starts)}] new best {hms(elapsed)} from {key[0]}@{key[1]}")
    print(f"best: {hms(best_elapsed)} from {best_key}, {time.time()-t0:.0f}s wall")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": 5000,
                                        "notes": f"multi-greedy best from {best_key[0]}@{best_key[1]}"},
                               "path": [[node_stop[n], node_t[n]] for n in best_path]}))
    print(f"wrote {len(best_path)} nodes -> {out}")
    return 0


def cmd_greedy(args) -> int:
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    start = _pick_start(G, args.start, args.after)
    print(f"start: {args.start} @ t={args.after}+ (node {start})")
    path = greedy_construct(G, si, start)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": 5000, "notes": f"greedy from {args.start}"},
                               "path": path}))
    print(f"wrote {len(path)} nodes -> {out}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Construct Subway Challenge routes.")
    sub = p.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("greedy", help="Greedy nearest-unvisited construction.")
    g.add_argument("--start", default="101S", help="Start platform stop id.")
    g.add_argument("--after", type=int, default=5 * 3600, help="Earliest start t (sec in week).")
    g.add_argument("--out", default="solutions/greedy.json")
    g.set_defaults(func=cmd_greedy)
    m = sub.add_parser("multi", help="Multi-start/multi-time greedy; keep the best.")
    m.add_argument("--times", default="16200,18000,21600,25200",
                   help="Comma start times in sec (default 04:30,05:00,06:00,07:00).")
    m.add_argument("--out", default="solutions/multi.json")
    m.set_defaults(func=cmd_multi)
    gr = sub.add_parser("grasp", help="Randomized restarts (GRASP); keep the best.")
    gr.add_argument("--k", type=int, default=3, help="Pick among k nearest unvisited.")
    gr.add_argument("--iters", type=int, default=300)
    gr.add_argument("--time-budget", type=float, default=180, help="Seconds.")
    gr.add_argument("--seed", type=int, default=0)
    gr.add_argument("--out", default="solutions/grasp.json")
    gr.set_defaults(func=cmd_grasp)
    o = sub.add_parser("optimize", help="Or-opt/2-opt local search on anchor order.")
    o.add_argument("--start", default="A03S")
    o.add_argument("--after", type=int, default=21600)
    o.add_argument("--maxseg", type=int, default=8, help="Max segment length for reversals.")
    o.add_argument("--time-budget", type=float, default=180)
    o.add_argument("--seed", type=int, default=0)
    o.add_argument("--out", default="solutions/optimized.json")
    o.set_defaults(func=cmd_optimize)
    ot = sub.add_parser("optimize-tail", help="Freeze prefix; local-search the tail anchors.")
    ot.add_argument("--start", default="A03S")
    ot.add_argument("--after", type=int, default=21600)
    ot.add_argument("--split", type=float, default=0.6, help="Fraction of anchors to freeze.")
    ot.add_argument("--maxseg", type=int, default=6)
    ot.add_argument("--time-budget", type=float, default=180)
    ot.add_argument("--t0", type=float, default=600, help="SA start temperature (sec).")
    ot.add_argument("--tend", type=float, default=20, help="SA end temperature (sec).")
    ot.add_argument("--seed", type=int, default=0)
    ot.add_argument("--out", default="solutions/tail.json")
    ot.set_defaults(func=cmd_optimize_tail)
    ts = sub.add_parser("tsp", help="Global TSP-order construction (static metric) + realize.")
    ts.add_argument("--start", default="A03S")
    ts.add_argument("--after", type=int, default=21600)
    ts.add_argument("--passes", type=int, default=4)
    ts.add_argument("--out", default="solutions/tsp.json")
    ts.set_defaults(func=cmd_tsp)
    sw = sub.add_parser("sweep", help="Windowed schedule-aware sweep along TSP order.")
    sw.add_argument("--starts", default="A03S,A02S,257N,701S,257S")
    sw.add_argument("--after", type=int, default=21600)
    sw.add_argument("--windows", default="3,6,10,16,25")
    sw.add_argument("--passes", type=int, default=5)
    sw.add_argument("--out", default="solutions/sweep.json")
    sw.set_defaults(func=cmd_sweep)
    pf = sub.add_parser("portfolio", help="Tail-SA from several starts; keep the best.")
    pf.add_argument("--starts", default="A03S,A02S,257N,257S,101S,701S")
    pf.add_argument("--after", type=int, default=21600)
    pf.add_argument("--split", type=float, default=0.5)
    pf.add_argument("--maxseg", type=int, default=8)
    pf.add_argument("--time-budget", type=float, default=300)
    pf.add_argument("--t0", type=float, default=800)
    pf.add_argument("--tend", type=float, default=15)
    pf.add_argument("--seed", type=int, default=1)
    pf.add_argument("--out", default="solutions/portfolio.json")
    pf.set_defaults(func=cmd_portfolio)
    pm = sub.add_parser("postman", help="Rural-Postman (eulerized per component) construction.")
    pm.add_argument("--after", type=int, default=21600)
    pm.add_argument("--out", default="solutions/postman.json")
    pm.set_defaults(func=cmd_postman)
    ln = sub.add_parser("lns", help="Time-dependent LNS sweep: ruin + regret-recreate + SA.")
    ln.add_argument("--start", default="A02S,A03S,257N,257S,101S", help="Comma start stops.")
    ln.add_argument("--splits", default="0.3,0.45,0.6", help="Comma freeze fractions.")
    ln.add_argument("--seeds", type=int, default=2, help="Seeds per (start,split).")
    ln.add_argument("--after", type=int, default=21600)
    ln.add_argument("--seed-from", default=None, help="Seed anchors from a solution JSON (iterated LNS).")
    ln.add_argument("--use-runs", action="store_true", help="Allow run shortcuts from ALL stations.")
    ln.add_argument("--terminal-runs", action="store_true", help="Allow runs only FROM dead-end terminals (high-value).")
    ln.add_argument("--cluster-ruin", action="store_true", help="Also ruin whole geographic regions, not just order slices.")
    ln.add_argument("--run-radius", type=float, default=1500, help="Max run distance (m) for runs.")
    ln.add_argument("--wmin", type=int, default=3)
    ln.add_argument("--wmax", type=int, default=18)
    ln.add_argument("--time-budget", type=float, default=600)
    ln.add_argument("--t0", type=float, default=600)
    ln.add_argument("--tend", type=float, default=10)
    ln.add_argument("--out", default="solutions/lns.json")
    ln.set_defaults(func=cmd_lns)
    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
