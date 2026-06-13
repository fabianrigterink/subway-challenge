"""Route optimizer for the Subway Challenge: time-dependent Large-Neighborhood
Search (LNS) with out-of-system "runs" between dead-end terminals.

The route is represented as an ordered list of target stations (anchors).
``realize_from`` prices an anchor order on the real timetable (earliest-arrival
per leg, optionally using runs). LNS then repeatedly **ruins** a window of the
order and **recreates** it with regret-2 insertion, accepting via simulated
annealing. The recreate is guided by a static station-to-station metric that is
made *run-aware* (terminal runs folded in). Iterate with ``--seed-from``.

This is the method that produced ``solutions/best.json`` (24:45:00). Run:

    python -m subway_challenge.search lns --seed-from solutions/best.json --terminal-runs
"""
from __future__ import annotations

import argparse
import collections
import heapq
import json
import math
import pickle
import random
import sys
import time
from pathlib import Path

from .solver import GRAPH_PKL, hms
from .stations import StationIndex

INF = float("inf")


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
    return path, anchors, node_t[current] - node_t[start_node]


def realize_from(adj, tables, current, visited0, anchors, cap=300000, runs=None):
    """Realize an anchor order on the schedule from node ``current`` (with
    ``visited0`` already covered). Anchors covered incidentally are skipped.
    Returns (suffix_path_nodes, end_node, visited) or (None, None, None)."""
    node_off = tables[0]
    visited = set(visited0)
    suffix = []
    for a in anchors:
        if a in visited:
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
            tend_temp, seed, jitter, anchors, runs):
    """One LNS run: freeze a prefix of the order, then ruin+regret-recreate+SA on
    the tail. Seeds from ``anchors`` if given, else greedy. Returns (path, elapsed)."""
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
        T = t0_temp * (tend_temp / t0_temp) ** ((time.time() - t0) / budget)
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


def cmd_lns(args) -> int:
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables

    runs, extra_edges = (None, [])
    if args.terminal_runs:
        runs, extra_edges = _terminal_run_layer(G, si, node_off, args.run_radius)
    D, _ = static_station_metric(G, node_off, extra_edges)

    splits = [float(x) for x in args.splits.split(",")]
    seed_anchors = None
    if args.seed_from:                                # iterated LNS: seed from a solution
        spath = json.loads(Path(args.seed_from).read_text())["path"]
        starts = [spath[0][0]]
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
        full, e = lns_run(G._adj, tables, si.canonical_stations, D, snode, sp, args.wmin,
                          args.wmax, per, args.t0, args.tend, sd, args.jitter, seed_anchors, runs)
        if full and e < best[0]:
            best = (e, full, (st, sp, sd))
            print(f"  new best {hms(e)} from start={st} split={sp} seed={sd}")

    e, full, key = best
    print(f"LNS sweep best: {hms(e)} {key}")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": {"radius_m": 5000, "notes": f"LNS {key}"},
                               "path": [[node_stop[n], node_t[n]] for n in full]}))
    print(f"wrote {len(full)} nodes -> {out}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Subway Challenge optimizer (LNS + terminal runs).")
    sub = p.add_subparsers(dest="cmd", required=True)
    ln = sub.add_parser("lns", help="Time-dependent LNS: ruin + regret recreate + SA.")
    ln.add_argument("--start", default="A02S", help="Comma start stops (greedy seed).")
    ln.add_argument("--seed-from", default=None, help="Seed anchor order from a solution JSON (iterate).")
    ln.add_argument("--splits", default="0.1,0.15,0.2,0.3", help="Comma freeze fractions.")
    ln.add_argument("--seeds", type=int, default=4, help="SA seeds per (start, split).")
    ln.add_argument("--after", type=int, default=21600, help="Earliest start t (sec in week).")
    ln.add_argument("--terminal-runs", action="store_true", help="Enable dead-end terminal runs.")
    ln.add_argument("--run-radius", type=float, default=2500, help="Run layer radius (m).")
    ln.add_argument("--wmin", type=int, default=3)
    ln.add_argument("--wmax", type=int, default=22)
    ln.add_argument("--jitter", type=float, default=0.15)
    ln.add_argument("--time-budget", type=float, default=600)
    ln.add_argument("--t0", type=float, default=600, help="SA start temperature (sec).")
    ln.add_argument("--tend", type=float, default=10, help="SA end temperature (sec).")
    ln.add_argument("--out", default="solutions/lns.json")
    ln.set_defaults(func=cmd_lns)
    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
