"""Seed-column extraction for column-generation experiments.

Columns are schedule-realized path slices cut out of existing route JSONs. They
are not an optimizer by themselves; they are the initial restricted master input
for later set-covering/path formulations.
"""
from __future__ import annotations

import argparse
import bisect
import heapq
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

from ortools.linear_solver import pywraplp
from ortools.sat.python import cp_model
from ortools.constraint_solver import pywrapcp, routing_enums_pb2

from .build_graph import WEEK
from .run_layer import RunLayer
from .search import (
    _elapsed,
    _dijkstra_to_station,
    _first_event_at_or_after,
    _node_tables,
    _stop_index,
    _terminal_stops,
    _terminal_run_layer,
    realize_from,
    static_station_metric,
)
from .solver import GRAPH_PKL, Problem, hms, load_problem, result_line, validate
from .stations import StationIndex

BIG_COST = 10**9
INF = float("inf")
DEFAULT_EXACT_ARC_CACHE = Path("reports/optimization_runs/exact_connector_cache.jsonl")


class _NoRunLayer:
    platform_index = {}

    def run_successors(self, _node):
        return ()


def _load_route(path: Path):
    data = json.loads(path.read_text())
    route = data.get("path", data)
    meta = data.get("meta", {})
    return route, meta


def _nodes_for_path(prob, path):
    nodes = []
    for stop, t, *_rest in path:
        node = prob.node_of(stop, int(t))
        if node is None:
            raise ValueError(f"node {stop}@{int(t)} not in graph")
        nodes.append(node)
    return nodes


def _station(prob, node):
    return prob.stations.resolve(prob.G.nodes[node]["stop"])


def _stop_time(prob, node):
    d = prob.G.nodes[node]
    return d["stop"], int(d["t"])


def _transition_series(prob, nodes):
    series = []
    for u, v in zip(nodes, nodes[1:]):
        weight, mode = prob.transition(u, v)
        if weight is None:
            us, ut = _stop_time(prob, u)
            vs, vt = _stop_time(prob, v)
            raise ValueError(f"illegal transition {us}@{ut} -> {vs}@{vt}")
        series.append((int(weight), mode))
    return series


def _split_points(prob, nodes, transitions, max_new_stations, max_elapsed_s, split_on_run):
    points = [0]
    seen_in_segment = {_station(prob, nodes[0])}
    elapsed = 0
    for i, (weight, mode) in enumerate(transitions, start=1):
        elapsed += weight
        seen_in_segment.add(_station(prob, nodes[i]))
        should_split = False
        if split_on_run and mode == "run":
            should_split = True
        if max_new_stations and len(seen_in_segment) >= max_new_stations:
            should_split = True
        if max_elapsed_s and elapsed >= max_elapsed_s:
            should_split = True
        if should_split and i < len(nodes) - 1:
            if i != points[-1]:
                points.append(i)
            seen_in_segment = {_station(prob, nodes[i])}
            elapsed = 0
    if points[-1] != len(nodes) - 1:
        points.append(len(nodes) - 1)
    return points


def _column_record(prob, source, column_idx, path, nodes, transitions, a, b, include_path):
    seg_nodes = nodes[a:b + 1]
    seg_transitions = transitions[a:b]
    start_stop, start_t = _stop_time(prob, seg_nodes[0])
    end_stop, end_t = _stop_time(prob, seg_nodes[-1])
    covered = {_station(prob, n) for n in seg_nodes} & prob.canonical
    modes = Counter(mode for _weight, mode in seg_transitions)
    elapsed = sum(weight for weight, _mode in seg_transitions)
    record = {
        "column_id": f"{source.stem}:{column_idx:03d}",
        "source": str(source),
        "start": {
            "stop": start_stop,
            "station": _station(prob, seg_nodes[0]),
            "time": start_t,
        },
        "end": {
            "stop": end_stop,
            "station": _station(prob, seg_nodes[-1]),
            "time": end_t,
        },
        "elapsed_s": int(elapsed),
        "elapsed": hms(elapsed),
        "steps": len(seg_nodes),
        "modes": dict(modes),
        "covered_count": len(covered),
        "covered_stations": sorted(covered),
    }
    if include_path:
        record["path"] = path[a:b + 1]
    return record


def _column_from_nodes(prob, source_label, column_idx, nodes, include_path=True):
    path = [[prob.G.nodes[n]["stop"], int(prob.G.nodes[n]["t"])] for n in nodes]
    transitions = _transition_series(prob, nodes)
    start_stop, start_t = _stop_time(prob, nodes[0])
    end_stop, end_t = _stop_time(prob, nodes[-1])
    covered = {_station(prob, n) for n in nodes} & prob.canonical
    modes = Counter(mode for _weight, mode in transitions)
    elapsed = sum(weight for weight, _mode in transitions)
    record = {
        "column_id": f"{source_label}:{column_idx:05d}",
        "source": source_label,
        "start": {
            "stop": start_stop,
            "station": _station(prob, nodes[0]),
            "time": start_t,
        },
        "end": {
            "stop": end_stop,
            "station": _station(prob, nodes[-1]),
            "time": end_t,
        },
        "elapsed_s": int(elapsed),
        "elapsed": hms(elapsed),
        "steps": len(nodes),
        "modes": dict(modes),
        "covered_count": len(covered),
        "covered_stations": sorted(covered),
    }
    if include_path:
        record["path"] = path
    return record


def extract_columns(prob, route_path: Path, max_new_stations: int,
                    max_elapsed_s: int, split_on_run: bool, include_path: bool):
    path, _meta = _load_route(route_path)
    nodes = _nodes_for_path(prob, path)
    transitions = _transition_series(prob, nodes)
    points = _split_points(
        prob, nodes, transitions, max_new_stations, max_elapsed_s, split_on_run)
    records = []
    for j, (a, b) in enumerate(zip(points, points[1:]), start=1):
        if b <= a:
            continue
        records.append(_column_record(
            prob, route_path, j, path, nodes, transitions, a, b, include_path))
    return records


def cmd_extract(args) -> int:
    route_paths = [Path(p) for p in args.routes]
    if not route_paths:
        route_paths = sorted(Path("solutions").glob("*.json"))
    prob = load_problem(radius_m=args.radius)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    skipped = []
    with out.open("w") as f:
        for route_path in route_paths:
            try:
                records = extract_columns(
                    prob,
                    route_path,
                    max_new_stations=args.max_new_stations,
                    max_elapsed_s=int(args.max_segment_minutes * 60),
                    split_on_run=args.split_on_run,
                    include_path=not args.no_path,
                )
            except Exception as exc:  # keep a portfolio run moving
                skipped.append((str(route_path), str(exc)))
                continue
            for record in records:
                f.write(json.dumps(record, sort_keys=True) + "\n")
            total += len(records)

    print(f"wrote {total} columns -> {out}")
    if skipped:
        print(f"skipped {len(skipped)} route(s):")
        for route_path, reason in skipped[:8]:
            print(f"  {route_path}: {reason}")
    return 0 if not skipped else 2


def cmd_stats(args) -> int:
    rows = [json.loads(line) for line in Path(args.file).read_text().splitlines()
            if line.strip()]
    if not rows:
        print("columns=0")
        return 0
    coverage = set()
    modes = Counter()
    for row in rows:
        coverage.update(row.get("covered_stations", []))
        modes.update(row.get("modes", {}))
    durations = sorted(row["elapsed_s"] for row in rows)
    print(f"columns={len(rows)}")
    print(f"coverage_union={len(coverage)}")
    print(f"duration_min={hms(durations[0])}")
    print(f"duration_median={hms(durations[len(durations) // 2])}")
    print(f"duration_max={hms(durations[-1])}")
    print(f"mode_edges={dict(modes)}")
    return 0


def cmd_cover(args) -> int:
    rows = [json.loads(line) for line in Path(args.file).read_text().splitlines()
            if line.strip()]
    if not rows:
        raise SystemExit("no columns to solve")

    canonical = set(StationIndex.load().canonical_stations)
    by_station = {station: [] for station in canonical}
    for i, row in enumerate(rows):
        for station in row.get("covered_stations", []):
            if station in by_station:
                by_station[station].append(i)
    missing = sorted(st for st, cols in by_station.items() if not cols)
    if missing:
        raise SystemExit(f"column pool misses {len(missing)} station(s), e.g. {missing[:8]}")

    model = cp_model.CpModel()
    x = [model.NewBoolVar(f"col_{i}") for i in range(len(rows))]
    for station, cols in by_station.items():
        model.Add(sum(x[i] for i in cols) >= 1).WithName(f"cover_{station}")
    model.Minimize(sum((int(row["elapsed_s"]) + args.column_penalty) * x[i]
                       for i, row in enumerate(rows)))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.time_limit
    solver.parameters.num_search_workers = args.workers
    status = solver.Solve(model)
    selected = [i for i, var in enumerate(x) if solver.Value(var)]
    selected_rows = [rows[i] for i in selected]

    coverage = set()
    modes = Counter()
    for row in selected_rows:
        coverage.update(row.get("covered_stations", []))
        modes.update(row.get("modes", {}))
    objective = sum(int(row["elapsed_s"]) for row in selected_rows)
    source_counts = Counter(row["source"] for row in selected_rows)
    result = {
        "status": solver.StatusName(status),
        "objective_s_excluding_penalty": int(objective),
        "objective": hms(objective),
        "solver_objective_with_penalty": int(solver.ObjectiveValue()),
        "selected_count": len(selected_rows),
        "covered_count": len(coverage & canonical),
        "mode_edges": dict(modes),
        "source_counts": dict(source_counts),
        "selected_columns": selected_rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, sort_keys=True))
    print(f"status={result['status']}")
    print(f"selected={result['selected_count']} covered={result['covered_count']}/472 "
          f"objective={result['objective']}")
    print(f"mode_edges={result['mode_edges']}")
    print(f"wrote {out}")
    return 0


def _column_anchors(row, si):
    if "path" not in row:
        raise ValueError(f"column {row.get('column_id')} has no path slice")
    anchors = []
    for stop, _t, *_rest in row["path"]:
        station = si.resolve(stop)
        if station not in anchors:
            anchors.append(station)
    return anchors


def _static_connector_costs(rows, D):
    costs = {}
    for i, row_i in enumerate(rows):
        a = row_i["end"]["station"]
        for j, row_j in enumerate(rows):
            if i == j:
                continue
            b = row_j["start"]["station"]
            costs[i, j] = 0 if a == b else int(D.get(a, {}).get(b, BIG_COST))
    return costs


def _candidate_connector_arcs(rows, D, top_k):
    arcs = []
    for i, row_i in enumerate(rows):
        a = row_i["end"]["station"]
        candidates = []
        for j, row_j in enumerate(rows):
            if i == j:
                continue
            b = row_j["start"]["station"]
            cost = 0 if a == b else int(D.get(a, {}).get(b, BIG_COST))
            if cost < BIG_COST:
                candidates.append((cost, i, j))
        candidates.sort()
        arcs.extend(candidates[:top_k])
    return arcs


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
        runs = RunLayer.from_graph(G, radius_m=radius_m)
        return runs, _run_metric_edges(runs, si)
    raise ValueError(f"unknown run mode: {mode!r}")


def _exact_event_node(stop_events, node_t, stop, t):
    times, nodes = stop_events.get(stop, ((), ()))
    i = bisect.bisect_left(times, int(t))
    if i < len(times) and times[i] == int(t):
        return nodes[i]
    return None


def _row_endpoint_nodes(stop_events, node_t, row):
    start = _exact_event_node(stop_events, node_t, row["start"]["stop"], row["start"]["time"])
    end = _exact_event_node(stop_events, node_t, row["end"]["stop"], row["end"]["time"])
    if start is None:
        raise ValueError(f"missing start event for {row.get('column_id')}")
    if end is None:
        raise ValueError(f"missing end event for {row.get('column_id')}")
    return start, end


def _row_endpoint_nodes_from_tables(node_stop, node_t, row):
    """Resolve one row's endpoints without building a full stop-event index."""
    start_stop = row["start"]["stop"]
    start_t = int(row["start"]["time"])
    end_stop = row["end"]["stop"]
    end_t = int(row["end"]["time"])
    start = end = None
    for nid, stop in enumerate(node_stop):
        if stop == start_stop and int(node_t[nid]) == start_t:
            start = nid
        if stop == end_stop and int(node_t[nid]) == end_t:
            end = nid
        if start is not None and end is not None:
            break
    if start is None:
        raise ValueError(f"missing start event for {row.get('column_id')}")
    if end is None:
        raise ValueError(f"missing end event for {row.get('column_id')}")
    return start, end


def _row_path_nodes(stop_events, node_t, row):
    nodes = []
    for stop, t, *_rest in row["path"]:
        node = _exact_event_node(stop_events, node_t, stop, int(t))
        if node is None:
            raise ValueError(f"missing path event {stop}@{int(t)} for {row.get('column_id')}")
        nodes.append(node)
    return nodes


def _dijkstra_to_node(adj, src, target, cap, runs=None, max_cost=None):
    dist = {src: 0.0}
    prev = {}
    pq = [(0.0, src)]
    popped = 0
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, INF):
            continue
        if max_cost is not None and d > max_cost:
            return None, prev
        popped += 1
        if popped > cap:
            return None, prev
        if u == target:
            return u, prev
        for v, ed in adj[u].items():
            nd = d + ed["weight"]
            if max_cost is not None and nd > max_cost:
                continue
            if nd < dist.get(v, INF):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
        if runs is not None:
            for v, w, _info in runs.run_successors(u):
                nd = d + w
                if max_cost is not None and nd > max_cost:
                    continue
                if nd < dist.get(v, INF):
                    dist[v] = nd
                    prev[v] = u
                    heapq.heappush(pq, (nd, v))
    return None, prev


def _dijkstra_to_station_bounded(adj, node_off, src, target, cap, runs=None, max_cost=None):
    dist = {src: 0.0}
    prev = {}
    pq = [(0.0, src)]
    popped = 0
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, INF):
            continue
        if max_cost is not None and d > max_cost:
            return None, prev
        popped += 1
        if popped > cap:
            return None, prev
        if u != src and node_off[u] == target:
            return u, prev
        for v, ed in adj[u].items():
            nd = d + ed["weight"]
            if max_cost is not None and nd > max_cost:
                continue
            if nd < dist.get(v, INF):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
        if runs is not None:
            for v, w, _info in runs.run_successors(u):
                nd = d + w
                if max_cost is not None and nd > max_cost:
                    continue
                if nd < dist.get(v, INF):
                    dist[v] = nd
                    prev[v] = u
                    heapq.heappush(pq, (nd, v))
    return None, prev


def _path_from_prev(prev, src, target):
    if target == src:
        return []
    seg = []
    x = target
    while x != src:
        seg.append(x)
        if x not in prev:
            return None
        x = prev[x]
    return list(reversed(seg))


def _dynamic_connector_costs(G, tables, rows, runs, cap):
    node_off, node_stop, node_t = tables
    stop_events = _stop_index(node_stop, node_t)
    costs = {}
    for i, row_i in enumerate(rows):
        src = _exact_event_node(
            stop_events, node_t, row_i["end"]["stop"], row_i["end"]["time"])
        if src is None:
            raise ValueError(f"missing end event for {row_i.get('column_id')}")
        a = row_i["end"]["station"]
        for j, row_j in enumerate(rows):
            if i == j:
                continue
            b = row_j["start"]["station"]
            if a == b:
                costs[i, j] = 0
                continue
            tgt, _prev = _dijkstra_to_station(G._adj, node_off, src, b, cap, runs=runs)
            costs[i, j] = BIG_COST if tgt is None else _elapsed(node_t[src], node_t[tgt])
        if (i + 1) % 5 == 0 or i + 1 == len(rows):
            print(f"dynamic connectors: {i + 1}/{len(rows)} rows", flush=True)
    return costs


def _load_column_rows(path: Path):
    text = path.read_text().strip()
    if not text:
        return []
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) == 1 and text[0] == "{":
        data = json.loads(text)
        if "selected_columns" in data:
            return data["selected_columns"]
        if "columns" in data:
            return data["columns"]
        raise ValueError(f"{path} does not contain selected_columns")
    return [json.loads(line) for line in lines]


def _exclude_column_rows(rows, column_ids=None, column_prefixes=None):
    excluded_ids = {
        column_id.strip()
        for column_id in str(column_ids or "").split(",")
        if column_id.strip()
    }
    excluded_prefixes = [
        prefix.strip()
        for prefix in str(column_prefixes or "").split(",")
        if prefix.strip()
    ]
    if not excluded_ids and not excluded_prefixes:
        return rows, []
    kept = []
    excluded = []
    for row in rows:
        column_id = str(row.get("column_id", ""))
        if column_id in excluded_ids or any(
            column_id.startswith(prefix) for prefix in excluded_prefixes
        ):
            excluded.append(column_id)
        else:
            kept.append(row)
    return kept, excluded


def _phase_proxy_gap(source_end_t, target_start_t, static_cost):
    gap = (int(target_start_t) - int(source_end_t)) % WEEK
    if gap >= static_cost:
        return gap
    return gap + WEEK


def _direction(stop):
    return stop[-1] if stop and stop[-1] in {"N", "S"} else ""


def _load_exact_arc_cache(path):
    if not path:
        return {}
    cache_path = Path(path)
    if not cache_path.exists():
        return {}
    cache = {}
    with cache_path.open() as f:
        for line in f:
            if not line.strip():
                continue
            try:
                item = json.loads(line)
            except json.JSONDecodeError:
                continue
            key = tuple(item.get("key", []))
            if key:
                cache[key] = item.get("cost_s")
    return cache


def _append_exact_arc_cache(path, records):
    if not path or not records:
        return
    cache_path = Path(path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("a") as f:
        for key, cost_s in records:
            f.write(json.dumps({"key": list(key), "cost_s": cost_s}, sort_keys=True) + "\n")


def _exact_arc_cache_key(node_stop, node_t, src, dst, run_mode, run_radius, cap, max_cost):
    return (
        node_stop[src],
        int(node_t[src]),
        node_stop[dst],
        int(node_t[dst]),
        str(run_mode),
        int(round(float(run_radius) * 1000)),
        int(cap),
        int(max_cost),
    )


def _candidate_exact_arc_proxy(rows, D, i, j,
                               min_same_station_gap_s=0,
                               min_opposite_direction_gap_s=0):
    if i == j:
        return None
    row_i = rows[i]
    row_j = rows[j]
    a = row_i["end"]["station"]
    end_t = int(row_i["end"]["time"])
    end_stop = row_i["end"]["stop"]
    b = row_j["start"]["station"]
    static = 0 if a == b else int(D.get(a, {}).get(b, BIG_COST))
    if static >= BIG_COST:
        return None
    raw_gap = (int(row_j["start"]["time"]) - end_t) % WEEK
    start_stop = row_j["start"]["stop"]
    if a == b and end_stop != start_stop and raw_gap < min_same_station_gap_s:
        return None
    if (_direction(end_stop) and _direction(start_stop)
            and _direction(end_stop) != _direction(start_stop)
            and raw_gap < min_opposite_direction_gap_s):
        return None
    return _phase_proxy_gap(end_t, row_j["start"]["time"], static)


def _candidate_exact_arcs(rows, D, top_k, max_proxy_s,
                          min_same_station_gap_s=0,
                          min_opposite_direction_gap_s=0):
    arcs = []
    for i, _row_i in enumerate(rows):
        candidates = []
        for j, _row_j in enumerate(rows):
            proxy = _candidate_exact_arc_proxy(
                rows,
                D,
                i,
                j,
                min_same_station_gap_s=min_same_station_gap_s,
                min_opposite_direction_gap_s=min_opposite_direction_gap_s,
            )
            if proxy is None or proxy > max_proxy_s:
                continue
            candidates.append((proxy, i, j))
        candidates.sort()
        arcs.extend(candidates[:top_k])
    return arcs


def _price_exact_arc_candidates(G, tables, rows, runs, candidates,
                                max_connector_s, cap, run_mode, run_radius,
                                cache_path=None, progress_label="exact connectors"):
    node_off, node_stop, node_t = tables
    stop_events = _stop_index(node_stop, node_t)
    arcs = []
    cache = _load_exact_arc_cache(cache_path)
    cache_hits = 0
    cache_misses = 0
    cache_new = []
    for k, (_proxy, i, j) in enumerate(candidates, start=1):
        _start_i, end_i = _row_endpoint_nodes(stop_events, node_t, rows[i])
        start_j, _end_j = _row_endpoint_nodes(stop_events, node_t, rows[j])
        target_elapsed = _elapsed(node_t[end_i], node_t[start_j])
        if end_i == start_j:
            cost = 0
        else:
            max_cost = min(max_connector_s, target_elapsed)
            key = _exact_arc_cache_key(
                node_stop,
                node_t,
                end_i,
                start_j,
                run_mode,
                run_radius,
                cap,
                max_cost,
            )
            if key in cache:
                cache_hits += 1
                cached_cost = cache[key]
                cost = BIG_COST if cached_cost is None else int(cached_cost)
            else:
                cache_misses += 1
                tgt, _prev = _dijkstra_to_node(
                    G._adj, end_i, start_j, cap, runs=runs, max_cost=max_cost)
                cost = BIG_COST if tgt is None else target_elapsed
                cache[key] = None if cost >= BIG_COST else int(cost)
                cache_new.append((key, cache[key]))
        if cost <= max_connector_s:
            arcs.append((int(cost), i, j))
        if k % 500 == 0 or k == len(candidates):
            print(f"{progress_label}: {k}/{len(candidates)} candidates "
                  f"kept={len(arcs)}", flush=True)
    _append_exact_arc_cache(cache_path, cache_new)
    if cache_path:
        print(f"{progress_label} cache: hits={cache_hits} misses={cache_misses} "
              f"new={len(cache_new)} entries={len(cache)} path={cache_path}",
              flush=True)
    return arcs


def _exact_connector_arcs(G, tables, rows, runs, D, top_k, max_proxy_s,
                          max_connector_s, cap, run_mode, run_radius,
                          cache_path=None,
                          min_same_station_gap_s=0,
                          min_opposite_direction_gap_s=0):
    candidates = _candidate_exact_arcs(
        rows, D, top_k, max_proxy_s,
        min_same_station_gap_s=min_same_station_gap_s,
        min_opposite_direction_gap_s=min_opposite_direction_gap_s)
    return _price_exact_arc_candidates(
        G,
        tables,
        rows,
        runs,
        candidates,
        max_connector_s=max_connector_s,
        cap=cap,
        run_mode=run_mode,
        run_radius=run_radius,
        cache_path=cache_path,
    )


def _required_neighbor_exact_arcs(G, tables, rows, runs, D, required_row_indices,
                                  top_k, max_proxy_s, max_connector_s, cap,
                                  run_mode, run_radius, cache_path=None,
                                  min_same_station_gap_s=0,
                                  min_opposite_direction_gap_s=0):
    required = sorted(set(required_row_indices or []))
    if not required or top_k <= 0:
        return []
    by_pair = {}
    for j in required:
        incoming = []
        for i in range(len(rows)):
            proxy = _candidate_exact_arc_proxy(
                rows,
                D,
                i,
                j,
                min_same_station_gap_s=min_same_station_gap_s,
                min_opposite_direction_gap_s=min_opposite_direction_gap_s,
            )
            if proxy is None or proxy > max_proxy_s:
                continue
            incoming.append((proxy, i, j))
        for proxy, i, j in sorted(incoming)[:top_k]:
            by_pair[(i, j)] = (proxy, i, j)

    for i in required:
        outgoing = []
        for j in range(len(rows)):
            proxy = _candidate_exact_arc_proxy(
                rows,
                D,
                i,
                j,
                min_same_station_gap_s=min_same_station_gap_s,
                min_opposite_direction_gap_s=min_opposite_direction_gap_s,
            )
            if proxy is None or proxy > max_proxy_s:
                continue
            outgoing.append((proxy, i, j))
        for proxy, i, j in sorted(outgoing)[:top_k]:
            by_pair[(i, j)] = (proxy, i, j)

    candidates = sorted(by_pair.values())
    print(f"required-neighbor arc candidates: required={len(required)} "
          f"top_k={top_k} candidates={len(candidates)}", flush=True)
    return _price_exact_arc_candidates(
        G,
        tables,
        rows,
        runs,
        candidates,
        max_connector_s=max_connector_s,
        cap=cap,
        run_mode=run_mode,
        run_radius=run_radius,
        cache_path=cache_path,
        progress_label="required exact connectors",
    )


def _merge_arcs_min_cost(arcs, extra_arcs):
    by_pair = {}
    for cost, i, j in list(arcs) + list(extra_arcs):
        pair = (i, j)
        if pair not in by_pair or cost < by_pair[pair]:
            by_pair[pair] = int(cost)
    return [(cost, i, j) for (i, j), cost in by_pair.items()]


def _forced_solution_arcs(G, tables, rows, runs, solution_files, cap, max_connector_s):
    if not solution_files:
        return []
    by_id = {row.get("column_id"): i for i, row in enumerate(rows)}
    node_t = tables[2]
    stop_events = _stop_index(tables[1], node_t)
    forced = []
    missing = 0
    infeasible = 0
    missing_pairs = []
    infeasible_pairs = []
    for solution_file in solution_files:
        solution_ids = _selected_column_ids_from_solution(Path(solution_file))
        for column_i, column_j in zip(solution_ids, solution_ids[1:]):
            i = by_id.get(column_i)
            j = by_id.get(column_j)
            if i is None or j is None:
                missing += 1
                if len(missing_pairs) < 8:
                    missing_pairs.append((column_i, column_j))
                continue
            _start_i, end_i = _row_endpoint_nodes(stop_events, node_t, rows[i])
            start_j, _end_j = _row_endpoint_nodes(stop_events, node_t, rows[j])
            gap = _elapsed(node_t[end_i], node_t[start_j])
            if end_i == start_j:
                forced.append((0, i, j))
                continue
            max_cost = min(max_connector_s, gap)
            tgt, _prev = _dijkstra_to_node(
                G._adj, end_i, start_j, cap, runs=runs, max_cost=max_cost)
            if tgt is None:
                infeasible += 1
                if len(infeasible_pairs) < 8:
                    infeasible_pairs.append((column_i, column_j))
                continue
            forced.append((gap, i, j))
    print(f"forced solution arcs: kept={len(forced)} missing={missing} "
          f"infeasible={infeasible}", flush=True)
    if missing_pairs:
        print(f"forced solution missing arcs e.g. {missing_pairs}", flush=True)
    if infeasible_pairs:
        print(f"forced solution infeasible arcs e.g. {infeasible_pairs}", flush=True)
    return forced


def _exact_connector_cost(G, tables, rows, runs, i, j, cap):
    node_t = tables[2]
    stop_events = _stop_index(tables[1], node_t)
    _start_i, end_i = _row_endpoint_nodes(stop_events, node_t, rows[i])
    start_j, _end_j = _row_endpoint_nodes(stop_events, node_t, rows[j])
    if end_i == start_j:
        return 0
    max_cost = _elapsed(node_t[end_i], node_t[start_j])
    tgt, _prev = _dijkstra_to_node(
        G._adj, end_i, start_j, cap, runs=runs, max_cost=max_cost)
    if tgt is None:
        return None
    return max_cost


def _replace_arc_cost(arcs, pair, cost):
    out = []
    replaced = False
    for old_cost, i, j in arcs:
        if (i, j) == pair:
            out.append((int(cost), i, j))
            replaced = True
        else:
            out.append((old_cost, i, j))
    return out, replaced


def _realize_exact_column_path(G, tables, rows, order, runs, cap):
    node_stop, node_t = tables[1], tables[2]
    stop_events = _stop_index(node_stop, node_t)
    route_nodes = []
    for pos, row_idx in enumerate(order):
        row_nodes = _row_path_nodes(stop_events, node_t, rows[row_idx])
        if not route_nodes:
            route_nodes.extend(row_nodes)
            continue
        src = route_nodes[-1]
        target = row_nodes[0]
        if src != target:
            max_cost = _elapsed(node_t[src], node_t[target])
            tgt, prev = _dijkstra_to_node(
                G._adj, src, target, cap, runs=runs, max_cost=max_cost)
            if tgt is None:
                src_stop, src_t = node_stop[src], node_t[src]
                dst_stop, dst_t = node_stop[target], node_t[target]
                print("exact replay failed: "
                      f"{rows[order[pos - 1]].get('column_id')} "
                      f"{src_stop}@{src_t} -> "
                      f"{rows[row_idx].get('column_id')} {dst_stop}@{dst_t}",
                      flush=True)
                return None, (order[pos - 1], row_idx)
            connector = _path_from_prev(prev, src, target)
            if connector is None:
                print("exact replay failed: predecessor chain missing", flush=True)
                return None, (order[pos - 1], row_idx)
            route_nodes.extend(connector)
        route_nodes.extend(row_nodes[1:])
    return [[node_stop[n], node_t[n]] for n in route_nodes], None


def _solve_column_order(rows, connector_costs, time_limit_s, log_search=False):
    n = len(rows)
    start_dummy = n
    end_dummy = n + 1
    manager = pywrapcp.RoutingIndexManager(n + 2, 1, [start_dummy], [end_dummy])
    routing = pywrapcp.RoutingModel(manager)

    def connector(i, j):
        return int(connector_costs.get((i, j), BIG_COST))

    def transit(from_index, to_index):
        a = manager.IndexToNode(from_index)
        b = manager.IndexToNode(to_index)
        if a == end_dummy:
            return BIG_COST
        if a == start_dummy:
            return 0 if b < n else BIG_COST
        if b == end_dummy:
            return int(rows[a]["elapsed_s"])
        if b == start_dummy:
            return BIG_COST
        return int(rows[a]["elapsed_s"]) + connector(a, b)

    transit_idx = routing.RegisterTransitCallback(transit)
    routing.SetArcCostEvaluatorOfAllVehicles(transit_idx)

    params = pywrapcp.DefaultRoutingSearchParameters()
    params.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
    params.local_search_metaheuristic = routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
    params.time_limit.FromSeconds(max(1, int(time_limit_s)))
    params.log_search = bool(log_search)

    solution = routing.SolveWithParameters(params)
    if solution is None:
        return None, None

    order = []
    idx = routing.Start(0)
    while not routing.IsEnd(idx):
        node = manager.IndexToNode(idx)
        if node < n:
            order.append(node)
        idx = solution.Value(routing.NextVar(idx))
    return order, int(solution.ObjectiveValue())


def _realize_column_order(G, si, rows, order, start_stop, target_t,
                          terminal_runs, run_radius):
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    stop_events = _stop_index(node_stop, node_t)
    start_node = _first_event_at_or_after(stop_events, start_stop, target_t)
    if start_node is None:
        raise ValueError(f"no event found for start stop {start_stop!r}")

    runs = None
    if terminal_runs:
        runs, _extra_edges = _terminal_run_layer(G, si, node_off, run_radius)

    anchors = []
    seen = set()
    for idx in order:
        for station in _column_anchors(rows[idx], si):
            if station not in seen:
                seen.add(station)
                anchors.append(station)
    start_station = node_off[start_node]
    target_anchors = anchors[1:] if anchors and anchors[0] == start_station else anchors
    suffix, end, visited = realize_from(
        G._adj, tables, start_node, {start_station}, target_anchors, runs=runs)
    if suffix is None:
        return None, None, None, None

    nodes = [start_node] + suffix
    path = [[node_stop[n], node_t[n]] for n in nodes]
    elapsed = _elapsed(node_t[start_node], node_t[end])
    covered = len(visited & si.canonical_stations)
    return path, elapsed, covered, anchors


def cmd_sequence(args) -> int:
    import pickle

    cover = json.loads(Path(args.cover_solution).read_text())
    rows = cover.get("selected_columns", [])
    if not rows:
        raise SystemExit("cover solution has no selected_columns")
    if args.max_columns:
        rows = rows[:args.max_columns]

    si = StationIndex.load()
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off = tables[0]
    extra_edges = []
    runs = None
    if args.terminal_runs:
        runs, extra_edges = _terminal_run_layer(G, si, node_off, args.run_radius)
    if args.connector_mode == "dynamic":
        print(f"building dynamic connector matrix for {len(rows)} columns", flush=True)
        connector_costs = _dynamic_connector_costs(
            G, tables, rows, runs, int(args.connector_cap))
    else:
        print(f"building static connector metric for {len(rows)} columns", flush=True)
        D, _ = static_station_metric(G, node_off, extra_edges)
        connector_costs = _static_connector_costs(rows, D)

    order, static_obj = _solve_column_order(
        rows, connector_costs, args.time_limit, args.log_search)
    if not order:
        raise SystemExit("OR-Tools failed to sequence columns")

    first = rows[order[0]]
    start_stop = args.start or first["start"]["stop"]
    target_t = int(args.time) if args.time is not None else int(first["start"]["time"])
    print(f"column objective={hms(static_obj)} columns={len(order)} "
          f"start={start_stop}@{target_t}", flush=True)

    path, elapsed, covered, anchors = _realize_column_order(
        G, si, rows, order, start_stop, target_t, args.terminal_runs, args.run_radius)
    if path is None:
        raise SystemExit("failed to realize column sequence on schedule")

    ordered_rows = [rows[i] for i in order]
    meta = {
        "radius_m": args.run_radius if args.terminal_runs else 5000,
        "elapsed_s": int(elapsed),
        "notes": ("Column cover sequencing prototype "
                  f"cover={args.cover_solution} connector={args.connector_mode} "
                  f"sequence_obj={static_obj}"),
        "sequence_objective_s": int(static_obj),
        "sequence_connector_mode": args.connector_mode,
        "sequence_columns": [row["column_id"] for row in ordered_rows],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": meta, "path": path}, sort_keys=True))
    print(f"realized: elapsed={hms(elapsed)} covered={covered}/472 "
          f"anchors={len(anchors)} steps={len(path)} -> {out}", flush=True)

    if args.validate:
        full_runs = RunLayer.from_graph(G, radius_m=meta["radius_m"])
        result = validate(Problem(G, full_runs, si), path)
        print(result_line(result), flush=True)
    return 0


def _solve_path_cover(rows, arcs, time_limit_s, workers, column_penalty,
                      require_pricing_kind=None, min_required_pricing=0,
                      required_row_indices=None,
                      required_index_groups=None,
                      allowed_start_rows=None, allowed_end_rows=None,
                      uncovered_penalty_s=0,
                      min_covered_count=None,
                      max_total_elapsed_s=None,
                      hint_order_indices=None,
                      protected_stations=None,
                      protected_station_groups=None,
                      uncovered_penalty_groups=None,
                      ordered_station_groups=None,
                      strict_order_station_groups=False,
                      first_hit_order_station_groups=False,
                      stop_after_first_solution=False,
                      forbidden_row_indices=None):
    canonical = set(StationIndex.load().canonical_stations)
    protected_stations = set(protected_stations or []) & canonical
    protected_station_groups = list(protected_station_groups or [])
    uncovered_penalty_groups = list(uncovered_penalty_groups or [])
    ordered_station_groups = list(ordered_station_groups or [])
    extra_uncovered_penalty_by_station = defaultdict(int)
    for group in uncovered_penalty_groups:
        for station in group["stations"]:
            if station in canonical:
                extra_uncovered_penalty_by_station[station] += int(group["penalty_s"])
    soft_coverage = bool(uncovered_penalty_s or extra_uncovered_penalty_by_station)
    by_station = {station: [] for station in canonical}
    for i, row in enumerate(rows):
        for station in row.get("covered_stations", []):
            if station in by_station:
                by_station[station].append(i)
    missing = sorted(st for st, cols in by_station.items() if not cols)
    if missing and not soft_coverage:
        raise ValueError(f"column pool misses {len(missing)} station(s), e.g. {missing[:8]}")
    missing_protected = sorted(st for st in protected_stations if not by_station.get(st))
    if missing_protected:
        raise ValueError(
            f"column pool misses protected station(s), e.g. {missing_protected[:8]}")
    for group in protected_station_groups:
        min_hits = int(group["min_hits"])
        if min_hits < 0 or min_hits > len(group["stations"]):
            raise ValueError(
                f"protected group {group['name']!r} has impossible min_hits "
                f"{min_hits} for {len(group['stations'])} station(s)")
    ordered_group_candidates = []
    for group in ordered_station_groups:
        min_hits = int(group["min_hits"])
        stations = set(group["stations"])
        candidates = [
            i
            for i, row in enumerate(rows)
            if len(stations.intersection(row.get("covered_stations", []))) >= min_hits
        ]
        if not candidates:
            raise ValueError(
                f"ordered group {group['label']!r} has no row covering "
                f"{min_hits}/{len(stations)} station(s)")
        ordered_group_candidates.append(candidates)

    model = cp_model.CpModel()
    n = len(rows)
    x = [model.NewBoolVar(f"col_{i}") for i in range(n)]
    start = [model.NewBoolVar(f"start_{i}") for i in range(n)]
    end = [model.NewBoolVar(f"end_{i}") for i in range(n)]
    order = [model.NewIntVar(0, n, f"order_{i}") for i in range(n)]
    y = {(i, j): model.NewBoolVar(f"arc_{i}_{j}") for _cost, i, j in arcs}

    incoming = {i: [] for i in range(n)}
    outgoing = {i: [] for i in range(n)}
    for _cost, i, j in arcs:
        outgoing[i].append(y[i, j])
        incoming[j].append(y[i, j])

    cover_vars = {}
    for station, cols in by_station.items():
        if station in protected_stations:
            model.Add(sum(x[i] for i in cols) >= 1).WithName(f"protect_{station}")
        elif soft_coverage:
            covered = model.NewBoolVar(f"covered_{station}")
            cover_vars[station] = covered
            if cols:
                model.Add(sum(x[i] for i in cols) >= covered).WithName(
                    f"covered_if_selected_{station}")
            else:
                model.Add(covered == 0).WithName(f"uncoverable_{station}")
        else:
            model.Add(sum(x[i] for i in cols) >= 1).WithName(f"cover_{station}")
    if protected_station_groups:
        if not soft_coverage:
            print("protected station groups are redundant without relaxed coverage",
                  flush=True)
        for group_index, group in enumerate(protected_station_groups):
            fixed_hits = 0
            terms = []
            for station in group["stations"]:
                if station in protected_stations:
                    fixed_hits += 1
                elif soft_coverage:
                    terms.append(cover_vars[station])
            if soft_coverage:
                model.Add(sum(terms) + fixed_hits >= int(group["min_hits"])).WithName(
                    f"protect_group_{group_index}_{group['name']}")
    if min_covered_count is not None:
        min_covered_count = int(min_covered_count)
        if min_covered_count < 0 or min_covered_count > len(canonical):
            raise ValueError(
                f"min covered count must be between 0 and {len(canonical)}, "
                f"got {min_covered_count}")
        if soft_coverage:
            model.Add(
                sum(cover_vars.values()) + len(protected_stations) >= min_covered_count
            ).WithName("min_covered_count")
    if require_pricing_kind and min_required_pricing:
        required = [i for i, row in enumerate(rows) if require_pricing_kind in row]
        if len(required) < min_required_pricing:
            raise ValueError(
                f"only {len(required)} columns contain {require_pricing_kind!r}, "
                f"cannot require {min_required_pricing}"
            )
        model.Add(sum(x[i] for i in required) >= min_required_pricing).WithName(
            f"require_{require_pricing_kind}")
    for i in sorted(set(required_row_indices or [])):
        model.Add(x[i] == 1).WithName(f"require_column_{i}")
    for i in sorted(set(forbidden_row_indices or [])):
        model.Add(x[i] == 0).WithName(f"forbid_column_{i}")
    for name, indices, min_count in required_index_groups or []:
        model.Add(sum(x[i] for i in indices) >= int(min_count)).WithName(name)
    ordered_group_vars = []
    for group_index, (group, candidates) in enumerate(
        zip(ordered_station_groups, ordered_group_candidates)
    ):
        witness = []
        pos = model.NewIntVar(0, n, f"ordered_group_pos_{group_index}")
        for i in candidates:
            z = model.NewBoolVar(f"ordered_group_{group_index}_row_{i}")
            model.Add(z <= x[i])
            model.Add(pos == order[i]).OnlyEnforceIf(z)
            if first_hit_order_station_groups:
                model.Add(pos <= order[i]).OnlyEnforceIf(x[i])
            witness.append((i, z))
        model.Add(sum(z for _i, z in witness) == 1).WithName(
            f"ordered_group_witness_{group_index}_{group['name']}")
        ordered_group_vars.append({
            "group": group,
            "candidates": candidates,
            "position": pos,
            "witness": witness,
        })
    for left, right in zip(ordered_group_vars, ordered_group_vars[1:]):
        if strict_order_station_groups:
            model.Add(left["position"] + 1 <= right["position"]).WithName(
                f"ordered_group_strict_before_{left['group']['name']}_{right['group']['name']}")
        else:
            model.Add(left["position"] <= right["position"]).WithName(
                f"ordered_group_before_{left['group']['name']}_{right['group']['name']}")
    model.Add(sum(start) == 1)
    model.Add(sum(end) == 1)
    if allowed_start_rows is not None:
        allowed = set(allowed_start_rows)
        if not allowed:
            raise ValueError("start time window leaves no eligible start columns")
        for i in range(n):
            if i not in allowed:
                model.Add(start[i] == 0)
    if allowed_end_rows is not None:
        allowed = set(allowed_end_rows)
        if not allowed:
            raise ValueError("end time window leaves no eligible end columns")
        for i in range(n):
            if i not in allowed:
                model.Add(end[i] == 0)
    for i in range(n):
        model.Add(start[i] + sum(incoming[i]) == x[i])
        model.Add(end[i] + sum(outgoing[i]) == x[i])
        model.Add(order[i] <= n * x[i])
        model.Add(order[i] >= x[i])

    for _cost, i, j in arcs:
        model.Add(order[j] >= order[i] + 1 - n * (1 - y[i, j]))

    arc_cost = {(i, j): cost for cost, i, j in arcs}
    elapsed_terms = (
        sum(int(row["elapsed_s"]) * x[i] for i, row in enumerate(rows))
        + sum(arc_cost[i, j] * var for (i, j), var in y.items())
    )
    max_possible_elapsed = (
        sum(int(row["elapsed_s"]) for row in rows)
        + sum(int(cost) for cost, _i, _j in arcs)
    )
    elapsed_var = model.NewIntVar(0, max_possible_elapsed, "total_elapsed")
    model.Add(elapsed_var == elapsed_terms).WithName("bind_total_elapsed")
    if max_total_elapsed_s is not None:
        model.Add(elapsed_var <= int(max_total_elapsed_s)).WithName("max_total_elapsed")
    objective = elapsed_var
    if column_penalty:
        objective += sum(int(column_penalty) * x[i] for i in range(n))
    if soft_coverage:
        objective += sum(int(uncovered_penalty_s) * (1 - var)
                         for var in cover_vars.values())
        objective += sum(
            int(extra_uncovered_penalty_by_station.get(station, 0)) * (1 - var)
            for station, var in cover_vars.items()
        )
    model.Minimize(objective)

    hint_order_indices = list(hint_order_indices or [])
    if hint_order_indices:
        hint_set = set(hint_order_indices)
        for i in range(n):
            model.AddHint(x[i], 1 if i in hint_set else 0)
            model.AddHint(start[i], 1 if i == hint_order_indices[0] else 0)
            model.AddHint(end[i], 1 if i == hint_order_indices[-1] else 0)
        for pos, i in enumerate(hint_order_indices, start=1):
            model.AddHint(order[i], pos)
        hint_arcs = set(zip(hint_order_indices, hint_order_indices[1:]))
        for pair, var in y.items():
            model.AddHint(var, 1 if pair in hint_arcs else 0)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_search_workers = workers
    solver.parameters.stop_after_first_solution = bool(stop_after_first_solution)
    status = solver.Solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return solver, status, [], {}

    selected = [i for i, var in enumerate(x) if solver.Value(var)]
    start_nodes = [i for i, var in enumerate(start) if solver.Value(var)]
    next_by = {i: j for (i, j), var in y.items() if solver.Value(var)}
    ordered = []
    current = start_nodes[0] if start_nodes else None
    while current is not None and current not in ordered:
        ordered.append(current)
        current = next_by.get(current)
    if set(ordered) != set(selected):
        ordered = sorted(selected, key=lambda i: solver.Value(order[i]))

    return solver, status, ordered, {
        "x": x,
        "y": y,
        "arc_cost": arc_cost,
        "covered": cover_vars,
        "elapsed_var": elapsed_var,
        "ordered_group_vars": ordered_group_vars,
    }


def cmd_path_cover(args) -> int:
    import pickle

    rows = [json.loads(line) for line in Path(args.columns_file).read_text().splitlines()
            if line.strip()]
    if args.max_columns:
        rows = rows[:args.max_columns]
    if not rows:
        raise SystemExit("no columns to solve")

    si = StationIndex.load()
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off = tables[0]
    extra_edges = []
    if args.terminal_runs:
        _runs, extra_edges = _terminal_run_layer(G, si, node_off, args.run_radius)
    print(f"building static connector arcs for {len(rows)} columns", flush=True)
    D, _ = static_station_metric(G, node_off, extra_edges)
    arcs = _candidate_connector_arcs(rows, D, args.top_k)
    print(f"path-cover model: columns={len(rows)} arcs={len(arcs)} "
          f"top_k={args.top_k}", flush=True)

    solver, status, order, _vars = _solve_path_cover(
        rows, arcs, args.time_limit, args.workers, args.column_penalty)
    status_name = solver.StatusName(status)
    if not order:
        print(f"status={status_name}")
        return 1

    ordered_rows = [rows[i] for i in order]
    coverage = set()
    modes = Counter()
    for row in ordered_rows:
        coverage.update(row.get("covered_stations", []))
        modes.update(row.get("modes", {}))
    raw_elapsed = sum(int(row["elapsed_s"]) for row in ordered_rows)
    result = {
        "status": status_name,
        "objective": int(solver.ObjectiveValue()),
        "selected_count": len(ordered_rows),
        "covered_count": len(coverage & si.canonical_stations),
        "raw_column_elapsed_s": raw_elapsed,
        "raw_column_elapsed": hms(raw_elapsed),
        "mode_edges": dict(modes),
        "selected_columns": ordered_rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, sort_keys=True))
    print(f"status={status_name} selected={len(ordered_rows)} "
          f"covered={result['covered_count']}/472 raw={hms(raw_elapsed)} "
          f"objective={hms(int(solver.ObjectiveValue()))}")
    print(f"wrote {out}")

    if args.route_out:
        first = ordered_rows[0]
        start_stop = args.start or first["start"]["stop"]
        target_t = int(args.time) if args.time is not None else int(first["start"]["time"])
        path, elapsed, covered, anchors = _realize_column_order(
            G, si, ordered_rows, list(range(len(ordered_rows))), start_stop, target_t,
            args.terminal_runs, args.run_radius)
        if path is None:
            raise SystemExit("failed to realize path-cover order on schedule")
        meta = {
            "radius_m": args.run_radius if args.terminal_runs else 5000,
            "elapsed_s": int(elapsed),
            "notes": ("Path-cover column master "
                      f"columns={args.columns_file} top_k={args.top_k} "
                      f"objective={int(solver.ObjectiveValue())}"),
            "path_cover_objective_s": int(solver.ObjectiveValue()),
            "path_cover_columns": [row["column_id"] for row in ordered_rows],
        }
        route_out = Path(args.route_out)
        route_out.parent.mkdir(parents=True, exist_ok=True)
        route_out.write_text(json.dumps({"meta": meta, "path": path}, sort_keys=True))
        print(f"realized: elapsed={hms(elapsed)} covered={covered}/472 "
              f"anchors={len(anchors)} steps={len(path)} -> {route_out}")
        if args.validate:
            full_runs = RunLayer.from_graph(G, radius_m=meta["radius_m"])
            result = validate(Problem(G, full_runs, si), path)
            print(result_line(result), flush=True)
    return 0


def _selected_column_ids_from_solution(path: Path) -> list[str]:
    data = json.loads(path.read_text())
    if isinstance(data.get("selected_order"), list):
        return [str(column_id) for column_id in data["selected_order"]]
    meta = data.get("meta")
    if isinstance(meta, dict) and isinstance(meta.get("selected_order"), list):
        return [str(column_id) for column_id in meta["selected_order"]]
    if isinstance(meta, dict) and isinstance(meta.get("exact_cover_columns"), list):
        return [str(column_id) for column_id in meta["exact_cover_columns"]]
    if isinstance(data.get("selected_columns"), list):
        return [
            str(row["column_id"])
            for row in data["selected_columns"]
            if isinstance(row, dict) and row.get("column_id")
        ]
    raise ValueError(f"{path} has no selected_order or selected_columns")


def _load_station_json(path: Path) -> list[str]:
    data = json.loads(path.read_text())
    if isinstance(data, list):
        return [str(station) for station in data]
    for key in ("stations", "protected_stations", "relaxed_uncovered_stations"):
        if isinstance(data, dict) and key in data:
            return [str(station) for station in data[key]]
    raise ValueError(
        f"{path} must be a JSON list or contain stations/protected_stations/"
        "relaxed_uncovered_stations")


def _load_exact_cover_station_set(args, si) -> list[str]:
    stations = []
    if args.protect_stations:
        stations.extend(
            station.strip()
            for station in str(args.protect_stations).split(",")
            if station.strip()
        )
    for path in args.protect_stations_file or []:
        stations.extend(_load_station_json(Path(path)))
    bad = [station for station in stations if station not in si.canonical_stations]
    if bad:
        raise SystemExit(f"--protect-stations contains non-canonical ids, e.g. {bad[:8]}")
    seen = set()
    return [station for station in stations if not (station in seen or seen.add(station))]


def _load_exact_cover_station_groups(args, si) -> list[dict]:
    groups = []
    for raw_group in str(args.protect_station_groups or "").split(";"):
        raw_group = raw_group.strip()
        if not raw_group:
            continue
        try:
            name, raw_min_hits, raw_stations = raw_group.split(":", 2)
        except ValueError as exc:
            raise SystemExit(
                "--protect-station-groups entries must be NAME:MIN:station,station"
            ) from exc
        name = name.strip() or f"group{len(groups) + 1}"
        safe_name = "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in name)
        try:
            min_hits = int(raw_min_hits)
        except ValueError as exc:
            raise SystemExit(
                f"--protect-station-groups has non-integer min hits for {name!r}"
            ) from exc
        stations = [
            station.strip()
            for station in raw_stations.split(",")
            if station.strip()
        ]
        bad = [station for station in stations if station not in si.canonical_stations]
        if bad:
            raise SystemExit(
                f"--protect-station-groups {name!r} has non-canonical ids, e.g. {bad[:8]}")
        seen = set()
        stations = [
            station for station in stations
            if not (station in seen or seen.add(station))
        ]
        if not stations:
            raise SystemExit(f"--protect-station-groups {name!r} has no stations")
        if min_hits < 0 or min_hits > len(stations):
            raise SystemExit(
                f"--protect-station-groups {name!r} min {min_hits} "
                f"must be between 0 and {len(stations)}")
        groups.append({
            "name": safe_name,
            "label": name,
            "min_hits": min_hits,
            "stations": stations,
        })
    return groups


def _load_exact_cover_penalty_groups(args, si) -> list[dict]:
    groups = []
    for raw_group in str(args.uncovered_penalty_groups or "").split(";"):
        raw_group = raw_group.strip()
        if not raw_group:
            continue
        try:
            name, raw_penalty, raw_stations = raw_group.split(":", 2)
        except ValueError as exc:
            raise SystemExit(
                "--uncovered-penalty-groups entries must be "
                "NAME:PENALTY_SECONDS:station,station"
            ) from exc
        name = name.strip() or f"group{len(groups) + 1}"
        safe_name = "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in name)
        try:
            penalty_s = int(raw_penalty)
        except ValueError as exc:
            raise SystemExit(
                f"--uncovered-penalty-groups has non-integer penalty for {name!r}"
            ) from exc
        if penalty_s < 0:
            raise SystemExit(
                f"--uncovered-penalty-groups {name!r} penalty must be non-negative")
        stations = [
            station.strip()
            for station in raw_stations.split(",")
            if station.strip()
        ]
        bad = [station for station in stations if station not in si.canonical_stations]
        if bad:
            raise SystemExit(
                f"--uncovered-penalty-groups {name!r} has non-canonical ids, "
                f"e.g. {bad[:8]}")
        seen = set()
        stations = [
            station for station in stations
            if not (station in seen or seen.add(station))
        ]
        if not stations:
            raise SystemExit(f"--uncovered-penalty-groups {name!r} has no stations")
        groups.append({
            "name": safe_name,
            "label": name,
            "penalty_s": penalty_s,
            "stations": stations,
        })
    return groups


def _load_ordered_station_groups(raw_value, si) -> list[dict]:
    groups = []
    for raw_group in str(raw_value or "").split(";"):
        raw_group = raw_group.strip()
        if not raw_group:
            continue
        try:
            name, raw_min_hits, raw_stations = raw_group.split(":", 2)
        except ValueError as exc:
            raise SystemExit(
                "--order-station-groups entries must be NAME:MIN:station,station"
            ) from exc
        name = name.strip() or f"ordered_group_{len(groups) + 1}"
        safe_name = "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in name)
        try:
            min_hits = int(raw_min_hits)
        except ValueError as exc:
            raise SystemExit(
                f"--order-station-groups has non-integer min hits for {name!r}"
            ) from exc
        stations = [
            station.strip()
            for station in raw_stations.split(",")
            if station.strip()
        ]
        bad = [station for station in stations if station not in si.canonical_stations]
        if bad:
            raise SystemExit(
                f"--order-station-groups {name!r} has non-canonical ids, e.g. {bad[:8]}")
        seen = set()
        stations = [
            station for station in stations
            if not (station in seen or seen.add(station))
        ]
        if not stations:
            raise SystemExit(f"--order-station-groups {name!r} has no stations")
        if min_hits < 1 or min_hits > len(stations):
            raise SystemExit(
                f"--order-station-groups {name!r} min {min_hits} "
                f"must be between 1 and {len(stations)}")
        groups.append({
            "name": safe_name,
            "label": name,
            "min_hits": min_hits,
            "stations": stations,
        })
    return groups


def _load_prefix_requirement_groups(raw_value):
    groups = []
    for raw_group in str(raw_value or "").split(";"):
        raw_group = raw_group.strip()
        if not raw_group:
            continue
        parts = raw_group.split(":", 2)
        if len(parts) == 2:
            label = f"prefix_group_{len(groups) + 1}"
            raw_min, raw_prefixes = parts
        elif len(parts) == 3:
            label, raw_min, raw_prefixes = parts
            label = label.strip() or f"prefix_group_{len(groups) + 1}"
        else:
            raise SystemExit(
                "--require-column-id-prefix-group entries must be "
                "[LABEL:]MIN:PREFIX|PREFIX")
        try:
            min_count = int(raw_min)
        except ValueError as exc:
            raise SystemExit(
                f"--require-column-id-prefix-group {label!r} has "
                f"non-integer min count {raw_min!r}") from exc
        if min_count < 0:
            raise SystemExit(
                f"--require-column-id-prefix-group {label!r} min count "
                "must be non-negative")
        prefixes = [
            prefix.strip()
            for prefix in raw_prefixes.replace(",", "|").split("|")
            if prefix.strip()
        ]
        if not prefixes:
            raise SystemExit(
                f"--require-column-id-prefix-group {label!r} has no prefixes")
        groups.append({
            "label": label,
            "min_count": min_count,
            "prefixes": prefixes,
        })
    return groups


def cmd_exact_cover(args) -> int:
    import pickle

    rows = _load_column_rows(Path(args.columns_source))
    if args.max_columns:
        rows = rows[:args.max_columns]
    rows, excluded_rows = _exclude_column_rows(
        rows,
        column_ids=args.exclude_column_id,
        column_prefixes=args.exclude_column_id_prefix,
    )
    if excluded_rows:
        print(f"excluded columns: {len(excluded_rows)} e.g. {excluded_rows[:8]}",
              flush=True)
    if not rows:
        raise SystemExit("no columns to solve")
    missing_path = [row.get("column_id", str(i)) for i, row in enumerate(rows) if "path" not in row]
    if missing_path:
        raise SystemExit(f"{len(missing_path)} columns have no path slices, e.g. {missing_path[:8]}")

    si = StationIndex.load()
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off = tables[0]
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)

    print(f"building {args.arc_pricing} exact-time connector arcs for {len(rows)} columns",
          flush=True)
    D, _ = static_station_metric(G, node_off, extra_edges)
    if args.arc_pricing == "exact":
        arcs = _exact_connector_arcs(
            G,
            tables,
            rows,
            runs,
            D,
            top_k=args.top_k,
            max_proxy_s=int(args.max_proxy_hours * 3600),
            max_connector_s=int(args.max_connector_hours * 3600),
            cap=args.connector_cap,
            run_mode=run_mode,
            run_radius=args.run_radius,
            cache_path=None if args.no_exact_arc_cache else args.exact_arc_cache,
            min_same_station_gap_s=args.min_same_station_gap,
            min_opposite_direction_gap_s=args.min_opposite_direction_gap,
        )
    else:
        arcs = _candidate_exact_arcs(
            rows,
            D,
            top_k=args.top_k,
            max_proxy_s=int(args.max_proxy_hours * 3600),
            min_same_station_gap_s=args.min_same_station_gap,
            min_opposite_direction_gap_s=args.min_opposite_direction_gap,
        )
    if args.force_solution_arcs:
        forced_arcs = _forced_solution_arcs(
            G,
            tables,
            rows,
            runs,
            args.force_solution_arcs,
            cap=args.connector_cap,
            max_connector_s=int(args.max_connector_hours * 3600),
        )
        before = len(arcs)
        arcs = _merge_arcs_min_cost(arcs, forced_arcs)
        print(f"merged forced arcs: {before} + {len(forced_arcs)} -> {len(arcs)}",
              flush=True)
    print(f"exact-cover model: columns={len(rows)} arcs={len(arcs)} "
          f"top_k={args.top_k}", flush=True)
    if not arcs:
        raise SystemExit("no exact connector arcs survived filtering")

    start_window = _parse_week_time_window(args.start_time_window)
    end_window = _parse_week_time_window(args.end_time_window)
    max_total_elapsed_s = _parse_duration_seconds(args.max_total_elapsed)
    allowed_start_rows = _rows_in_time_window(rows, "start", start_window)
    allowed_end_rows = _rows_in_time_window(rows, "end", end_window)
    protected_stations = _load_exact_cover_station_set(args, si)
    protected_station_groups = _load_exact_cover_station_groups(args, si)
    uncovered_penalty_groups = _load_exact_cover_penalty_groups(args, si)
    ordered_station_groups = _load_ordered_station_groups(
        args.order_station_groups, si)
    protected_row_start_stops = {
        stop.strip()
        for stop in str(args.protected_row_start_stop or "").split(",")
        if stop.strip()
    }
    protected_row_prefixes = _column_id_prefixes(args.protected_row_column_id_prefix)
    protected_filter_stations = set(protected_stations)
    for group in protected_station_groups:
        protected_filter_stations.update(group["stations"])
    protected_row_allowed_indices = []
    protected_row_forbidden_indices = []
    protected_row_filter_active = bool(protected_row_start_stops or protected_row_prefixes)
    if protected_row_filter_active and not protected_filter_stations:
        print("--protected-row-* filters ignored; no protected stations/groups",
              flush=True)
    elif protected_row_filter_active:
        protected_cover_indices = [
            i
            for i, row in enumerate(rows)
            if protected_filter_stations.intersection(row.get("covered_stations", []))
        ]

        def _allowed_protected_row(row):
            if protected_row_start_stops:
                if str(row.get("start", {}).get("stop", "")) not in protected_row_start_stops:
                    return False
            if protected_row_prefixes:
                column_id = str(row.get("column_id", ""))
                if not any(column_id.startswith(prefix) for prefix in protected_row_prefixes):
                    return False
            return True

        protected_row_allowed_indices = [
            i for i in protected_cover_indices if _allowed_protected_row(rows[i])
        ]
        if not protected_row_allowed_indices:
            raise SystemExit(
                "--protected-row-* filters exclude all rows covering protected stations")
        allowed_set = set(protected_row_allowed_indices)
        protected_row_forbidden_indices = [
            i for i in protected_cover_indices if i not in allowed_set
        ]
        print("protected row filters: "
              f"cover_rows={len(protected_cover_indices)} "
              f"allowed={len(protected_row_allowed_indices)} "
              f"forbidden={len(protected_row_forbidden_indices)}",
              flush=True)
    row_by_id = {row.get("column_id"): i for i, row in enumerate(rows)}
    required_column_ids = [
        column_id.strip()
        for column_id in str(args.require_column_id or "").split(",")
        if column_id.strip()
    ]
    required_row_indices = []
    required_column_id_prefixes = [
        prefix.strip()
        for prefix in str(args.require_column_id_prefix or "").split(",")
        if prefix.strip()
    ]
    required_prefix_groups = _load_prefix_requirement_groups(
        args.require_column_id_prefix_group)
    required_pricing_indices = []
    required_index_groups = []
    required_neighbor_indices = []
    if required_column_ids:
        missing_required = [
            column_id for column_id in required_column_ids if column_id not in row_by_id
        ]
        if missing_required:
            raise SystemExit(
                f"--require-column-id missing ids: {missing_required[:8]}")
        required_row_indices = [row_by_id[column_id] for column_id in required_column_ids]
        required_neighbor_indices.extend(required_row_indices)
        print(f"required column ids: {len(required_row_indices)}", flush=True)
    for prefix in required_column_id_prefixes:
        group = [
            i for i, row in enumerate(rows)
            if str(row.get("column_id", "")).startswith(prefix)
        ]
        min_count = int(args.min_required_column_prefix)
        if len(group) < min_count:
            raise SystemExit(
                f"--require-column-id-prefix {prefix!r} matched {len(group)} "
                f"column(s), cannot require {min_count}")
        required_index_groups.append((
            f"require_prefix_{len(required_index_groups)}",
            group,
            min_count,
        ))
        required_neighbor_indices.extend(group)
        print(f"required column prefix {prefix!r}: matched={len(group)} "
              f"min={min_count}", flush=True)
    for group_spec in required_prefix_groups:
        group = [
            i for i, row in enumerate(rows)
            if any(
                str(row.get("column_id", "")).startswith(prefix)
                for prefix in group_spec["prefixes"]
            )
        ]
        min_count = int(group_spec["min_count"])
        if len(group) < min_count:
            raise SystemExit(
                f"--require-column-id-prefix-group {group_spec['label']!r} "
                f"matched {len(group)} column(s), cannot require {min_count}")
        required_index_groups.append((
            f"require_prefix_group_{len(required_index_groups)}",
            group,
            min_count,
        ))
        required_neighbor_indices.extend(group)
        print(
            f"required column prefix group {group_spec['label']!r}: "
            f"prefixes={len(group_spec['prefixes'])} matched={len(group)} "
            f"min={min_count}",
            flush=True,
        )
    if args.require_pricing_kind and args.min_required_pricing:
        required_pricing_indices = [
            i for i, row in enumerate(rows)
            if args.require_pricing_kind in row
        ]
        if len(required_pricing_indices) < int(args.min_required_pricing):
            raise SystemExit(
                f"--require-pricing-kind {args.require_pricing_kind!r} matched "
                f"{len(required_pricing_indices)} column(s), cannot require "
                f"{args.min_required_pricing}")
        required_neighbor_indices.extend(required_pricing_indices)
        print(f"required pricing kind {args.require_pricing_kind!r}: "
              f"matched={len(required_pricing_indices)} "
              f"min={args.min_required_pricing}", flush=True)
    protected_neighbor_indices = []
    if args.protected_arc_top_k:
        protected_arc_stations = set(protected_stations)
        for group in protected_station_groups:
            protected_arc_stations.update(group["stations"])
        if protected_arc_stations:
            allowed_protected_set = set(protected_row_allowed_indices)
            protected_neighbor_indices = sorted({
                i
                for i, row in enumerate(rows)
                if protected_arc_stations.intersection(row.get("covered_stations", []))
                and (
                    not protected_row_filter_active
                    or i in allowed_protected_set
                )
            })
            print(f"protected-neighbor stations={len(protected_arc_stations)} "
                  f"matched_rows={len(protected_neighbor_indices)} "
                  f"top_k={args.protected_arc_top_k}", flush=True)
        else:
            print("--protected-arc-top-k ignored; no protected stations/groups",
                  flush=True)
    if required_neighbor_indices and args.required_arc_top_k:
        required_arcs = _required_neighbor_exact_arcs(
            G,
            tables,
            rows,
            runs,
            D,
            required_neighbor_indices,
            top_k=args.required_arc_top_k,
            max_proxy_s=int(args.max_proxy_hours * 3600),
            max_connector_s=int(args.max_connector_hours * 3600),
            cap=args.connector_cap,
            run_mode=run_mode,
            run_radius=args.run_radius,
            cache_path=None if args.no_exact_arc_cache else args.exact_arc_cache,
            min_same_station_gap_s=args.min_same_station_gap,
            min_opposite_direction_gap_s=args.min_opposite_direction_gap,
        )
        before = len(arcs)
        arcs = _merge_arcs_min_cost(arcs, required_arcs)
        print(f"merged required-neighbor arcs: {before} + "
              f"{len(required_arcs)} -> {len(arcs)}", flush=True)
    if protected_neighbor_indices and args.protected_arc_top_k:
        protected_arcs = _required_neighbor_exact_arcs(
            G,
            tables,
            rows,
            runs,
            D,
            protected_neighbor_indices,
            top_k=args.protected_arc_top_k,
            max_proxy_s=int(args.max_proxy_hours * 3600),
            max_connector_s=int(args.max_connector_hours * 3600),
            cap=args.connector_cap,
            run_mode=run_mode,
            run_radius=args.run_radius,
            cache_path=None if args.no_exact_arc_cache else args.exact_arc_cache,
            min_same_station_gap_s=args.min_same_station_gap,
            min_opposite_direction_gap_s=args.min_opposite_direction_gap,
        )
        before = len(arcs)
        arcs = _merge_arcs_min_cost(arcs, protected_arcs)
        print(f"merged protected-neighbor arcs: {before} + "
              f"{len(protected_arcs)} -> {len(arcs)}", flush=True)
    hint_order_indices = []
    ignored_hint_solutions = []
    if args.hint_solution:
        hint_paths = [Path(path) for path in args.hint_solution]
        existing_hint_paths = [path for path in hint_paths if path.exists()]
        missing_hint_files = [str(path) for path in hint_paths if not path.exists()]
        if not existing_hint_paths:
            print("hint solution: no existing hint file found; continuing without hint",
                  flush=True)
            if missing_hint_files:
                print(f"hint files missing e.g. {missing_hint_files[:3]}", flush=True)
        else:
            hint_path = existing_hint_paths[0]
            ignored_hint_solutions = [
                str(path)
                for path in hint_paths
                if path != hint_path
            ]
            if ignored_hint_solutions:
                print(f"hint solution: using {hint_path}; ignoring "
                      f"{len(ignored_hint_solutions)} extra/missing hint file(s)",
                      flush=True)
            if missing_hint_files:
                print(f"hint files missing e.g. {missing_hint_files[:3]}", flush=True)
            hint_column_ids = _selected_column_ids_from_solution(hint_path)
            seen_hint_ids = set()
            hint_column_ids = [
                column_id
                for column_id in hint_column_ids
                if not (column_id in seen_hint_ids or seen_hint_ids.add(column_id))
            ]
            missing_hint = [
                column_id for column_id in hint_column_ids if column_id not in row_by_id
            ]
            hint_order_indices = [
                row_by_id[column_id]
                for column_id in hint_column_ids
                if column_id in row_by_id
            ]
            print(f"hint solution columns: {len(hint_order_indices)} "
                  f"missing={len(missing_hint)}", flush=True)
            if missing_hint:
                print(f"hint missing e.g. {missing_hint[:8]}", flush=True)
    if start_window is not None:
        print(f"start time window {start_window}: {len(allowed_start_rows)} eligible columns",
              flush=True)
    if end_window is not None:
        print(f"end time window {end_window}: {len(allowed_end_rows)} eligible columns",
              flush=True)
    if max_total_elapsed_s is not None:
        print(f"max modeled elapsed: {hms(max_total_elapsed_s)}", flush=True)
    if args.min_covered_count is not None:
        print(f"min covered count: {args.min_covered_count}", flush=True)
    if protected_stations:
        print(f"protected stations: {len(protected_stations)}", flush=True)
    if protected_station_groups:
        print("protected station groups: "
              + ", ".join(
                  f"{group['label']}={group['min_hits']}/{len(group['stations'])}"
                  for group in protected_station_groups),
              flush=True)
    if uncovered_penalty_groups:
        print("uncovered penalty groups: "
              + ", ".join(
                  f"{group['label']}=+{group['penalty_s']}s/"
                  f"{len(group['stations'])}"
                  for group in uncovered_penalty_groups),
              flush=True)
    if ordered_station_groups:
        print("ordered station groups: "
              + " -> ".join(
                  f"{group['label']}={group['min_hits']}/{len(group['stations'])}"
                  for group in ordered_station_groups),
              flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    banned = set()
    exact_cost_cache = {}
    max_attempts = 1
    if args.arc_pricing == "proxy":
        max_attempts += max(0, int(args.repair_cuts))

    for attempt in range(1, max_attempts + 1):
        active_arcs = [arc for arc in arcs if (arc[1], arc[2]) not in banned]
        if not active_arcs:
            raise SystemExit("all exact-cover arcs were cut")
        if banned:
            print(f"exact-cover attempt {attempt}: active_arcs={len(active_arcs)} "
                  f"cuts={len(banned)}", flush=True)

        solver, status, order, solve_vars = _solve_path_cover(
            rows,
            active_arcs,
            args.time_limit,
            args.workers,
            args.column_penalty,
            require_pricing_kind=args.require_pricing_kind,
            min_required_pricing=args.min_required_pricing,
            required_row_indices=required_row_indices,
            required_index_groups=required_index_groups,
            allowed_start_rows=allowed_start_rows,
            allowed_end_rows=allowed_end_rows,
            uncovered_penalty_s=args.uncovered_penalty_s,
            min_covered_count=args.min_covered_count,
            max_total_elapsed_s=max_total_elapsed_s,
            hint_order_indices=hint_order_indices,
            protected_stations=protected_stations,
            protected_station_groups=protected_station_groups,
            uncovered_penalty_groups=uncovered_penalty_groups,
            ordered_station_groups=ordered_station_groups,
            strict_order_station_groups=args.strict_order_station_groups,
            first_hit_order_station_groups=args.first_hit_order_station_groups,
            stop_after_first_solution=args.stop_after_first_solution,
            forbidden_row_indices=protected_row_forbidden_indices,
        )
        status_name = solver.StatusName(status)
        if not order:
            print(f"status={status_name}")
            return 1

        ordered_rows = [rows[i] for i in order]
        coverage = set()
        modes = Counter()
        for row in ordered_rows:
            coverage.update(row.get("covered_stations", []))
            modes.update(row.get("modes", {}))
        relaxed_uncovered = []
        if args.uncovered_penalty_s or uncovered_penalty_groups:
            relaxed_uncovered = sorted(
                station
                for station, var in solve_vars.get("covered", {}).items()
                if not solver.Value(var)
            )
        raw_elapsed = sum(int(row["elapsed_s"]) for row in ordered_rows)
        modeled_elapsed = None
        if "elapsed_var" in solve_vars:
            modeled_elapsed = int(solver.Value(solve_vars["elapsed_var"]))
        active_arc_cost = {}
        for cost, i, j in active_arcs:
            pair = (i, j)
            if pair not in active_arc_cost or int(cost) < active_arc_cost[pair]:
                active_arc_cost[pair] = int(cost)
        selected_arcs = []
        for i, j in zip(order, order[1:]):
            cost = active_arc_cost.get((i, j))
            selected_arcs.append({
                "from_index": int(i),
                "to_index": int(j),
                "from_column": rows[i].get("column_id"),
                "to_column": rows[j].get("column_id"),
                "connector_elapsed_s": cost,
                "connector_elapsed": hms(cost) if cost is not None else None,
            })
        ordered_group_results = []
        for group_var in solve_vars.get("ordered_group_vars", []):
            witness_rows = []
            for i, var in group_var["witness"]:
                if solver.Value(var):
                    witness_rows.append({
                        "index": int(i),
                        "column_id": rows[i].get("column_id"),
                        "order": int(solver.Value(group_var["position"])),
                        "start": rows[i].get("start"),
                        "end": rows[i].get("end"),
                        "hits": sorted(
                            set(rows[i].get("covered_stations", []))
                            & set(group_var["group"]["stations"])
                        ),
                    })
            ordered_group_results.append({
                "label": group_var["group"]["label"],
                "min_hits": int(group_var["group"]["min_hits"]),
                "stations": group_var["group"]["stations"],
                "matched_row_count": len(group_var["candidates"]),
                "witness_rows": witness_rows,
            })
        result = {
            "status": status_name,
            "objective_s": int(solver.ObjectiveValue()),
            "objective": hms(int(solver.ObjectiveValue())),
            "selected_count": len(ordered_rows),
            "covered_count": len(coverage & si.canonical_stations),
            "raw_column_elapsed_s": raw_elapsed,
            "raw_column_elapsed": hms(raw_elapsed),
            "modeled_elapsed_s": modeled_elapsed,
            "modeled_elapsed": hms(modeled_elapsed) if modeled_elapsed is not None else None,
            "max_total_elapsed_s": max_total_elapsed_s,
            "max_total_elapsed": (
                hms(max_total_elapsed_s) if max_total_elapsed_s is not None else None
            ),
            "mode_edges": dict(modes),
            "proxy_cuts": [[rows[i]["column_id"], rows[j]["column_id"]]
                           for i, j in sorted(banned)],
            "uncovered_penalty_s": int(args.uncovered_penalty_s),
            "min_covered_count": args.min_covered_count,
            "relaxed_uncovered_stations": relaxed_uncovered,
            "protected_stations": protected_stations,
            "protected_station_groups": protected_station_groups,
            "uncovered_penalty_groups": uncovered_penalty_groups,
            "ordered_station_groups": ordered_group_results,
            "strict_order_station_groups": bool(args.strict_order_station_groups),
            "first_hit_order_station_groups": bool(args.first_hit_order_station_groups),
            "protected_row_start_stops": sorted(protected_row_start_stops),
            "protected_row_column_id_prefixes": protected_row_prefixes,
            "protected_row_allowed_count": len(protected_row_allowed_indices),
            "protected_row_forbidden_count": len(protected_row_forbidden_indices),
            "start_time_window": list(start_window) if start_window is not None else None,
            "end_time_window": list(end_window) if end_window is not None else None,
            "required_column_ids": required_column_ids,
            "required_column_id_prefixes": required_column_id_prefixes,
            "min_required_column_prefix": int(args.min_required_column_prefix),
            "excluded_column_count": len(excluded_rows),
            "excluded_columns": excluded_rows,
            "hint_solution": list(args.hint_solution or []),
            "ignored_hint_solutions": ignored_hint_solutions,
            "selected_order": [row.get("column_id") for row in ordered_rows],
            "selected_arcs": selected_arcs,
            "selected_columns": ordered_rows,
        }
        out.write_text(json.dumps(result, sort_keys=True))
        print(f"status={status_name} selected={len(ordered_rows)} "
              f"covered={result['covered_count']}/472 raw={hms(raw_elapsed)} "
              f"modeled={result['modeled_elapsed']} objective={result['objective']}")
        if relaxed_uncovered:
            print(f"relaxed uncovered={len(relaxed_uncovered)} "
                  f"e.g. {relaxed_uncovered[:12]}")
        print(f"wrote {out}")

        if not args.route_out:
            return 0

        if args.refine_selected_arcs and args.arc_pricing == "proxy":
            changed = False
            for i, j in zip(order, order[1:]):
                pair = (i, j)
                if pair not in exact_cost_cache:
                    exact_cost_cache[pair] = _exact_connector_cost(
                        G, tables, rows, runs, i, j, args.connector_cap)
                exact_cost = exact_cost_cache[pair]
                if exact_cost is None:
                    banned.add(pair)
                    print(f"refined arc infeasible: {rows[i]['column_id']} -> "
                          f"{rows[j]['column_id']}", flush=True)
                    changed = True
                    continue
                current_cost = next((c for c, a, b in active_arcs if (a, b) == pair), None)
                if current_cost is not None and exact_cost > current_cost:
                    arcs, replaced = _replace_arc_cost(arcs, pair, exact_cost)
                    if replaced:
                        print(f"refined arc cost: {rows[i]['column_id']} -> "
                              f"{rows[j]['column_id']} {hms(current_cost)} -> "
                              f"{hms(exact_cost)}", flush=True)
                        changed = True
            if changed and attempt < max_attempts:
                continue

        path, failed = _realize_exact_column_path(
            G, tables, rows, order, runs, args.connector_cap)
        if path is None:
            if failed is not None and attempt < max_attempts and args.arc_pricing == "proxy":
                banned.add(failed)
                i, j = failed
                print(f"added proxy cut: {rows[i]['column_id']} -> {rows[j]['column_id']}",
                      flush=True)
                continue
            raise SystemExit("failed to realize exact column path")
        elapsed = _elapsed(path[0][1], path[-1][1])
        meta = {
            "radius_m": args.validation_radius,
            "elapsed_s": int(elapsed),
            "notes": ("Exact-time column path-cover "
                      f"source={args.columns_source} top_k={args.top_k} "
                      f"run_mode={run_mode} "
                      f"pricing={args.arc_pricing} "
                      f"objective={int(solver.ObjectiveValue())}"),
            "exact_cover_objective_s": int(solver.ObjectiveValue()),
            "exact_cover_arc_pricing": args.arc_pricing,
            "exact_cover_run_mode": run_mode,
            "exact_cover_columns": [row["column_id"] for row in ordered_rows],
        }
        route_out = Path(args.route_out)
        route_out.parent.mkdir(parents=True, exist_ok=True)
        route_out.write_text(json.dumps({"meta": meta, "path": path}, sort_keys=True))
        print(f"realized exact path: steps={len(path)} -> {route_out}", flush=True)
        if args.validate:
            full_runs = RunLayer.from_graph(G, radius_m=meta["radius_m"])
            validation = validate(Problem(G, full_runs, si), path)
            print(result_line(validation), flush=True)
        return 0
    return 1


def _source_order_elapsed_s(data, rows, order):
    if data.get("modeled_elapsed_s") is not None:
        return int(data["modeled_elapsed_s"])
    meta = data.get("meta")
    if isinstance(meta, dict) and meta.get("elapsed_s") is not None:
        return int(meta["elapsed_s"])
    total = sum(int(rows[i]["elapsed_s"]) for i in order)
    for i, j in zip(order, order[1:]):
        total += _elapsed(rows[i]["end"]["time"], rows[j]["start"]["time"])
    return int(total)


def _replacement_route_is_better(validation, best_validation,
                                 min_covered_count, max_total_elapsed_s):
    if validation is None:
        return False
    if validation["covered"] < min_covered_count:
        return False
    if validation["elapsed_s"] > max_total_elapsed_s:
        return False
    if best_validation is None:
        return True
    return (
        int(validation["covered"]),
        -int(validation["elapsed_s"]),
    ) > (
        int(best_validation["covered"]),
        -int(best_validation["elapsed_s"]),
    )


def cmd_block_replace(args) -> int:
    import heapq
    import pickle

    rows = _load_column_rows(Path(args.columns_source))
    if not rows:
        raise SystemExit("no columns to search")
    row_by_id = {row.get("column_id"): i for i, row in enumerate(rows)}
    source_data = json.loads(Path(args.source_order).read_text())
    selected_ids = _selected_column_ids_from_solution(Path(args.source_order))
    missing_selected = [column_id for column_id in selected_ids if column_id not in row_by_id]
    if missing_selected:
        raise SystemExit(f"source order has ids not in pool, e.g. {missing_selected[:8]}")
    selected = [row_by_id[column_id] for column_id in selected_ids]
    if len(selected) < 2:
        raise SystemExit("source order must contain at least two columns")
    selected_set = set(selected)

    si = StationIndex.load()
    canonical = set(si.canonical_stations)
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off = tables[0]
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    D, _ = static_station_metric(G, node_off, extra_edges)
    validation_runs = RunLayer.from_graph(G, radius_m=args.validation_radius)
    prob = Problem(G, validation_runs, si)

    row_cov = [set(row.get("covered_stations", [])) & canonical for row in rows]
    selected_cov = [row_cov[i] for i in selected]
    base_coverage = set().union(*selected_cov)
    target_stations = set(
        _parse_station_list(args.target_stations, si, "--target-stations")
    ) if args.target_stations else (
        canonical - base_coverage
    )
    if not target_stations:
        raise SystemExit("no target stations to improve")
    base_elapsed_s = _source_order_elapsed_s(source_data, rows, selected)
    max_total_elapsed_s = _parse_duration_seconds(args.max_total_elapsed)
    if max_total_elapsed_s is None:
        max_total_elapsed_s = base_elapsed_s
    min_covered_count = (
        len(base_coverage) + 1
        if args.min_covered_count is None
        else int(args.min_covered_count)
    )

    n = len(selected)
    row_elapsed = [int(rows[i]["elapsed_s"]) for i in selected]
    arc_elapsed = [
        _elapsed(rows[i]["end"]["time"], rows[j]["start"]["time"])
        for i, j in zip(selected, selected[1:])
    ]
    prefix_rows = [0]
    for value in row_elapsed:
        prefix_rows.append(prefix_rows[-1] + value)
    prefix_arcs = [0]
    for value in arc_elapsed:
        prefix_arcs.append(prefix_arcs[-1] + value)
    prefix_cov = []
    seen = set()
    for cov in selected_cov:
        seen = seen | cov
        prefix_cov.append(set(seen))
    suffix_cov = [set() for _ in range(n)]
    seen = set()
    for idx in range(n - 1, -1, -1):
        seen = seen | selected_cov[idx]
        suffix_cov[idx] = set(seen)

    candidates = []
    checked = 0
    for k, row in enumerate(rows):
        if k in selected_set:
            continue
        hits = target_stations & row_cov[k]
        if not hits:
            continue
        for p in range(n):
            prev = rows[selected[p]]
            static_in = (
                0 if prev["end"]["station"] == row["start"]["station"]
                else int(D.get(prev["end"]["station"], {}).get(row["start"]["station"], BIG_COST))
            )
            if static_in >= BIG_COST:
                continue
            conn_in = _phase_proxy_gap(prev["end"]["time"], row["start"]["time"], static_in)
            if conn_in >= WEEK:
                continue
            for q in range(p + 1, n + 1):
                conn_out = 0
                if q < n:
                    nxt = rows[selected[q]]
                    static_out = (
                        0 if row["end"]["station"] == nxt["start"]["station"]
                        else int(D.get(row["end"]["station"], {}).get(
                            nxt["start"]["station"], BIG_COST))
                    )
                    if static_out >= BIG_COST:
                        continue
                    conn_out = _phase_proxy_gap(
                        row["end"]["time"], nxt["start"]["time"], static_out)
                    if conn_out >= WEEK:
                        continue
                checked += 1
                arc_prefix_end = q if q < n else n - 1
                old_span = (
                    prefix_rows[q] - prefix_rows[p + 1]
                    + prefix_arcs[arc_prefix_end] - prefix_arcs[p]
                )
                new_span = conn_in + int(row["elapsed_s"]) + conn_out
                new_elapsed_s = base_elapsed_s - old_span + new_span
                if new_elapsed_s > max_total_elapsed_s:
                    continue
                coverage = prefix_cov[p] | (suffix_cov[q] if q < n else set()) | row_cov[k]
                covered_count = len(coverage)
                if covered_count < min_covered_count:
                    continue
                candidates.append((
                    -covered_count,
                    int(new_elapsed_s),
                    int(new_span - old_span),
                    q - p - 1,
                    k,
                    p,
                    q,
                    sorted(hits),
                    sorted(canonical - coverage),
                ))

    candidates = heapq.nsmallest(max(0, int(args.max_candidates)), candidates)
    print(f"block-replace base covered={len(base_coverage)}/472 "
          f"elapsed={hms(base_elapsed_s)} target_stations={len(target_stations)}",
          flush=True)
    print(f"proxy checked={checked} candidates={len(candidates)} "
          f"min_covered={min_covered_count} max_elapsed={hms(max_total_elapsed_s)}",
          flush=True)

    summary = {
        "source_order": args.source_order,
        "columns_source": args.columns_source,
        "base_covered_count": len(base_coverage),
        "base_elapsed_s": base_elapsed_s,
        "base_elapsed": hms(base_elapsed_s),
        "target_stations": sorted(target_stations),
        "min_covered_count": min_covered_count,
        "max_total_elapsed_s": max_total_elapsed_s,
        "max_total_elapsed": hms(max_total_elapsed_s),
        "run_mode": run_mode,
        "run_radius": args.run_radius,
        "validation_radius": args.validation_radius,
        "candidates": [],
    }

    out_route = Path(args.route_out) if args.route_out else None
    best_route_validation = None
    for rank, rec in enumerate(candidates, start=1):
        neg_covered, proxy_elapsed_s, proxy_delta_s, removed_count, k, p, q, hits, missing = rec
        order = selected[:p + 1] + [k] + selected[q:]
        path, failed = _realize_exact_column_path(
            G, tables, rows, order, validation_runs, args.connector_cap)
        validation = None
        route_path = None
        if path is not None:
            validation = validate(prob, path)
            print(f"candidate {rank}: {result_line(validation)} "
                  f"replacement={rows[k]['column_id']}", flush=True)
            if (
                out_route is not None
                and _replacement_route_is_better(
                    validation,
                    best_route_validation,
                    min_covered_count,
                    max_total_elapsed_s,
                )
            ):
                out_route.parent.mkdir(parents=True, exist_ok=True)
                meta = {
                    "radius_m": args.validation_radius,
                    "elapsed_s": int(_elapsed(path[0][1], path[-1][1])),
                    "notes": "One-column block replacement from selected column order",
                    "source_order": args.source_order,
                    "columns_source": args.columns_source,
                    "replacement_column": rows[k]["column_id"],
                    "replaced_after": rows[selected[p]]["column_id"],
                    "replaced_before": (
                        rows[selected[q]]["column_id"] if q < n else None
                    ),
                    "removed_columns": [rows[idx]["column_id"] for idx in selected[p + 1:q]],
                    "proxy_covered_count": -neg_covered,
                    "proxy_missing_stations": missing,
                    "proxy_elapsed_s": proxy_elapsed_s,
                    "proxy_elapsed": hms(proxy_elapsed_s),
                    "selected_order": [rows[idx]["column_id"] for idx in order],
                }
                out_route.write_text(json.dumps({"meta": meta, "path": path}, sort_keys=True))
                route_path = str(out_route)
                best_route_validation = validation
                print(f"wrote {out_route}", flush=True)
                if args.stop_after_first:
                    summary["candidates"].append({
                        "rank": rank,
                        "replacement_column": rows[k]["column_id"],
                        "proxy_covered_count": -neg_covered,
                        "proxy_elapsed_s": proxy_elapsed_s,
                        "proxy_delta_s": proxy_delta_s,
                        "removed_count": removed_count,
                        "target_hits": hits,
                        "failed_arc": failed,
                        "validation": validation,
                        "route_out": route_path,
                    })
                    break
        else:
            print(f"candidate {rank}: failed replacement={rows[k]['column_id']} "
                  f"failed_arc={failed}", flush=True)
        summary["candidates"].append({
            "rank": rank,
            "replacement_column": rows[k]["column_id"],
            "replaced_after": rows[selected[p]]["column_id"],
            "replaced_before": rows[selected[q]]["column_id"] if q < n else None,
            "removed_columns": [rows[idx]["column_id"] for idx in selected[p + 1:q]],
            "proxy_covered_count": -neg_covered,
            "proxy_elapsed_s": proxy_elapsed_s,
            "proxy_elapsed": hms(proxy_elapsed_s),
            "proxy_delta_s": proxy_delta_s,
            "removed_count": removed_count,
            "target_hits": hits,
            "proxy_missing_stations": missing,
            "failed_arc": list(failed) if failed is not None else None,
            "validation": validation,
            "route_out": route_path,
        })

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(summary, sort_keys=True))
        print(f"wrote {out}", flush=True)
    return 0


def _column_id_prefixes(value):
    return [prefix.strip() for prefix in str(value or "").split(",") if prefix.strip()]


def _matches_any_prefix(row, prefixes):
    if not prefixes:
        return True
    column_id = str(row.get("column_id", ""))
    return any(column_id.startswith(prefix) for prefix in prefixes)


def _replacement_source_state(args):
    rows = _load_column_rows(Path(args.columns_source))
    if not rows:
        raise SystemExit("no columns to search")
    row_by_id = {row.get("column_id"): i for i, row in enumerate(rows)}
    source_data = json.loads(Path(args.source_order).read_text())
    selected_ids = _selected_column_ids_from_solution(Path(args.source_order))
    missing_selected = [column_id for column_id in selected_ids if column_id not in row_by_id]
    if missing_selected:
        raise SystemExit(f"source order has ids not in pool, e.g. {missing_selected[:8]}")
    selected = [row_by_id[column_id] for column_id in selected_ids]
    if len(selected) < 2:
        raise SystemExit("source order must contain at least two columns")
    missing_path = [
        row.get("column_id", str(i))
        for i, row in enumerate(rows)
        if "path" not in row
    ]
    if missing_path:
        raise SystemExit(f"{len(missing_path)} columns have no path slices, e.g. {missing_path[:8]}")
    return rows, source_data, selected


def _selected_route_series(rows, selected):
    row_elapsed = [int(rows[i]["elapsed_s"]) for i in selected]
    arc_elapsed = [
        _elapsed(rows[i]["end"]["time"], rows[j]["start"]["time"])
        for i, j in zip(selected, selected[1:])
    ]
    prefix_rows = [0]
    for value in row_elapsed:
        prefix_rows.append(prefix_rows[-1] + value)
    prefix_arcs = [0]
    for value in arc_elapsed:
        prefix_arcs.append(prefix_arcs[-1] + value)
    return prefix_rows, prefix_arcs


def _old_selected_span_elapsed(prefix_rows, prefix_arcs, n, start, end):
    p = start - 1
    arc_prefix_end = end if end < n else n - 1
    return (
        prefix_rows[end] - prefix_rows[start]
        + prefix_arcs[arc_prefix_end] - prefix_arcs[p]
    )


def cmd_pair_replace(args) -> int:
    import heapq
    import pickle

    rows, source_data, selected = _replacement_source_state(args)
    selected_set = set(selected)

    si = StationIndex.load()
    canonical = set(si.canonical_stations)
    row_cov = [set(row.get("covered_stations", [])) & canonical for row in rows]
    selected_cov = [row_cov[i] for i in selected]
    base_coverage = set().union(*selected_cov)
    target_stations = set(
        _parse_station_list(args.target_stations, si, "--target-stations")
    ) if args.target_stations else (
        canonical - base_coverage
    )
    if not target_stations:
        raise SystemExit("no target stations to improve")
    base_elapsed_s = _source_order_elapsed_s(source_data, rows, selected)
    max_total_elapsed_s = _parse_duration_seconds(args.max_total_elapsed)
    if max_total_elapsed_s is None:
        max_total_elapsed_s = base_elapsed_s
    min_covered_count = (
        len(base_coverage) + 1
        if args.min_covered_count is None
        else int(args.min_covered_count)
    )

    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off = tables[0]
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    D, _ = static_station_metric(G, node_off, extra_edges)
    validation_runs = RunLayer.from_graph(G, radius_m=args.validation_radius)
    prob = Problem(G, validation_runs, si)

    n = len(selected)
    prefix_rows, prefix_arcs = _selected_route_series(rows, selected)
    prefix_cov = []
    seen = set()
    for cov in selected_cov:
        seen = seen | cov
        prefix_cov.append(set(seen))
    suffix_cov = [set() for _ in range(n)]
    seen = set()
    for idx in range(n - 1, -1, -1):
        seen = seen | selected_cov[idx]
        suffix_cov[idx] = set(seen)

    shared_prefixes = _column_id_prefixes(args.column_id_prefix)
    first_prefixes = _column_id_prefixes(args.first_column_id_prefix) or shared_prefixes
    second_prefixes = _column_id_prefixes(args.second_column_id_prefix) or shared_prefixes

    def candidate_rows(prefixes, max_rows):
        candidates = []
        for i, row in enumerate(rows):
            if i in selected_set or not _matches_any_prefix(row, prefixes):
                continue
            hits = row_cov[i] & target_stations
            if not hits and not args.include_non_target_candidates:
                continue
            candidates.append((
                -len(hits),
                -len(row_cov[i]),
                int(row["elapsed_s"]),
                str(row.get("column_id", "")),
                i,
            ))
        candidates.sort()
        if max_rows:
            candidates = candidates[:max_rows]
        return [i for *_sort, i in candidates]

    first_candidates = candidate_rows(first_prefixes, args.max_first_rows)
    second_candidates = candidate_rows(second_prefixes, args.max_second_rows)
    if not first_candidates or not second_candidates:
        raise SystemExit(
            f"no candidate rows after filtering: first={len(first_candidates)} "
            f"second={len(second_candidates)}")

    min_start = max(1, int(args.min_replace_start_index))
    max_start = int(args.max_replace_start_index or (n - 1))
    max_start = min(max_start, n - 1)
    min_end_arg = int(args.min_replace_end_index or 0)
    max_end = int(args.max_replace_end_index or n)
    max_end = min(max_end, n)
    if min_start > max_start:
        raise SystemExit("replacement start-index window is empty")

    proxy_cache = {}

    def proxy(i, j):
        key = (i, j)
        if key not in proxy_cache:
            proxy_cache[key] = _candidate_exact_arc_proxy(
                rows,
                D,
                i,
                j,
                min_same_station_gap_s=args.min_same_station_gap,
                min_opposite_direction_gap_s=args.min_opposite_direction_gap,
            )
        return proxy_cache[key]

    candidates = []
    checked = 0
    proxy_feasible = 0
    elapsed_feasible = 0
    for start in range(min_start, max_start + 1):
        end_min = max(start + int(args.min_removed_count), min_end_arg or start + 1)
        end_max = max_end
        if args.max_removed_count:
            end_max = min(end_max, start + int(args.max_removed_count))
        if end_min > end_max:
            continue
        p = start - 1
        prev_idx = selected[p]
        for end in range(end_min, end_max + 1):
            if end <= start or end > n:
                continue
            q_idx = selected[end] if end < n else None
            old_span = _old_selected_span_elapsed(prefix_rows, prefix_arcs, n, start, end)
            outside_cov = prefix_cov[p] | (suffix_cov[end] if end < n else set())
            for first_idx in first_candidates:
                conn_in = proxy(prev_idx, first_idx)
                if conn_in is None:
                    continue
                for second_idx in second_candidates:
                    if second_idx == first_idx:
                        continue
                    checked += 1
                    conn_mid = proxy(first_idx, second_idx)
                    if conn_mid is None:
                        continue
                    conn_out = 0
                    if q_idx is not None:
                        conn_out = proxy(second_idx, q_idx)
                        if conn_out is None:
                            continue
                    proxy_feasible += 1
                    new_span = (
                        int(conn_in)
                        + int(rows[first_idx]["elapsed_s"])
                        + int(conn_mid)
                        + int(rows[second_idx]["elapsed_s"])
                        + int(conn_out)
                    )
                    new_elapsed_s = base_elapsed_s - old_span + new_span
                    if new_elapsed_s > max_total_elapsed_s:
                        continue
                    elapsed_feasible += 1
                    coverage = outside_cov | row_cov[first_idx] | row_cov[second_idx]
                    covered_count = len(coverage)
                    if covered_count < min_covered_count:
                        continue
                    hits = sorted((row_cov[first_idx] | row_cov[second_idx]) & target_stations)
                    candidates.append((
                        -covered_count,
                        int(new_elapsed_s),
                        int(new_span - old_span),
                        end - start,
                        first_idx,
                        second_idx,
                        start,
                        end,
                        hits,
                        sorted(canonical - coverage),
                    ))
        if args.progress_every_starts and (
            (start - min_start + 1) % int(args.progress_every_starts) == 0
            or start == max_start
        ):
            print(f"pair-replace progress start={start}/{max_start} "
                  f"checked={checked} candidates={len(candidates)}",
                  flush=True)

    candidates = heapq.nsmallest(max(0, int(args.max_candidates)), candidates)
    print(f"pair-replace base covered={len(base_coverage)}/472 "
          f"elapsed={hms(base_elapsed_s)} target_stations={len(target_stations)}",
          flush=True)
    print(f"candidate rows first={len(first_candidates)} second={len(second_candidates)} "
          f"start_index={min_start}-{max_start} end_max={max_end}", flush=True)
    print(f"proxy checked={checked} feasible={proxy_feasible} "
          f"elapsed_feasible={elapsed_feasible} candidates={len(candidates)} "
          f"min_covered={min_covered_count} max_elapsed={hms(max_total_elapsed_s)}",
          flush=True)

    summary = {
        "source_order": args.source_order,
        "columns_source": args.columns_source,
        "base_covered_count": len(base_coverage),
        "base_elapsed_s": base_elapsed_s,
        "base_elapsed": hms(base_elapsed_s),
        "target_stations": sorted(target_stations),
        "min_covered_count": min_covered_count,
        "max_total_elapsed_s": max_total_elapsed_s,
        "max_total_elapsed": hms(max_total_elapsed_s),
        "run_mode": run_mode,
        "run_radius": args.run_radius,
        "validation_radius": args.validation_radius,
        "first_candidate_count": len(first_candidates),
        "second_candidate_count": len(second_candidates),
        "checked": checked,
        "proxy_feasible": proxy_feasible,
        "elapsed_feasible": elapsed_feasible,
        "candidates": [],
    }

    out_route = Path(args.route_out) if args.route_out else None
    best_route_validation = None
    for rank, rec in enumerate(candidates, start=1):
        neg_covered, proxy_elapsed_s, proxy_delta_s, removed_count, first_idx, second_idx, start, end, hits, missing = rec
        order = selected[:start] + [first_idx, second_idx] + selected[end:]
        path, failed = _realize_exact_column_path(
            G, tables, rows, order, validation_runs, args.connector_cap)
        validation = None
        route_path = None
        if path is not None:
            validation = validate(prob, path)
            print(f"candidate {rank}: {result_line(validation)} "
                  f"replacement={rows[first_idx]['column_id']} + "
                  f"{rows[second_idx]['column_id']}", flush=True)
            if (
                out_route is not None
                and _replacement_route_is_better(
                    validation,
                    best_route_validation,
                    min_covered_count,
                    max_total_elapsed_s,
                )
            ):
                out_route.parent.mkdir(parents=True, exist_ok=True)
                meta = {
                    "radius_m": args.validation_radius,
                    "elapsed_s": int(_elapsed(path[0][1], path[-1][1])),
                    "notes": "Two-column block replacement from selected column order",
                    "source_order": args.source_order,
                    "columns_source": args.columns_source,
                    "replacement_columns": [
                        rows[first_idx]["column_id"],
                        rows[second_idx]["column_id"],
                    ],
                    "replaced_after": rows[selected[start - 1]]["column_id"],
                    "replaced_before": rows[selected[end]]["column_id"] if end < n else None,
                    "removed_columns": [rows[idx]["column_id"] for idx in selected[start:end]],
                    "proxy_covered_count": -neg_covered,
                    "proxy_missing_stations": missing,
                    "proxy_elapsed_s": proxy_elapsed_s,
                    "proxy_elapsed": hms(proxy_elapsed_s),
                    "selected_order": [rows[idx]["column_id"] for idx in order],
                }
                out_route.write_text(json.dumps({"meta": meta, "path": path}, sort_keys=True))
                route_path = str(out_route)
                best_route_validation = validation
                print(f"wrote {out_route}", flush=True)
                if args.stop_after_first:
                    summary["candidates"].append({
                        "rank": rank,
                        "replacement_columns": [
                            rows[first_idx]["column_id"],
                            rows[second_idx]["column_id"],
                        ],
                        "proxy_covered_count": -neg_covered,
                        "proxy_elapsed_s": proxy_elapsed_s,
                        "proxy_delta_s": proxy_delta_s,
                        "removed_count": removed_count,
                        "target_hits": hits,
                        "failed_arc": failed,
                        "validation": validation,
                        "route_out": route_path,
                    })
                    break
        else:
            print(f"candidate {rank}: failed replacement={rows[first_idx]['column_id']} "
                  f"+ {rows[second_idx]['column_id']} failed_arc={failed}",
                  flush=True)
        summary["candidates"].append({
            "rank": rank,
            "replacement_columns": [
                rows[first_idx]["column_id"],
                rows[second_idx]["column_id"],
            ],
            "replaced_after": rows[selected[start - 1]]["column_id"],
            "replaced_before": rows[selected[end]]["column_id"] if end < n else None,
            "removed_columns": [rows[idx]["column_id"] for idx in selected[start:end]],
            "proxy_covered_count": -neg_covered,
            "proxy_elapsed_s": proxy_elapsed_s,
            "proxy_elapsed": hms(proxy_elapsed_s),
            "proxy_delta_s": proxy_delta_s,
            "removed_count": removed_count,
            "target_hits": hits,
            "proxy_missing_stations": missing,
            "failed_arc": list(failed) if failed is not None else None,
            "validation": validation,
            "route_out": route_path,
        })

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(summary, sort_keys=True))
        print(f"wrote {out}", flush=True)
    return 0


def cmd_chain_replace(args) -> int:
    import heapq
    import pickle

    replacement_count = int(args.replacement_count)
    if replacement_count < 1:
        raise SystemExit("--replacement-count must be positive")

    rows, source_data, selected = _replacement_source_state(args)
    selected_set = set(selected)

    si = StationIndex.load()
    canonical = set(si.canonical_stations)
    row_cov = [set(row.get("covered_stations", [])) & canonical for row in rows]
    selected_cov = [row_cov[i] for i in selected]
    base_coverage = set().union(*selected_cov)
    target_stations = set(
        _parse_station_list(args.target_stations, si, "--target-stations")
    ) if args.target_stations else (
        canonical - base_coverage
    )
    if not target_stations:
        raise SystemExit("no target stations to improve")
    base_elapsed_s = _source_order_elapsed_s(source_data, rows, selected)
    max_total_elapsed_s = _parse_duration_seconds(args.max_total_elapsed)
    if max_total_elapsed_s is None:
        max_total_elapsed_s = base_elapsed_s
    min_covered_count = (
        len(base_coverage) + 1
        if args.min_covered_count is None
        else int(args.min_covered_count)
    )

    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off = tables[0]
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    D, _ = static_station_metric(G, node_off, extra_edges)
    validation_runs = RunLayer.from_graph(G, radius_m=args.validation_radius)
    prob = Problem(G, validation_runs, si)

    n = len(selected)
    prefix_rows, prefix_arcs = _selected_route_series(rows, selected)
    prefix_cov = []
    seen = set()
    for cov in selected_cov:
        seen = seen | cov
        prefix_cov.append(frozenset(seen))
    suffix_cov = [frozenset() for _ in range(n)]
    seen = set()
    for idx in range(n - 1, -1, -1):
        seen = seen | selected_cov[idx]
        suffix_cov[idx] = frozenset(seen)

    candidate_prefixes = _column_id_prefixes(args.column_id_prefix)
    position_prefixes = []
    if args.chain_column_id_prefixes:
        position_prefixes = [
            _column_id_prefixes(part)
            for part in str(args.chain_column_id_prefixes).split(";")
        ]
        if len(position_prefixes) != replacement_count:
            raise SystemExit(
                "--chain-column-id-prefixes must contain one semicolon-separated "
                "prefix group per replacement position")

    def ranked_candidates(prefixes):
        candidate_rows = []
        for i, row in enumerate(rows):
            if i in selected_set or not _matches_any_prefix(row, prefixes):
                continue
            hits = row_cov[i] & target_stations
            if not hits and not args.include_non_target_candidates:
                continue
            candidate_rows.append((
                -len(hits),
                -len(row_cov[i]),
                int(row["elapsed_s"]),
                str(row.get("column_id", "")),
                i,
            ))
        by_index = {}
        for rec in sorted(candidate_rows)[:max(0, int(args.max_candidate_rows))]:
            by_index[rec[-1]] = rec
        if args.extra_fast_candidate_rows:
            for rec in sorted(
                candidate_rows,
                key=lambda item: (item[2], item[0], item[1], item[3]),
            )[:int(args.extra_fast_candidate_rows)]:
                by_index[rec[-1]] = rec
        return [i for *_sort, i in sorted(by_index.values())]

    if position_prefixes:
        candidates_by_depth = [
            ranked_candidates(prefixes)
            for prefixes in position_prefixes
        ]
        empty_depths = [
            idx + 1 for idx, indices in enumerate(candidates_by_depth) if not indices
        ]
        if empty_depths:
            raise SystemExit(
                "--chain-column-id-prefixes produced no candidates at "
                f"position(s) {empty_depths}")
        candidates_idx = sorted(set().union(*map(set, candidates_by_depth)))
    else:
        candidates_idx = ranked_candidates(candidate_prefixes)
        candidates_by_depth = [candidates_idx for _ in range(replacement_count)]
        if len(candidates_idx) < replacement_count:
            raise SystemExit(
                f"need at least {replacement_count} candidate rows after filtering; "
                f"got {len(candidates_idx)}")

    min_start = max(1, int(args.min_replace_start_index))
    max_start = int(args.max_replace_start_index or (n - 1))
    max_start = min(max_start, n - 1)
    min_end_arg = int(args.min_replace_end_index or 0)
    max_end = int(args.max_replace_end_index or n)
    max_end = min(max_end, n)
    if min_start > max_start:
        raise SystemExit("replacement start-index window is empty")

    row_elapsed = [int(row["elapsed_s"]) for row in rows]
    proxy_cache = {}

    def proxy(i, j):
        key = (i, j)
        if key not in proxy_cache:
            proxy_cache[key] = _candidate_exact_arc_proxy(
                rows,
                D,
                i,
                j,
                min_same_station_gap_s=args.min_same_station_gap,
                min_opposite_direction_gap_s=args.min_opposite_direction_gap,
            )
        return proxy_cache[key]

    def state_key(span, cov, chain):
        hits = len(cov & target_stations)
        return (
            int(span)
            - int(args.target_reward_s) * hits
            - int(args.cover_reward_s) * len(cov),
            int(span),
            -hits,
            -len(cov),
            tuple(rows[i].get("column_id", "") for i in chain),
        )

    candidate_heap = []
    checked_extensions = 0
    final_checked = 0
    elapsed_feasible = 0
    blocks_scanned = 0

    for start in range(min_start, max_start + 1):
        end_min = max(start + int(args.min_removed_count), min_end_arg or start + 1)
        end_max = max_end
        if args.max_removed_count:
            end_max = min(end_max, start + int(args.max_removed_count))
        if end_min > end_max:
            continue
        p = start - 1
        prev_left = selected[p]
        for end in range(end_min, end_max + 1):
            if end <= start or end > n:
                continue
            blocks_scanned += 1
            q_idx = selected[end] if end < n else None
            old_span = _old_selected_span_elapsed(prefix_rows, prefix_arcs, n, start, end)
            outside_cov = prefix_cov[p] | (suffix_cov[end] if end < n else frozenset())
            states = [(0, (), outside_cov)]
            for _depth in range(replacement_count):
                next_states = []
                seen_state_keys = set()
                depth_candidates = candidates_by_depth[_depth]
                for span, chain, cov in states:
                    prev_idx = prev_left if not chain else chain[-1]
                    used = set(chain)
                    for row_idx in depth_candidates:
                        if row_idx in used:
                            continue
                        checked_extensions += 1
                        conn = proxy(prev_idx, row_idx)
                        if conn is None:
                            continue
                        new_span = int(span) + int(conn) + row_elapsed[row_idx]
                        if base_elapsed_s - old_span + new_span > max_total_elapsed_s:
                            continue
                        new_chain = chain + (row_idx,)
                        new_cov = cov | frozenset(row_cov[row_idx])
                        dedupe_key = (
                            row_idx,
                            tuple(sorted(new_cov & target_stations)),
                            int(new_span // max(1, int(args.time_bucket_s))),
                        )
                        if dedupe_key in seen_state_keys:
                            continue
                        seen_state_keys.add(dedupe_key)
                        next_states.append((new_span, new_chain, new_cov))
                next_states.sort(key=lambda item: state_key(item[0], item[2], item[1]))
                states = next_states[:int(args.beam_size)]
                if not states:
                    break
            for span, chain, cov in states:
                final_checked += 1
                conn_out = 0
                if q_idx is not None:
                    conn_out = proxy(chain[-1], q_idx)
                    if conn_out is None:
                        continue
                chain_span = int(span) + int(conn_out)
                new_elapsed_s = base_elapsed_s - old_span + chain_span
                if new_elapsed_s > max_total_elapsed_s:
                    continue
                elapsed_feasible += 1
                covered_count = len(cov)
                if covered_count < min_covered_count:
                    continue
                hits = sorted(set().union(*(row_cov[i] for i in chain)) & target_stations)
                rec = (
                    -covered_count,
                    int(new_elapsed_s),
                    int(chain_span - old_span),
                    end - start,
                    chain,
                    start,
                    end,
                    hits,
                    sorted(canonical - set(cov)),
                )
                heapq.heappush(candidate_heap, rec)
                if len(candidate_heap) > max(1, int(args.keep_proxy_candidates)):
                    candidate_heap = heapq.nsmallest(
                        int(args.keep_proxy_candidates), candidate_heap)
                    heapq.heapify(candidate_heap)
        if args.progress_every_starts and (
            (start - min_start + 1) % int(args.progress_every_starts) == 0
            or start == max_start
        ):
            print(f"chain-replace progress start={start}/{max_start} "
                  f"blocks={blocks_scanned} extensions={checked_extensions} "
                  f"proxy_candidates={len(candidate_heap)}",
                  flush=True)

    proxy_candidates = heapq.nsmallest(max(0, int(args.max_candidates)), candidate_heap)
    print(f"chain-replace base covered={len(base_coverage)}/472 "
          f"elapsed={hms(base_elapsed_s)} target_stations={len(target_stations)}",
          flush=True)
    print(f"candidate rows={len(candidates_idx)} replacement_count={replacement_count} "
          f"beam={args.beam_size} blocks={blocks_scanned}", flush=True)
    print(f"extensions={checked_extensions} final_states={final_checked} "
          f"elapsed_feasible={elapsed_feasible} candidates={len(proxy_candidates)} "
          f"min_covered={min_covered_count} max_elapsed={hms(max_total_elapsed_s)}",
          flush=True)

    summary = {
        "source_order": args.source_order,
        "columns_source": args.columns_source,
        "base_covered_count": len(base_coverage),
        "base_elapsed_s": base_elapsed_s,
        "base_elapsed": hms(base_elapsed_s),
        "target_stations": sorted(target_stations),
        "min_covered_count": min_covered_count,
        "max_total_elapsed_s": max_total_elapsed_s,
        "max_total_elapsed": hms(max_total_elapsed_s),
        "run_mode": run_mode,
        "run_radius": args.run_radius,
        "validation_radius": args.validation_radius,
        "candidate_count": len(candidates_idx),
        "position_candidate_counts": [len(indices) for indices in candidates_by_depth],
        "chain_column_id_prefixes": args.chain_column_id_prefixes,
        "replacement_count": replacement_count,
        "beam_size": int(args.beam_size),
        "blocks_scanned": blocks_scanned,
        "checked_extensions": checked_extensions,
        "final_states": final_checked,
        "elapsed_feasible": elapsed_feasible,
        "candidates": [],
    }

    out_route = Path(args.route_out) if args.route_out else None
    best_route_validation = None
    for rank, rec in enumerate(proxy_candidates, start=1):
        neg_covered, proxy_elapsed_s, proxy_delta_s, removed_count, chain, start, end, hits, missing = rec
        order = selected[:start] + list(chain) + selected[end:]
        path, failed = _realize_exact_column_path(
            G, tables, rows, order, validation_runs, args.connector_cap)
        validation = None
        route_path = None
        replacement_columns = [rows[i]["column_id"] for i in chain]
        if path is not None:
            validation = validate(prob, path)
            print(f"candidate {rank}: {result_line(validation)} "
                  f"replacement={replacement_columns}", flush=True)
            if (
                out_route is not None
                and _replacement_route_is_better(
                    validation,
                    best_route_validation,
                    min_covered_count,
                    max_total_elapsed_s,
                )
            ):
                out_route.parent.mkdir(parents=True, exist_ok=True)
                meta = {
                    "radius_m": args.validation_radius,
                    "elapsed_s": int(_elapsed(path[0][1], path[-1][1])),
                    "notes": "Beam-screened multi-column block replacement",
                    "source_order": args.source_order,
                    "columns_source": args.columns_source,
                    "replacement_columns": replacement_columns,
                    "replaced_after": rows[selected[start - 1]]["column_id"],
                    "replaced_before": rows[selected[end]]["column_id"] if end < n else None,
                    "removed_columns": [rows[idx]["column_id"] for idx in selected[start:end]],
                    "proxy_covered_count": -neg_covered,
                    "proxy_missing_stations": missing,
                    "proxy_elapsed_s": proxy_elapsed_s,
                    "proxy_elapsed": hms(proxy_elapsed_s),
                    "selected_order": [rows[idx]["column_id"] for idx in order],
                }
                out_route.write_text(json.dumps({"meta": meta, "path": path}, sort_keys=True))
                route_path = str(out_route)
                best_route_validation = validation
                print(f"wrote {out_route}", flush=True)
                if args.stop_after_first:
                    summary["candidates"].append({
                        "rank": rank,
                        "replacement_columns": replacement_columns,
                        "proxy_covered_count": -neg_covered,
                        "proxy_elapsed_s": proxy_elapsed_s,
                        "proxy_delta_s": proxy_delta_s,
                        "removed_count": removed_count,
                        "target_hits": hits,
                        "failed_arc": failed,
                        "validation": validation,
                        "route_out": route_path,
                    })
                    break
        else:
            print(f"candidate {rank}: failed replacement={replacement_columns} "
                  f"failed_arc={failed}", flush=True)
        summary["candidates"].append({
            "rank": rank,
            "replacement_columns": replacement_columns,
            "replaced_after": rows[selected[start - 1]]["column_id"],
            "replaced_before": rows[selected[end]]["column_id"] if end < n else None,
            "removed_columns": [rows[idx]["column_id"] for idx in selected[start:end]],
            "proxy_covered_count": -neg_covered,
            "proxy_elapsed_s": proxy_elapsed_s,
            "proxy_elapsed": hms(proxy_elapsed_s),
            "proxy_delta_s": proxy_delta_s,
            "removed_count": removed_count,
            "target_hits": hits,
            "proxy_missing_stations": missing,
            "failed_arc": list(failed) if failed is not None else None,
            "validation": validation,
            "route_out": route_path,
        })

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(summary, sort_keys=True))
        print(f"wrote {out}", flush=True)
    return 0


def cmd_price_terminals(args) -> int:
    import pickle

    si = StationIndex.load()
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, _node_stop, node_t = tables
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)

    seed_path, _meta = _load_route(Path(args.seed))
    seed_nodes = _nodes_for_path(prob, seed_path)
    source_nodes = seed_nodes[::max(1, args.source_stride)]
    if args.max_sources:
        source_nodes = source_nodes[:args.max_sources]

    if args.targets == "all":
        targets = sorted(si.canonical_stations)
    else:
        targets = sorted({si.resolve(stop) for stop in _terminal_stops()
                          if si.resolve(stop) in si.canonical_stations})
    if args.max_targets:
        targets = targets[:args.max_targets]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    generated = []
    seen = set()
    t0 = time.time()
    max_elapsed_s = int(args.max_elapsed_minutes * 60)

    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")
        for sidx, src in enumerate(source_nodes, start=1):
            for target in targets:
                if node_off[src] == target:
                    continue
                tgt, prev = _dijkstra_to_station(
                    G._adj, node_off, src, target, args.connector_cap, runs=runs)
                if tgt is None:
                    continue
                elapsed = _elapsed(node_t[src], node_t[tgt])
                if elapsed <= 0 or elapsed > max_elapsed_s:
                    continue
                seg = _path_from_prev(prev, src, tgt)
                if seg is None:
                    continue
                nodes = [src] + seg
                covered = {node_off[n] for n in nodes} & si.canonical_stations
                if len(covered) < args.min_covered:
                    continue
                key = (node_t[src], node_off[src], node_t[tgt], node_off[tgt])
                if key in seen:
                    continue
                seen.add(key)
                record = _column_from_nodes(
                    prob, args.label, len(base_rows) + len(generated) + 1, nodes)
                generated.append(record)
                f.write(json.dumps(record, sort_keys=True) + "\n")
            if sidx % 10 == 0 or sidx == len(source_nodes):
                print(f"priced sources {sidx}/{len(source_nodes)} "
                      f"generated={len(generated)} wall={time.time() - t0:.1f}s",
                      flush=True)

    print(f"wrote {len(base_rows)} base + {len(generated)} priced columns -> {out}")
    return 0


def _parse_ints(value):
    return [int(x.strip()) for x in str(value).split(",") if x.strip()]


def _linear_status_name(status):
    names = {
        pywraplp.Solver.OPTIMAL: "OPTIMAL",
        pywraplp.Solver.FEASIBLE: "FEASIBLE",
        pywraplp.Solver.INFEASIBLE: "INFEASIBLE",
        pywraplp.Solver.UNBOUNDED: "UNBOUNDED",
        pywraplp.Solver.ABNORMAL: "ABNORMAL",
        pywraplp.Solver.NOT_SOLVED: "NOT_SOLVED",
    }
    return names.get(status, str(status))


def cmd_dual_cover(args) -> int:
    rows = _load_column_rows(Path(args.columns_file))
    if not rows:
        raise SystemExit("no columns to solve")

    canonical = set(StationIndex.load().canonical_stations)
    by_station = {station: [] for station in canonical}
    row_coverage = []
    for i, row in enumerate(rows):
        covered = sorted({station for station in row.get("covered_stations", [])
                          if station in canonical})
        row_coverage.append(covered)
        for station in covered:
            by_station[station].append(i)
    missing = sorted(st for st, cols in by_station.items() if not cols)
    if missing:
        raise SystemExit(f"column pool misses {len(missing)} station(s), e.g. {missing[:8]}")

    solver = pywraplp.Solver.CreateSolver(args.solver)
    if solver is None:
        raise SystemExit(f"linear solver is not available: {args.solver}")
    if args.time_limit > 0:
        solver.SetTimeLimit(int(args.time_limit * 1000))

    upper = 1.0 if args.unit_upper_bound else solver.infinity()
    x = [solver.NumVar(0.0, upper, f"col_{i}") for i in range(len(rows))]
    constraints = {}
    for station, cols in by_station.items():
        ct = solver.Constraint(1.0, solver.infinity(), f"cover_{station}")
        for i in cols:
            ct.SetCoefficient(x[i], 1.0)
        constraints[station] = ct

    objective = solver.Objective()
    for i, row in enumerate(rows):
        objective.SetCoefficient(x[i], float(int(row["elapsed_s"]) + args.column_penalty))
    objective.SetMinimization()

    status = solver.Solve()
    status_name = _linear_status_name(status)
    if status not in (pywraplp.Solver.OPTIMAL, pywraplp.Solver.FEASIBLE):
        raise SystemExit(f"dual-cover failed: status={status_name}")

    station_duals = {}
    for station, ct in constraints.items():
        dual = float(ct.dual_value())
        station_duals[station] = {
            "dual_s": dual,
            "reward_s": max(0.0, dual),
            "covering_columns": len(by_station[station]),
        }

    positive = [(st, info["reward_s"], info["covering_columns"])
                for st, info in station_duals.items() if info["reward_s"] > 1e-7]
    positive.sort(key=lambda item: (-item[1], item[0]))

    selected = []
    for i, var in enumerate(x):
        value = float(var.solution_value())
        if value < args.min_x:
            continue
        reward = sum(station_duals[st]["reward_s"] for st in row_coverage[i])
        cost = int(rows[i]["elapsed_s"]) + args.column_penalty
        selected.append({
            "row": i,
            "column_id": rows[i].get("column_id"),
            "x": value,
            "cost_s": cost,
            "dual_reward_s": reward,
            "station_reduced_cost_s": cost - reward,
            "covered_count": len(row_coverage[i]),
            "start": rows[i].get("start"),
            "end": rows[i].get("end"),
        })
    selected.sort(key=lambda item: (-item["x"], item["station_reduced_cost_s"]))

    result = {
        "source": str(args.columns_file),
        "solver": args.solver,
        "status": status_name,
        "objective_s": float(objective.Value()),
        "objective": hms(int(round(objective.Value()))),
        "columns": len(rows),
        "stations": len(canonical),
        "positive_dual_count": len(positive),
        "positive_dual_sum_s": sum(reward for _st, reward, _count in positive),
        "column_penalty_s": args.column_penalty,
        "unit_upper_bound": bool(args.unit_upper_bound),
        "station_duals": station_duals,
        "top_duals": [
            {"station": st, "reward_s": reward, "covering_columns": count}
            for st, reward, count in positive[:args.top_duals]
        ],
        "active_columns": selected[:args.max_active_columns],
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, sort_keys=True))

    print(f"status={status_name} objective={result['objective']} "
          f"positive_duals={len(positive)}/{len(canonical)}")
    print("top dual stations:")
    for item in result["top_duals"][:12]:
        print(f"  {item['station']}: {hms(int(round(item['reward_s'])))} "
              f"columns={item['covering_columns']}")
    print(f"active_columns>={args.min_x:g}: {len(selected)}")
    print(f"wrote {out}")
    return 0


def _load_dual_rewards(path: Path):
    data = json.loads(path.read_text())
    rewards = {}
    for station, info in data.get("station_duals", {}).items():
        reward = float(info.get("reward_s", max(0.0, float(info.get("dual_s", 0.0)))))
        if reward > 0:
            rewards[station] = reward
    return data, rewards


def _row_dual_reward(row, rewards):
    return sum(rewards.get(station, 0.0) for station in row.get("covered_stations", []))


def _row_pricing_reduced_cost(row):
    values = []
    for key in ("pricing_dual", "pricing_dual_beam"):
        if key in row and "reduced_cost_s" in row[key]:
            values.append(float(row[key]["reduced_cost_s"]))
    if "pricing_connector" in row and "score_s" in row["pricing_connector"]:
        values.append(float(row["pricing_connector"]["score_s"]))
    if "pricing_window_event" in row and "score_s" in row["pricing_window_event"]:
        values.append(float(row["pricing_window_event"]["score_s"]))
    if "pricing_phase_variant" in row and "score_s" in row["pricing_phase_variant"]:
        values.append(float(row["pricing_phase_variant"]["score_s"]))
    if "pricing_phase_window" in row and "score_s" in row["pricing_phase_window"]:
        values.append(float(row["pricing_phase_window"]["score_s"]))
    if "pricing_cluster_corridor" in row and "score_s" in row["pricing_cluster_corridor"]:
        values.append(float(row["pricing_cluster_corridor"]["score_s"]))
    if "pricing_late_tail_beam" in row and "score_s" in row["pricing_late_tail_beam"]:
        values.append(float(row["pricing_late_tail_beam"]["score_s"]))
    if "pricing_resource_chain" in row and "score_s" in row["pricing_resource_chain"]:
        values.append(float(row["pricing_resource_chain"]["score_s"]))
    if "pricing_stage_chain" in row and "score_s" in row["pricing_stage_chain"]:
        values.append(float(row["pricing_stage_chain"]["score_s"]))
    return min(values) if values else None


def _active_pool_score(row, rewards, column_penalty):
    priced = _row_pricing_reduced_cost(row)
    if priced is not None:
        return priced
    return int(row["elapsed_s"]) + column_penalty - _row_dual_reward(row, rewards)


def cmd_active_pool(args) -> int:
    rows = []
    by_id = {}
    duplicates = 0
    for source in args.columns_files:
        for row in _load_column_rows(Path(source)):
            column_id = row.get("column_id")
            if not column_id:
                raise SystemExit(f"column without column_id in {source}")
            if column_id in by_id:
                duplicates += 1
                continue
            rows.append(row)
            by_id[column_id] = row
    if not rows:
        raise SystemExit("no columns to select")

    _dual_data, rewards = _load_dual_rewards(Path(args.duals)) if args.duals else ({}, {})
    canonical = set(StationIndex.load().canonical_stations)
    selected = set()
    reasons = defaultdict(set)
    missing_solution_ids = []

    for solution_path in args.include_solution or []:
        solution_ids = _selected_column_ids_from_solution(Path(solution_path))
        label = f"solution:{Path(solution_path).name}"
        for column_id in solution_ids:
            if column_id in by_id:
                selected.add(column_id)
                reasons[column_id].add(label)
            else:
                missing_solution_ids.append((column_id, solution_path))

    scored = sorted(
        rows,
        key=lambda row: (
            _active_pool_score(row, rewards, args.column_penalty),
            int(row["elapsed_s"]),
            -int(row.get("covered_count", len(row.get("covered_stations", [])))),
            row["column_id"],
        ),
    )

    if args.top_priced:
        priced_rows = [row for row in rows if _row_pricing_reduced_cost(row) is not None]
        priced_rows.sort(key=lambda row: (
            _row_pricing_reduced_cost(row),
            int(row["elapsed_s"]),
            row["column_id"],
        ))
        for row in priced_rows[:args.top_priced]:
            selected.add(row["column_id"])
            reasons[row["column_id"]].add("priced")

    for row in scored:
        if len(selected) >= args.top_columns:
            break
        selected.add(row["column_id"])
        reasons[row["column_id"]].add("score")

    def selected_coverage():
        covered = set()
        for column_id in selected:
            covered.update(set(by_id[column_id].get("covered_stations", [])) & canonical)
        return covered

    if not args.no_ensure_coverage:
        covered = selected_coverage()
        while len(covered) < len(canonical):
            missing = canonical - covered
            best = None
            for row in rows:
                column_id = row["column_id"]
                if column_id in selected:
                    continue
                newly = len(set(row.get("covered_stations", [])) & missing)
                if newly <= 0:
                    continue
                score = _active_pool_score(row, rewards, args.column_penalty)
                candidate = (score / newly, score, -newly, int(row["elapsed_s"]), column_id)
                if best is None or candidate < best[0]:
                    best = (candidate, row)
            if best is None:
                break
            row = best[1]
            selected.add(row["column_id"])
            reasons[row["column_id"]].add("coverage")
            covered = selected_coverage()

    selected_rows = [row for row in rows if row["column_id"] in selected]
    coverage = set()
    modes = Counter()
    reason_counts = Counter()
    for row in selected_rows:
        coverage.update(set(row.get("covered_stations", [])) & canonical)
        modes.update(row.get("modes", {}))
        for reason in reasons[row["column_id"]]:
            reason_counts[reason] += 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        for row in selected_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")

    print(f"loaded={len(rows)} duplicates_skipped={duplicates}")
    print(f"selected={len(selected_rows)} covered={len(coverage)}/{len(canonical)} "
          f"reasons={dict(reason_counts)} modes={dict(modes)}")
    if missing_solution_ids:
        print(f"missing selected ids from solution files: {len(missing_solution_ids)} "
              f"e.g. {missing_solution_ids[:6]}")
    if len(coverage) < len(canonical):
        missing = sorted(canonical - coverage)
        raise SystemExit(f"active pool misses {len(missing)} station(s), e.g. {missing[:8]}")
    print(f"wrote {out}")
    return 0


def _parse_minutes_to_seconds(value):
    return [int(round(float(x.strip()) * 60))
            for x in str(value).split(",") if x.strip()]


def _parse_duration_seconds(value):
    if value is None:
        return None
    raw = str(value).strip()
    if not raw:
        return None
    if ":" not in raw:
        return int(round(float(raw)))
    parts = raw.split(":")
    if len(parts) > 3:
        raise ValueError(f"duration must be seconds, MM:SS, or HH:MM:SS, got {value!r}")
    try:
        values = [int(part) for part in parts]
    except ValueError as exc:
        raise ValueError(f"duration contains a non-integer field: {value!r}") from exc
    if any(part < 0 for part in values):
        raise ValueError(f"duration fields must be non-negative: {value!r}")
    total = 0
    for part in values:
        total = total * 60 + part
    return total


def _parse_week_time_window(value):
    if value is None:
        return None
    raw = str(value).strip()
    if not raw:
        return None
    sep = ".." if ".." in raw else "-"
    if sep not in raw:
        raise ValueError(f"time window must be START{sep}END, got {value!r}")
    left, right = [part.strip() for part in raw.split(sep, 1)]
    if not left or not right:
        raise ValueError(f"time window must have both endpoints, got {value!r}")
    return int(left), int(right)


def _rows_in_time_window(rows, endpoint, window):
    if window is None:
        return None
    lo, hi = window
    key = "start" if endpoint == "start" else "end"
    if lo <= hi:
        return [i for i, row in enumerate(rows) if lo <= int(row[key]["time"]) <= hi]
    return [i for i, row in enumerate(rows)
            if int(row[key]["time"]) >= lo or int(row[key]["time"]) <= hi]


def _time_in_window(t, window):
    if window is None:
        return True
    lo, hi = window
    t = int(t)
    if lo <= hi:
        return lo <= t <= hi
    return t >= lo or t <= hi


def _row_anchor_sequence(row, si, mode):
    if "path" not in row:
        return []
    if mode not in {"first", "transitions", "revisit"}:
        raise ValueError(f"unknown anchor mode: {mode!r}")
    anchors = []
    seen = set()
    last = None
    for stop, _t, *_rest in row["path"]:
        station = si.resolve(stop)
        if station == last:
            continue
        last = station
        if mode == "transitions":
            anchors.append(station)
        elif station not in seen:
            seen.add(station)
            anchors.append(station)
        elif mode == "revisit":
            anchors.append(station)
    return anchors


def _station_platforms(node_off, node_stop):
    by_station = defaultdict(list)
    seen = defaultdict(set)
    for station, stop in zip(node_off, node_stop):
        if not station or not stop or stop in seen[station]:
            continue
        seen[station].add(stop)
        by_station[station].append(stop)
    return {station: sorted(stops) for station, stops in by_station.items()}


def _phase_source_rows(rows, rewards, max_rows, column_penalty):
    def key(row):
        covered_count = int(row.get("covered_count", len(row.get("covered_stations", []))))
        if rewards:
            return (_active_pool_score(row, rewards, column_penalty),
                    -covered_count,
                    int(row["elapsed_s"]),
                    row.get("column_id", ""))
        return (-covered_count, int(row["elapsed_s"]), row.get("column_id", ""))

    candidates = [row for row in rows if "path" in row and int(row.get("elapsed_s", 0)) > 0]
    candidates.sort(key=key)
    return candidates[:max_rows] if max_rows else candidates


def cmd_price_phase_variants(args) -> int:
    import pickle

    source_rows = _load_column_rows(Path(args.source_columns))
    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    _dual_data, rewards = _load_dual_rewards(Path(args.duals)) if args.duals else ({}, {})

    si = StationIndex.load()
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    stop_events = _stop_index(node_stop, node_t)
    station_platforms = _station_platforms(node_off, node_stop)
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)

    offsets = _parse_minutes_to_seconds(args.offset_minutes)
    if not offsets:
        raise SystemExit("--offset-minutes must include at least one value")
    selected_sources = _phase_source_rows(
        source_rows, rewards, args.max_source_rows, args.column_penalty)
    if not selected_sources:
        raise SystemExit("no source columns with path slices selected")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    generated = []
    seen = set()
    for row in base_rows:
        seen.add((
            row["start"]["time"],
            row["start"]["station"],
            row["end"]["time"],
            row["end"]["station"],
            tuple(row.get("covered_stations", [])),
        ))

    max_extra_s = int(args.max_extra_minutes * 60)
    max_elapsed_s = int(args.max_elapsed_minutes * 60)
    max_score = INF if args.max_score_s is None else float(args.max_score_s)
    t0 = time.time()

    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")

        for ridx, row in enumerate(selected_sources, start=1):
            anchors = _row_anchor_sequence(row, si, args.anchor_mode)
            if len(anchors) < args.min_anchors:
                continue
            original_covered = {station for station in row.get("covered_stations", [])
                                if station in si.canonical_stations}
            if len(original_covered) < args.min_covered:
                continue

            start_station = row["start"]["station"]
            target_anchors = anchors[1:] if anchors and anchors[0] == start_station else anchors
            if not target_anchors:
                continue

            if args.source_stop_mode == "same-platform":
                source_stops = [row["start"]["stop"]]
            else:
                source_stops = station_platforms.get(start_station, [row["start"]["stop"]])

            for source_stop in source_stops:
                for offset in offsets:
                    target_t = int(row["start"]["time"]) + offset
                    src = _first_event_at_or_after(stop_events, source_stop, target_t)
                    if src is None:
                        continue
                    if (not args.include_original
                            and source_stop == row["start"]["stop"]
                            and int(node_t[src]) == int(row["start"]["time"])):
                        continue
                    if node_off[src] != start_station:
                        continue

                    suffix, end, _visited = realize_from(
                        G._adj,
                        tables,
                        src,
                        {start_station},
                        target_anchors,
                        cap=args.connector_cap,
                        runs=runs,
                        skip_visited=args.skip_visited_anchors,
                    )
                    if suffix is None:
                        continue
                    nodes = [src] + suffix
                    elapsed = sum(weight for weight, _mode in _transition_series(prob, nodes))
                    if elapsed <= 0 or elapsed > max_elapsed_s:
                        continue
                    if elapsed > int(row["elapsed_s"]) + max_extra_s:
                        continue
                    covered = {node_off[n] for n in nodes} & si.canonical_stations
                    if len(covered) < args.min_covered:
                        continue
                    overlap = len(covered & original_covered) / max(1, len(original_covered))
                    if overlap < args.min_original_overlap:
                        continue
                    reward = (_row_dual_reward({"covered_stations": covered}, rewards)
                              if rewards else args.cover_reward_s * len(covered))
                    score = elapsed + args.column_penalty - reward
                    if score > max_score:
                        continue
                    key = (node_t[src], node_off[src], node_t[end], node_off[end],
                           tuple(sorted(covered)), tuple(target_anchors))
                    if key in seen:
                        continue
                    seen.add(key)
                    record = _column_from_nodes(
                        prob, args.label, len(base_rows) + len(generated) + 1, nodes)
                    record["pricing_phase_variant"] = {
                        "source_column": row.get("column_id"),
                        "source_start_time": int(row["start"]["time"]),
                        "source_end_time": int(row["end"]["time"]),
                        "source_elapsed_s": int(row["elapsed_s"]),
                        "source_stop": row["start"]["stop"],
                        "variant_stop": source_stop,
                        "offset_s": int(offset),
                        "offset_minutes": round(offset / 60, 3),
                        "anchor_mode": args.anchor_mode,
                        "anchor_count": len(anchors),
                        "target_anchor_count": len(target_anchors),
                        "original_covered_count": len(original_covered),
                        "overlap": overlap,
                        "reward_s": reward,
                        "score_s": score,
                        "elapsed_delta_s": int(elapsed) - int(row["elapsed_s"]),
                        "run_mode": run_mode,
                    }
                    generated.append(record)
                    f.write(json.dumps(record, sort_keys=True) + "\n")
                    if args.max_generated and len(generated) >= args.max_generated:
                        break
                if args.max_generated and len(generated) >= args.max_generated:
                    break
            if ridx % 10 == 0 or ridx == len(selected_sources):
                print(f"phase variants {ridx}/{len(selected_sources)} "
                      f"generated={len(generated)} wall={time.time() - t0:.1f}s",
                      flush=True)
            if args.max_generated and len(generated) >= args.max_generated:
                break

    print(f"phase source rows={len(selected_sources)} offsets={len(offsets)}")
    print(f"wrote {len(base_rows)} base + {len(generated)} phase-variant columns -> {out}")
    return 0


def _combined_path(rows):
    path = []
    for row in rows:
        row_path = row.get("path", [])
        if not row_path:
            continue
        if path and row_path and path[-1][:2] == row_path[0][:2]:
            path.extend(row_path[1:])
        else:
            path.extend(row_path)
    return path


def _anchor_sequence_from_path(path, si, mode):
    if mode not in {"first", "transitions", "revisit"}:
        raise ValueError(f"unknown anchor mode: {mode!r}")
    anchors = []
    seen = set()
    last = None
    for stop, _t, *_rest in path:
        station = si.resolve(stop)
        if station == last:
            continue
        last = station
        if mode == "transitions":
            anchors.append(station)
        elif station not in seen:
            seen.add(station)
            anchors.append(station)
        elif mode == "revisit":
            anchors.append(station)
    return anchors


def cmd_price_phase_windows(args) -> int:
    import pickle

    solution_rows = _load_column_rows(Path(args.source_solution))
    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    _dual_data, rewards = _load_dual_rewards(Path(args.duals)) if args.duals else ({}, {})
    widths = _parse_ints(args.widths)
    if not widths:
        raise SystemExit("--widths must include at least one positive integer")

    si = StationIndex.load()
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    stop_events = _stop_index(node_stop, node_t)
    station_platforms = _station_platforms(node_off, node_stop)
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)
    offsets = _parse_minutes_to_seconds(args.offset_minutes)
    if not offsets:
        raise SystemExit("--offset-minutes must include at least one value")
    variant_start_window = _parse_week_time_window(args.variant_start_time_window)
    variant_end_window = _parse_week_time_window(args.variant_end_time_window)

    jobs = []
    for start in range(0, len(solution_rows), max(1, args.window_stride)):
        for width in widths:
            end = start + width - 1
            if end >= len(solution_rows):
                continue
            window_rows = solution_rows[start:end + 1]
            if any("path" not in row for row in window_rows):
                continue
            path = _combined_path(window_rows)
            anchors = _anchor_sequence_from_path(path, si, args.anchor_mode)
            if len(anchors) < args.min_anchors:
                continue
            start_row, end_row = window_rows[0], window_rows[-1]
            original_elapsed = _phase_proxy_gap(
                start_row["start"]["time"], end_row["end"]["time"], 0)
            if original_elapsed <= 0 or original_elapsed > int(args.max_original_minutes * 60):
                continue
            original_covered = set()
            for row in window_rows:
                original_covered.update(
                    station for station in row.get("covered_stations", [])
                    if station in si.canonical_stations)
            if len(original_covered) < args.min_covered:
                continue
            reward = (_row_dual_reward({"covered_stations": original_covered}, rewards)
                      if rewards else args.cover_reward_s * len(original_covered))
            score = original_elapsed - reward
            jobs.append({
                "score": score,
                "start": start,
                "end": end,
                "width": width,
                "rows": window_rows,
                "anchors": anchors,
                "original_elapsed_s": original_elapsed,
                "original_covered": original_covered,
                "original_reward_s": reward,
                "window_columns": [row["column_id"] for row in window_rows],
            })

    jobs.sort(key=lambda job: (-job["score"], -job["original_elapsed_s"], job["start"]))
    if args.max_windows:
        jobs = jobs[:args.max_windows]
    if not jobs:
        raise SystemExit("no phase-window jobs selected")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    generated = []
    seen = set()
    for row in base_rows:
        seen.add((
            row["start"]["time"],
            row["start"]["station"],
            row["end"]["time"],
            row["end"]["station"],
            tuple(row.get("covered_stations", [])),
        ))

    max_extra_s = int(args.max_extra_minutes * 60)
    max_elapsed_s = int(args.max_elapsed_minutes * 60)
    max_score = INF if args.max_score_s is None else float(args.max_score_s)
    t0 = time.time()

    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")

        for jidx, job in enumerate(jobs, start=1):
            first = job["rows"][0]
            start_station = first["start"]["station"]
            anchors = job["anchors"]
            target_anchors = anchors[1:] if anchors and anchors[0] == start_station else anchors
            if not target_anchors:
                continue
            if args.source_stop_mode == "same-platform":
                source_stops = [first["start"]["stop"]]
            else:
                source_stops = station_platforms.get(start_station, [first["start"]["stop"]])

            for source_stop in source_stops:
                for offset in offsets:
                    target_t = int(first["start"]["time"]) + offset
                    src = _first_event_at_or_after(stop_events, source_stop, target_t)
                    if src is None:
                        continue
                    if not _time_in_window(node_t[src], variant_start_window):
                        continue
                    if (not args.include_original
                            and source_stop == first["start"]["stop"]
                            and int(node_t[src]) == int(first["start"]["time"])):
                        continue
                    if node_off[src] != start_station:
                        continue
                    suffix, end, _visited = realize_from(
                        G._adj,
                        tables,
                        src,
                        {start_station},
                        target_anchors,
                        cap=args.connector_cap,
                        runs=runs,
                        skip_visited=args.skip_visited_anchors,
                    )
                    if suffix is None:
                        continue
                    nodes = [src] + suffix
                    elapsed = sum(weight for weight, _mode in _transition_series(prob, nodes))
                    if elapsed <= 0 or elapsed > max_elapsed_s:
                        continue
                    if not _time_in_window(node_t[end], variant_end_window):
                        continue
                    if elapsed > job["original_elapsed_s"] + max_extra_s:
                        continue
                    covered = {node_off[n] for n in nodes} & si.canonical_stations
                    if len(covered) < args.min_covered:
                        continue
                    overlap = len(covered & job["original_covered"]) / max(
                        1, len(job["original_covered"]))
                    if overlap < args.min_original_overlap:
                        continue
                    reward = (_row_dual_reward({"covered_stations": covered}, rewards)
                              if rewards else args.cover_reward_s * len(covered))
                    score = elapsed + args.column_penalty - reward
                    if args.window_credit:
                        score -= args.window_credit * job["original_elapsed_s"]
                    if max_score < INF and score > max_score:
                        continue
                    key = (node_t[src], node_off[src], node_t[end], node_off[end],
                           tuple(sorted(covered)), tuple(target_anchors))
                    if key in seen:
                        continue
                    seen.add(key)
                    record = _column_from_nodes(
                        prob, args.label, len(base_rows) + len(generated) + 1, nodes)
                    record["pricing_phase_window"] = {
                        "source_solution": args.source_solution,
                        "start_index": job["start"],
                        "end_index": job["end"],
                        "width": job["width"],
                        "window_columns": job["window_columns"],
                        "source_start_time": int(first["start"]["time"]),
                        "source_end_time": int(job["rows"][-1]["end"]["time"]),
                        "source_elapsed_s": int(job["original_elapsed_s"]),
                        "source_stop": first["start"]["stop"],
                        "variant_stop": source_stop,
                        "offset_s": int(offset),
                        "offset_minutes": round(offset / 60, 3),
                        "anchor_mode": args.anchor_mode,
                        "anchor_count": len(anchors),
                        "target_anchor_count": len(target_anchors),
                        "original_covered_count": len(job["original_covered"]),
                        "overlap": overlap,
                        "reward_s": reward,
                        "score_s": score,
                        "elapsed_delta_s": int(elapsed) - int(job["original_elapsed_s"]),
                        "run_mode": run_mode,
                    }
                    generated.append(record)
                    f.write(json.dumps(record, sort_keys=True) + "\n")
                    if args.max_generated and len(generated) >= args.max_generated:
                        break
                if args.max_generated and len(generated) >= args.max_generated:
                    break
            if jidx % 10 == 0 or jidx == len(jobs):
                print(f"phase windows {jidx}/{len(jobs)} generated={len(generated)} "
                      f"wall={time.time() - t0:.1f}s", flush=True)
            if args.max_generated and len(generated) >= args.max_generated:
                break

    print(f"phase-window jobs={len(jobs)} offsets={len(offsets)} "
          f"variant_start_window={variant_start_window} "
          f"variant_end_window={variant_end_window}")
    print(f"wrote {len(base_rows)} base + {len(generated)} phase-window columns -> {out}")
    return 0


def cmd_price_duals(args) -> int:
    import pickle

    _dual_data, rewards = _load_dual_rewards(Path(args.duals))
    if not rewards:
        raise SystemExit("dual file has no positive station rewards")

    si = StationIndex.load()
    candidate_targets = [
        (station, reward)
        for station, reward in rewards.items()
        if station in si.canonical_stations and reward >= args.min_dual_s
    ]
    candidate_targets.sort(key=lambda item: (-item[1], item[0]))
    if args.max_targets:
        candidate_targets = candidate_targets[:args.max_targets]
    targets = [station for station, _reward in candidate_targets]
    if not targets:
        raise SystemExit("no dual targets selected")

    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, _node_stop, node_t = tables
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)

    seed_specs = args.seed or ["solutions/best.json"]
    source_events = []
    for seed_spec in seed_specs:
        seed_path, _meta = _load_route(Path(seed_spec))
        seed_nodes = _nodes_for_path(prob, seed_path)
        for local_idx, node in enumerate(seed_nodes[::max(1, args.source_stride)]):
            source_events.append((str(seed_spec), local_idx, node))
    if args.max_sources:
        source_events = source_events[:args.max_sources]
    if not source_events:
        raise SystemExit("no seed source events selected")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    generated = []
    seen = set()
    t0 = time.time()
    max_elapsed_s = int(args.max_elapsed_minutes * 60)

    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")
        for sidx, (seed_spec, local_idx, src) in enumerate(source_events, start=1):
            for target in targets:
                if node_off[src] == target:
                    continue
                tgt, prev = _dijkstra_to_station(
                    G._adj, node_off, src, target, args.connector_cap, runs=runs)
                if tgt is None:
                    continue
                elapsed = _elapsed(node_t[src], node_t[tgt])
                if elapsed <= 0 or elapsed > max_elapsed_s:
                    continue
                seg = _path_from_prev(prev, src, tgt)
                if seg is None:
                    continue
                nodes = [src] + seg
                covered = {node_off[n] for n in nodes} & si.canonical_stations
                if len(covered) < args.min_covered:
                    continue
                reward = sum(rewards.get(station, 0.0) for station in covered)
                reduced_cost = elapsed + args.column_penalty - reward
                if reward < args.min_reward_s:
                    continue
                if reduced_cost > args.max_reduced_cost_s:
                    continue
                key = (node_t[src], node_off[src], node_t[tgt], node_off[tgt],
                       tuple(sorted(covered)))
                if key in seen:
                    continue
                seen.add(key)
                record = _column_from_nodes(
                    prob, args.label, len(base_rows) + len(generated) + 1, nodes)
                record["pricing_dual"] = {
                    "dual_file": args.duals,
                    "seed": seed_spec,
                    "source_index": local_idx,
                    "source_station": node_off[src],
                    "source_time": int(node_t[src]),
                    "target_station": target,
                    "target_reward_s": rewards.get(target, 0.0),
                    "path_reward_s": reward,
                    "reduced_cost_s": reduced_cost,
                    "run_mode": run_mode,
                }
                generated.append(record)
                f.write(json.dumps(record, sort_keys=True) + "\n")

            if sidx % 10 == 0 or sidx == len(source_events):
                print(f"priced dual sources {sidx}/{len(source_events)} "
                      f"generated={len(generated)} wall={time.time() - t0:.1f}s",
                      flush=True)

    print(f"dual targets={len(targets)} min_dual={args.min_dual_s:g}s")
    print(f"wrote {len(base_rows)} base + {len(generated)} dual-priced columns -> {out}")
    return 0


def _label_path_nodes(labels, label_id):
    nodes = []
    while label_id >= 0:
        node, parent, _elapsed_s, _reward_s, _mask = labels[label_id]
        nodes.append(node)
        label_id = parent
    return list(reversed(nodes))


def _beam_label_dominated(labels, label_ids, elapsed_s, reward_s, mask):
    for other_id in label_ids:
        _node, _parent, other_elapsed, other_reward, other_mask = labels[other_id]
        if (other_elapsed <= elapsed_s
                and other_reward >= reward_s
                and (other_mask | mask) == other_mask):
            return True
    return False


def cmd_price_dual_beam(args) -> int:
    import pickle

    _dual_data, all_rewards = _load_dual_rewards(Path(args.duals))
    if not all_rewards:
        raise SystemExit("dual file has no positive station rewards")

    si = StationIndex.load()
    target_rewards = [
        (station, reward)
        for station, reward in all_rewards.items()
        if station in si.canonical_stations and reward >= args.min_dual_s
    ]
    target_rewards.sort(key=lambda item: (-item[1], item[0]))
    if args.top_reward_stations:
        target_rewards = target_rewards[:args.top_reward_stations]
    if not target_rewards:
        raise SystemExit("no reward stations selected")
    search_rewards = {station: reward for station, reward in target_rewards}
    reward_bit = {station: i for i, (station, _reward) in enumerate(target_rewards)}

    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, _node_stop, node_t = tables
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)

    seed_specs = args.seed or ["solutions/best.json"]
    source_events = []
    for seed_spec in seed_specs:
        seed_path, _meta = _load_route(Path(seed_spec))
        seed_nodes = _nodes_for_path(prob, seed_path)
        for local_idx, node in enumerate(seed_nodes[::max(1, args.source_stride)]):
            source_events.append((str(seed_spec), local_idx, node))
    if args.max_sources:
        source_events = source_events[:args.max_sources]
    if not source_events:
        raise SystemExit("no seed source events selected")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    generated = []
    seen = set()
    max_elapsed_s = int(args.max_elapsed_minutes * 60)
    t0 = time.time()

    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")

        for sidx, (seed_spec, local_idx, src) in enumerate(source_events, start=1):
            start_station = node_off[src]
            start_mask = 0
            start_reward = 0.0
            bit = reward_bit.get(start_station)
            if bit is not None:
                start_mask |= 1 << bit
                start_reward += search_rewards[start_station]

            labels = [(src, -1, 0, start_reward, start_mask)]
            labels_by_node = defaultdict(list)
            labels_by_node[src].append(0)
            heap = [(args.column_penalty - start_reward, 0, 0)]
            candidate_ids = []
            expansions = 0

            while heap and expansions < args.max_expansions_per_source:
                _priority, _elapsed_key, label_id = heapq.heappop(heap)
                node, parent, elapsed_s, reward_s, mask = labels[label_id]
                if label_id not in labels_by_node.get(node, ()):
                    continue
                expansions += 1

                reduced = elapsed_s + args.column_penalty - reward_s
                if (node != src
                        and reduced <= args.max_reduced_cost_s
                        and reward_s >= args.min_reward_s
                        and mask.bit_count() >= args.min_reward_stations):
                    candidate_ids.append(label_id)
                    if (args.candidate_pool_per_source
                            and len(candidate_ids) > args.candidate_pool_per_source):
                        candidate_ids.sort(
                            key=lambda lid: (
                                labels[lid][2] + args.column_penalty - labels[lid][3],
                                -labels[lid][3],
                                labels[lid][2],
                            ))
                        del candidate_ids[args.candidate_pool_per_source:]

                successors = [(v, ed["weight"]) for v, ed in G._adj[node].items()]
                if runs is not None:
                    successors.extend((v, w) for v, w, _info in runs.run_successors(node))

                for nxt, weight in successors:
                    next_elapsed = elapsed_s + int(weight)
                    if next_elapsed <= 0 or next_elapsed > max_elapsed_s:
                        continue
                    next_station = node_off[nxt]
                    next_reward = reward_s
                    next_mask = mask
                    bit = reward_bit.get(next_station)
                    if bit is not None and not (next_mask & (1 << bit)):
                        next_mask |= 1 << bit
                        next_reward += search_rewards[next_station]
                    if _beam_label_dominated(
                            labels, labels_by_node.get(nxt, ()),
                            next_elapsed, next_reward, next_mask):
                        continue

                    next_id = len(labels)
                    labels.append((nxt, label_id, next_elapsed, next_reward, next_mask))
                    bucket = labels_by_node[nxt]
                    bucket.append(next_id)
                    bucket.sort(key=lambda lid: (
                        labels[lid][2] + args.column_penalty - labels[lid][3],
                        labels[lid][2],
                        -labels[lid][3],
                    ))
                    if len(bucket) > args.max_labels_per_node:
                        del bucket[args.max_labels_per_node:]
                    if next_id in bucket:
                        heapq.heappush(
                            heap,
                            (next_elapsed + args.column_penalty - next_reward,
                             next_elapsed, next_id),
                        )

            candidate_ids = sorted(set(candidate_ids), key=lambda lid: (
                labels[lid][2] + args.column_penalty - labels[lid][3],
                -labels[lid][3],
                labels[lid][2],
            ))
            emitted = 0
            for label_id in candidate_ids:
                if emitted >= args.emit_per_source:
                    break
                nodes = _label_path_nodes(labels, label_id)
                if len(nodes) < 2:
                    continue
                covered = {node_off[n] for n in nodes} & si.canonical_stations
                if len(covered) < args.min_covered:
                    continue
                full_reward = sum(all_rewards.get(station, 0.0) for station in covered)
                elapsed = _elapsed(node_t[nodes[0]], node_t[nodes[-1]])
                reduced = elapsed + args.column_penalty - full_reward
                if full_reward < args.min_reward_s or reduced > args.max_reduced_cost_s:
                    continue
                key = (node_t[nodes[0]], node_off[nodes[0]],
                       node_t[nodes[-1]], node_off[nodes[-1]],
                       tuple(sorted(covered)))
                if key in seen:
                    continue
                seen.add(key)
                record = _column_from_nodes(
                    prob, args.label, len(base_rows) + len(generated) + 1, nodes)
                record["pricing_dual_beam"] = {
                    "dual_file": args.duals,
                    "seed": seed_spec,
                    "source_index": local_idx,
                    "source_station": start_station,
                    "source_time": int(node_t[src]),
                    "search_reward_s": labels[label_id][3],
                    "path_reward_s": full_reward,
                    "reduced_cost_s": reduced,
                    "expansions": expansions,
                    "reward_stations": args.top_reward_stations,
                    "run_mode": run_mode,
                }
                generated.append(record)
                f.write(json.dumps(record, sort_keys=True) + "\n")
                emitted += 1

            if sidx % 5 == 0 or sidx == len(source_events):
                print(f"beam-priced sources {sidx}/{len(source_events)} "
                      f"generated={len(generated)} last_expansions={expansions} "
                      f"wall={time.time() - t0:.1f}s",
                      flush=True)

    print(f"beam reward stations={len(target_rewards)}")
    print(f"wrote {len(base_rows)} base + {len(generated)} beam-priced columns -> {out}")
    return 0


def _selected_connector_rows(G, tables, rows, runs, cap):
    node_t = tables[2]
    stop_events = _stop_index(tables[1], node_t)
    out = []
    for pos, (row_i, row_j) in enumerate(zip(rows, rows[1:])):
        start_i, end_i = _row_endpoint_nodes(stop_events, node_t, row_i)
        start_j, _end_j = _row_endpoint_nodes(stop_events, node_t, row_j)
        gap = _elapsed(node_t[end_i], node_t[start_j])
        if end_i == start_j:
            cost = 0
        else:
            tgt, _prev = _dijkstra_to_node(
                G._adj, end_i, start_j, cap, runs=runs, max_cost=gap)
            cost = gap if tgt is not None else None
        out.append({
            "pos": pos,
            "from": row_i,
            "to": row_j,
            "from_start_node": start_i,
            "from_end_node": end_i,
            "to_start_node": start_j,
            "connector_s": cost,
            "phase_gap_s": gap,
        })
    return out


def cmd_price_connectors(args) -> int:
    import pickle

    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    _dual_data, rewards = _load_dual_rewards(Path(args.duals)) if args.duals else ({}, {})

    si = StationIndex.load()
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, _node_stop, node_t = tables
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)

    connector_jobs = []
    for solution_file in args.solution_files:
        rows = _load_column_rows(Path(solution_file))
        connectors = _selected_connector_rows(
            G, tables, rows, runs, args.connector_cap)
        connectors = [
            conn for conn in connectors
            if conn["connector_s"] is not None
            and conn["connector_s"] >= int(args.min_connector_minutes * 60)
        ]
        connectors.sort(key=lambda conn: (-conn["connector_s"], conn["pos"]))
        if args.top_connectors:
            connectors = connectors[:args.top_connectors]
        for conn in connectors:
            conn["solution_file"] = solution_file
            conn["solution_rows"] = rows
            connector_jobs.append(conn)

    if not connector_jobs:
        raise SystemExit("no connector pricing jobs selected")

    source_points = {item.strip() for item in str(args.source_points).split(",") if item.strip()}
    bad_source_points = sorted(source_points - {"start", "end"})
    if bad_source_points:
        raise SystemExit(f"--source-points must contain only start/end, got {bad_source_points}")
    if not source_points:
        raise SystemExit("--source-points must include at least one of start,end")

    max_elapsed_s = int(args.max_elapsed_minutes * 60)
    max_score = INF if args.max_score_s is None else float(args.max_score_s)
    generated = []
    t0 = time.time()
    seen = set()
    for row in base_rows:
        key = (
            row["start"]["time"],
            row["start"]["station"],
            row["end"]["time"],
            row["end"]["station"],
            tuple(row.get("covered_stations", [])),
        )
        seen.add(key)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")

        for cidx, conn in enumerate(connector_jobs, start=1):
            rows = conn["solution_rows"]
            start = conn["pos"] + 1
            end = min(len(rows), start + args.lookahead_columns)
            downstream = rows[start:end]
            target_scores = {}
            for row in downstream:
                for station in [row["start"]["station"], row["end"]["station"]]:
                    target_scores[station] = max(target_scores.get(station, 0.0),
                                                 rewards.get(station, 0.0))
                for station in row.get("covered_stations", []):
                    target_scores[station] = max(target_scores.get(station, 0.0),
                                                 rewards.get(station, 0.0))
            if args.include_global_duals and rewards:
                for station, reward in sorted(rewards.items(), key=lambda item: -item[1]):
                    if station in si.canonical_stations:
                        target_scores[station] = max(target_scores.get(station, 0.0), reward)
                    if len(target_scores) >= args.max_targets_per_connector:
                        break

            targets = sorted(target_scores, key=lambda st: (-target_scores[st], st))
            if args.max_targets_per_connector:
                targets = targets[:args.max_targets_per_connector]

            sources = []
            if "end" in source_points:
                sources.append(("end", conn["from_end_node"]))
            if "start" in source_points:
                sources.append(("start", conn["from_start_node"]))

            for source_point, src in sources:
                for target in targets:
                    if node_off[src] == target:
                        continue
                    tgt, prev = _dijkstra_to_station_bounded(
                        G._adj,
                        node_off,
                        src,
                        target,
                        args.connector_cap,
                        runs=runs,
                        max_cost=max_elapsed_s,
                    )
                    if tgt is None:
                        continue
                    elapsed = _elapsed(node_t[src], node_t[tgt])
                    if elapsed <= 0 or elapsed > max_elapsed_s:
                        continue
                    seg = _path_from_prev(prev, src, tgt)
                    if seg is None:
                        continue
                    nodes = [src] + seg
                    covered = {node_off[n] for n in nodes} & si.canonical_stations
                    if len(covered) < args.min_covered:
                        continue
                    reward = (_row_dual_reward({"covered_stations": covered}, rewards)
                              if rewards else args.cover_reward_s * len(covered))
                    score = elapsed - reward - args.connector_credit * conn["connector_s"]
                    if score > max_score:
                        continue
                    key = (node_t[src], node_off[src], node_t[tgt], node_off[tgt],
                           tuple(sorted(covered)))
                    if key in seen:
                        continue
                    seen.add(key)
                    record = _column_from_nodes(
                        prob, args.label, len(base_rows) + len(generated) + 1, nodes)
                    record["pricing_connector"] = {
                        "solution_file": conn["solution_file"],
                        "position": conn["pos"],
                        "source_point": source_point,
                        "from_column": conn["from"]["column_id"],
                        "to_column": conn["to"]["column_id"],
                        "connector_s": conn["connector_s"],
                        "phase_gap_s": conn["phase_gap_s"],
                        "target_station": target,
                        "path_reward_s": reward,
                        "score_s": score,
                        "run_mode": run_mode,
                        "lookahead_columns": args.lookahead_columns,
                        "downstream_columns": [row["column_id"] for row in downstream],
                    }
                    generated.append(record)
                    f.write(json.dumps(record, sort_keys=True) + "\n")

            if cidx % 5 == 0 or cidx == len(connector_jobs):
                print(f"priced connectors {cidx}/{len(connector_jobs)} "
                      f"generated={len(generated)} wall={time.time() - t0:.1f}s",
                      flush=True)

    print(f"connector jobs={len(connector_jobs)}")
    print(f"wrote {len(base_rows)} base + {len(generated)} connector-priced columns -> {out}")
    return 0


def _parse_choice_set(value, allowed, name):
    choices = {item.strip() for item in str(value).split(",") if item.strip()}
    bad = sorted(choices - set(allowed))
    if bad:
        raise SystemExit(f"{name} must contain only {sorted(allowed)}, got {bad}")
    if not choices:
        raise SystemExit(f"{name} must include at least one value")
    return choices


def cmd_price_window_events(args) -> int:
    import pickle

    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    _dual_data, rewards = _load_dual_rewards(Path(args.duals)) if args.duals else ({}, {})
    widths = _parse_ints(args.widths)
    if not widths:
        raise SystemExit("--widths must include at least one positive integer")
    widths = sorted({w for w in widths if w > 0})
    source_points = _parse_choice_set(args.source_points, {"start", "end"}, "--source-points")
    target_points = _parse_choice_set(
        args.target_points, {"start", "end", "next-start"}, "--target-points")

    si = StationIndex.load()
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    stop_events = _stop_index(node_stop, node_t)
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)

    jobs = []
    for solution_file in args.solution_files:
        rows = _load_column_rows(Path(solution_file))
        endpoints = [_row_endpoint_nodes(stop_events, node_t, row) for row in rows]
        for start in range(0, len(rows), max(1, args.window_stride)):
            for width in widths:
                end = start + width - 1
                if end >= len(rows):
                    continue
                window_rows = rows[start:end + 1]
                window_stations = set()
                for row in window_rows:
                    window_stations.update(
                        station for station in row.get("covered_stations", [])
                        if station in si.canonical_stations)
                window_reward = (_row_dual_reward(
                    {"covered_stations": window_stations}, rewards)
                    if rewards else args.cover_reward_s * len(window_stations))
                for source_point in source_points:
                    src = endpoints[start][0] if source_point == "start" else endpoints[start][1]
                    for target_point in target_points:
                        if target_point == "start":
                            target_index = end
                            target = endpoints[target_index][0]
                        elif target_point == "end":
                            target_index = end
                            target = endpoints[target_index][1]
                        else:
                            target_index = end + 1
                            if target_index >= len(rows):
                                continue
                            target = endpoints[target_index][0]
                        elapsed = _elapsed(node_t[src], node_t[target])
                        if elapsed <= 0 or elapsed > int(args.max_window_minutes * 60):
                            continue
                        density_score = elapsed - window_reward
                        jobs.append({
                            "score": density_score,
                            "solution_file": solution_file,
                            "rows": rows,
                            "start": start,
                            "end": end,
                            "source_point": source_point,
                            "target_point": target_point,
                            "target_index": target_index,
                            "src": src,
                            "target": target,
                            "elapsed_s": elapsed,
                            "window_reward_s": window_reward,
                            "window_stations": window_stations,
                            "window_columns": [row["column_id"] for row in window_rows],
                        })
    jobs.sort(key=lambda job: (job["score"], -len(job["window_stations"])))
    if args.high_score_windows:
        jobs.sort(key=lambda job: (-job["score"], job["elapsed_s"]))
    if args.max_windows:
        jobs = jobs[:args.max_windows]
    if not jobs:
        raise SystemExit("no window-event pricing jobs selected")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    generated = []
    seen = set()
    for row in base_rows:
        seen.add((
            row["start"]["time"],
            row["start"]["station"],
            row["end"]["time"],
            row["end"]["station"],
            tuple(row.get("covered_stations", [])),
        ))
    max_score = INF if args.max_score_s is None else float(args.max_score_s)
    max_cost_margin = int(args.max_cost_margin_minutes * 60)
    t0 = time.time()

    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")

        for jidx, job in enumerate(jobs, start=1):
            src = job["src"]
            target = job["target"]
            if src == target:
                continue
            max_cost = job["elapsed_s"] + max_cost_margin
            tgt, prev = _dijkstra_to_node(
                G._adj, src, target, args.connector_cap, runs=runs, max_cost=max_cost)
            if tgt is None:
                continue
            seg = _path_from_prev(prev, src, target)
            if seg is None:
                continue
            nodes = [src] + seg
            elapsed = sum(weight for weight, _mode in _transition_series(prob, nodes))
            if elapsed <= 0 or elapsed > int(args.max_elapsed_minutes * 60):
                continue
            covered = {node_off[n] for n in nodes} & si.canonical_stations
            if len(covered) < args.min_covered:
                continue
            reward = (_row_dual_reward({"covered_stations": covered}, rewards)
                      if rewards else args.cover_reward_s * len(covered))
            score = elapsed - reward - args.window_credit * job["elapsed_s"]
            if score > max_score:
                continue
            key = (node_t[src], node_off[src], node_t[target], node_off[target],
                   tuple(sorted(covered)))
            if key in seen:
                continue
            seen.add(key)
            record = _column_from_nodes(
                prob, args.label, len(base_rows) + len(generated) + 1, nodes)
            record["pricing_window_event"] = {
                "solution_file": job["solution_file"],
                "start_index": job["start"],
                "end_index": job["end"],
                "source_point": job["source_point"],
                "target_point": job["target_point"],
                "target_index": job["target_index"],
                "window_elapsed_s": job["elapsed_s"],
                "window_elapsed": hms(job["elapsed_s"]),
                "window_reward_s": job["window_reward_s"],
                "path_reward_s": reward,
                "score_s": score,
                "run_mode": run_mode,
                "window_columns": job["window_columns"],
            }
            generated.append(record)
            f.write(json.dumps(record, sort_keys=True) + "\n")

            if jidx % 25 == 0 or jidx == len(jobs):
                print(f"priced window-events {jidx}/{len(jobs)} "
                      f"generated={len(generated)} wall={time.time() - t0:.1f}s",
                      flush=True)

    print(f"window-event jobs={len(jobs)}")
    print(f"wrote {len(base_rows)} base + {len(generated)} window-event columns -> {out}")
    return 0


def _load_target_stations(args, si):
    targets = set()
    if args.target_stations:
        targets.update(
            station.strip()
            for station in str(args.target_stations).split(",")
            if station.strip())
    if args.targets_file:
        data = json.loads(Path(args.targets_file).read_text())
        if isinstance(data, list):
            targets.update(str(station) for station in data)
        elif "relaxed_uncovered_stations" in data:
            targets.update(str(station) for station in data["relaxed_uncovered_stations"])
        elif "stations" in data:
            targets.update(str(station) for station in data["stations"])
        else:
            raise SystemExit(
                f"{args.targets_file} must be a list or contain relaxed_uncovered_stations")
    targets = {station for station in targets if station in si.canonical_stations}
    if not targets:
        raise SystemExit("no canonical target stations supplied")
    return targets


def _endpoint_candidates(rows, stop_events, node_t, endpoint_points, time_window,
                         max_candidates):
    candidates = []
    seen = set()
    for ridx, row in enumerate(rows):
        start, end = _row_endpoint_nodes(stop_events, node_t, row)
        for point, node in (("start", start), ("end", end)):
            if point not in endpoint_points:
                continue
            if not _time_in_window(node_t[node], time_window):
                continue
            key = (node, point)
            if key in seen:
                continue
            seen.add(key)
            candidates.append({
                "row_index": ridx,
                "row": row,
                "point": point,
                "node": node,
                "time": int(node_t[node]),
            })
    candidates.sort(key=lambda item: (
        item["time"],
        item["row"].get("column_id", ""),
        item["point"],
    ))
    if max_candidates:
        candidates = candidates[:max_candidates]
    return candidates


def cmd_price_cluster_corridors(args) -> int:
    import pickle

    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    destination_rows = (
        _load_column_rows(Path(args.destination_columns))
        if args.destination_columns else base_rows)
    source_rows = _load_column_rows(Path(args.source_solution))

    si = StationIndex.load()
    targets = _load_target_stations(args, si)
    via_stations = sorted(targets)
    if args.via_stations:
        via_stations = [
            station.strip()
            for station in str(args.via_stations).split(",")
            if station.strip() and station.strip() in targets
        ]
    if args.max_vias:
        via_stations = via_stations[:args.max_vias]
    if not via_stations:
        raise SystemExit("no canonical via stations selected")

    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    stop_events = _stop_index(node_stop, node_t)
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)

    source_points = _parse_choice_set(
        args.source_points, {"start", "end"}, "--source-points")
    destination_points = _parse_choice_set(
        args.destination_points, {"start", "end"}, "--destination-points")
    source_window = _parse_week_time_window(args.source_time_window)
    destination_window = _parse_week_time_window(args.destination_time_window)
    sources = _endpoint_candidates(
        source_rows, stop_events, node_t, source_points, source_window, args.max_sources)
    destinations = _endpoint_candidates(
        destination_rows,
        stop_events,
        node_t,
        destination_points,
        destination_window,
        args.max_destinations,
    )
    if args.destination_mode == "station":
        by_station = {}
        for item in destinations:
            station = node_off[item["node"]]
            if station not in by_station:
                by_station[station] = dict(item, station=station)
        destinations = sorted(
            by_station.values(),
            key=lambda item: (item["time"], item["station"], item["row"].get("column_id", "")),
        )
        if args.max_destinations:
            destinations = destinations[:args.max_destinations]
    if not sources:
        raise SystemExit("no source endpoint candidates selected")
    if not destinations:
        raise SystemExit("no destination endpoint candidates selected")

    max_elapsed_s = int(args.max_elapsed_minutes * 60)
    max_leg_s = int(args.max_leg_minutes * 60) if args.max_leg_minutes else max_elapsed_s
    max_score = INF if args.max_score_s is None else float(args.max_score_s)
    generated = []
    seen = set()
    first_leg_cache = {}
    second_leg_cache = {}
    for row in base_rows:
        seen.add((
            row["start"]["time"],
            row["start"]["station"],
            row["end"]["time"],
            row["end"]["station"],
            tuple(row.get("covered_stations", [])),
        ))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    jobs = 0
    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")

        for source in sources:
            src = source["node"]
            for destination in destinations:
                if args.destination_mode == "exact":
                    initial_dst = destination["node"]
                    total_gap = _elapsed(node_t[src], node_t[initial_dst])
                    if total_gap <= 0 or total_gap > max_elapsed_s:
                        continue
                else:
                    initial_dst = None
                for via in via_stations:
                    if node_off[src] == via:
                        via_node = src
                        seg1 = []
                    else:
                        first_key = (src, via)
                        if first_key not in first_leg_cache:
                            hit, prev1 = _dijkstra_to_station_bounded(
                                G._adj,
                                node_off,
                                src,
                                via,
                                args.connector_cap,
                                runs=runs,
                                max_cost=max_leg_s,
                            )
                            if hit is None:
                                first_leg_cache[first_key] = None
                            else:
                                first_leg_cache[first_key] = (
                                    hit,
                                    _path_from_prev(prev1, src, hit),
                                )
                        cached_first = first_leg_cache[first_key]
                        if cached_first is None:
                            continue
                        via_node, seg1 = cached_first
                        if seg1 is None:
                            continue
                    if args.destination_mode == "exact":
                        dst = initial_dst
                        remaining = _elapsed(node_t[via_node], node_t[dst])
                        if remaining > max_elapsed_s:
                            continue
                        max_second_leg = min(max_leg_s, remaining, max_elapsed_s)
                    else:
                        dst = None
                        max_second_leg = max_leg_s
                    if dst is not None and via_node == dst:
                        seg2 = []
                    else:
                        second_key = (
                            args.destination_mode,
                            via_node,
                            destination.get("station", dst),
                            max_second_leg,
                        )
                        if second_key not in second_leg_cache:
                            if args.destination_mode == "exact":
                                hit, prev2 = _dijkstra_to_node(
                                    G._adj,
                                    via_node,
                                    dst,
                                    args.connector_cap,
                                    runs=runs,
                                    max_cost=max_second_leg,
                                )
                            else:
                                hit, prev2 = _dijkstra_to_station_bounded(
                                    G._adj,
                                    node_off,
                                    via_node,
                                    destination["station"],
                                    args.connector_cap,
                                    runs=runs,
                                    max_cost=max_second_leg,
                                )
                            if hit is None:
                                second_leg_cache[second_key] = None
                            else:
                                second_leg_cache[second_key] = (
                                    hit,
                                    _path_from_prev(prev2, via_node, hit),
                                )
                        cached_second = second_leg_cache[second_key]
                        if cached_second is None:
                            continue
                        hit, seg2 = cached_second
                        dst = hit
                        if seg2 is None:
                            continue
                    nodes = [src] + seg1
                    if seg2:
                        if nodes and nodes[-1] == via_node:
                            nodes.extend(seg2)
                        else:
                            nodes.append(via_node)
                            nodes.extend(seg2)
                    elif not nodes or nodes[-1] != via_node:
                        nodes.append(via_node)
                    if nodes[-1] != dst:
                        continue
                    elapsed = sum(weight for weight, _mode in _transition_series(prob, nodes))
                    if elapsed <= 0 or elapsed > max_elapsed_s:
                        continue
                    covered = {node_off[n] for n in nodes} & si.canonical_stations
                    target_hits = covered & targets
                    if len(target_hits) < args.min_target_hits:
                        continue
                    if len(covered) < args.min_covered:
                        continue
                    reward = (
                        args.target_reward_s * len(target_hits)
                        + args.cover_reward_s * len(covered)
                    )
                    score = elapsed + args.column_penalty - reward
                    if max_score < INF and score > max_score:
                        continue
                    key = (
                        node_t[src],
                        node_off[src],
                        node_t[dst],
                        node_off[dst],
                        tuple(sorted(covered)),
                        tuple(sorted(target_hits)),
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    record = _column_from_nodes(
                        prob, args.label, len(base_rows) + len(generated) + 1, nodes)
                    record["pricing_cluster_corridor"] = {
                        "source_solution": args.source_solution,
                        "destination_columns": args.destination_columns,
                        "source_column": source["row"].get("column_id"),
                        "source_point": source["point"],
                        "destination_column": destination["row"].get("column_id"),
                        "destination_point": destination["point"],
                        "destination_mode": args.destination_mode,
                        "destination_station": destination.get("station", node_off[dst]),
                        "via_station": via,
                        "target_stations": sorted(targets),
                        "target_hits": sorted(target_hits),
                        "target_hit_count": len(target_hits),
                        "reward_s": reward,
                        "score_s": score,
                        "run_mode": run_mode,
                    }
                    generated.append(record)
                    f.write(json.dumps(record, sort_keys=True) + "\n")
                    if args.max_generated and len(generated) >= args.max_generated:
                        break
                jobs += 1
                if jobs % 10 == 0:
                    print(f"cluster corridors jobs={jobs} generated={len(generated)} "
                          f"wall={time.time() - t0:.1f}s", flush=True)
                if args.max_generated and len(generated) >= args.max_generated:
                    break
            if args.max_generated and len(generated) >= args.max_generated:
                break

    print(f"cluster-corridor sources={len(sources)} destinations={len(destinations)} "
          f"vias={len(via_stations)} jobs={jobs}")
    print(f"wrote {len(base_rows)} base + {len(generated)} cluster-corridor columns -> {out}")
    return 0


def _parse_station_list(value, si, name):
    stations = [
        station.strip()
        for station in str(value or "").split(",")
        if station.strip()
    ]
    bad = [station for station in stations if station not in si.canonical_stations]
    if bad:
        raise SystemExit(f"{name} contains non-canonical station ids, e.g. {bad[:8]}")
    return stations


def _load_anchor_sequences(args, si):
    sequences = []
    if args.anchor_sequences:
        for raw_sequence in str(args.anchor_sequences).split(";"):
            raw_sequence = raw_sequence.strip()
            if not raw_sequence:
                continue
            sequences.append(
                _parse_station_list(raw_sequence, si, "--anchor-sequences")
            )
    if args.anchor_sequences_file:
        data = json.loads(Path(args.anchor_sequences_file).read_text())
        if not isinstance(data, list):
            raise SystemExit("--anchor-sequences-file must contain a JSON list")
        for item in data:
            if isinstance(item, str):
                sequences.append(_parse_station_list(item, si, "--anchor-sequences-file"))
            elif isinstance(item, list):
                sequence = [str(station) for station in item]
                bad = [station for station in sequence if station not in si.canonical_stations]
                if bad:
                    raise SystemExit(
                        "--anchor-sequences-file contains non-canonical station ids, "
                        f"e.g. {bad[:8]}")
                sequences.append(sequence)
            else:
                raise SystemExit(
                    "--anchor-sequences-file entries must be strings or station-id lists")
    if not sequences:
        raise SystemExit("provide --anchor-sequences or --anchor-sequences-file")
    seen = set()
    unique = []
    for sequence in sequences:
        key = tuple(sequence)
        if key in seen:
            continue
        seen.add(key)
        unique.append(sequence)
    return unique


def cmd_price_anchor_sequences(args) -> int:
    import pickle

    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    source_rows = _load_column_rows(Path(args.source_columns))

    si = StationIndex.load()
    sequences = _load_anchor_sequences(args, si)
    targets = set(_parse_station_list(args.target_stations, si, "--target-stations"))

    source_row = next(
        (row for row in source_rows if row.get("column_id") == args.source_column_id),
        None,
    )
    if source_row is None:
        raise SystemExit(f"source column not found: {args.source_column_id}")

    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, runs or _NoRunLayer(), si)

    source_start, source_end = _row_endpoint_nodes_from_tables(
        node_stop, node_t, source_row)
    src = source_start if args.source_point == "start" else source_end
    max_elapsed_s = int(args.max_elapsed_minutes * 60)
    max_leg_s = int(args.max_leg_minutes * 60) if args.max_leg_minutes else max_elapsed_s
    max_end_time = int(args.max_end_time) if args.max_end_time is not None else None

    seen = set()
    for row in base_rows:
        seen.add((
            row["start"]["time"],
            row["start"]["station"],
            row["end"]["time"],
            row["end"]["station"],
            tuple(row.get("covered_stations", [])),
        ))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    generated = []
    skipped = Counter()
    t0 = time.time()

    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")

        for sidx, sequence in enumerate(sequences, start=1):
            nodes = [src]
            current = src
            elapsed_so_far = 0
            failed = None
            for anchor in sequence:
                if node_off[current] == anchor:
                    continue
                remaining_elapsed = max_elapsed_s - elapsed_so_far
                if remaining_elapsed <= 0:
                    failed = (anchor, "elapsed_cap")
                    break
                remaining_time = None
                if max_end_time is not None:
                    remaining_time = max_end_time - int(node_t[current])
                    if remaining_time < 0:
                        failed = (anchor, "end_time_cap")
                        break
                max_cost = min(max_leg_s, remaining_elapsed)
                if remaining_time is not None:
                    max_cost = min(max_cost, remaining_time)
                hit, prev = _dijkstra_to_station_bounded(
                    G._adj,
                    node_off,
                    current,
                    anchor,
                    args.connector_cap,
                    runs=runs,
                    max_cost=max_cost,
                )
                if hit is None:
                    failed = (anchor, "unreachable")
                    break
                seg = _path_from_prev(prev, current, hit)
                if seg is None:
                    failed = (anchor, "missing_path")
                    break
                elapsed_so_far += sum(
                    weight for weight, _mode in _transition_series(prob, [current] + seg))
                nodes.extend(seg)
                current = hit

            if failed is not None:
                skipped[failed[1]] += 1
                if args.verbose:
                    print(f"skip sequence {sidx}: failed {failed} {sequence}", flush=True)
                continue
            elapsed = elapsed_so_far
            if elapsed <= 0 or elapsed > max_elapsed_s:
                skipped["too_long"] += 1
                continue
            if max_end_time is not None and int(node_t[nodes[-1]]) > max_end_time:
                skipped["too_late"] += 1
                continue
            covered = {node_off[n] for n in nodes} & si.canonical_stations
            if len(covered) < args.min_covered:
                skipped["too_few_covered"] += 1
                continue
            target_hits = covered & targets if targets else set()
            if args.min_target_hits and len(target_hits) < args.min_target_hits:
                skipped["too_few_target_hits"] += 1
                continue
            key = (
                int(node_t[nodes[0]]),
                node_off[nodes[0]],
                int(node_t[nodes[-1]]),
                node_off[nodes[-1]],
                tuple(sorted(covered)),
            )
            if key in seen:
                skipped["duplicate"] += 1
                continue
            seen.add(key)
            record = _column_from_nodes(
                prob, args.label, len(base_rows) + len(generated) + 1, nodes)
            record["pricing_anchor_sequence"] = {
                "source_columns": args.source_columns,
                "source_column": args.source_column_id,
                "source_point": args.source_point,
                "anchor_sequence": sequence,
                "target_stations": sorted(targets),
                "target_hits": sorted(target_hits),
                "target_hit_count": len(target_hits),
                "run_mode": run_mode,
                "max_end_time": max_end_time,
            }
            generated.append(record)
            f.write(json.dumps(record, sort_keys=True) + "\n")

    print(f"anchor-sequences source={args.source_column_id} sequences={len(sequences)} "
          f"generated={len(generated)} skipped={dict(skipped)} "
          f"wall={time.time() - t0:.1f}s")
    print(f"wrote {len(base_rows)} base + {len(generated)} anchor-sequence columns -> {out}")
    return 0


def _load_stage_groups(args, si):
    groups = []
    if args.stage_groups:
        for raw_group in str(args.stage_groups).split(";"):
            group = _parse_station_list(raw_group, si, "--stage-groups")
            if group:
                groups.append(group)
    if args.stage_groups_file:
        data = json.loads(Path(args.stage_groups_file).read_text())
        if not isinstance(data, list):
            raise SystemExit("--stage-groups-file must contain a JSON list")
        for item in data:
            if isinstance(item, str):
                group = _parse_station_list(item, si, "--stage-groups-file")
            elif isinstance(item, list):
                group = [str(station) for station in item]
                bad = [station for station in group if station not in si.canonical_stations]
                if bad:
                    raise SystemExit(
                        "--stage-groups-file contains non-canonical station ids, "
                        f"e.g. {bad[:8]}")
            else:
                raise SystemExit(
                    "--stage-groups-file entries must be strings or station-id lists")
            if group:
                groups.append(group)
    if not groups:
        raise SystemExit("provide --stage-groups or --stage-groups-file")
    unique_groups = []
    for group in groups:
        seen = set()
        unique_groups.append([
            station for station in group if not (station in seen or seen.add(station))
        ])
    return unique_groups


def _parse_stage_min_hits(value, stage_count):
    if value is None:
        return [1] * stage_count
    values = [
        int(part.strip())
        for part in str(value).split(",")
        if part.strip()
    ]
    if not values:
        return [1] * stage_count
    if len(values) == 1:
        return values * stage_count
    if len(values) != stage_count:
        raise SystemExit(
            f"--stage-min-hits must contain one value or {stage_count} values")
    if any(value < 0 for value in values):
        raise SystemExit("--stage-min-hits values must be non-negative")
    return values


def _load_resource_groups(args, si):
    groups = []

    def add_group(name, raw_stations, source):
        if isinstance(raw_stations, str):
            stations = _parse_station_list(raw_stations, si, source)
        elif isinstance(raw_stations, list):
            stations = [str(station) for station in raw_stations]
            bad = [station for station in stations if station not in si.canonical_stations]
            if bad:
                raise SystemExit(
                    f"{source} contains non-canonical station ids, e.g. {bad[:8]}")
        else:
            raise SystemExit(f"{source} station group must be a string or list")
        seen = set()
        stations = [
            station for station in stations
            if not (station in seen or seen.add(station))
        ]
        if not stations:
            return
        groups.append({
            "name": str(name or f"resource{len(groups) + 1}"),
            "stations": stations,
        })

    if args.resource_groups:
        for raw_group in str(args.resource_groups).split(";"):
            raw_group = raw_group.strip()
            if not raw_group:
                continue
            if ":" in raw_group:
                name, raw_stations = raw_group.split(":", 1)
            elif "=" in raw_group:
                name, raw_stations = raw_group.split("=", 1)
            else:
                name, raw_stations = f"resource{len(groups) + 1}", raw_group
            add_group(name.strip(), raw_stations, "--resource-groups")

    if args.resource_groups_file:
        data = json.loads(Path(args.resource_groups_file).read_text())
        if isinstance(data, dict):
            iterable = data.items()
        elif isinstance(data, list):
            iterable = []
            for item in data:
                if isinstance(item, dict):
                    name = item.get("name", f"resource{len(groups) + len(iterable) + 1}")
                    stations = item.get("stations", item.get("group"))
                    iterable.append((name, stations))
                elif isinstance(item, str):
                    if ":" in item:
                        name, stations = item.split(":", 1)
                    elif "=" in item:
                        name, stations = item.split("=", 1)
                    else:
                        name, stations = (
                            f"resource{len(groups) + len(iterable) + 1}",
                            item,
                        )
                    iterable.append((name, stations))
                else:
                    raise SystemExit(
                        "--resource-groups-file list entries must be objects or strings")
        else:
            raise SystemExit("--resource-groups-file must contain a JSON object or list")
        for name, stations in iterable:
            add_group(name, stations, "--resource-groups-file")

    if not groups:
        raise SystemExit("provide --resource-groups or --resource-groups-file")

    seen_names = Counter(group["name"] for group in groups)
    if any(count > 1 for count in seen_names.values()):
        renamed = []
        used = Counter()
        for group in groups:
            used[group["name"]] += 1
            if seen_names[group["name"]] == 1:
                renamed.append(group)
            else:
                renamed.append({
                    **group,
                    "name": f"{group['name']}_{used[group['name']]}",
                })
        groups = renamed
    return groups


def _parse_resource_min_hits(value, resource_count):
    if value is None:
        return [1] * resource_count
    values = [
        int(part.strip())
        for part in str(value).split(",")
        if part.strip()
    ]
    if not values:
        return [1] * resource_count
    if len(values) == 1:
        values = values * resource_count
    if len(values) != resource_count:
        raise SystemExit(
            f"--resource-min-hits must contain one value or {resource_count} values")
    if any(value < 0 for value in values):
        raise SystemExit("--resource-min-hits values must be non-negative")
    return values


def cmd_price_resource_chains(args) -> int:
    import pickle

    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    source_rows = _load_column_rows(Path(args.source_columns))

    si = StationIndex.load()
    resources = _load_resource_groups(args, si)
    resource_min_hits = _parse_resource_min_hits(
        args.resource_min_hits,
        len(resources),
    )
    resource_sets = [set(resource["stations"]) for resource in resources]
    resource_targets = {station for group in resource_sets for station in group}
    targets = set(resource_targets)
    if args.target_stations or args.targets_file:
        targets.update(_load_target_stations(args, si))
    target_order = []
    for resource in resources:
        for station in resource["stations"]:
            if station not in target_order:
                target_order.append(station)
    for station in _parse_station_list(args.target_stations, si, "--target-stations"):
        if station not in target_order:
            target_order.append(station)
    for station in sorted(targets):
        if station not in target_order:
            target_order.append(station)
    target_order_index = {
        station: idx
        for idx, station in enumerate(target_order)
    }

    finals = _parse_station_list(args.final_stations, si, "--final-stations")
    if not finals:
        raise SystemExit("--final-stations must include at least one canonical station")
    min_resource_count = int(args.min_resource_count)
    if min_resource_count < 0 or min_resource_count > len(resources):
        raise SystemExit("--min-resource-count must be between 0 and the resource count")
    frontier_min_resource_count = (
        min_resource_count
        if args.frontier_min_resource_count is None
        else int(args.frontier_min_resource_count)
    )
    if frontier_min_resource_count < 0 or frontier_min_resource_count > len(resources):
        raise SystemExit(
            "--frontier-min-resource-count must be between 0 and the resource count")
    frontier_min_target_hits = (
        int(args.min_target_hits)
        if args.frontier_min_target_hits is None
        else int(args.frontier_min_target_hits)
    )
    frontier_min_covered = (
        int(args.min_covered)
        if args.frontier_min_covered is None
        else int(args.frontier_min_covered)
    )
    resource_name_to_index = {resource["name"]: idx for idx, resource in enumerate(resources)}
    required_resource_names = [
        name.strip()
        for name in str(args.require_resources or "").split(",")
        if name.strip()
    ]
    missing_required_resources = [
        name for name in required_resource_names if name not in resource_name_to_index
    ]
    if missing_required_resources:
        raise SystemExit(
            "--require-resources contains unknown names, "
            f"e.g. {missing_required_resources[:8]}")
    required_resource_mask = 0
    for name in required_resource_names:
        required_resource_mask |= 1 << resource_name_to_index[name]

    source_row = next(
        (row for row in source_rows if row.get("column_id") == args.source_column_id),
        None,
    )
    if source_row is None:
        raise SystemExit(f"source column not found: {args.source_column_id}")

    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, runs or _NoRunLayer(), si)

    source_start, source_end = _row_endpoint_nodes_from_tables(
        node_stop, node_t, source_row)
    src = source_start if args.source_point == "start" else source_end
    max_end_time = int(args.max_end_time) if args.max_end_time is not None else None
    if max_end_time is not None and int(node_t[src]) > max_end_time:
        raise SystemExit("source starts after --max-end-time")
    max_elapsed_s = int(args.max_elapsed_minutes * 60)
    max_leg_s = int(args.max_leg_minutes * 60) if args.max_leg_minutes else max_elapsed_s
    max_score = INF if args.max_score_s is None else float(args.max_score_s)
    time_bucket_s = max(1, int(args.time_bucket_s))
    max_expand_targets = max(0, int(args.max_expand_targets))

    leg_cache = {}
    final_cache = {}
    persistent_cache = {}
    persistent_cache_new = []
    persistent_cache_hits = 0
    persistent_cache_misses = 0
    persistent_cache_path = (
        None
        if args.no_resource_chain_cache
        else Path(args.resource_chain_cache)
    )
    if persistent_cache_path is not None and persistent_cache_path.exists():
        with persistent_cache_path.open() as f:
            for line in f:
                if not line.strip():
                    continue
                try:
                    item = json.loads(line)
                except json.JSONDecodeError:
                    continue
                key = tuple(item.get("key", []))
                if not key:
                    continue
                if item.get("reachable"):
                    persistent_cache[key] = (
                        int(item["hit"]),
                        [int(node) for node in item.get("seg", [])],
                    )
                else:
                    persistent_cache[key] = None
    t0 = time.time()
    persistent_cache_written = 0

    def flush_persistent_cache():
        nonlocal persistent_cache_new, persistent_cache_written
        if persistent_cache_path is None or not persistent_cache_new:
            return
        persistent_cache_path.parent.mkdir(parents=True, exist_ok=True)
        with persistent_cache_path.open("a") as f:
            for item in persistent_cache_new:
                f.write(json.dumps(item, sort_keys=True) + "\n")
        persistent_cache_written += len(persistent_cache_new)
        persistent_cache_new = []

    def cached_to_station(cur, station, cap_s, cache):
        nonlocal persistent_cache_hits, persistent_cache_misses
        key = (cur, station, int(cap_s))
        if key not in cache:
            nonlocal_key = (
                int(cur),
                str(station),
                int(cap_s),
                run_mode,
                int(round(float(args.run_radius))),
                int(args.connector_cap),
            )
            if persistent_cache_path is not None and nonlocal_key in persistent_cache:
                persistent_cache_hits += 1
                cache[key] = persistent_cache[nonlocal_key]
                return cache[key]
            if persistent_cache_path is not None:
                persistent_cache_misses += 1
            hit, prev = _dijkstra_to_station_bounded(
                G._adj,
                node_off,
                cur,
                station,
                args.connector_cap,
                runs=runs,
                max_cost=cap_s,
            )
            if hit is None:
                cache[key] = None
                if persistent_cache_path is not None:
                    persistent_cache[nonlocal_key] = None
                    persistent_cache_new.append({
                        "key": list(nonlocal_key),
                        "reachable": False,
                    })
            else:
                seg = _path_from_prev(prev, cur, hit)
                cache[key] = (hit, seg)
                if persistent_cache_path is not None:
                    persistent_cache[nonlocal_key] = (int(hit), [int(node) for node in seg])
                    persistent_cache_new.append({
                        "key": list(nonlocal_key),
                        "reachable": True,
                        "hit": int(hit),
                        "seg": [int(node) for node in seg],
                    })
            if len(persistent_cache_new) >= 1000:
                flush_persistent_cache()
        return cache[key]

    def resource_hits_for(covered):
        return tuple(frozenset(covered & group) for group in resource_sets)

    def resource_mask_for(resource_hits):
        mask = 0
        for idx, hits in enumerate(resource_hits):
            if len(hits) >= resource_min_hits[idx]:
                mask |= 1 << idx
        return mask

    def score_state(state):
        resource_hit_total = sum(len(hits) for hits in state["resource_hits"])
        return (
            int(state["elapsed_s"])
            + int(args.column_penalty)
            - float(args.target_reward_s) * len(state["hits"])
            - float(args.cover_reward_s) * len(state["covered"])
            - float(args.resource_reward_s) * int(state["resource_mask"].bit_count())
            - float(args.resource_hit_reward_s) * resource_hit_total
        )

    def extend_to_station(state, station, cache):
        cur = state["current"]
        remaining_elapsed = max_elapsed_s - int(state["elapsed_s"])
        if remaining_elapsed < 0:
            return None
        if max_end_time is None:
            remaining_time = remaining_elapsed
        else:
            remaining_time = max_end_time - int(node_t[cur])
            if remaining_time < 0:
                return None
        cap_s = min(max_leg_s, remaining_elapsed, remaining_time)
        if node_off[cur] == station:
            hit = cur
            seg = []
        else:
            cached = cached_to_station(cur, station, cap_s, cache)
            if cached is None:
                return None
            hit, seg = cached
            if seg is None:
                return None
        if max_end_time is not None and int(node_t[hit]) > max_end_time:
            return None
        nodes = state["nodes"] + tuple(seg)
        add_elapsed = (
            sum(weight for weight, _mode in _transition_series(prob, [cur] + list(seg)))
            if seg else 0
        )
        elapsed = int(state["elapsed_s"]) + int(add_elapsed)
        if elapsed > max_elapsed_s:
            return None
        covered = set(state["covered"])
        covered.update(node_off[n] for n in seg if node_off[n] in si.canonical_stations)
        covered.add(node_off[hit])
        resource_hits = resource_hits_for(covered)
        resource_mask = resource_mask_for(resource_hits)
        return {
            "nodes": nodes,
            "current": hit,
            "covered": frozenset(covered),
            "hits": frozenset(covered & targets),
            "resource_hits": resource_hits,
            "resource_mask": resource_mask,
            "anchors": state["anchors"] + (station,),
            "elapsed_s": elapsed,
        }

    def expansion_targets(state):
        candidates = []
        for station in target_order:
            if station in state["hits"]:
                continue
            resource_need = 0
            for idx, group in enumerate(resource_sets):
                if station not in group or station in state["covered"]:
                    continue
                if state["resource_mask"] & (1 << idx):
                    continue
                resource_need += 1
            candidates.append((
                -resource_need,
                target_order_index.get(station, len(target_order_index)),
                station,
            ))
        candidates.sort()
        stations = [station for _priority, _order, station in candidates]
        if max_expand_targets:
            stations = stations[:max_expand_targets]
        return stations

    initial_covered = {node_off[src]} & si.canonical_stations
    initial_resource_hits = resource_hits_for(initial_covered)
    beam = [{
        "nodes": (src,),
        "current": src,
        "covered": frozenset(initial_covered),
        "hits": frozenset(initial_covered & targets),
        "resource_hits": initial_resource_hits,
        "resource_mask": resource_mask_for(initial_resource_hits),
        "anchors": (),
        "elapsed_s": 0,
    }]
    best_by_key = {}
    final_candidates = []
    seen_final = set()
    frontier_pool = {}

    def record_frontier_state(state):
        if not args.emit_frontier_top:
            return
        if (state["resource_mask"] & required_resource_mask) != required_resource_mask:
            return
        if int(state["resource_mask"].bit_count()) < frontier_min_resource_count:
            return
        if len(state["hits"]) < frontier_min_target_hits:
            return
        if len(state["covered"]) < frontier_min_covered:
            return
        nodes = state["nodes"]
        key = (
            int(node_t[nodes[0]]),
            node_off[nodes[0]],
            int(node_t[nodes[-1]]),
            node_off[nodes[-1]],
            tuple(sorted(state["covered"])),
            tuple(sorted(state["hits"])),
            state["resource_mask"],
        )
        candidate = {
            **state,
            "final_station": None,
            "score_s": float(score_state(state)),
            "resource_chain_frontier": True,
        }
        old = frontier_pool.get(key)
        if old is None or candidate["score_s"] < old["score_s"]:
            frontier_pool[key] = candidate

    def add_final_candidate(state, final_station):
        final_state = extend_to_station(state, final_station, final_cache)
        if final_state is None:
            return
        if (final_state["resource_mask"] & required_resource_mask) != required_resource_mask:
            return
        if int(final_state["resource_mask"].bit_count()) < min_resource_count:
            return
        if len(final_state["hits"]) < args.min_target_hits:
            return
        if len(final_state["covered"]) < args.min_covered:
            return
        score = score_state(final_state)
        if max_score < INF and score > max_score:
            return
        nodes = final_state["nodes"]
        key = (
            int(node_t[nodes[0]]),
            node_off[nodes[0]],
            int(node_t[nodes[-1]]),
            node_off[nodes[-1]],
            tuple(sorted(final_state["covered"])),
            tuple(sorted(final_state["hits"])),
            final_state["resource_mask"],
        )
        if key in seen_final:
            return
        seen_final.add(key)
        final_candidates.append({
            **final_state,
            "final_station": final_station,
            "score_s": float(score),
        })

    for depth in range(args.max_depth + 1):
        for state in beam:
            record_frontier_state(state)
            if not args.frontier_only:
                for final_station in finals:
                    add_final_candidate(state, final_station)
        if depth == args.max_depth:
            break

        expansion_jobs = []
        for state_index, state in enumerate(beam):
            for station in expansion_targets(state):
                expansion_jobs.append((
                    score_state(state),
                    state["elapsed_s"],
                    int(node_t[state["current"]]),
                    state_index,
                    station,
                    state,
                ))
        expansion_jobs.sort(key=lambda item: item[:5])
        skipped_expansions = 0
        if args.max_expansions_per_depth and len(expansion_jobs) > args.max_expansions_per_depth:
            skipped_expansions = len(expansion_jobs) - int(args.max_expansions_per_depth)
            expansion_jobs = expansion_jobs[:int(args.max_expansions_per_depth)]

        next_states = []
        for _score, _elapsed, _time, _state_index, station, state in expansion_jobs:
            extended = extend_to_station(state, station, leg_cache)
            if extended is None:
                continue
            if (
                len(extended["hits"]) <= len(state["hits"])
                and extended["resource_mask"] == state["resource_mask"]
            ):
                continue
            key = (
                extended["current"],
                extended["resource_mask"],
                extended["hits"],
                int(extended["elapsed_s"] // time_bucket_s),
            )
            old = best_by_key.get(key)
            if old is not None and score_state(old) <= score_state(extended):
                continue
            best_by_key[key] = extended
            next_states.append(extended)

        dedup = {}
        for state in next_states:
            key = (
                state["current"],
                state["resource_mask"],
                state["hits"],
                int(state["elapsed_s"] // time_bucket_s),
            )
            old = dedup.get(key)
            if old is None or score_state(state) < score_state(old):
                dedup[key] = state
        beam = sorted(
            dedup.values(),
            key=lambda state: (
                score_state(state),
                -int(state["resource_mask"].bit_count()),
                -len(state["hits"]),
                -len(state["covered"]),
                state["elapsed_s"],
                int(node_t[state["current"]]),
                state["anchors"],
            ),
        )[:args.beam_size]
        mask_counts = Counter(state["resource_mask"] for state in beam)
        print(f"resource-chain depth={depth + 1} states={len(beam)} "
              f"finals={len(final_candidates)} masks={len(mask_counts)} "
              f"expansions={len(expansion_jobs)} skipped={skipped_expansions} "
              f"leg_cache={len(leg_cache)} wall={time.time() - t0:.1f}s",
              flush=True)
        flush_persistent_cache()
        if not beam:
            break

    final_candidates.sort(
        key=lambda state: (
            state["score_s"],
            -int(state["resource_mask"].bit_count()),
            -len(state["hits"]),
            -len(state["covered"]),
            state["elapsed_s"],
            int(node_t[state["nodes"][-1]]),
            state["anchors"],
        )
    )
    if args.emit_top:
        final_candidates = final_candidates[:args.emit_top]
    if args.max_generated:
        final_candidates = final_candidates[:args.max_generated]

    frontier_candidates = []
    if args.emit_frontier_top:
        frontier_candidates = list(frontier_pool.values())
        frontier_candidates.sort(
            key=lambda state: (
                state["score_s"],
                -int(state["resource_mask"].bit_count()),
                -len(state["hits"]),
                -len(state["covered"]),
                state["elapsed_s"],
                int(node_t[state["nodes"][-1]]),
                state["anchors"],
            )
        )
        frontier_candidates = frontier_candidates[:int(args.emit_frontier_top)]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    generated = []
    seen_base = set()
    for row in base_rows:
        seen_base.add((
            row["start"]["time"],
            row["start"]["station"],
            row["end"]["time"],
            row["end"]["station"],
            tuple(row.get("covered_stations", [])),
        ))

    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")
        for state in list(final_candidates) + frontier_candidates:
            nodes = state["nodes"]
            key = (
                int(node_t[nodes[0]]),
                node_off[nodes[0]],
                int(node_t[nodes[-1]]),
                node_off[nodes[-1]],
                tuple(sorted(state["covered"])),
            )
            if key in seen_base:
                continue
            seen_base.add(key)
            record = _column_from_nodes(
                prob, args.label, len(base_rows) + len(generated) + 1, nodes)
            record["pricing_resource_chain"] = {
                "source_columns": args.source_columns,
                "source_column": args.source_column_id,
                "source_point": args.source_point,
                "resources": [
                    {
                        "name": resource["name"],
                        "stations": resource["stations"],
                        "min_hits": resource_min_hits[idx],
                        "hits": sorted(state["resource_hits"][idx]),
                    }
                    for idx, resource in enumerate(resources)
                ],
                "resource_mask": int(state["resource_mask"]),
                "resource_count": int(state["resource_mask"].bit_count()),
                "required_resources": required_resource_names,
                "anchors": list(state["anchors"]),
                "final_station": state["final_station"],
                "frontier": bool(state.get("resource_chain_frontier", False)),
                "target_stations": sorted(targets),
                "target_hits": sorted(state["hits"]),
                "target_hit_count": len(state["hits"]),
                "score_s": state["score_s"],
                "max_end_time": max_end_time,
                "run_mode": run_mode,
            }
            generated.append(record)
            f.write(json.dumps(record, sort_keys=True) + "\n")

    flush_persistent_cache()

    print(f"resource-chain source={args.source_column_id} resources={len(resources)} "
          f"final_candidates={len(final_candidates)} "
          f"frontier_candidates={len(frontier_candidates)} "
          f"generated={len(generated)} "
          f"leg_cache={len(leg_cache)} final_cache={len(final_cache)} "
          f"persistent_cache_hits={persistent_cache_hits} "
          f"persistent_cache_misses={persistent_cache_misses} "
          f"persistent_cache_new={len(persistent_cache_new)} "
          f"persistent_cache_written={persistent_cache_written} "
          f"wall={time.time() - t0:.1f}s")
    if persistent_cache_path is not None:
        print(f"resource-chain cache path={persistent_cache_path}")
    print(f"wrote {len(base_rows)} base + {len(generated)} resource-chain columns -> {out}")
    return 0


def cmd_price_stage_chains(args) -> int:
    import pickle

    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    source_rows = _load_column_rows(Path(args.source_columns))

    si = StationIndex.load()
    stage_groups = _load_stage_groups(args, si)
    stage_min_hits = _parse_stage_min_hits(args.stage_min_hits, len(stage_groups))
    staged_targets = {station for group in stage_groups for station in group}
    targets = set(staged_targets)
    if args.target_stations or args.targets_file:
        targets.update(_load_target_stations(args, si))
    finals = _parse_station_list(args.final_stations, si, "--final-stations")
    if not finals:
        raise SystemExit("--final-stations must include at least one canonical station")

    source_row = next(
        (row for row in source_rows if row.get("column_id") == args.source_column_id),
        None,
    )
    if source_row is None:
        raise SystemExit(f"source column not found: {args.source_column_id}")

    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, runs or _NoRunLayer(), si)

    source_start, source_end = _row_endpoint_nodes_from_tables(
        node_stop, node_t, source_row)
    src = source_start if args.source_point == "start" else source_end
    max_end_time = int(args.max_end_time) if args.max_end_time is not None else None
    if max_end_time is not None and int(node_t[src]) > max_end_time:
        raise SystemExit("source starts after --max-end-time")
    max_elapsed_s = int(args.max_elapsed_minutes * 60)
    max_leg_s = int(args.max_leg_minutes * 60) if args.max_leg_minutes else max_elapsed_s
    max_score = INF if args.max_score_s is None else float(args.max_score_s)
    time_bucket_s = max(1, int(args.time_bucket_s))

    leg_cache = {}

    def cached_to_station(cur, station, cap_s):
        key = (cur, station, int(cap_s))
        if key not in leg_cache:
            hit, prev = _dijkstra_to_station_bounded(
                G._adj,
                node_off,
                cur,
                station,
                args.connector_cap,
                runs=runs,
                max_cost=cap_s,
            )
            if hit is None:
                leg_cache[key] = None
            else:
                leg_cache[key] = (hit, _path_from_prev(prev, cur, hit))
        return leg_cache[key]

    def score_state(state):
        return (
            state["elapsed_s"]
            + args.column_penalty
            - args.target_reward_s * len(state["hits"])
            - args.cover_reward_s * len(state["covered"])
            - args.stage_reward_s * len(state["stage_hits"])
        )

    def extend_to_station(state, station, stage_index, anchor_prefix, record_stage=True):
        cur = state["current"]
        remaining_elapsed = max_elapsed_s - int(state["elapsed_s"])
        if remaining_elapsed < 0:
            return None
        if max_end_time is None:
            remaining_time = remaining_elapsed
        else:
            remaining_time = max_end_time - int(node_t[cur])
            if remaining_time < 0:
                return None
        cap_s = min(max_leg_s, remaining_elapsed, remaining_time)
        if node_off[cur] == station:
            hit = cur
            seg = ()
        else:
            cached = cached_to_station(cur, station, cap_s)
            if cached is None:
                return None
            hit, seg = cached
            if seg is None:
                return None
        if max_end_time is not None and int(node_t[hit]) > max_end_time:
            return None
        add_elapsed = _elapsed(node_t[cur], node_t[hit])
        elapsed = int(state["elapsed_s"]) + int(add_elapsed)
        if elapsed > max_elapsed_s:
            return None
        covered = set(state["covered"])
        covered.update(node_off[n] for n in seg if node_off[n] in si.canonical_stations)
        covered.add(node_off[hit])
        hits = frozenset(covered & targets)
        stage_hits = (
            tuple(list(state["stage_hits"]) + [station])
            if record_stage else state["stage_hits"]
        )
        return {
            "nodes": state["nodes"] + tuple(seg),
            "current": hit,
            "covered": frozenset(covered),
            "hits": hits,
            "stage_hits": stage_hits,
            "anchors": state["anchors"] + (f"{anchor_prefix}:{station}",),
            "elapsed_s": int(elapsed),
        }

    initial_covered = {node_off[src]} & si.canonical_stations
    beam = [{
        "nodes": (src,),
        "current": src,
        "covered": frozenset(initial_covered),
        "hits": frozenset(initial_covered & targets),
        "stage_hits": (),
        "anchors": (),
        "elapsed_s": 0,
    }]
    t0 = time.time()
    frontier_pool = {}

    def record_frontier_state(state, stage_index):
        if not args.emit_frontier_top:
            return
        if len(state["hits"]) < args.frontier_min_target_hits:
            return
        if len(state["covered"]) < args.frontier_min_covered:
            return
        nodes = state["nodes"]
        key = (
            int(node_t[nodes[0]]),
            node_off[nodes[0]],
            int(node_t[nodes[-1]]),
            node_off[nodes[-1]],
            tuple(sorted(state["covered"])),
            tuple(sorted(state["hits"])),
            state["anchors"],
            int(stage_index),
        )
        candidate = {
            **state,
            "final_station": None,
            "score_s": score_state(state),
            "stage_chain_frontier": True,
            "frontier_stage_index": int(stage_index),
        }
        old = frontier_pool.get(key)
        if old is None or candidate["score_s"] < old["score_s"]:
            frontier_pool[key] = candidate

    def mark_stage_satisfied(state, stage_index, stage_set, required_hits):
        if required_hits <= 0:
            chosen = []
        else:
            chosen = sorted(state["covered"] & stage_set)[:required_hits]
        old_hits = list(state["stage_hits"])
        old_set = set(old_hits)
        additions = [station for station in chosen if station not in old_set]
        marker = ",".join(chosen) if chosen else "none"
        return {
            **state,
            "stage_hits": tuple(old_hits + additions),
            "anchors": state["anchors"] + (f"stage{stage_index}:covered:{marker}",),
        }

    for stage_index, group in enumerate(stage_groups, start=1):
        stage_set = set(group)
        required_hits = min(int(stage_min_hits[stage_index - 1]), len(stage_set))
        stage_targets = group[:args.max_stage_targets] if args.max_stage_targets else group
        max_stage_depth = (
            int(args.max_stage_depth)
            if args.max_stage_depth
            else max(1, required_hits)
        )
        stage_frontier = beam
        satisfied = []
        for depth in range(max_stage_depth + 1):
            next_frontier = []
            for state in stage_frontier:
                covered_in_stage = state["covered"] & stage_set
                if (
                    len(covered_in_stage) >= required_hits
                    and (not args.force_stage_move or depth > 0)
                ):
                    satisfied.append(
                        mark_stage_satisfied(
                            state,
                            stage_index,
                            stage_set,
                            required_hits,
                        )
                    )
                if depth >= max_stage_depth:
                    continue
                for station in stage_targets:
                    if station in covered_in_stage:
                        continue
                    extended = extend_to_station(
                        state, station, stage_index, f"stage{stage_index}")
                    if extended is not None:
                        next_frontier.append(extended)
            dedup = {}
            for state in next_frontier:
                key = (
                    state["current"],
                    state["hits"],
                    int(state["elapsed_s"] // time_bucket_s),
                    len(state["covered"] & stage_set),
                )
                old = dedup.get(key)
                if old is None or score_state(state) < score_state(old):
                    dedup[key] = state
            stage_frontier = sorted(
                dedup.values(),
                key=lambda state: (
                    score_state(state),
                    -len(state["hits"]),
                    -len(state["covered"] & stage_set),
                    state["elapsed_s"],
                    int(node_t[state["current"]]),
                    state["anchors"],
                ),
            )[:args.beam_size]
            if not stage_frontier:
                break

        dedup = {}
        for state in satisfied:
            key = (
                state["current"],
                state["hits"],
                int(state["elapsed_s"] // time_bucket_s),
                tuple(sorted(state["covered"] & stage_set)),
            )
            old = dedup.get(key)
            if old is None or score_state(state) < score_state(old):
                dedup[key] = state
        beam = sorted(
            dedup.values(),
            key=lambda state: (
                score_state(state),
                -len(state["hits"]),
                -len(state["covered"] & stage_set),
                state["elapsed_s"],
                int(node_t[state["current"]]),
                state["anchors"],
            ),
        )[:args.beam_size]
        print(f"stage-chain stage={stage_index}/{len(stage_groups)} "
              f"min_hits={required_hits} states={len(beam)} "
              f"leg_cache={len(leg_cache)} "
              f"wall={time.time() - t0:.1f}s", flush=True)
        for state in beam:
            record_frontier_state(state, stage_index)
        if not beam:
            break

    final_candidates = []
    seen_final = set()
    if not args.frontier_only:
        for state in beam:
            for final_station in finals:
                final_state = extend_to_station(
                    state,
                    final_station,
                    len(stage_groups) + 1,
                    "final",
                    record_stage=False,
                )
                if final_state is None:
                    continue
                if len(final_state["hits"]) < args.min_target_hits:
                    continue
                if len(final_state["covered"]) < args.min_covered:
                    continue
                score = score_state(final_state)
                if max_score < INF and score > max_score:
                    continue
                key = (
                    int(node_t[final_state["nodes"][0]]),
                    node_off[final_state["nodes"][0]],
                    int(node_t[final_state["nodes"][-1]]),
                    node_off[final_state["nodes"][-1]],
                    tuple(sorted(final_state["covered"])),
                    tuple(sorted(final_state["hits"])),
                    final_state["anchors"],
                )
                if key in seen_final:
                    continue
                seen_final.add(key)
                final_candidates.append({
                    **final_state,
                    "final_station": final_station,
                    "score_s": score,
                    "stage_chain_frontier": False,
                    "frontier_stage_index": None,
                })

    frontier_candidates = []
    if args.emit_frontier_top:
        frontier_candidates = list(frontier_pool.values())
        frontier_candidates.sort(
            key=lambda state: (
                state["score_s"],
                -len(state["hits"]),
                -len(state["covered"]),
                state["elapsed_s"],
                int(node_t[state["nodes"][-1]]),
                state["anchors"],
            )
        )
        frontier_candidates = frontier_candidates[:int(args.emit_frontier_top)]

    final_candidates.sort(
        key=lambda state: (
            state["score_s"],
            -len(state["hits"]),
            -len(state["covered"]),
            state["elapsed_s"],
            int(node_t[state["nodes"][-1]]),
            state["anchors"],
        )
    )
    if args.emit_top:
        final_candidates = final_candidates[:args.emit_top]
    if args.max_generated:
        final_candidates = final_candidates[:args.max_generated]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    generated = []
    seen_base = set()
    for row in base_rows:
        seen_base.add((
            row["start"]["time"],
            row["start"]["station"],
            row["end"]["time"],
            row["end"]["station"],
            tuple(row.get("covered_stations", [])),
        ))

    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")
        for state in list(final_candidates) + frontier_candidates:
            nodes = state["nodes"]
            key = (
                int(node_t[nodes[0]]),
                node_off[nodes[0]],
                int(node_t[nodes[-1]]),
                node_off[nodes[-1]],
                tuple(sorted(state["covered"])),
            )
            if key in seen_base:
                continue
            seen_base.add(key)
            record = _column_from_nodes(
                prob, args.label, len(base_rows) + len(generated) + 1, nodes)
            record["pricing_stage_chain"] = {
                "source_columns": args.source_columns,
                "source_column": args.source_column_id,
                "source_point": args.source_point,
                "stage_groups": stage_groups,
                "stage_min_hits": stage_min_hits,
                "stage_hits": list(state["stage_hits"]),
                "anchors": list(state["anchors"]),
                "final_station": state["final_station"],
                "frontier": bool(state.get("stage_chain_frontier", False)),
                "frontier_stage_index": state.get("frontier_stage_index"),
                "target_stations": sorted(targets),
                "target_hits": sorted(state["hits"]),
                "target_hit_count": len(state["hits"]),
                "score_s": state["score_s"],
                "max_end_time": max_end_time,
                "run_mode": run_mode,
            }
            generated.append(record)
            f.write(json.dumps(record, sort_keys=True) + "\n")

    print(f"stage-chain source={args.source_column_id} stages={len(stage_groups)} "
          f"final_candidates={len(final_candidates)} "
          f"frontier_candidates={len(frontier_candidates)} "
          f"generated={len(generated)} "
          f"leg_cache={len(leg_cache)} wall={time.time() - t0:.1f}s")
    print(f"wrote {len(base_rows)} base + {len(generated)} stage-chain columns -> {out}")
    return 0


def _state_score(state, target_reward_s, cover_reward_s):
    return (
        target_reward_s * len(state["hits"])
        + cover_reward_s * len(state["covered"])
        - state["elapsed_s"]
    )


def cmd_price_late_tail_beam(args) -> int:
    import pickle

    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    source_rows = _load_column_rows(Path(args.source_columns))

    si = StationIndex.load()
    targets = _load_target_stations(args, si)
    target_order = _parse_station_list(args.target_stations, si, "--target-stations")
    for station in sorted(targets):
        if station not in target_order:
            target_order.append(station)
    finals = _parse_station_list(args.final_stations, si, "--final-stations")
    if not finals:
        raise SystemExit("--final-stations must include at least one canonical station")

    source_row = next(
        (row for row in source_rows if row.get("column_id") == args.source_column_id),
        None,
    )
    if source_row is None:
        raise SystemExit(f"source column not found: {args.source_column_id}")

    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, runs or _NoRunLayer(), si)

    source_start, source_end = _row_endpoint_nodes_from_tables(
        node_stop, node_t, source_row)
    src = source_start if args.source_point == "start" else source_end
    if node_t[src] > args.max_end_time:
        raise SystemExit("source starts after --max-end-time")

    max_elapsed_s = int(args.max_elapsed_minutes * 60)
    max_leg_s = int(args.max_leg_minutes * 60) if args.max_leg_minutes else max_elapsed_s
    max_score = INF if args.max_score_s is None else float(args.max_score_s)
    initial_covered = {node_off[src]} & si.canonical_stations
    initial_hits = initial_covered & targets
    beam = [{
        "nodes": (src,),
        "current": src,
        "covered": frozenset(initial_covered),
        "hits": frozenset(initial_hits),
        "anchors": (),
        "elapsed_s": 0,
    }]
    best_by_key = {}
    leg_cache = {}
    final_cache = {}
    final_candidates = []
    seen_final = set()
    t0 = time.time()

    def cached_to_station(cur, station, cap_s, cache):
        key = (cur, station, cap_s)
        if key not in cache:
            hit, prev = _dijkstra_to_station_bounded(
                G._adj,
                node_off,
                cur,
                station,
                args.connector_cap,
                runs=runs,
                max_cost=cap_s,
            )
            if hit is None:
                cache[key] = None
            else:
                cache[key] = (hit, _path_from_prev(prev, cur, hit))
        return cache[key]

    def add_final_candidate(state, final_station):
        cur = state["current"]
        remaining_time = int(args.max_end_time) - int(node_t[cur])
        if remaining_time < 0:
            return
        if node_off[cur] == final_station:
            hit = cur
            seg = []
        else:
            cached = cached_to_station(
                cur,
                final_station,
                min(max_leg_s, remaining_time, max_elapsed_s),
                final_cache,
            )
            if cached is None:
                return
            hit, seg = cached
            if seg is None:
                return
        nodes = state["nodes"] + tuple(seg)
        if int(node_t[nodes[-1]]) > int(args.max_end_time):
            return
        elapsed = sum(weight for weight, _mode in _transition_series(prob, nodes))
        if elapsed <= 0 or elapsed > max_elapsed_s:
            return
        covered = set(state["covered"])
        covered.update(node_off[n] for n in seg if node_off[n] in si.canonical_stations)
        hits = frozenset(covered & targets)
        if len(hits) < args.min_target_hits or len(covered) < args.min_covered:
            return
        score = elapsed + args.column_penalty - (
            args.target_reward_s * len(hits)
            + args.cover_reward_s * len(covered)
        )
        if max_score < INF and score > max_score:
            return
        key = (
            node_t[nodes[0]],
            node_off[nodes[0]],
            node_t[nodes[-1]],
            node_off[nodes[-1]],
            tuple(sorted(covered)),
            tuple(sorted(hits)),
        )
        if key in seen_final:
            return
        seen_final.add(key)
        final_candidates.append({
            "nodes": nodes,
            "covered": frozenset(covered),
            "hits": hits,
            "anchors": state["anchors"] + (f"final:{final_station}",),
            "elapsed_s": int(elapsed),
            "score_s": float(score),
            "final_station": final_station,
        })

    for depth in range(args.max_depth + 1):
        for state in beam:
            for final_station in finals:
                add_final_candidate(state, final_station)
        if depth == args.max_depth:
            break

        next_states = []
        for state in beam:
            cur = state["current"]
            remaining_time = int(args.max_end_time) - int(node_t[cur])
            if remaining_time <= 0:
                continue
            available = [station for station in target_order if station not in state["hits"]]
            if args.max_expand_targets:
                available = available[:args.max_expand_targets]
            for target in available:
                cached = cached_to_station(
                    cur,
                    target,
                    min(max_leg_s, remaining_time, max_elapsed_s),
                    leg_cache,
                )
                if cached is None:
                    continue
                hit, seg = cached
                if seg is None or not seg:
                    continue
                if int(node_t[hit]) > int(args.max_end_time):
                    continue
                nodes = state["nodes"] + tuple(seg)
                elapsed = sum(weight for weight, _mode in _transition_series(prob, nodes))
                if elapsed <= 0 or elapsed > max_elapsed_s:
                    continue
                covered = set(state["covered"])
                covered.update(node_off[n] for n in seg if node_off[n] in si.canonical_stations)
                hits = frozenset(covered & targets)
                if len(hits) <= len(state["hits"]):
                    continue
                new_state = {
                    "nodes": nodes,
                    "current": hit,
                    "covered": frozenset(covered),
                    "hits": hits,
                    "anchors": state["anchors"] + (target,),
                    "elapsed_s": int(elapsed),
                }
                dkey = (hit, hits)
                old = best_by_key.get(dkey)
                new_score = _state_score(new_state, args.target_reward_s, args.cover_reward_s)
                if old is not None:
                    old_score = _state_score(old, args.target_reward_s, args.cover_reward_s)
                    if (old["elapsed_s"], -old_score) <= (new_state["elapsed_s"], -new_score):
                        continue
                best_by_key[dkey] = new_state
                next_states.append(new_state)
        next_states.sort(
            key=lambda state: (
                -_state_score(state, args.target_reward_s, args.cover_reward_s),
                state["elapsed_s"],
                int(node_t[state["current"]]),
                state["anchors"],
            )
        )
        beam = next_states[:args.beam_size]
        print(f"late-tail depth={depth + 1} states={len(beam)} "
              f"finals={len(final_candidates)} leg_cache={len(leg_cache)} "
              f"wall={time.time() - t0:.1f}s", flush=True)
        if not beam:
            break

    final_candidates.sort(
        key=lambda item: (
            item["score_s"],
            -len(item["hits"]),
            item["elapsed_s"],
            int(node_t[item["nodes"][-1]]),
        )
    )
    if args.emit_top:
        final_candidates = final_candidates[:args.emit_top]
    if args.max_generated:
        final_candidates = final_candidates[:args.max_generated]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    generated = []
    seen_base = set()
    for row in base_rows:
        seen_base.add((
            row["start"]["time"],
            row["start"]["station"],
            row["end"]["time"],
            row["end"]["station"],
            tuple(row.get("covered_stations", [])),
        ))
    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")
        for item in final_candidates:
            nodes = item["nodes"]
            covered = item["covered"]
            key = (
                node_t[nodes[0]],
                node_off[nodes[0]],
                node_t[nodes[-1]],
                node_off[nodes[-1]],
                tuple(sorted(covered)),
            )
            if key in seen_base:
                continue
            seen_base.add(key)
            record = _column_from_nodes(
                prob, args.label, len(base_rows) + len(generated) + 1, nodes)
            record["pricing_late_tail_beam"] = {
                "source_columns": args.source_columns,
                "source_column": args.source_column_id,
                "source_point": args.source_point,
                "anchors": list(item["anchors"]),
                "final_station": item["final_station"],
                "target_stations": sorted(targets),
                "target_hits": sorted(item["hits"]),
                "target_hit_count": len(item["hits"]),
                "reward_s": (
                    args.target_reward_s * len(item["hits"])
                    + args.cover_reward_s * len(covered)
                ),
                "score_s": item["score_s"],
                "max_end_time": int(args.max_end_time),
                "run_mode": run_mode,
            }
            generated.append(record)
            f.write(json.dumps(record, sort_keys=True) + "\n")

    print(f"late-tail beam generated={len(generated)} finals_seen={len(final_candidates)} "
          f"beam_size={args.beam_size} max_depth={args.max_depth}")
    print(f"wrote {len(base_rows)} base + {len(generated)} late-tail columns -> {out}")
    return 0


def cmd_split_columns(args) -> int:
    import pickle

    si = StationIndex.load()
    split_stations = _parse_station_list(args.split_stations, si, "--split-stations")
    if not split_stations:
        raise SystemExit("--split-stations must include at least one canonical station")
    split_set = set(split_stations)
    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    column_id_prefixes = [
        prefix.strip()
        for prefix in str(args.column_id_prefix or "").split(",")
        if prefix.strip()
    ]
    source_rows = []
    for source_file in args.source_columns:
        for row in _load_column_rows(Path(source_file)):
            if args.pricing_kind and args.pricing_kind not in row:
                continue
            if column_id_prefixes and not any(
                str(row.get("column_id", "")).startswith(prefix)
                for prefix in column_id_prefixes
            ):
                continue
            source_rows.append((source_file, row))

    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off = tables[0]
    stop_events = _stop_index(tables[1], tables[2])
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)

    seen = set()
    for row in base_rows:
        seen.add((
            row["start"]["time"],
            row["start"]["station"],
            row["end"]["time"],
            row["end"]["station"],
            tuple(row.get("covered_stations", [])),
        ))

    generated = []
    skipped = Counter()
    max_elapsed_s = int(args.max_elapsed_minutes * 60) if args.max_elapsed_minutes else 0
    for source_file, row in source_rows:
        try:
            nodes = _row_path_nodes(stop_events, tables[2], row)
        except ValueError as exc:
            skipped["missing_path_event"] += 1
            if args.verbose:
                print(f"skip {row.get('column_id')}: {exc}", flush=True)
            continue
        split_points = [0]
        for idx, node in enumerate(nodes[1:-1], start=1):
            if node_off[node] in split_set:
                split_points.append(idx)
        split_points.append(len(nodes) - 1)
        split_points = sorted(set(split_points))
        if len(split_points) < 2:
            skipped["no_split_points"] += 1
            continue

        for left_pos, left in enumerate(split_points[:-1]):
            right_limit = len(split_points)
            if args.max_split_gap:
                right_limit = min(right_limit, left_pos + args.max_split_gap + 1)
            for right in split_points[left_pos + 1:right_limit]:
                if right <= left:
                    continue
                seg_nodes = nodes[left:right + 1]
                if len(seg_nodes) < args.min_steps:
                    skipped["too_few_steps"] += 1
                    continue
                transitions = _transition_series(prob, seg_nodes)
                elapsed = sum(weight for weight, _mode in transitions)
                if elapsed <= 0:
                    skipped["zero_elapsed"] += 1
                    continue
                if max_elapsed_s and elapsed > max_elapsed_s:
                    skipped["too_long"] += 1
                    continue
                covered = {node_off[n] for n in seg_nodes} & si.canonical_stations
                if len(covered) < args.min_covered:
                    skipped["too_few_covered"] += 1
                    continue
                if args.min_target_hits:
                    hits = covered & split_set
                    if len(hits) < args.min_target_hits:
                        skipped["too_few_target_hits"] += 1
                        continue
                key = (
                    int(tables[2][seg_nodes[0]]),
                    node_off[seg_nodes[0]],
                    int(tables[2][seg_nodes[-1]]),
                    node_off[seg_nodes[-1]],
                    tuple(sorted(covered)),
                )
                if key in seen:
                    skipped["duplicate"] += 1
                    continue
                seen.add(key)
                record = _column_from_nodes(
                    prob, args.label, len(base_rows) + len(generated) + 1, seg_nodes)
                record["pricing_split"] = {
                    "source_columns": source_file,
                    "source_column": row.get("column_id"),
                    "source_pricing_kind": args.pricing_kind,
                    "source_start": row.get("start"),
                    "source_end": row.get("end"),
                    "source_anchors": (
                        row.get(args.pricing_kind, {}).get("anchors")
                        if args.pricing_kind else None
                    ),
                    "split_start_index": left,
                    "split_end_index": right,
                    "split_start_station": node_off[seg_nodes[0]],
                    "split_end_station": node_off[seg_nodes[-1]],
                    "split_stations": split_stations,
                }
                generated.append(record)
                if args.max_generated and len(generated) >= args.max_generated:
                    break
            if args.max_generated and len(generated) >= args.max_generated:
                break
        if args.max_generated and len(generated) >= args.max_generated:
            break

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")
        for row in generated:
            f.write(json.dumps(row, sort_keys=True) + "\n")
    print(f"split source rows={len(source_rows)} generated={len(generated)} "
          f"skipped={dict(skipped)}")
    print(f"wrote {len(base_rows)} base + {len(generated)} split columns -> {out}")
    return 0


def cmd_price_windows(args) -> int:
    import pickle

    si = StationIndex.load()
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, _node_stop, node_t = tables
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)

    seed_path, _meta = _load_route(Path(args.seed))
    seed_nodes = _nodes_for_path(prob, seed_path)
    widths = _parse_ints(args.widths)
    if not widths:
        raise SystemExit("--widths must contain at least one positive integer")

    windows = []
    for start in range(0, len(seed_nodes) - 2, max(1, args.window_stride)):
        for width in widths:
            end = min(len(seed_nodes) - 1, start + width)
            if end <= start:
                continue
            elapsed = _elapsed(node_t[seed_nodes[start]], node_t[seed_nodes[end]])
            stations = {node_off[n] for n in seed_nodes[start:end + 1]} & si.canonical_stations
            score = elapsed - args.cover_reward_s * len(stations)
            windows.append((score, elapsed, start, end, stations))
    windows.sort(reverse=True)
    if args.max_windows:
        windows = windows[:args.max_windows]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    base_rows = _load_column_rows(Path(args.base_columns)) if args.base_columns else []
    generated = []
    seen = set()
    t0 = time.time()
    max_elapsed_s = int(args.max_elapsed_minutes * 60)

    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")
        for widx, (_score, win_elapsed, start, end, _stations) in enumerate(windows, start=1):
            src = seed_nodes[start]
            target_end = min(len(seed_nodes) - 1, end + args.lookahead)
            target_nodes = seed_nodes[start + 1:target_end + 1]
            targets = []
            seen_targets = set()
            for node in target_nodes:
                station = node_off[node]
                if station in si.canonical_stations and station not in seen_targets:
                    seen_targets.add(station)
                    targets.append(station)
            if args.max_targets_per_window:
                targets = targets[:args.max_targets_per_window]

            for target in targets:
                if node_off[src] == target:
                    continue
                tgt, prev = _dijkstra_to_station(
                    G._adj, node_off, src, target, args.connector_cap, runs=runs)
                if tgt is None:
                    continue
                elapsed = _elapsed(node_t[src], node_t[tgt])
                if elapsed <= 0 or elapsed > max_elapsed_s:
                    continue
                seg = _path_from_prev(prev, src, tgt)
                if seg is None:
                    continue
                nodes = [src] + seg
                covered = {node_off[n] for n in nodes} & si.canonical_stations
                if len(covered) < args.min_covered:
                    continue
                key = (node_t[src], node_off[src], node_t[tgt], node_off[tgt])
                if key in seen:
                    continue
                seen.add(key)
                record = _column_from_nodes(
                    prob, args.label, len(base_rows) + len(generated) + 1, nodes)
                record["pricing_window"] = {
                    "seed": args.seed,
                    "start_index": start,
                    "end_index": end,
                    "elapsed_s": int(win_elapsed),
                    "elapsed": hms(win_elapsed),
                }
                generated.append(record)
                f.write(json.dumps(record, sort_keys=True) + "\n")

            if widx % 10 == 0 or widx == len(windows):
                print(f"priced windows {widx}/{len(windows)} generated={len(generated)} "
                      f"wall={time.time() - t0:.1f}s", flush=True)

    print(f"wrote {len(base_rows)} base + {len(generated)} window-priced columns -> {out}")
    return 0


def _load_proxy_cuts(paths):
    cuts = []
    for path in paths:
        data = json.loads(Path(path).read_text())
        for a, b in data.get("proxy_cuts", []):
            cuts.append((a, b, str(path)))
    return cuts


def cmd_price_cuts(args) -> int:
    import pickle

    if not args.cut_files:
        raise SystemExit("provide at least one --cut-file")
    base_rows = _load_column_rows(Path(args.base_columns))
    by_id = {row["column_id"]: row for row in base_rows}
    cuts = _load_proxy_cuts(args.cut_files)
    if args.max_cuts:
        cuts = cuts[:args.max_cuts]

    si = StationIndex.load()
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    tables = _node_tables(G, si)
    node_off, node_stop, node_t = tables
    stop_events = _stop_index(node_stop, node_t)
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    generated = []
    missing = []
    seen = set()
    max_elapsed_s = int(args.max_elapsed_minutes * 60)
    t0 = time.time()

    with out.open("w") as f:
        for row in base_rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")
        for cidx, (src_id, dst_id, cut_source) in enumerate(cuts, start=1):
            src_row = by_id.get(src_id)
            dst_row = by_id.get(dst_id)
            if src_row is None or dst_row is None:
                missing.append((src_id, dst_id))
                continue
            _src_start, src_end = _row_endpoint_nodes(stop_events, node_t, src_row)
            targets = []
            for station in [dst_row["start"]["station"], dst_row["end"]["station"]]:
                if station not in targets:
                    targets.append(station)
            if args.include_covered:
                for station in dst_row.get("covered_stations", []):
                    if station not in targets:
                        targets.append(station)
                    if args.max_targets_per_cut and len(targets) >= args.max_targets_per_cut:
                        break
            if args.max_targets_per_cut:
                targets = targets[:args.max_targets_per_cut]

            for target in targets:
                if node_off[src_end] == target:
                    continue
                tgt, prev = _dijkstra_to_station(
                    G._adj, node_off, src_end, target, args.connector_cap, runs=runs)
                if tgt is None:
                    continue
                elapsed = _elapsed(node_t[src_end], node_t[tgt])
                if elapsed <= 0 or elapsed > max_elapsed_s:
                    continue
                seg = _path_from_prev(prev, src_end, tgt)
                if seg is None:
                    continue
                nodes = [src_end] + seg
                covered = {node_off[n] for n in nodes} & si.canonical_stations
                if len(covered) < args.min_covered:
                    continue
                key = (src_id, dst_id, node_t[src_end], node_off[src_end],
                       node_t[tgt], node_off[tgt])
                if key in seen:
                    continue
                seen.add(key)
                record = _column_from_nodes(
                    prob, args.label, len(base_rows) + len(generated) + 1, nodes)
                record["pricing_cut"] = {
                    "source": cut_source,
                    "from_column": src_id,
                    "to_column": dst_id,
                }
                generated.append(record)
                f.write(json.dumps(record, sort_keys=True) + "\n")

            if cidx % 10 == 0 or cidx == len(cuts):
                print(f"priced cuts {cidx}/{len(cuts)} generated={len(generated)} "
                      f"wall={time.time() - t0:.1f}s", flush=True)

    if missing:
        print(f"missing {len(missing)} cut rows, e.g. {missing[:5]}")
    print(f"wrote {len(base_rows)} base + {len(generated)} cut-priced columns -> {out}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Extract route columns for optimization.")
    sub = p.add_subparsers(dest="cmd", required=True)

    ex = sub.add_parser("extract", help="Extract schedule-realized columns from route JSONs.")
    ex.add_argument("routes", nargs="*", help="Route JSON files; defaults to solutions/*.json.")
    ex.add_argument("--radius", type=int, default=5000)
    ex.add_argument("--max-new-stations", type=int, default=24)
    ex.add_argument("--max-segment-minutes", type=float, default=120)
    ex.add_argument("--split-on-run", action="store_true")
    ex.add_argument("--no-path", action="store_true", help="Do not store path slices.")
    ex.add_argument("--out", default="reports/optimization_runs/seed_columns.jsonl")
    ex.set_defaults(func=cmd_extract)

    st = sub.add_parser("stats", help="Summarize an extracted columns JSONL file.")
    st.add_argument("file")
    st.set_defaults(func=cmd_stats)

    cv = sub.add_parser("cover", help="Solve a sequencing-free set-cover master.")
    cv.add_argument("file")
    cv.add_argument("--time-limit", type=float, default=20)
    cv.add_argument("--workers", type=int, default=8)
    cv.add_argument("--column-penalty", type=int, default=0,
                    help="Extra seconds per selected column to discourage fragmentation.")
    cv.add_argument("--out", default="reports/optimization_runs/seed_cover_solution.json")
    cv.set_defaults(func=cmd_cover)

    dc = sub.add_parser(
        "dual-cover",
        help="Solve the LP set-cover relaxation and export station dual rewards.")
    dc.add_argument("columns_file")
    dc.add_argument("--solver", default="GLOP")
    dc.add_argument("--time-limit", type=float, default=60)
    dc.add_argument("--column-penalty", type=int, default=0,
                    help="Extra seconds per selected column in the LP objective.")
    dc.add_argument("--unit-upper-bound", action="store_true",
                    help="Use x<=1 bounds for the set-cover relaxation.")
    dc.add_argument("--min-x", type=float, default=1e-6,
                    help="Minimum LP value for reporting active columns.")
    dc.add_argument("--top-duals", type=int, default=80)
    dc.add_argument("--max-active-columns", type=int, default=200)
    dc.add_argument("--out", default="reports/optimization_runs/column_duals.json")
    dc.set_defaults(func=cmd_dual_cover)

    ap = sub.add_parser(
        "active-pool",
        help="Select a compact, coverage-complete column pool for exact connector pricing.")
    ap.add_argument("columns_files", nargs="+",
                    help="Column JSONL or selected-columns JSON files to merge.")
    ap.add_argument("--duals", default=None,
                    help="Optional dual-cover JSON used for reduced-cost ranking.")
    ap.add_argument("--include-solution", action="append",
                    help="Selected-columns JSON whose column ids should be forced in.")
    ap.add_argument("--top-columns", type=int, default=800,
                    help="Target active-pool size before coverage repair.")
    ap.add_argument("--top-priced", type=int, default=0,
                    help="Also force in this many rows with explicit pricing reduced costs.")
    ap.add_argument("--column-penalty", type=int, default=0)
    ap.add_argument("--no-ensure-coverage", action="store_true",
                    help="Do not greedily add columns for missing stations.")
    ap.add_argument("--out", default="reports/optimization_runs/active_columns.jsonl")
    ap.set_defaults(func=cmd_active_pool)

    seq = sub.add_parser("sequence", help="Sequence selected columns and realize a route.")
    seq.add_argument("cover_solution", help="JSON output from the cover command.")
    seq.add_argument("--terminal-runs", action="store_true")
    seq.add_argument("--run-radius", type=float, default=2500)
    seq.add_argument("--start", default=None, help="Optional platform stop id.")
    seq.add_argument("--time", type=int, default=None,
                     help="Optional absolute week-second target start time.")
    seq.add_argument("--max-columns", type=int, default=0)
    seq.add_argument("--time-limit", type=float, default=30)
    seq.add_argument("--connector-mode", choices=("static", "dynamic"), default="static")
    seq.add_argument("--connector-cap", type=int, default=300000,
                     help="Dijkstra pop cap for dynamic column connector pricing.")
    seq.add_argument("--validate", action="store_true")
    seq.add_argument("--log-search", action="store_true")
    seq.add_argument("--out", default="reports/optimization_runs/column_sequence_route.json")
    seq.set_defaults(func=cmd_sequence)

    pc = sub.add_parser("path-cover", help="Select and sequence columns in one CP-SAT model.")
    pc.add_argument("columns_file")
    pc.add_argument("--terminal-runs", action="store_true")
    pc.add_argument("--run-radius", type=float, default=2500)
    pc.add_argument("--top-k", type=int, default=12,
                    help="Candidate outgoing connector arcs retained per column.")
    pc.add_argument("--column-penalty", type=int, default=0)
    pc.add_argument("--max-columns", type=int, default=0)
    pc.add_argument("--time-limit", type=float, default=60)
    pc.add_argument("--workers", type=int, default=8)
    pc.add_argument("--start", default=None, help="Optional platform stop id.")
    pc.add_argument("--time", type=int, default=None,
                    help="Optional absolute week-second target start time.")
    pc.add_argument("--route-out", default=None,
                    help="Also realize the ordered columns to a route JSON.")
    pc.add_argument("--validate", action="store_true")
    pc.add_argument("--out", default="reports/optimization_runs/path_cover_solution.json")
    pc.set_defaults(func=cmd_path_cover)

    ec = sub.add_parser(
        "exact-cover",
        help="Select and chain exact timed column slices with legal connectors.")
    ec.add_argument("columns_source",
                    help="Columns JSONL or a JSON file containing selected_columns.")
    ec.add_argument("--terminal-runs", action="store_true")
    ec.add_argument("--run-mode", choices=("none", "terminal", "all"), default="none",
                    help="Connector run policy; --terminal-runs aliases terminal.")
    ec.add_argument("--run-radius", type=float, default=2500)
    ec.add_argument("--validation-radius", type=float, default=5000,
                    help="Run radius used when validating/writing exact routes.")
    ec.add_argument("--top-k", type=int, default=12,
                    help="Candidate outgoing exact-time connectors per column.")
    ec.add_argument("--arc-pricing", choices=("proxy", "exact"), default="proxy",
                    help="Use phase proxy arc costs, or prove every candidate by Dijkstra.")
    ec.add_argument("--max-proxy-hours", type=float, default=6.0,
                    help="Filter candidate arcs by static+phase proxy gap.")
    ec.add_argument("--max-connector-hours", type=float, default=6.0,
                    help="Drop exact connector arcs above this duration.")
    ec.add_argument("--connector-cap", type=int, default=300000,
                    help="Dijkstra pop cap for exact connector pricing/replay.")
    ec.add_argument("--exact-arc-cache", default=str(DEFAULT_EXACT_ARC_CACHE),
                    help="JSONL cache for exact event-to-event connector feasibility.")
    ec.add_argument("--no-exact-arc-cache", action="store_true",
                    help="Disable the exact connector cache for this run.")
    ec.add_argument("--repair-cuts", type=int, default=0,
                    help="For proxy pricing, iteratively cut failed replay arcs.")
    ec.add_argument("--refine-selected-arcs", action="store_true",
                    help="Exact-price selected proxy arcs and re-solve if needed.")
    ec.add_argument("--min-same-station-gap", type=int, default=120,
                    help="Drop proxy arcs with too little same-station platform-change slack.")
    ec.add_argument("--min-opposite-direction-gap", type=int, default=180,
                    help="Drop proxy arcs with too little opposite-direction slack.")
    ec.add_argument("--force-solution-arcs", nargs="*", default=[],
                    help="Selected-columns JSON files whose adjacent arcs must be included.")
    ec.add_argument("--hint-solution", nargs="*", default=[],
                    help=("Selected-columns JSON file to use as a CP-SAT solution hint. "
                          "If multiple files are supplied, only the first is used."))
    ec.add_argument("--require-pricing-kind", default=None,
                    help="Require selected columns to include this pricing metadata key.")
    ec.add_argument("--min-required-pricing", type=int, default=0,
                    help="Minimum selected columns with --require-pricing-kind.")
    ec.add_argument("--require-column-id", default=None,
                    help="Comma-separated exact column ids that must be selected.")
    ec.add_argument("--require-column-id-prefix", default=None,
                    help=("Comma-separated column-id prefixes; at least "
                          "--min-required-column-prefix matching columns must be selected."))
    ec.add_argument("--min-required-column-prefix", type=int, default=1,
                    help="Minimum selected columns per --require-column-id-prefix group.")
    ec.add_argument("--require-column-id-prefix-group", default=None,
                    help=("Semicolon-separated [LABEL:]MIN:PREFIX|PREFIX groups. "
                          "At least MIN columns matching any prefix in each group "
                          "must be selected."))
    ec.add_argument("--exclude-column-id", default=None,
                    help="Comma-separated exact column ids to remove from the pool.")
    ec.add_argument("--exclude-column-id-prefix", default=None,
                    help="Comma-separated column-id prefixes to remove from the pool.")
    ec.add_argument("--required-arc-top-k", type=int, default=0,
                    help=("When columns are required by id, prefix, or pricing kind, "
                          "exact-price this many "
                          "extra incoming and outgoing proxy-nearest arcs around "
                          "each matching row. Useful for "
                          "forced-column diagnostics without increasing global --top-k."))
    ec.add_argument("--protected-arc-top-k", type=int, default=0,
                    help=("When stations or station groups are protected, exact-price "
                          "this many extra incoming and outgoing proxy-nearest arcs "
                          "around rows that cover those stations. Useful for small "
                          "residual-repair diagnostics without forcing specific rows."))
    ec.add_argument("--protected-row-start-stop", default=None,
                    help=("Comma-separated start stop ids. When protection is active, "
                          "rows that cover protected stations but start elsewhere are "
                          "forbidden. Useful for partitioning broad protected families."))
    ec.add_argument("--protected-row-column-id-prefix", default=None,
                    help=("Comma-separated column-id prefixes. When protection is "
                          "active, rows that cover protected stations but do not match "
                          "one of these prefixes are forbidden. Combines with "
                          "--protected-row-start-stop."))
    ec.add_argument("--start-time-window", default=None,
                    help="Restrict the route start column to week-second range START-END.")
    ec.add_argument("--end-time-window", default=None,
                    help="Restrict the route end column to week-second range START-END.")
    ec.add_argument("--max-total-elapsed", default=None,
                    help=("Hard cap on modeled column+connector elapsed. Accepts seconds, "
                          "MM:SS, or HH:MM:SS."))
    ec.add_argument("--column-penalty", type=int, default=0)
    ec.add_argument("--uncovered-penalty-s", type=int, default=0,
                    help=("If positive, allow uncovered stations with this penalty "
                          "instead of enforcing hard coverage. Diagnostic only."))
    ec.add_argument("--uncovered-penalty-groups", default=None,
                    help=("Semicolon-separated NAME:PENALTY_SECONDS:station,station "
                          "groups that add extra objective penalty for leaving "
                          "listed stations uncovered. Also enables relaxed coverage."))
    ec.add_argument("--min-covered-count", type=int, default=None,
                    help=("With relaxed coverage, require at least this many canonical "
                          "stations to be covered."))
    ec.add_argument("--protect-stations", default=None,
                    help=("Comma-separated canonical station ids that must be covered "
                          "even when --uncovered-penalty-s is used."))
    ec.add_argument("--protect-stations-file", nargs="*", default=[],
                    help=("JSON list/object of canonical station ids to protect. "
                          "Object keys: stations, protected_stations, or "
                          "relaxed_uncovered_stations."))
    ec.add_argument("--protect-station-groups", default=None,
                    help=("Semicolon-separated NAME:MIN:station,station groups. "
                          "With relaxed coverage, each group must cover at "
                          "least MIN stations without requiring every station."))
    ec.add_argument("--order-station-groups", default=None,
                    help=("Semicolon-separated NAME:MIN:station,station groups. "
                          "Each group chooses one selected witness row covering "
                          "at least MIN listed stations, and witness rows must "
                          "appear in the given order. Diagnostic packet-state "
                          "constraint."))
    ec.add_argument("--strict-order-station-groups", action="store_true",
                    help=("Require ordered station-group witnesses to appear at "
                          "strictly increasing row positions, so one selected "
                          "row cannot satisfy adjacent ordered groups."))
    ec.add_argument("--first-hit-order-station-groups", action="store_true",
                    help=("Interpret ordered station-group positions as the "
                          "earliest selected row covering each group's minimum "
                          "hit threshold. This prevents incidental earlier "
                          "packet coverage from being ignored."))
    ec.add_argument("--max-columns", type=int, default=0)
    ec.add_argument("--time-limit", type=float, default=60)
    ec.add_argument("--workers", type=int, default=8)
    ec.add_argument("--stop-after-first-solution", action="store_true",
                    help=("Stop CP-SAT after the first feasible selected order. "
                          "Useful for large relaxed diagnostics where any "
                          "record-capped basin is more informative than proof "
                          "of optimality."))
    ec.add_argument("--route-out", default=None,
                    help="Also write a route JSON made from exact column slices.")
    ec.add_argument("--validate", action="store_true")
    ec.add_argument("--out", default="reports/optimization_runs/exact_cover_solution.json")
    ec.set_defaults(func=cmd_exact_cover)

    br = sub.add_parser(
        "block-replace",
        help="Try one-column replacements of contiguous blocks in a selected column order.")
    br.add_argument("columns_source",
                    help="Column JSONL or selected-columns JSON containing replacement rows.")
    br.add_argument("source_order",
                    help="Selected-columns JSON or route JSON with meta.selected_order.")
    br.add_argument("--target-stations", default=None,
                    help="Comma-separated canonical stations to improve; defaults to uncovered.")
    br.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal",
                    help="Run policy used for proxy/static candidate screening.")
    br.add_argument("--terminal-runs", action="store_true",
                    help="Alias for --run-mode terminal.")
    br.add_argument("--run-radius", type=float, default=2500)
    br.add_argument("--validation-radius", type=float, default=5000)
    br.add_argument("--connector-cap", type=int, default=300000)
    br.add_argument("--min-covered-count", type=int, default=None,
                    help="Minimum validated station count to write as an improvement.")
    br.add_argument("--max-total-elapsed", default=None,
                    help="Maximum route elapsed; defaults to the source order elapsed.")
    br.add_argument("--max-candidates", type=int, default=20,
                    help="Exact-realize this many best proxy replacement candidates.")
    br.add_argument("--stop-after-first", action="store_true",
                    help="Stop once a candidate meeting the thresholds is written.")
    br.add_argument("--route-out",
                    default="reports/optimization_runs/block_replace_route.json")
    br.add_argument("--out",
                    default="reports/optimization_runs/block_replace_summary.json")
    br.set_defaults(func=cmd_block_replace)

    pair = sub.add_parser(
        "pair-replace",
        help="Try two-column replacements of contiguous blocks in a selected column order.")
    pair.add_argument("columns_source",
                      help="Column JSONL or selected-columns JSON containing replacement rows.")
    pair.add_argument("source_order",
                      help="Selected-columns JSON or route JSON with meta.selected_order.")
    pair.add_argument("--target-stations", default=None,
                      help="Comma-separated canonical stations to improve; defaults to uncovered.")
    pair.add_argument("--column-id-prefix", default=None,
                      help="Comma-separated prefixes allowed for both replacement columns.")
    pair.add_argument("--first-column-id-prefix", default=None,
                      help="Comma-separated prefixes allowed for the first replacement column.")
    pair.add_argument("--second-column-id-prefix", default=None,
                      help="Comma-separated prefixes allowed for the second replacement column.")
    pair.add_argument("--include-non-target-candidates", action="store_true",
                      help="Allow replacement columns that do not hit current target stations.")
    pair.add_argument("--max-first-rows", type=int, default=900,
                      help="Keep this many best first-column candidates after filtering.")
    pair.add_argument("--max-second-rows", type=int, default=900,
                      help="Keep this many best second-column candidates after filtering.")
    pair.add_argument("--min-replace-start-index", type=int, default=1,
                      help=("First selected-column index that may be removed. "
                            "Index 0 is preserved as the left route anchor."))
    pair.add_argument("--max-replace-start-index", type=int, default=0,
                      help="Last selected-column start index that may be removed; 0 means last.")
    pair.add_argument("--min-replace-end-index", type=int, default=0,
                      help="Minimum half-open selected-column end index; 0 means start+1.")
    pair.add_argument("--max-replace-end-index", type=int, default=0,
                      help="Maximum half-open selected-column end index; 0 means route end.")
    pair.add_argument("--min-removed-count", type=int, default=1)
    pair.add_argument("--max-removed-count", type=int, default=0,
                      help="Maximum removed selected columns; 0 means no limit.")
    pair.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal",
                      help="Run policy used for proxy/static candidate screening.")
    pair.add_argument("--terminal-runs", action="store_true",
                      help="Alias for --run-mode terminal.")
    pair.add_argument("--run-radius", type=float, default=2500)
    pair.add_argument("--validation-radius", type=float, default=5000)
    pair.add_argument("--connector-cap", type=int, default=300000)
    pair.add_argument("--min-same-station-gap", type=int, default=120,
                      help="Drop proxy arcs with too little same-station platform-change slack.")
    pair.add_argument("--min-opposite-direction-gap", type=int, default=180,
                      help="Drop proxy arcs with too little opposite-direction slack.")
    pair.add_argument("--min-covered-count", type=int, default=None,
                      help="Minimum validated station count to write as an improvement.")
    pair.add_argument("--max-total-elapsed", default=None,
                      help="Maximum route elapsed; defaults to the source order elapsed.")
    pair.add_argument("--max-candidates", type=int, default=20,
                      help="Exact-realize this many best proxy replacement candidates.")
    pair.add_argument("--progress-every-starts", type=int, default=0,
                      help="Print scan progress after this many replacement start indexes.")
    pair.add_argument("--stop-after-first", action="store_true",
                      help="Stop once a candidate meeting the thresholds is written.")
    pair.add_argument("--route-out",
                      default="reports/optimization_runs/pair_replace_route.json")
    pair.add_argument("--out",
                      default="reports/optimization_runs/pair_replace_summary.json")
    pair.set_defaults(func=cmd_pair_replace)

    chain = sub.add_parser(
        "chain-replace",
        help="Beam-screen multi-column replacements of contiguous selected-column blocks.")
    chain.add_argument("columns_source",
                       help="Column JSONL or selected-columns JSON containing replacement rows.")
    chain.add_argument("source_order",
                       help="Selected-columns JSON or route JSON with meta.selected_order.")
    chain.add_argument("--replacement-count", type=int, default=3,
                       help="Number of ordered replacement columns to insert.")
    chain.add_argument("--target-stations", default=None,
                       help="Comma-separated canonical stations to improve; defaults to uncovered.")
    chain.add_argument("--column-id-prefix", default=None,
                       help="Comma-separated prefixes allowed for replacement columns.")
    chain.add_argument("--chain-column-id-prefixes", default=None,
                       help=("Semicolon-separated prefix groups for each replacement "
                             "position. Each group is a comma-separated list."))
    chain.add_argument("--include-non-target-candidates", action="store_true",
                       help="Allow replacement columns that do not hit current target stations.")
    chain.add_argument("--max-candidate-rows", type=int, default=180,
                       help="Keep this many best target-ranked rows after filtering.")
    chain.add_argument("--extra-fast-candidate-rows", type=int, default=0,
                       help="Also include this many fastest filtered rows.")
    chain.add_argument("--beam-size", type=int, default=1200,
                       help="Partial replacement chains retained per block and depth.")
    chain.add_argument("--keep-proxy-candidates", type=int, default=200,
                       help="Keep this many best proxy-complete chains before exact replay.")
    chain.add_argument("--target-reward-s", type=int, default=1800,
                       help="Beam reward per currently targeted station covered.")
    chain.add_argument("--cover-reward-s", type=int, default=30,
                       help="Beam reward per canonical station covered.")
    chain.add_argument("--time-bucket-s", type=int, default=300,
                       help="Deduplicate partial states by this elapsed-time bucket.")
    chain.add_argument("--min-replace-start-index", type=int, default=1,
                       help=("First selected-column index that may be removed. "
                             "Index 0 is preserved as the left route anchor."))
    chain.add_argument("--max-replace-start-index", type=int, default=0,
                       help="Last selected-column start index that may be removed; 0 means last.")
    chain.add_argument("--min-replace-end-index", type=int, default=0,
                       help="Minimum half-open selected-column end index; 0 means start+1.")
    chain.add_argument("--max-replace-end-index", type=int, default=0,
                       help="Maximum half-open selected-column end index; 0 means route end.")
    chain.add_argument("--min-removed-count", type=int, default=1)
    chain.add_argument("--max-removed-count", type=int, default=0,
                       help="Maximum removed selected columns; 0 means no limit.")
    chain.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal",
                       help="Run policy used for proxy/static candidate screening.")
    chain.add_argument("--terminal-runs", action="store_true",
                       help="Alias for --run-mode terminal.")
    chain.add_argument("--run-radius", type=float, default=2500)
    chain.add_argument("--validation-radius", type=float, default=5000)
    chain.add_argument("--connector-cap", type=int, default=300000)
    chain.add_argument("--min-same-station-gap", type=int, default=120,
                       help="Drop proxy arcs with too little same-station platform-change slack.")
    chain.add_argument("--min-opposite-direction-gap", type=int, default=180,
                       help="Drop proxy arcs with too little opposite-direction slack.")
    chain.add_argument("--min-covered-count", type=int, default=None,
                       help="Minimum validated station count to write as an improvement.")
    chain.add_argument("--max-total-elapsed", default=None,
                       help="Maximum route elapsed; defaults to the source order elapsed.")
    chain.add_argument("--max-candidates", type=int, default=20,
                       help="Exact-realize this many best proxy replacement candidates.")
    chain.add_argument("--progress-every-starts", type=int, default=0,
                       help="Print scan progress after this many replacement start indexes.")
    chain.add_argument("--stop-after-first", action="store_true",
                       help="Stop once a candidate meeting the thresholds is written.")
    chain.add_argument("--route-out",
                       default="reports/optimization_runs/chain_replace_route.json")
    chain.add_argument("--out",
                       default="reports/optimization_runs/chain_replace_summary.json")
    chain.set_defaults(func=cmd_chain_replace)

    pr = sub.add_parser(
        "price-terminals",
        help="Generate new exact columns from seed events to terminal/all targets.")
    pr.add_argument("--seed", default="solutions/best.json")
    pr.add_argument("--base-columns", default=None,
                    help="Optional existing column pool to copy before priced columns.")
    pr.add_argument("--targets", choices=("terminals", "all"), default="terminals")
    pr.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    pr.add_argument("--terminal-runs", action="store_true",
                    help="Alias for --run-mode terminal.")
    pr.add_argument("--run-radius", type=float, default=2500)
    pr.add_argument("--validation-radius", type=float, default=5000)
    pr.add_argument("--source-stride", type=int, default=12)
    pr.add_argument("--max-sources", type=int, default=0)
    pr.add_argument("--max-targets", type=int, default=0)
    pr.add_argument("--max-elapsed-minutes", type=float, default=90)
    pr.add_argument("--min-covered", type=int, default=4)
    pr.add_argument("--connector-cap", type=int, default=300000)
    pr.add_argument("--label", default="priced_terminal")
    pr.add_argument("--out", default="reports/optimization_runs/seed_columns_priced.jsonl")
    pr.set_defaults(func=cmd_price_terminals)

    ppv = sub.add_parser(
        "price-phase-variants",
        help="Replay useful column anchor patterns from alternate departure phases.")
    ppv.add_argument("--source-columns", required=True,
                     help="Column JSONL/selected-columns JSON whose path patterns are replayed.")
    ppv.add_argument("--base-columns", default=None,
                     help="Optional existing column pool to copy before priced columns.")
    ppv.add_argument("--duals", default=None,
                     help="Optional dual-cover JSON used to score covered stations.")
    ppv.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    ppv.add_argument("--terminal-runs", action="store_true",
                     help="Alias for --run-mode terminal.")
    ppv.add_argument("--run-radius", type=float, default=2500)
    ppv.add_argument("--validation-radius", type=float, default=5000)
    ppv.add_argument("--anchor-mode", choices=("first", "transitions", "revisit"),
                     default="transitions")
    ppv.add_argument("--skip-visited-anchors", dest="skip_visited_anchors",
                     action="store_true", default=False,
                     help="Skip anchors already covered incidentally while replaying a fragment.")
    ppv.add_argument("--no-skip-visited-anchors", dest="skip_visited_anchors",
                     action="store_false")
    ppv.add_argument("--source-stop-mode", choices=("same-platform", "station-platforms"),
                     default="same-platform")
    ppv.add_argument("--offset-minutes", default="-120,-90,-60,-30,30,60,90,120",
                     help="Comma offsets from each source column's start event.")
    ppv.add_argument("--include-original", action="store_true")
    ppv.add_argument("--max-source-rows", type=int, default=200)
    ppv.add_argument("--max-generated", type=int, default=0)
    ppv.add_argument("--min-anchors", type=int, default=3)
    ppv.add_argument("--min-covered", type=int, default=2)
    ppv.add_argument("--min-original-overlap", type=float, default=0.6,
                     help="Minimum fraction of the source column's stations retained.")
    ppv.add_argument("--max-elapsed-minutes", type=float, default=180)
    ppv.add_argument("--max-extra-minutes", type=float, default=45,
                     help="Drop variants this much slower than their source column.")
    ppv.add_argument("--cover-reward-s", type=float, default=180.0)
    ppv.add_argument("--column-penalty", type=int, default=0)
    ppv.add_argument("--max-score-s", type=float, default=None,
                     help="Keep variants with elapsed + penalty - reward <= this.")
    ppv.add_argument("--connector-cap", type=int, default=300000)
    ppv.add_argument("--label", default="priced_phase")
    ppv.add_argument("--out", default="reports/optimization_runs/seed_columns_phase_variants.jsonl")
    ppv.set_defaults(func=cmd_price_phase_variants)

    ppw = sub.add_parser(
        "price-phase-windows",
        help="Replay adjacent selected-column windows from alternate departure phases.")
    ppw.add_argument("--source-solution", required=True,
                     help="Selected-columns JSON whose ordered windows are replayed.")
    ppw.add_argument("--base-columns", default=None,
                     help="Optional existing column pool to copy before priced columns.")
    ppw.add_argument("--duals", default=None,
                     help="Optional dual-cover JSON used to score covered stations.")
    ppw.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    ppw.add_argument("--terminal-runs", action="store_true",
                     help="Alias for --run-mode terminal.")
    ppw.add_argument("--run-radius", type=float, default=2500)
    ppw.add_argument("--validation-radius", type=float, default=5000)
    ppw.add_argument("--anchor-mode", choices=("first", "transitions", "revisit"),
                     default="transitions")
    ppw.add_argument("--skip-visited-anchors", dest="skip_visited_anchors",
                     action="store_true", default=False,
                     help="Skip anchors already covered incidentally while replaying a window.")
    ppw.add_argument("--no-skip-visited-anchors", dest="skip_visited_anchors",
                     action="store_false")
    ppw.add_argument("--source-stop-mode", choices=("same-platform", "station-platforms"),
                     default="same-platform")
    ppw.add_argument("--widths", default="2,3,4,6",
                     help="Comma-separated selected-column window widths.")
    ppw.add_argument("--window-stride", type=int, default=1)
    ppw.add_argument("--max-windows", type=int, default=200)
    ppw.add_argument("--offset-minutes", default="-180,-120,-90,-60,-30,30,60,90,120,180",
                     help="Comma offsets from each window's first start event.")
    ppw.add_argument("--variant-start-time-window", default=None,
                     help="Keep only generated windows starting in week-second range START-END.")
    ppw.add_argument("--variant-end-time-window", default=None,
                     help="Keep only generated windows ending in week-second range START-END.")
    ppw.add_argument("--include-original", action="store_true")
    ppw.add_argument("--max-generated", type=int, default=0)
    ppw.add_argument("--min-anchors", type=int, default=4)
    ppw.add_argument("--min-covered", type=int, default=4)
    ppw.add_argument("--min-original-overlap", type=float, default=0.6,
                     help="Minimum fraction of the source window's stations retained.")
    ppw.add_argument("--max-original-minutes", type=float, default=360)
    ppw.add_argument("--max-elapsed-minutes", type=float, default=360)
    ppw.add_argument("--max-extra-minutes", type=float, default=45,
                     help="Drop variants this much slower than their source window.")
    ppw.add_argument("--cover-reward-s", type=float, default=180.0)
    ppw.add_argument("--window-credit", type=float, default=0.0,
                     help="Extra score credit per second of source window elapsed.")
    ppw.add_argument("--column-penalty", type=int, default=0)
    ppw.add_argument("--max-score-s", type=float, default=None,
                     help="Keep variants with adjusted score <= this.")
    ppw.add_argument("--connector-cap", type=int, default=300000)
    ppw.add_argument("--label", default="priced_phase_window")
    ppw.add_argument("--out", default="reports/optimization_runs/seed_columns_phase_windows.jsonl")
    ppw.set_defaults(func=cmd_price_phase_windows)

    pd = sub.add_parser(
        "price-duals",
        help="Generate exact columns from seed events to high-dual station targets.")
    pd.add_argument("--duals", required=True,
                    help="JSON output from the dual-cover command.")
    pd.add_argument("--seed", action="append",
                    help="Seed route JSON. Can be passed more than once; defaults to best.json.")
    pd.add_argument("--base-columns", default=None,
                    help="Optional existing column pool to copy before priced columns.")
    pd.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    pd.add_argument("--terminal-runs", action="store_true",
                    help="Alias for --run-mode terminal.")
    pd.add_argument("--run-radius", type=float, default=2500)
    pd.add_argument("--validation-radius", type=float, default=5000)
    pd.add_argument("--source-stride", type=int, default=12)
    pd.add_argument("--max-sources", type=int, default=0)
    pd.add_argument("--max-targets", type=int, default=80)
    pd.add_argument("--min-dual-s", type=float, default=0.0)
    pd.add_argument("--max-elapsed-minutes", type=float, default=120)
    pd.add_argument("--min-covered", type=int, default=2)
    pd.add_argument("--min-reward-s", type=float, default=0.0)
    pd.add_argument("--max-reduced-cost-s", type=float, default=0.0)
    pd.add_argument("--column-penalty", type=int, default=0)
    pd.add_argument("--connector-cap", type=int, default=300000)
    pd.add_argument("--label", default="priced_dual")
    pd.add_argument("--out", default="reports/optimization_runs/seed_columns_dual_priced.jsonl")
    pd.set_defaults(func=cmd_price_duals)

    pdb = sub.add_parser(
        "price-dual-beam",
        help="Beam-search exact columns using dual station rewards as path prizes.")
    pdb.add_argument("--duals", required=True,
                     help="JSON output from the dual-cover command.")
    pdb.add_argument("--seed", action="append",
                     help="Seed route JSON. Can be passed more than once; defaults to best.json.")
    pdb.add_argument("--base-columns", default=None,
                     help="Optional existing column pool to copy before priced columns.")
    pdb.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    pdb.add_argument("--terminal-runs", action="store_true",
                     help="Alias for --run-mode terminal.")
    pdb.add_argument("--run-radius", type=float, default=2500)
    pdb.add_argument("--validation-radius", type=float, default=5000)
    pdb.add_argument("--source-stride", type=int, default=24)
    pdb.add_argument("--max-sources", type=int, default=0)
    pdb.add_argument("--top-reward-stations", type=int, default=48)
    pdb.add_argument("--min-dual-s", type=float, default=0.0)
    pdb.add_argument("--max-elapsed-minutes", type=float, default=90)
    pdb.add_argument("--min-covered", type=int, default=4)
    pdb.add_argument("--min-reward-stations", type=int, default=2)
    pdb.add_argument("--min-reward-s", type=float, default=0.0)
    pdb.add_argument("--max-reduced-cost-s", type=float, default=0.0)
    pdb.add_argument("--column-penalty", type=int, default=0)
    pdb.add_argument("--max-expansions-per-source", type=int, default=50000)
    pdb.add_argument("--max-labels-per-node", type=int, default=3)
    pdb.add_argument("--candidate-pool-per-source", type=int, default=200)
    pdb.add_argument("--emit-per-source", type=int, default=12)
    pdb.add_argument("--label", default="priced_dual_beam")
    pdb.add_argument("--out", default="reports/optimization_runs/seed_columns_dual_beam_priced.jsonl")
    pdb.set_defaults(func=cmd_price_dual_beam)

    pconn = sub.add_parser(
        "price-connectors",
        help="Generate columns from expensive exact-cover connector handoffs.")
    pconn.add_argument("--solution-files", nargs="+", required=True,
                       help="Exact-cover selected-columns JSON files to mine for connectors.")
    pconn.add_argument("--base-columns", default=None,
                       help="Optional existing column pool to copy before priced columns.")
    pconn.add_argument("--duals", default=None,
                       help="Optional dual-cover JSON used to score covered stations.")
    pconn.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    pconn.add_argument("--terminal-runs", action="store_true",
                       help="Alias for --run-mode terminal.")
    pconn.add_argument("--run-radius", type=float, default=2500)
    pconn.add_argument("--validation-radius", type=float, default=5000)
    pconn.add_argument("--top-connectors", type=int, default=24)
    pconn.add_argument("--min-connector-minutes", type=float, default=4)
    pconn.add_argument("--lookahead-columns", type=int, default=4)
    pconn.add_argument("--source-points", default="end",
                       help="Comma-separated source points to price from: end,start.")
    pconn.add_argument("--max-targets-per-connector", type=int, default=48)
    pconn.add_argument("--include-global-duals", action="store_true",
                       help="Also include globally high-dual stations as targets.")
    pconn.add_argument("--max-elapsed-minutes", type=float, default=90)
    pconn.add_argument("--min-covered", type=int, default=2)
    pconn.add_argument("--cover-reward-s", type=float, default=180.0)
    pconn.add_argument("--connector-credit", type=float, default=0.5,
                       help="Seconds of score credit per second of expensive connector.")
    pconn.add_argument("--max-score-s", type=float, default=0.0,
                       help="Keep columns with elapsed - reward - connector_credit*connector <= this.")
    pconn.add_argument("--connector-cap", type=int, default=300000)
    pconn.add_argument("--label", default="priced_connector")
    pconn.add_argument("--out", default="reports/optimization_runs/seed_columns_connector_priced.jsonl")
    pconn.set_defaults(func=cmd_price_connectors)

    pwe = sub.add_parser(
        "price-window-events",
        help="Generate exact-event replacement columns over multi-column solution windows.")
    pwe.add_argument("--solution-files", nargs="+", required=True,
                     help="Exact-cover selected-columns JSON files to mine for windows.")
    pwe.add_argument("--base-columns", default=None,
                     help="Optional existing column pool to copy before priced columns.")
    pwe.add_argument("--duals", default=None,
                     help="Optional dual-cover JSON used to score covered stations.")
    pwe.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    pwe.add_argument("--terminal-runs", action="store_true",
                     help="Alias for --run-mode terminal.")
    pwe.add_argument("--run-radius", type=float, default=2500)
    pwe.add_argument("--validation-radius", type=float, default=5000)
    pwe.add_argument("--widths", default="2,3,4,6,8",
                     help="Comma-separated selected-column window widths.")
    pwe.add_argument("--window-stride", type=int, default=1)
    pwe.add_argument("--source-points", default="start,end",
                     help="Comma-separated source event anchors: start,end.")
    pwe.add_argument("--target-points", default="end,next-start",
                     help="Comma-separated target anchors: start,end,next-start.")
    pwe.add_argument("--max-windows", type=int, default=500)
    pwe.add_argument("--high-score-windows", action="store_true",
                     help="Price high elapsed-minus-reward windows first instead of low-score ones.")
    pwe.add_argument("--max-window-minutes", type=float, default=180)
    pwe.add_argument("--max-elapsed-minutes", type=float, default=180)
    pwe.add_argument("--max-cost-margin-minutes", type=float, default=0)
    pwe.add_argument("--min-covered", type=int, default=2)
    pwe.add_argument("--cover-reward-s", type=float, default=180.0)
    pwe.add_argument("--window-credit", type=float, default=0.5,
                     help="Seconds of score credit per second in the replaced event window.")
    pwe.add_argument("--max-score-s", type=float, default=0.0,
                     help="Keep columns with elapsed - reward - window_credit*window <= this.")
    pwe.add_argument("--connector-cap", type=int, default=300000)
    pwe.add_argument("--label", default="priced_window_event")
    pwe.add_argument("--out", default="reports/optimization_runs/seed_columns_window_event_priced.jsonl")
    pwe.set_defaults(func=cmd_price_window_events)

    pcc = sub.add_parser(
        "price-cluster-corridors",
        help="Generate exact corridors that pass through a stubborn target-station cluster.")
    pcc.add_argument("--source-solution", required=True,
                     help="Selected-columns JSON whose endpoint events are corridor sources.")
    pcc.add_argument("--destination-columns", default=None,
                     help="Columns JSONL/selected-columns JSON whose endpoint events are destinations.")
    pcc.add_argument("--base-columns", default=None,
                     help="Optional existing column pool to copy before generated corridors.")
    pcc.add_argument("--targets-file", default=None,
                     help="JSON list or relaxed-master JSON containing relaxed_uncovered_stations.")
    pcc.add_argument("--target-stations", default=None,
                     help="Comma-separated canonical station ids to reward and require.")
    pcc.add_argument("--via-stations", default=None,
                     help="Comma-separated target stations to force as one intermediate via.")
    pcc.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    pcc.add_argument("--terminal-runs", action="store_true",
                     help="Alias for --run-mode terminal.")
    pcc.add_argument("--run-radius", type=float, default=2500)
    pcc.add_argument("--validation-radius", type=float, default=5000)
    pcc.add_argument("--source-points", default="end",
                     help="Comma-separated source endpoint points: start,end.")
    pcc.add_argument("--destination-points", default="start",
                     help="Comma-separated destination endpoint points: start,end.")
    pcc.add_argument("--destination-mode", choices=("exact", "station"), default="exact",
                     help=("exact requires the corridor to end at the chosen endpoint event; "
                           "station ends at the first reachable event for that endpoint station."))
    pcc.add_argument("--source-time-window", default=None,
                     help="Keep only source endpoints in week-second range START-END.")
    pcc.add_argument("--destination-time-window", default=None,
                     help="Keep only destination endpoints in week-second range START-END.")
    pcc.add_argument("--max-sources", type=int, default=0)
    pcc.add_argument("--max-destinations", type=int, default=0)
    pcc.add_argument("--max-vias", type=int, default=0)
    pcc.add_argument("--max-generated", type=int, default=0)
    pcc.add_argument("--max-elapsed-minutes", type=float, default=240)
    pcc.add_argument("--max-leg-minutes", type=float, default=0,
                     help="Optional cap for each source-via and via-destination leg.")
    pcc.add_argument("--min-target-hits", type=int, default=2)
    pcc.add_argument("--min-covered", type=int, default=4)
    pcc.add_argument("--target-reward-s", type=float, default=1800.0)
    pcc.add_argument("--cover-reward-s", type=float, default=90.0)
    pcc.add_argument("--column-penalty", type=int, default=0)
    pcc.add_argument("--max-score-s", type=float, default=None)
    pcc.add_argument("--connector-cap", type=int, default=300000)
    pcc.add_argument("--label", default="priced_cluster_corridor")
    pcc.add_argument("--out", default="reports/optimization_runs/seed_columns_cluster_corridors.jsonl")
    pcc.set_defaults(func=cmd_price_cluster_corridors)

    pas = sub.add_parser(
        "price-anchor-sequences",
        help="Generate exact columns by forcing explicit ordered station-anchor sequences.")
    pas.add_argument("--source-columns", required=True,
                     help="Column JSONL/selected-columns JSON containing the source column.")
    pas.add_argument("--source-column-id", required=True)
    pas.add_argument("--source-point", choices=("start", "end"), default="end")
    pas.add_argument("--base-columns", default=None,
                     help="Optional existing column pool to copy before generated columns.")
    pas.add_argument("--anchor-sequences", default=None,
                     help="Semicolon-separated comma station-id sequences.")
    pas.add_argument("--anchor-sequences-file", default=None,
                     help="JSON list of comma strings or station-id lists.")
    pas.add_argument("--target-stations", default=None,
                     help="Optional canonical station ids used only for hit metadata/filtering.")
    pas.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    pas.add_argument("--terminal-runs", action="store_true",
                     help="Alias for --run-mode terminal.")
    pas.add_argument("--run-radius", type=float, default=2500)
    pas.add_argument("--validation-radius", type=float, default=5000)
    pas.add_argument("--max-end-time", type=int, default=None,
                     help="Optional absolute week-second cap for generated column end events.")
    pas.add_argument("--max-elapsed-minutes", type=float, default=240)
    pas.add_argument("--max-leg-minutes", type=float, default=0)
    pas.add_argument("--min-covered", type=int, default=4)
    pas.add_argument("--min-target-hits", type=int, default=0)
    pas.add_argument("--connector-cap", type=int, default=300000)
    pas.add_argument("--label", default="priced_anchor_sequence")
    pas.add_argument("--verbose", action="store_true")
    pas.add_argument("--out", default="reports/optimization_runs/seed_columns_anchor_sequences.jsonl")
    pas.set_defaults(func=cmd_price_anchor_sequences)

    psc = sub.add_parser(
        "price-stage-chains",
        help="Generate exact columns through ordered resource-stage station groups.")
    psc.add_argument("--source-columns", required=True,
                     help="Column JSONL/selected-columns JSON containing the source column.")
    psc.add_argument("--source-column-id", required=True)
    psc.add_argument("--source-point", choices=("start", "end"), default="end")
    psc.add_argument("--base-columns", default=None,
                     help="Optional existing column pool to copy before generated columns.")
    psc.add_argument("--stage-groups", default=None,
                     help="Semicolon-separated comma station-id groups to satisfy in order.")
    psc.add_argument("--stage-groups-file", default=None,
                     help="JSON list of comma strings or station-id lists.")
    psc.add_argument("--stage-min-hits", default=None,
                     help=("Comma-separated minimum covered stations per stage. "
                           "One value applies to every stage. Default is 1."))
    psc.add_argument("--target-stations", default=None,
                     help="Optional additional canonical station ids to reward.")
    psc.add_argument("--targets-file", default=None,
                     help="JSON list or relaxed-master JSON containing relaxed_uncovered_stations.")
    psc.add_argument("--final-stations", required=True,
                     help="Comma-separated canonical station ids allowed as generated endpoints.")
    psc.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    psc.add_argument("--terminal-runs", action="store_true",
                     help="Alias for --run-mode terminal.")
    psc.add_argument("--run-radius", type=float, default=2500)
    psc.add_argument("--validation-radius", type=float, default=5000)
    psc.add_argument("--max-end-time", type=int, default=None,
                     help="Optional absolute week-second cap for generated column end events.")
    psc.add_argument("--max-elapsed-minutes", type=float, default=420)
    psc.add_argument("--max-leg-minutes", type=float, default=150)
    psc.add_argument("--max-stage-targets", type=int, default=0,
                     help="Limit each stage to this many listed stations; 0 allows all.")
    psc.add_argument("--max-stage-depth", type=int, default=0,
                     help=("Maximum station-to-station expansions inside each stage. "
                           "0 uses that stage's minimum hit count."))
    psc.add_argument("--force-stage-move", action="store_true",
                     help="Do not let incidental prior coverage satisfy a later stage.")
    psc.add_argument("--beam-size", type=int, default=120)
    psc.add_argument("--emit-top", type=int, default=200)
    psc.add_argument("--emit-frontier-top", type=int, default=0,
                     help="Also emit this many best satisfied stage-boundary states.")
    psc.add_argument("--frontier-only", action="store_true",
                     help="Skip final-station extensions and emit only stage frontiers.")
    psc.add_argument("--max-generated", type=int, default=0)
    psc.add_argument("--min-target-hits", type=int, default=4)
    psc.add_argument("--min-covered", type=int, default=20)
    psc.add_argument("--frontier-min-target-hits", type=int, default=1,
                     help="Minimum target hits for emitted stage-frontier states.")
    psc.add_argument("--frontier-min-covered", type=int, default=1,
                     help="Minimum covered stations for emitted stage-frontier states.")
    psc.add_argument("--stage-reward-s", type=float, default=1200.0)
    psc.add_argument("--target-reward-s", type=float, default=1800.0)
    psc.add_argument("--cover-reward-s", type=float, default=45.0)
    psc.add_argument("--column-penalty", type=int, default=0)
    psc.add_argument("--time-bucket-s", type=int, default=300)
    psc.add_argument("--max-score-s", type=float, default=None)
    psc.add_argument("--connector-cap", type=int, default=300000)
    psc.add_argument("--label", default="priced_stage_chain")
    psc.add_argument("--out", default="reports/optimization_runs/seed_columns_stage_chains.jsonl")
    psc.set_defaults(func=cmd_price_stage_chains)

    prc = sub.add_parser(
        "price-resource-chains",
        help="Beam-search exact columns with unordered resource-group state.")
    prc.add_argument("--source-columns", required=True,
                     help="Column JSONL/selected-columns JSON containing the source column.")
    prc.add_argument("--source-column-id", required=True)
    prc.add_argument("--source-point", choices=("start", "end"), default="end")
    prc.add_argument("--base-columns", default=None,
                     help="Optional existing column pool to copy before generated columns.")
    prc.add_argument("--resource-groups", default=None,
                     help=("Semicolon-separated NAME:station,station resource groups. "
                           "NAME may also be omitted."))
    prc.add_argument("--resource-groups-file", default=None,
                     help=("JSON object name -> stations, or list of objects with "
                           "name/stations fields."))
    prc.add_argument("--resource-min-hits", default=None,
                     help=("Comma-separated minimum covered stations per resource. "
                           "One value applies to every resource. Default is 1."))
    prc.add_argument("--require-resources", default=None,
                     help="Comma-separated resource names that every emitted row must satisfy.")
    prc.add_argument("--target-stations", default=None,
                     help="Optional extra canonical stations to reward/expand.")
    prc.add_argument("--targets-file", default=None,
                     help="JSON list or relaxed-master JSON containing relaxed_uncovered_stations.")
    prc.add_argument("--final-stations", required=True,
                     help="Comma-separated canonical station ids allowed as generated endpoints.")
    prc.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    prc.add_argument("--terminal-runs", action="store_true",
                     help="Alias for --run-mode terminal.")
    prc.add_argument("--run-radius", type=float, default=2500)
    prc.add_argument("--validation-radius", type=float, default=5000)
    prc.add_argument("--max-end-time", type=int, default=None,
                     help="Optional absolute week-second cap for generated column end events.")
    prc.add_argument("--max-elapsed-minutes", type=float, default=420)
    prc.add_argument("--max-leg-minutes", type=float, default=150)
    prc.add_argument("--beam-size", type=int, default=160)
    prc.add_argument("--max-depth", type=int, default=8)
    prc.add_argument("--max-expand-targets", type=int, default=32)
    prc.add_argument("--max-expansions-per-depth", type=int, default=0,
                     help="Cap exact leg expansions per depth; 0 means no cap.")
    prc.add_argument("--emit-top", type=int, default=200)
    prc.add_argument("--emit-frontier-top", type=int, default=0,
                     help="Also emit this many best non-final beam states as exact columns.")
    prc.add_argument("--frontier-only", action="store_true",
                     help="Skip final-station extensions and emit only frontier states.")
    prc.add_argument("--max-generated", type=int, default=0)
    prc.add_argument("--min-resource-count", type=int, default=1)
    prc.add_argument("--min-target-hits", type=int, default=4)
    prc.add_argument("--min-covered", type=int, default=20)
    prc.add_argument("--frontier-min-resource-count", type=int, default=None,
                     help="Minimum satisfied resources for emitted frontier states.")
    prc.add_argument("--frontier-min-target-hits", type=int, default=None,
                     help="Minimum target hits for emitted frontier states.")
    prc.add_argument("--frontier-min-covered", type=int, default=None,
                     help="Minimum covered stations for emitted frontier states.")
    prc.add_argument("--resource-reward-s", type=float, default=2400.0)
    prc.add_argument("--resource-hit-reward-s", type=float, default=600.0)
    prc.add_argument("--target-reward-s", type=float, default=1200.0)
    prc.add_argument("--cover-reward-s", type=float, default=45.0)
    prc.add_argument("--column-penalty", type=int, default=0)
    prc.add_argument("--time-bucket-s", type=int, default=300)
    prc.add_argument("--max-score-s", type=float, default=None)
    prc.add_argument("--connector-cap", type=int, default=300000)
    prc.add_argument("--resource-chain-cache",
                     default="reports/optimization_runs/resource_chain_leg_cache.jsonl",
                     help="Append-only JSONL cache for exact resource-chain leg paths.")
    prc.add_argument("--no-resource-chain-cache", action="store_true",
                     help="Disable the resource-chain leg-path cache.")
    prc.add_argument("--label", default="priced_resource_chain")
    prc.add_argument("--out", default="reports/optimization_runs/seed_columns_resource_chains.jsonl")
    prc.set_defaults(func=cmd_price_resource_chains)

    plt = sub.add_parser(
        "price-late-tail-beam",
        help="Beam-search exact late-tail columns with station prizes and a hard end-time cap.")
    plt.add_argument("--source-columns", required=True,
                     help="Column JSONL/selected-columns JSON containing the source column.")
    plt.add_argument("--source-column-id", required=True)
    plt.add_argument("--source-point", choices=("start", "end"), default="end")
    plt.add_argument("--base-columns", default=None,
                     help="Optional existing column pool to copy before generated columns.")
    plt.add_argument("--targets-file", default=None,
                     help="JSON list or relaxed-master JSON containing relaxed_uncovered_stations.")
    plt.add_argument("--target-stations", default=None,
                     help="Comma-separated canonical station ids, in preferred expansion order.")
    plt.add_argument("--final-stations", required=True,
                     help="Comma-separated canonical station ids allowed as beam endpoints.")
    plt.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    plt.add_argument("--terminal-runs", action="store_true",
                     help="Alias for --run-mode terminal.")
    plt.add_argument("--run-radius", type=float, default=2500)
    plt.add_argument("--validation-radius", type=float, default=5000)
    plt.add_argument("--max-end-time", type=int, required=True,
                     help="Absolute week-second cap for generated column end events.")
    plt.add_argument("--max-elapsed-minutes", type=float, default=420)
    plt.add_argument("--max-leg-minutes", type=float, default=150)
    plt.add_argument("--beam-size", type=int, default=80)
    plt.add_argument("--max-depth", type=int, default=8)
    plt.add_argument("--max-expand-targets", type=int, default=24)
    plt.add_argument("--emit-top", type=int, default=200)
    plt.add_argument("--max-generated", type=int, default=0)
    plt.add_argument("--min-target-hits", type=int, default=8)
    plt.add_argument("--min-covered", type=int, default=40)
    plt.add_argument("--target-reward-s", type=float, default=2400.0)
    plt.add_argument("--cover-reward-s", type=float, default=60.0)
    plt.add_argument("--column-penalty", type=int, default=0)
    plt.add_argument("--max-score-s", type=float, default=None)
    plt.add_argument("--connector-cap", type=int, default=300000)
    plt.add_argument("--label", default="priced_late_tail_beam")
    plt.add_argument("--out", default="reports/optimization_runs/seed_columns_late_tail_beam.jsonl")
    plt.set_defaults(func=cmd_price_late_tail_beam)

    spl = sub.add_parser(
        "split-columns",
        help="Split exact priced columns into smaller exact subcolumns at gateway stations.")
    spl.add_argument("source_columns", nargs="+",
                     help="Column JSONL/selected-columns JSON files containing rows to split.")
    spl.add_argument("--base-columns", default=None,
                     help="Optional existing column pool to copy before split columns.")
    spl.add_argument("--pricing-kind", default=None,
                     help="Only split rows containing this pricing metadata key.")
    spl.add_argument("--column-id-prefix", default=None,
                     help="Comma-separated column-id prefixes to split.")
    spl.add_argument("--split-stations", required=True,
                     help="Comma-separated canonical station ids used as split gateways.")
    spl.add_argument("--min-covered", type=int, default=8)
    spl.add_argument("--min-target-hits", type=int, default=1,
                     help="Minimum split gateway stations covered by the subcolumn.")
    spl.add_argument("--min-steps", type=int, default=2)
    spl.add_argument("--max-elapsed-minutes", type=float, default=180)
    spl.add_argument("--max-split-gap", type=int, default=0,
                     help="Limit subcolumns to this many split intervals; 0 allows all pairs.")
    spl.add_argument("--max-generated", type=int, default=0)
    spl.add_argument("--validation-radius", type=float, default=5000)
    spl.add_argument("--label", default="priced_split")
    spl.add_argument("--verbose", action="store_true")
    spl.add_argument("--out", default="reports/optimization_runs/seed_columns_split.jsonl")
    spl.set_defaults(func=cmd_split_columns)

    pw = sub.add_parser(
        "price-windows",
        help="Generate columns around expensive seed-route windows.")
    pw.add_argument("--seed", default="solutions/best.json")
    pw.add_argument("--base-columns", default=None,
                    help="Optional existing column pool to copy before priced columns.")
    pw.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    pw.add_argument("--terminal-runs", action="store_true",
                    help="Alias for --run-mode terminal.")
    pw.add_argument("--run-radius", type=float, default=2500)
    pw.add_argument("--validation-radius", type=float, default=5000)
    pw.add_argument("--widths", default="24,36,48",
                    help="Comma-separated route-node window widths.")
    pw.add_argument("--window-stride", type=int, default=6)
    pw.add_argument("--lookahead", type=int, default=72)
    pw.add_argument("--max-windows", type=int, default=80)
    pw.add_argument("--max-targets-per-window", type=int, default=0)
    pw.add_argument("--max-elapsed-minutes", type=float, default=120)
    pw.add_argument("--min-covered", type=int, default=4)
    pw.add_argument("--cover-reward-s", type=int, default=180,
                    help="Seconds of reward per station when ranking expensive windows.")
    pw.add_argument("--connector-cap", type=int, default=300000)
    pw.add_argument("--label", default="priced_window")
    pw.add_argument("--out", default="reports/optimization_runs/seed_columns_window_priced.jsonl")
    pw.set_defaults(func=cmd_price_windows)

    pcuts = sub.add_parser(
        "price-cuts",
        help="Generate feasible bridge columns from failed proxy connector cuts.")
    pcuts.add_argument("--base-columns", required=True)
    pcuts.add_argument("--cut-files", nargs="+", required=True)
    pcuts.add_argument("--run-mode", choices=("none", "terminal", "all"), default="terminal")
    pcuts.add_argument("--terminal-runs", action="store_true",
                       help="Alias for --run-mode terminal.")
    pcuts.add_argument("--run-radius", type=float, default=2500)
    pcuts.add_argument("--validation-radius", type=float, default=5000)
    pcuts.add_argument("--max-cuts", type=int, default=0)
    pcuts.add_argument("--max-targets-per-cut", type=int, default=8)
    pcuts.add_argument("--include-covered", action="store_true")
    pcuts.add_argument("--max-elapsed-minutes", type=float, default=120)
    pcuts.add_argument("--min-covered", type=int, default=1)
    pcuts.add_argument("--connector-cap", type=int, default=300000)
    pcuts.add_argument("--label", default="priced_cut")
    pcuts.add_argument("--out", default="reports/optimization_runs/seed_columns_cut_priced.jsonl")
    pcuts.set_defaults(func=cmd_price_cuts)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
