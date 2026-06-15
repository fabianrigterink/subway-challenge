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


def _candidate_exact_arcs(rows, D, top_k, max_proxy_s,
                          min_same_station_gap_s=0,
                          min_opposite_direction_gap_s=0):
    arcs = []
    for i, row_i in enumerate(rows):
        a = row_i["end"]["station"]
        end_t = int(row_i["end"]["time"])
        end_stop = row_i["end"]["stop"]
        candidates = []
        for j, row_j in enumerate(rows):
            if i == j:
                continue
            b = row_j["start"]["station"]
            static = 0 if a == b else int(D.get(a, {}).get(b, BIG_COST))
            if static >= BIG_COST:
                continue
            raw_gap = (int(row_j["start"]["time"]) - end_t) % WEEK
            start_stop = row_j["start"]["stop"]
            if a == b and end_stop != start_stop and raw_gap < min_same_station_gap_s:
                continue
            if (_direction(end_stop) and _direction(start_stop)
                    and _direction(end_stop) != _direction(start_stop)
                    and raw_gap < min_opposite_direction_gap_s):
                continue
            proxy = _phase_proxy_gap(end_t, row_j["start"]["time"], static)
            if proxy <= max_proxy_s:
                candidates.append((proxy, i, j))
        candidates.sort()
        arcs.extend(candidates[:top_k])
    return arcs


def _exact_connector_arcs(G, tables, rows, runs, D, top_k, max_proxy_s,
                          max_connector_s, cap, run_mode, run_radius,
                          cache_path=None,
                          min_same_station_gap_s=0,
                          min_opposite_direction_gap_s=0):
    node_off, node_stop, node_t = tables
    stop_events = _stop_index(node_stop, node_t)
    candidates = _candidate_exact_arcs(
        rows, D, top_k, max_proxy_s,
        min_same_station_gap_s=min_same_station_gap_s,
        min_opposite_direction_gap_s=min_opposite_direction_gap_s)
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
            print(f"exact connectors: {k}/{len(candidates)} candidates "
                  f"kept={len(arcs)}", flush=True)
    _append_exact_arc_cache(cache_path, cache_new)
    if cache_path:
        print(f"exact connector cache: hits={cache_hits} misses={cache_misses} "
              f"new={len(cache_new)} entries={len(cache)} path={cache_path}",
              flush=True)
    return arcs


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
    for solution_file in solution_files:
        solution_rows = _load_column_rows(Path(solution_file))
        for row_i, row_j in zip(solution_rows, solution_rows[1:]):
            i = by_id.get(row_i.get("column_id"))
            j = by_id.get(row_j.get("column_id"))
            if i is None or j is None:
                missing += 1
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
                continue
            forced.append((gap, i, j))
    print(f"forced solution arcs: kept={len(forced)} missing={missing} "
          f"infeasible={infeasible}", flush=True)
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
                      allowed_start_rows=None, allowed_end_rows=None,
                      uncovered_penalty_s=0,
                      hint_order_indices=None):
    canonical = set(StationIndex.load().canonical_stations)
    by_station = {station: [] for station in canonical}
    for i, row in enumerate(rows):
        for station in row.get("covered_stations", []):
            if station in by_station:
                by_station[station].append(i)
    missing = sorted(st for st, cols in by_station.items() if not cols)
    if missing and not uncovered_penalty_s:
        raise ValueError(f"column pool misses {len(missing)} station(s), e.g. {missing[:8]}")

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
    if uncovered_penalty_s:
        for station, cols in by_station.items():
            covered = model.NewBoolVar(f"covered_{station}")
            cover_vars[station] = covered
            if cols:
                model.Add(sum(x[i] for i in cols) >= covered).WithName(
                    f"covered_if_selected_{station}")
            else:
                model.Add(covered == 0).WithName(f"uncoverable_{station}")
    else:
        for station, cols in by_station.items():
            model.Add(sum(x[i] for i in cols) >= 1).WithName(f"cover_{station}")
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
    objective = (
        sum((int(row["elapsed_s"]) + column_penalty) * x[i] for i, row in enumerate(rows))
        + sum(arc_cost[i, j] * var for (i, j), var in y.items())
    )
    if uncovered_penalty_s:
        objective += sum(int(uncovered_penalty_s) * (1 - var)
                         for var in cover_vars.values())
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
    if isinstance(data.get("selected_columns"), list):
        return [
            str(row["column_id"])
            for row in data["selected_columns"]
            if isinstance(row, dict) and row.get("column_id")
        ]
    raise ValueError(f"{path} has no selected_order or selected_columns")


def cmd_exact_cover(args) -> int:
    import pickle

    rows = _load_column_rows(Path(args.columns_source))
    if args.max_columns:
        rows = rows[:args.max_columns]
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
    allowed_start_rows = _rows_in_time_window(rows, "start", start_window)
    allowed_end_rows = _rows_in_time_window(rows, "end", end_window)
    row_by_id = {row.get("column_id"): i for i, row in enumerate(rows)}
    required_column_ids = [
        column_id.strip()
        for column_id in str(args.require_column_id or "").split(",")
        if column_id.strip()
    ]
    required_row_indices = []
    if required_column_ids:
        missing_required = [
            column_id for column_id in required_column_ids if column_id not in row_by_id
        ]
        if missing_required:
            raise SystemExit(
                f"--require-column-id missing ids: {missing_required[:8]}")
        required_row_indices = [row_by_id[column_id] for column_id in required_column_ids]
        print(f"required column ids: {len(required_row_indices)}", flush=True)
    hint_order_indices = []
    if args.hint_solution:
        hint_column_ids = []
        for hint_path in args.hint_solution:
            hint_column_ids.extend(_selected_column_ids_from_solution(Path(hint_path)))
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
            allowed_start_rows=allowed_start_rows,
            allowed_end_rows=allowed_end_rows,
            uncovered_penalty_s=args.uncovered_penalty_s,
            hint_order_indices=hint_order_indices,
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
        if args.uncovered_penalty_s:
            relaxed_uncovered = sorted(
                station
                for station, var in solve_vars.get("covered", {}).items()
                if not solver.Value(var)
            )
        raw_elapsed = sum(int(row["elapsed_s"]) for row in ordered_rows)
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
        result = {
            "status": status_name,
            "objective_s": int(solver.ObjectiveValue()),
            "objective": hms(int(solver.ObjectiveValue())),
            "selected_count": len(ordered_rows),
            "covered_count": len(coverage & si.canonical_stations),
            "raw_column_elapsed_s": raw_elapsed,
            "raw_column_elapsed": hms(raw_elapsed),
            "mode_edges": dict(modes),
            "proxy_cuts": [[rows[i]["column_id"], rows[j]["column_id"]]
                           for i, j in sorted(banned)],
            "uncovered_penalty_s": int(args.uncovered_penalty_s),
            "relaxed_uncovered_stations": relaxed_uncovered,
            "start_time_window": list(start_window) if start_window is not None else None,
            "end_time_window": list(end_window) if end_window is not None else None,
            "required_column_ids": required_column_ids,
            "hint_solution": list(args.hint_solution or []),
            "selected_order": [row.get("column_id") for row in ordered_rows],
            "selected_arcs": selected_arcs,
            "selected_columns": ordered_rows,
        }
        out.write_text(json.dumps(result, sort_keys=True))
        print(f"status={status_name} selected={len(ordered_rows)} "
              f"covered={result['covered_count']}/472 raw={hms(raw_elapsed)} "
              f"objective={result['objective']}")
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
        solution_rows = _load_column_rows(Path(solution_path))
        label = f"solution:{Path(solution_path).name}"
        for row in solution_rows:
            column_id = row.get("column_id")
            if not column_id:
                continue
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
    stop_events = _stop_index(node_stop, node_t)
    run_mode = "terminal" if args.terminal_runs else args.run_mode
    runs, _extra_edges = _run_layer_and_metric_edges(
        G, si, node_off, run_mode, args.run_radius)
    prob = Problem(G, RunLayer.from_graph(G, radius_m=args.validation_radius), si)

    source_start, source_end = _row_endpoint_nodes(stop_events, node_t, source_row)
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
    source_rows = []
    for source_file in args.source_columns:
        for row in _load_column_rows(Path(source_file)):
            if args.pricing_kind and args.pricing_kind not in row:
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
                    help="Selected-columns JSON files to use as CP-SAT solution hints.")
    ec.add_argument("--require-pricing-kind", default=None,
                    help="Require selected columns to include this pricing metadata key.")
    ec.add_argument("--min-required-pricing", type=int, default=0,
                    help="Minimum selected columns with --require-pricing-kind.")
    ec.add_argument("--require-column-id", default=None,
                    help="Comma-separated exact column ids that must be selected.")
    ec.add_argument("--start-time-window", default=None,
                    help="Restrict the route start column to week-second range START-END.")
    ec.add_argument("--end-time-window", default=None,
                    help="Restrict the route end column to week-second range START-END.")
    ec.add_argument("--column-penalty", type=int, default=0)
    ec.add_argument("--uncovered-penalty-s", type=int, default=0,
                    help=("If positive, allow uncovered stations with this penalty "
                          "instead of enforcing hard coverage. Diagnostic only."))
    ec.add_argument("--max-columns", type=int, default=0)
    ec.add_argument("--time-limit", type=float, default=60)
    ec.add_argument("--workers", type=int, default=8)
    ec.add_argument("--route-out", default=None,
                    help="Also write a route JSON made from exact column slices.")
    ec.add_argument("--validate", action="store_true")
    ec.add_argument("--out", default="reports/optimization_runs/exact_cover_solution.json")
    ec.set_defaults(func=cmd_exact_cover)

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
