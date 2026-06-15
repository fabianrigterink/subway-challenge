"""Solver-backed optimization experiments for the Subway Challenge.

The commands here deliberately keep a short path back to the canonical
validator: solve a compact/macro model, realize the resulting station order on
the real time-expanded graph, then emit a normal route JSON.
"""
from __future__ import annotations

import argparse
import json
import pickle
import time
from pathlib import Path

from ortools.constraint_solver import pywrapcp, routing_enums_pb2

from .run_layer import RunLayer
from .search import (
    _elapsed,
    _first_event_at_or_after,
    _load_seed_anchors,
    _node_tables,
    _parse_days,
    _parse_time_of_day,
    _stop_index,
    _terminal_run_layer,
    realize_from,
    static_station_metric,
)
from .solver import GRAPH_PKL, Problem, hms, result_line, validate
from .stations import StationIndex

BIG_COST = 10**9


def _station_subset(canonical, start_station, max_stations=0, seed_order=None):
    """Choose a deterministic station subset for smoke tests or full solves."""
    if max_stations and max_stations < 2:
        raise ValueError("--max-stations must be 0 or at least 2")
    if not max_stations:
        return sorted(canonical)

    ordered = []
    if seed_order:
        ordered.extend(s for s in seed_order if s in canonical)
    ordered.extend(sorted(canonical))

    seen = set()
    subset = []
    for station in [start_station] + ordered:
        if station in canonical and station not in seen:
            seen.add(station)
            subset.append(station)
        if len(subset) >= max_stations:
            break
    return subset


def _solve_ortools_path(stations, D, start_station, time_limit_s, seed_order=None,
                        log_search=False):
    """Solve an open ATSP path from start_station to a dummy zero-cost sink."""
    station_to_idx = {s: i for i, s in enumerate(stations)}
    if start_station not in station_to_idx:
        raise ValueError(f"start station {start_station!r} not in station set")

    dummy = len(stations)
    manager = pywrapcp.RoutingIndexManager(
        len(stations) + 1, 1, [station_to_idx[start_station]], [dummy])
    routing = pywrapcp.RoutingModel(manager)

    def transit(from_index, to_index):
        a = manager.IndexToNode(from_index)
        b = manager.IndexToNode(to_index)
        if b == dummy:
            return 0
        if a == dummy:
            return BIG_COST
        if a == b:
            return 0
        return int(D.get(stations[a], {}).get(stations[b], BIG_COST))

    transit_idx = routing.RegisterTransitCallback(transit)
    routing.SetArcCostEvaluatorOfAllVehicles(transit_idx)

    params = pywrapcp.DefaultRoutingSearchParameters()
    params.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
    params.local_search_metaheuristic = routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
    params.time_limit.FromSeconds(max(1, int(time_limit_s)))
    params.log_search = bool(log_search)

    assignment = None
    if seed_order:
        route = [station_to_idx[s] for s in seed_order
                 if s in station_to_idx and s != start_station]
        if route:
            assignment = routing.ReadAssignmentFromRoutes([route], True)

    if assignment is not None:
        solution = routing.SolveFromAssignmentWithParameters(assignment, params)
    else:
        solution = routing.SolveWithParameters(params)
    if solution is None:
        return None, None

    order = []
    idx = routing.Start(0)
    while not routing.IsEnd(idx):
        node = manager.IndexToNode(idx)
        if node != dummy:
            order.append(stations[node])
        idx = solution.Value(routing.NextVar(idx))
    return order, int(solution.ObjectiveValue())


def _route_path_json(G, tables, runs, start_node, station_order, canonical):
    node_off, node_stop, node_t = tables
    suffix, end, visited = realize_from(
        G._adj, tables, start_node, {node_off[start_node]}, station_order[1:], runs=runs)
    if suffix is None:
        return None, None, None
    nodes = [start_node] + suffix
    path = [[node_stop[n], node_t[n]] for n in nodes]
    elapsed = _elapsed(node_t[start_node], node_t[end])
    covered = len(visited & canonical)
    return path, elapsed, covered


def cmd_ortools_atsp(args) -> int:
    t0 = time.time()
    with open(GRAPH_PKL, "rb") as f:
        G = pickle.load(f)
    si = StationIndex.load()
    tables = _node_tables(G, si)
    node_off, _node_stop, node_t = tables
    stop_events = _stop_index(tables[1], node_t)

    target_t = _parse_days(args.day)[0] * 86400 + _parse_time_of_day(args.time)
    start_node = _first_event_at_or_after(stop_events, args.start, target_t)
    if start_node is None:
        raise SystemExit(f"no event found for start stop {args.start!r}")
    start_station = node_off[start_node]
    print(f"start event: {args.start}@{node_t[start_node]} station={start_station}", flush=True)

    runs, extra_edges = (None, [])
    if args.terminal_runs:
        runs, extra_edges = _terminal_run_layer(G, si, node_off, args.run_radius)

    print("building static station metric", flush=True)
    D, _ = static_station_metric(G, node_off, extra_edges)

    seed_order = None
    if args.seed_from:
        _spath, seed_order = _load_seed_anchors(args.seed_from, si)
    stations = _station_subset(
        si.canonical_stations, start_station, args.max_stations, seed_order)
    if start_station not in stations:
        stations.insert(0, start_station)
    print(f"OR-Tools ATSP path: stations={len(stations)} "
          f"time_limit={args.time_limit}s seed_from={args.seed_from or 'none'}", flush=True)

    order, static_obj = _solve_ortools_path(
        stations, D, start_station, args.time_limit, seed_order, args.log_search)
    if not order:
        raise SystemExit("OR-Tools failed to produce a route")
    print(f"OR-Tools static objective: {static_obj} ({hms(static_obj)})", flush=True)

    path, elapsed, covered = _route_path_json(
        G, tables, runs, start_node, order, si.canonical_stations)
    if path is None:
        raise SystemExit("failed to realize OR-Tools station order on schedule")

    meta = {
        "radius_m": args.run_radius if args.terminal_runs else 5000,
        "elapsed_s": int(elapsed),
        "notes": ("OR-Tools ATSP path relaxation "
                  f"start={args.start} day={args.day} time={args.time} "
                  f"stations={len(stations)} static_obj={static_obj}"),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"meta": meta, "path": path}))
    print(f"realized: elapsed={hms(elapsed)} covered={covered}/{len(si.canonical_stations)} "
          f"steps={len(path)} -> {out}", flush=True)

    if args.validate:
        full_runs = RunLayer.from_graph(G, radius_m=meta["radius_m"])
        result = validate(Problem(G, full_runs, si), path)
        print(result_line(result), flush=True)
    print(f"done in {time.time() - t0:.1f}s", flush=True)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Optimization-backed Subway Challenge experiments.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    atsp = sub.add_parser(
        "ortools-atsp",
        help="Solve a static ATSP path relaxation, then realize it on GTFS.")
    atsp.add_argument("--start", default="A02S", help="Platform stop id.")
    atsp.add_argument("--day", default="0", help="Service day, Mon=0 ... Sun=6.")
    atsp.add_argument("--time", default="06:00", help="HH:MM[:SS] time of day.")
    atsp.add_argument("--terminal-runs", action="store_true")
    atsp.add_argument("--run-radius", type=float, default=2500)
    atsp.add_argument("--seed-from", default=None, help="Optional route JSON initial order.")
    atsp.add_argument("--max-stations", type=int, default=0,
                      help="Smoke-test subset size; 0 means all canonical stations.")
    atsp.add_argument("--time-limit", type=float, default=60)
    atsp.add_argument("--validate", action="store_true")
    atsp.add_argument("--log-search", action="store_true")
    atsp.add_argument("--out", default="reports/optimization_runs/ortools_atsp.json")
    atsp.set_defaults(func=cmd_ortools_atsp)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
