"""Solution format, validator, and scorer for the Subway Challenge.

This is the harness a solver (or `/goal` loop) builds against. It does NOT solve
the problem; it *verifies and scores* a candidate route so the result lands in
the transcript for the `/goal` evaluator to read.

Solution file (JSON)::

    {
      "meta": {"radius_m": 5000, "notes": "..."},
      "path": [["127S", 28860], ["125S", 28950], ...]   # [stop_id, t_seconds_in_week]
    }

The ``path`` is the ordered list of time-expanded nodes the route passes
through. Every consecutive pair must be a legal transition -- either a graph edge
(train / wait / transfer) or a run from the on-demand :class:`RunLayer`. The
route is VALID iff it covers all 472 official stations; its score is the total
elapsed time = sum of transition weights.

CLI::

    python -m subway_challenge.solver validate solutions/candidate.json --record
    python -m subway_challenge.solver best
    python -m subway_challenge.solver selftest
"""
from __future__ import annotations

import argparse
import bisect
import json
import pickle
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from .build_graph import WEEK
from .run_layer import RunLayer
from .stations import StationIndex

GRAPH_PKL = Path("data/graph.pkl")
BEST = Path("solutions/best.json")
DEFAULT_RADIUS_M = 5000


def hms(seconds: int) -> str:
    seconds = int(seconds)
    return f"{seconds // 3600:02d}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


@dataclass
class Problem:
    G: object
    runs: RunLayer
    stations: StationIndex

    @property
    def canonical(self) -> frozenset:
        return self.stations.canonical_stations

    def node_of(self, stop: str, t: int):
        """Node id for (stop, t), or None if no such event exists."""
        idx = self.runs.platform_index.get(stop)
        if not idx:
            return None
        times, ids = idx
        i = bisect.bisect_left(times, int(t))
        return ids[i] if i < len(times) and times[i] == int(t) else None

    def transition(self, u, v):
        """(weight, mode) for a legal u->v transition, or (None, None)."""
        if self.G.has_edge(u, v):
            e = self.G[u][v]
            return e["weight"], e["mode"]
        best = None
        for rv, w, _info in self.runs.run_successors(u):
            if rv == v and (best is None or w < best):
                best = w
        return (best, "run") if best is not None else (None, None)


def load_problem(graph_pkl=GRAPH_PKL, radius_m=DEFAULT_RADIUS_M) -> Problem:
    with open(graph_pkl, "rb") as f:
        G = pickle.load(f)
    runs = RunLayer.from_graph(G, radius_m=radius_m)
    return Problem(G, runs, StationIndex.load())


def validate(prob: Problem, path: list) -> dict:
    """Validate and score a path of [stop, t] pairs."""
    if len(path) < 2:
        return {"valid": False, "reason": "path has fewer than 2 nodes"}

    nodes = []
    for entry in path:
        stop, t = entry[0], int(entry[1])
        nid = prob.node_of(stop, t)
        if nid is None:
            return {"valid": False, "reason": f"node {stop}@{t} not in graph"}
        nodes.append(nid)

    total = 0.0
    modes = Counter()
    for u, v in zip(nodes, nodes[1:]):
        w, mode = prob.transition(u, v)
        if w is None:
            du, dv = prob.G.nodes[u], prob.G.nodes[v]
            return {"valid": False,
                    "reason": f"illegal transition {du['stop']}@{int(du['t'])} -> "
                              f"{dv['stop']}@{int(dv['t'])}"}
        total += w
        modes[mode] += 1

    covered = {prob.stations.resolve(prob.G.nodes[n]["stop"]) for n in nodes} & prob.canonical
    missing = prob.canonical - covered
    return {
        "valid": not missing,
        "covered": len(covered),
        "total_stations": len(prob.canonical),
        "elapsed_s": int(total),
        "steps": len(nodes),
        "modes": dict(modes),
        "missing_sample": sorted(missing)[:8],
        "start": list(path[0]),
        "end": list(path[-1]),
    }


def result_line(r: dict) -> str:
    if "reason" in r and not r.get("valid"):
        return f"RESULT valid=false reason=\"{r['reason']}\""
    return (f"RESULT valid={str(r['valid']).lower()} "
            f"stations={r['covered']}/{r['total_stations']} "
            f"elapsed_s={r['elapsed_s']} elapsed={hms(r['elapsed_s'])} "
            f"steps={r['steps']} modes={r['modes']}")


def _load_solution(path_json: Path):
    data = json.loads(Path(path_json).read_text())
    return data.get("path", data), data.get("meta", {})


def cmd_validate(args) -> int:
    path, meta = _load_solution(args.file)
    prob = load_problem(radius_m=meta.get("radius_m", args.radius))
    r = validate(prob, path)
    print(result_line(r))
    if not r["valid"] and r.get("missing_sample"):
        print(f"  missing {r['total_stations']-r['covered']} stations, e.g. {r['missing_sample']}")

    if args.record and r["valid"]:
        prev = json.loads(BEST.read_text())["meta"]["elapsed_s"] if BEST.exists() else None
        if prev is None or r["elapsed_s"] < prev:
            BEST.parent.mkdir(parents=True, exist_ok=True)
            BEST.write_text(json.dumps({"meta": {**meta, "elapsed_s": r["elapsed_s"]},
                                        "path": path}))
            print(f"NEW BEST elapsed_s={r['elapsed_s']} ({hms(r['elapsed_s'])})"
                  + (f", previous {hms(prev)}" if prev else ""))
        else:
            print(f"not improved (best stays {hms(prev)})")
    return 0 if r["valid"] else 1


def cmd_best(args) -> int:
    if not BEST.exists():
        print("RESULT valid=false reason=\"no best solution yet\"")
        return 1
    path, meta = _load_solution(BEST)
    prob = load_problem(radius_m=meta.get("radius_m", args.radius))
    print(result_line(validate(prob, path)))
    return 0


def cmd_selftest(args) -> int:
    """Build a short legal path (train hops + a run) straight from the graph and
    validate it -- proves the harness end-to-end without a full tour."""
    prob = load_problem(radius_m=args.radius)
    G = prob.G
    start = next(n for n, d in G.nodes(data=True)
                 if 9 * 3600 <= d["t"] <= 9.05 * 3600
                 and any(G[n][v]["mode"] == "train" for v in G.successors(n)))
    nodes = [start]
    cur = start
    for _ in range(5):                       # follow train edges
        nxt = next((v for v in G.successors(cur) if G[cur][v]["mode"] == "train"), None)
        if nxt is None:
            break
        nodes.append(nxt)
        cur = nxt
    run = next(iter(prob.runs.run_successors(cur)), None)  # take one run
    if run:
        nodes.append(run[0])
    path = [[G.nodes[n]["stop"], int(G.nodes[n]["t"])] for n in nodes]
    r = validate(prob, path)
    print("selftest path:", " -> ".join(f"{s}@{t}" for s, t in path))
    print(result_line(r))
    ok = r.get("modes", {}) and "reason" not in r
    print("transitions all legal:", ok)
    return 0 if ok else 1


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Validate/score Subway Challenge solutions.")
    p.add_argument("--radius", type=int, default=DEFAULT_RADIUS_M,
                   help="Run-layer radius (m) for legality of run moves.")
    sub = p.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="Validate and score a solution file.")
    v.add_argument("file")
    v.add_argument("--record", action="store_true", help="Update solutions/best.json if better.")
    v.set_defaults(func=cmd_validate)
    sub.add_parser("best", help="Show the current best solution's score.").set_defaults(func=cmd_best)
    sub.add_parser("selftest", help="Self-check the validator on a synthetic path.").set_defaults(func=cmd_selftest)
    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
