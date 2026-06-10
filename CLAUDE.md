# Subway Challenge — solver brief

Goal: find a route that **visits all 472 official NYC subway stations in minimum
elapsed time**, riding trains and walking/running between nearby stations.

## Scope — what to touch

* **Your job**: write the solver. Put search/construction logic in a new
  `subway_challenge/search.py`; write routes to `solutions/`. That's it.
* **Read-only library** (already built & validated — do NOT rebuild or edit):
  `build_graph.py`, `walk_transfers.py`, `run_layer.py`, `stations.py`,
  `download_gtfs.py`, and `data/`. Call them; don't change them.
* `solver.py` is the **validator/scorer you build against** — read it, don't
  change its solution format.
* If something in the library seems wrong, note it and work around it — don't go
  refactoring prep. The data is done; the route is the work.

## How to work a `/goal` turn (read this first)

The `/goal` evaluator only reads the transcript — it cannot run code. So **every
turn must end by running the validator and surfacing its `RESULT` line.**

```bash
venv/bin/python -m subway_challenge.solver validate solutions/candidate.json --record
# -> RESULT valid=true stations=472/472 elapsed_s=83647 elapsed=23:14:07 steps=1843 modes={...}
# -> NEW BEST elapsed_s=83647 (23:14:07), previous 24:01:10
```

Workflow each turn: (1) produce or improve a solution file, (2) validate+record,
(3) print the `RESULT` line. Exit code is 0 if valid, 1 otherwise. `--record`
updates `solutions/best.json` only when the route is valid and faster.

## The problem, precisely

* **Visit** = the route passes through a node whose station resolves (via the
  official `Station ID`) to that station. All **472** must be covered.
* **Moves** (each has a `weight` = elapsed seconds; time only ever moves forward):
  `train` (ride a hop), `wait` (stay on a platform), `transfer` (in-system, incl.
  walk between linked complexes), `run` (out-of-system, on-demand — see below).
* **Objective**: minimize total elapsed time = sum of transition weights.
* Start/end station and start time-of-day are **free** unless a goal fixes them.

## Loading and solving (Python)

```python
from subway_challenge.solver import load_problem, validate
prob = load_problem(radius_m=5000)          # ~10s: loads graph.pkl + builds run layer
G, runs = prob.G, prob.runs

for v, weight, info in runs.neighbors(node): # graph edges + run options, uniform
    ...                                      # info["mode"]=="run" for runs
```

* `prob.G` — time-expanded `networkx.DiGraph`. Node attrs: `stop` (platform, e.g.
  `127S`), `t` (seconds in a cyclic week `[0,604800)`, Mon 00:00=0), `station`
  (GTFS parent). 1.58M nodes / 5.87M edges.
* `runs.neighbors(node)` — THE expansion relation for a solver (graph + runs).
* `prob.stations.resolve(stop)` → official Station ID; `prob.canonical` = the 472.
* `prob.node_of(stop, t)` → node id; `prob.transition(u, v)` → `(weight, mode)`.

## Solution format (`solutions/*.json`)

```json
{"meta": {"radius_m": 5000, "notes": "..."},
 "path": [["127S", 28860], ["125S", 28950], "..."]}
```

`path` = ordered `[stop_id, t]` nodes. **Build paths only via `runs.neighbors()`**
so every consecutive pair is a legal transition (the validator enforces this). A
run lands on the first boardable event at the target complex; later departures
are reached by appending `wait` steps.

## Invariants the solver MUST respect

* Use `weight` for cost — never differences of raw `t` (cyclic; would read as
  negative across the Sun→Mon wrap). Total time = sum of weights.
* Only `runs.neighbors()` transitions are legal. No teleporting.
* Coverage is by official Station ID (`prob.canonical`), i.e. 472 — the 3
  bi-level GTFS splits (145 St, W 4 St, Queensboro Plaza) collapse automatically.

## Run layer — high branching, prune it

At radius 5 km a node has ~170 successors, mostly runs (~64 target complexes).
The lookup is cheap but the **search must prune runs that don't beat staying on
transit**, or it explodes. A run is only useful if it reaches an unvisited
station (or a needed transfer) faster than the rail alternative. Consider a
smaller radius (1–2 km) early, widening only where it pays.

## Suggested approach (this is where agents help)

1. Construction: greedy / nearest-unvisited-station to get *a* valid 472 route
   (seeds `best.json`). Parallel agents from different start stations / seeds.
2. Improvement: local search (2-opt / or-opt on station order, reroute segments,
   swap in runs) to cut elapsed time. Parallel agents on disjoint segments.
3. Always validate+record; keep the best.

## Data layer (already built — reference only)

| Module | Purpose |
|---|---|
| `download_gtfs.py` | MTA subway GTFS → `data/gtfs/` |
| `build_graph.py` | time-expanded graph → `data/graph.pkl` (cyclic week, directed) |
| `stations.py` | `StationIndex`: GTFS parent → official Station ID; the 472 |
| `walk_transfers.py` | OSRM run distances (cached `data/walk/osrm_matrix.json`) |
| `run_layer.py` | `RunLayer`: on-demand run successors (zero static edges) |
| `solver.py` | **solution format + validator + scorer** (build against this) |

475 stations form one strongly-connected component (SIR excluded); train-only is
NOT connected (divisions don't interline) so transfers/runs are mandatory.

## Recommended `/goal` conditions

First, get any valid route:
```
/goal solutions/best.json holds a VALID route covering 472/472 stations, proven by running `python -m subway_challenge.solver best` and surfacing its RESULT line; or stop after 15 turns
```
Then improve:
```
/goal each turn, improve solutions/best.json and run the validator with --record; stop when elapsed drops below 79200 (22h) shown in a RESULT line, or after 30 turns with no improvement
```

Enable **auto mode** so each turn runs unattended. Allowed commands are
pre-approved in `.claude/settings.json`.
