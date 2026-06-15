# subway-challenge
Solving the NYC Subway Challenge using Agents.

## Results

Best automated route: **24:24:30**, visiting all **472** stations (valid, in
`solutions/best.json`; check with `python -m subway_challenge.solver best`).
Exact replay depends on the GTFS/graph cache; see
[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for hashes and rebuild notes.

| Milestone | Time |
|---|---|
| Greedy (nearest-unvisited) | 30:08:00 |
| Multi-start + tail simulated annealing | 28:10:30 |
| Time-dependent LNS (ruin + regret-2 recreate + SA, iterated) | 26:42:00 |
| **+ dead-end running between terminals** | **24:45:00** |
| Focused LNS / branch-order refinement | **24:24:30** |
| — lower bound (Rural Postman, this run model) | 22:01:18 |
| — current human world record (Kate Jones, 2023, 472 stns) | 22:14:10 |

**Run model, calibrated against the record.** Running between nearby dead-end
terminals is the key second-act lever (`--terminal-runs`). The committed model
(10.8 km/h, 60 s/end) has a lower bound of **22:01:18**, *below* the human record
(22:14) — i.e. consistent with what record-holders actually do (they run hard).
A stricter "realistic" model (9 km/h, 120 s/end, 1.5 km cap) was tested but
rejected: its lower bound (22:22) *exceeds* the record, proving it forbids a
route a human achieved. So the lenient model is the defensible one. (It does let
the optimizer use a few long runs; tune `DEFAULT_PACE_MPS` / `DEFAULT_RUN_MAX_METERS`
in `walk_transfers.py` to taste — the validator and optimizer share the setting.)

The lower bound comes from reducing the network to its 143 decision nodes / 554
track segments and solving the odd-degree (Rural-Postman) matching. Geometric-
optimal structures (TSP / Euler / postman) realize *worse* on the schedule
(32-48h); only schedule-aware search (greedy, LNS) works.

The big second-act gain came from modelling **out-of-system running between
nearby dead-end terminals** (finish a line, run to a different line's terminus,
continue) — exactly what human record-holders do. This is enabled with
`--terminal-runs` (runs allowed only from line termini, keeping branching low)
plus a run-aware recreate metric.

### Solver / search

The production route search is **time-dependent LNS** (ruin a window of the
visit order → regret-2 recreate → simulated-annealing accept), with
out-of-system **runs** between dead-end terminals and a run-aware recreate
metric. Iterate by seeding from the current best.

```bash
python -m subway_challenge.solver best                          # score the best route
python -m subway_challenge.search lns --seed-from solutions/best.json \
    --terminal-runs --run-radius 2500 --seeds 6                 # iterated LNS (the optimizer)
python -m subway_challenge.search start-grid --seed-from solutions/best.json \
    --terminal-runs --run-radius 2500                           # start/time/grid portfolio
python -m subway_challenge.solver validate <file> --record      # score/record any candidate
```

`search.py` exposes `lns` and `start-grid`. Both can use the default first-visit
anchor abstraction or stricter repeated-anchor modes (`--anchor-mode all` /
`--anchor-mode revisit`), and both support `--run-mode none|terminal|all`
(`--terminal-runs` is the high-value shorthand for terminal-only runs). Many
other heuristics — multi-start, GRASP, TSP-order, postman/Euler, etc. — were
tried during development and discarded; LNS + terminal-runs is what produced
`best.json`.

### Optimization experiments

The repo also contains solver-backed experimental tooling for attacking the
current `24:24:30` plateau. These tools are diagnostics and column-generation
infrastructure; every promoted route still has to pass `solver.py`.

```bash
python -m subway_challenge.optimize ortools-atsp --help
python -m subway_challenge.columns --help
python -m subway_challenge.repair --help
python -m subway_challenge.neos_client --help
```

The main research direction is documented in
[`OPTIMIZATION_STRATEGY.md`](OPTIMIZATION_STRATEGY.md) and
[`CODEX.md`](CODEX.md). Current findings:

| Route family | Status |
|---|---|
| Best full route | valid `472/472`, `24:24:30` |
| Best near-record partial | exact replay `454/472`, `22:16:30` |
| Current target gap | cover the remaining 18 stations without destroying late-middle coverage |

The strongest current partial route misses:

```text
108, 18, 193, 194, 195, 200, 201, 202, 203,
206, 207, 208, 209, 436, 437, 446, 7, 8
```

That evidence points away from single detours or one huge late suffix. The next
useful optimization model is a protected multi-corridor column-generation /
branch-and-price neighborhood: seed from the `454/472` exact partial route,
protect several existing coverage corridors, price replacement columns with
exact connector costs and time buckets, then reconstruct and validate any full
candidate.

### Route summary notebook

`notebooks/subway_route_summary.ipynb` summarizes the curated JSON routes in
`solutions/` for blog/reporting use: ranked durations, duration by mode, common
start/end stations, and common start/end times. It uses the local
`notebooks/blogstyle.py` copy so the notebook still runs when GitHub raw-file
SSL verification fails.

## Setup

```bash
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

## Pipeline

**1. Download GTFS** (MTA NYC Subway static feed → `data/gtfs/`):

```bash
python -m subway_challenge.download_gtfs            # skips if present
python -m subway_challenge.download_gtfs --force    # re-download
```

**2. Build the time-expanded graph** (`data/graph.pkl`):

```bash
python -m subway_challenge.build_graph                       # full cyclic week (Mon–Sun)
python -m subway_challenge.build_graph --days 0 --routes 1   # quick subset for testing
python -m subway_challenge.build_graph --no-transfers        # train + wait only
```

### Graph model

* **Node** = `(stop_id, t)` — a platform-level stop (e.g. `127N`) at time `t`
  seconds within a **cyclic week** `[0, 604800)`, Monday 00:00 = 0. After-midnight
  GTFS times and Sunday-night overflow wrap modulo one week.
* **Edges** (directed, `weight` in seconds, attribute `mode`):
  * `train` — scheduled ride between consecutive stops; carries `line`, `route`,
    `trip`, `direction`.
  * `wait` — stay on a platform until its next event (cyclic).
  * `transfer` — change platforms within a complex or walk between GTFS-linked
    complexes (`transfers.txt`), honoring `min_transfer_time`; `walk=True` for
    cross-complex.
* Map a platform to its station with `build_graph.station_of(stop_id)`.

Full week ≈ 1.58M nodes / 5.87M edges (~427 MB pickle, ~30 s build).

**Scope:** Staten Island Railway (`route_id "SI"`, 21 stations) is excluded by
default — it has no track/transfer link to the subway (ferry-only). The result is
**475 stations, one strongly-connected component** via train+transfer (matches the
~472 official challenge count). Note the network is *not* connected by train alone:
the IRT, B-division, L, 7, and 42 St shuttle are separate train components stitched
together only by in-system transfers, so a valid route must use transfers.

### Canonical stations (the official 472)

`stations.py` maps each GTFS parent to its official MTA `Station ID`, collapsing
the 3 bi-level stations that GTFS splits (145 St, W 4 St-Wash Sq, Queensboro
Plaza) so "visited every station" matches the Subway Challenge's **472** exactly.

```bash
python -m subway_challenge.stations          # prints the 472 count + the 3 collapses
```

```python
from subway_challenge.stations import StationIndex
idx = StationIndex.load()                     # downloads the MTA inventory if missing
idx.resolve("A12N")                           # -> official Station ID "151"
idx.canonical_stations                        # frozenset of all 472
```

`Station ID` (472) = one physical station; `Complex ID` (424) = one fare-controlled
complex (Times Sq etc.) — we use Station ID as identity, not Complex ID.

**4. Run distances** (`walk_transfers.py`) — computed automatically on first use

Out-of-system street distances between station complexes, for running between
stations. The full complex-to-complex matrix is fetched once from the public
**OSRM** demo server (~80 batched `table` calls) and cached to
`data/walk/osrm_matrix.json`; thereafter it's free. `complex_run_adjacency`
turns it into `{complex: [(other, seconds, meters)]}`, converting street distance
to time via a running pace + per-end access penalty (see the run-model note in
Results). This needs the MTA station/entrance datasets (auto-downloaded by
`stations.py` / on demand). You don't call this module directly — the run layer
builds on it.

**5. On-demand run layer** (`run_layer.py`)

The solver queries runs live rather than baking them into the graph. A run can
start at any time, so from node `(stop, t)` you may run to a nearby complex,
arrive at `t + run_seconds`, and board the first train there — **zero static
edges**, built in ~5 s from the cached matrix.

```python
from subway_challenge.run_layer import RunLayer
runs = RunLayer.from_graph(G, radius_m=5000)
for v, weight, info in runs.neighbors(node):        # graph edges + run options, uniform
    ...                                             # info["mode"] == "run" for runs
```

`weight` is total elapsed seconds (run + wait for the boarded train). The
optimizer restricts runs to *dead-end terminals* (`--terminal-runs`) to keep
branching low — that's the high-value subset (finish a line, run to another).
