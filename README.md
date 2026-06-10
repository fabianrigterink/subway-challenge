# subway-challenge
Solving the NYC Subway Challenge using Agents.

## Results

Best automated route: **24:45:00**, visiting all **472** stations (valid, in
`solutions/best.json`; check with `python -m subway_challenge.solver best`).

| Milestone | Time |
|---|---|
| Greedy (nearest-unvisited) | 30:08:00 |
| Multi-start + tail simulated annealing | 28:10:30 |
| Time-dependent LNS (ruin + regret-2 recreate + SA, iterated) | 26:42:00 |
| **+ dead-end "running" between terminals** (curated terminal-runs) | **24:45:00** |
| — proven lower bound (Rural Postman, runs allowed) | 22:01:18 |
| — current human world record (Kate Jones, 2023, 472 stns) | 22:14:10 |

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

```bash
python -m subway_challenge.solver best                          # score the best route
python -m subway_challenge.search greedy                        # baseline construction
python -m subway_challenge.search lns --seed-from solutions/best.json \
    --splits 0.2,0.35,0.5 --seeds 4                             # iterated LNS (the optimizer)
```

Other `search.py` commands: `multi`, `grasp`, `optimize-tail`, `tsp`, `sweep`,
`portfolio`, `postman`. Validate/record any candidate with
`solver.py validate <file> --record`.

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

**4. Out-of-system walk/run transfers** (`walk_transfers.py`)

Links between nearby *complexes* not connected by GTFS, for walking/running
between stations. Default backend is **hosted OSRM** (free, the public demo):
the full complex-to-complex street-distance matrix is fetched once (~80 batched
`table` calls, cached to `data/walk/osrm_matrix.json`), then any radius is
derived instantly. Street distance is symmetrized (one-way bias) and converted
to time with a **running pace** plus a per-end access penalty — we use OSRM's
*distance*, not the driving profile's speed.

```bash
python -m subway_challenge.walk_transfers --backend osrm --radius 5000   # -> data/walk/walk_transfers.csv
python -m subway_challenge.walk_transfers --backend osrm --radius 1000 --pace 3.0 --access-penalty 60
```

Optional precise short-hop backend — **Google**, entrance-to-entrance (MTA
entrances dataset, `computeRouteMatrix`/`WALK`, cached, ~$8 at radius 400 m):

```bash
python -m subway_challenge.walk_transfers --estimate                     # price first, no calls
export GOOGLE_MAPS_API_KEY=...
python -m subway_challenge.walk_transfers --backend google --radius 400 --k 3
```

Fold links into the graph as `walk=True` transfers:

```bash
python -m subway_challenge.build_graph --walk-transfers data/walk/run_1km.csv
```

**Graph-density caveat:** the builder materializes one transfer edge per source
train-event per link, so baking runs in scales edges fast (~22M extra at 1 km,
**~261M at 5 km** — infeasible). Prefer the on-demand layer below; only bake in
short walks (≤~1 km) if you need a self-contained static graph. Also note the
OSRM *driving* profile misses pedestrian-only passages (e.g. Park Place↔Fulton
St reads ~360 m vs the real ~0 m); the Google/entrances path captures those.

**5. On-demand run layer** (`run_layer.py`) — *recommended for the solver*

Rather than baking runs into the graph, the solver queries them live. A run can
start at any time, so from node `(stop, t)` you may run to any nearby complex,
arrive at `t + run_seconds`, and board the first train there. Full 5 km coverage,
**zero static edges**, built in ~5 s from the cached OSRM matrix.

```python
from subway_challenge.run_layer import RunLayer
runs = RunLayer.from_graph(G, radius_m=5000)        # pace/access_penalty configurable
for v, weight, info in runs.neighbors(node):        # graph edges + run options, uniform
    ...                                             # info["mode"] == "run" for runs
```

`weight` is total elapsed seconds (run + wait for the boarded train). Branching
is high at 5 km (~64 run targets/complex); the solver should prune runs that
don't beat the subway alternative.
