# Subway Challenge - Codex Optimization Charter

This is the operating brief for Codex in this repository. Treat the checkout as
a clean starting point for solving the NYC Subway Challenge with optimization
methods over the real MTA schedule data.

## Mission

Use mathematical optimization, solver-backed search, and parallel computation to
find and validate a route that visits all 472 official NYC subway stations in
strictly less than 22:14:10.

The current repository best is 24:24:30. Improving that baseline is useful, but
the mission is not complete unless the full validator proves:

```text
valid=true
stations=472/472
elapsed_s < 80050
elapsed < 22:14:10
```

A tie does not count. A heuristic estimate does not count. A solver objective on
a simplified model does not count. Only the full validator decides.

## Clean-Start Rule

Start each serious optimization pass from the current checkout and current data
cache, not from remembered experiments.

Begin by orienting from:

```bash
git status --short
git ls-files
venv/bin/python -m subway_challenge.solver best
```

Untracked files under `reports/` are scratch artifacts unless deliberately
promoted. Leave unrelated dirty files alone.

## Source Of Truth

Read these files before changing solver behavior:

| Path | Purpose |
|---|---|
| `README.md` | Project overview, baseline, data pipeline, and known commands. |
| `REPRODUCIBILITY.md` | GTFS/cache hashes for exact route replay. |
| `subway_challenge/solver.py` | Validator, scorer, legal moves, and solution format. |
| `subway_challenge/search.py` | LNS and start-grid search infrastructure. |
| `subway_challenge/build_graph.py` | GTFS-derived time-expanded graph construction. |
| `subway_challenge/run_layer.py` | Legal out-of-system run transitions. |
| `subway_challenge/walk_transfers.py` | Run-distance and OSRM model. |
| `subway_challenge/stations.py` | Official 472-station identity mapping. |
| `solutions/best.json` | Current best route, trusted only after validation. |

Do not weaken the validator, station identity mapping, legal move model, or
elapsed-time scoring to make a route pass. Improve the route, not the rules.

## Optimization Abstraction

The useful abstraction is:

> Minimum elapsed-time walk in a cyclic, time-expanded, multimodal graph,
> visiting every official station label at least once, with incidental station
> coverage along paths.

This is related to time-dependent GTSP/TSP, graphical TSP, prize/coverage path
problems, set-covering path models, and column generation. Static TSP/VRP
solutions are only useful when they can be realized and validated on the actual
schedule.

## Current Frontier

The multi-day optimization run did not find a record-beating full route. The
current state is:

```text
Best validated full route:      24:24:30, 472/472
World-record target:            22:14:10, 472/472
Best record-capped partial:     22:12:30, 456/472
Prior over-cap partial:         22:14:30, 456/472
```

The strongest world-record-capped partial route found during the ignored local
experiments was:

```text
reports/optimization_runs/exact_cover_anchor_g14_station1_rpacket_7tail_min456_refined_route.json
```

It validates/replays at 456/472 stations and 22:12:30, starting `A07N@22020`
and ending `227S@101970`. Its miss set is:

```text
18, 108, 139, 141, 142, 148, 193, 194, 195,
200, 201, 202, 203, 253, 436, 437
```

Those `reports/optimization_runs` files are gitignored scratch artifacts, not a
promoted solution portfolio. If the artifact directory is absent on a fresh
checkout, reproduce or regenerate the needed pool before using these exact
paths.

## Main Findings To Preserve

The useful negative result is that the current 456/472 record-capped frontier is
not missing one cheap station insertion. The remaining stations are coupled
packets:

| Packet | Stations / role |
|---|---|
| F39 / Neptune | `253` |
| Franklin / Myrtle | `139,141,142` |
| H / Rockaway | `193-203`, especially the `193-195` and `200-203` split |
| M branch | `108` plus nearby J/M structure |
| Harlem / north Manhattan | `148`, `436`, `437` |
| Earlier route packets | G/7/Astoria/Flushing and final-tail coverage that are easily destroyed |

Important closed branches:

- Forcing `253/F39` into B23, B10, F16, D28, G14, D38, and N07 style packets
  cannot preserve the current 456/472 frontier under the record cap. The first
  feasible touch-miss model drops to 455/472 by recovering `253` while losing
  `463/464`.
- Treat `193-209` as one coupled H/Rockaway resource. Protecting large pieces of
  it with M-branch and final-tail resources is infeasible in the tested active
  pools.
- B14/B23/F36/F26/G35/G24/G22/R03 suffix resource-chain probes cannot produce a
  strict-cap suffix that preserves both M branch and A/Rockaway-family coverage.
- R08 is early enough to generate exact M+A/Rockaway resource paths, but the
  opaque suffix rows are globally destructive. Splitting those rows improves the
  relaxed basin to roughly 450/472, still below the frontier.
- G24, G22, and R03 packet-stage generators can make large local
  F39/tail/H/Rockaway/M rows, but Flushing/Astoria cannot be appended before the
  record cap and the rows are infeasible in the global relaxed master.
- The A47 H-first split experiment generated useful Astoria and F39 fragments,
  but forcing a clean Astoria fragment plus a later F39 fragment was infeasible
  under 22:14:09 in the current active pool.

The practical implication: stop tuning one giant late suffix or one isolated
station repair. The next serious model has to move the macro order earlier than
the R/Astoria/R03 knot, or generate complete macro-route candidates with packet
resources and exact handoff timing inside the pricing problem.

## Useful Tooling Now In The Repo

`subway_challenge.columns` is the main research workbench. Durable capabilities
added during the optimization run include:

- hard elapsed caps for exact-cover via `--max-total-elapsed`;
- relaxed coverage with `--uncovered-penalty-s` and `--min-covered-count`;
- extra penalties for critical miss groups via `--uncovered-penalty-groups`;
- hard and partial protected station groups via `--protect-stations` and
  `--protect-station-groups`;
- required column ids, prefixes, prefix groups, pricing kinds, and excluded ids;
- targeted exact connector pricing around required/protected rows via
  `--required-arc-top-k` and `--protected-arc-top-k`;
- bounded replacement screens: `block-replace`, `pair-replace`, and
  `chain-replace`;
- exact stage/resource column generators: `price-stage-chains` and
  `price-resource-chains`;
- `split-columns` for cutting long schedule-realized rows into smaller
  subcolumns at gateway stations.

See `OPTIMIZATION_STRATEGY.md` for how these tools fit together. Keep generated
JSONL/JSON artifacts under `reports/optimization_runs/` unless they become
small, documented, and intentionally promoted.

## Available External Capabilities

Use `.env` for local secrets and machine-specific settings. It is gitignored.
Never print secret values.

Smoke-tested capabilities:

- `NEOS_EMAIL` works for NEOS XML-RPC submissions. A tiny `milp:Cbc:AMPL` job
  was submitted and solved successfully.
- `NVIDIA_API_KEY` works against NVIDIA endpoints, and the account can see the
  active `ai-nvidia-cuopt` function.
- `OPT_ARTIFACT_DIR` is writable and defaults to `reports/optimization_runs`.
- OR-Tools is installed locally.

Current caveats:

- Python's default SSL certificate discovery may fail on this machine. Use
  `certifi` explicitly for HTTPS/XML-RPC clients.
- Local AMPL, `amplpy`, Gurobi, CPLEX, and MOSEK are not currently installed or
  configured unless `.env` is updated.

## Next Continuation Goal

Build a protected packet-state branch-and-price or macro-route generator that:

1. starts before the R/Astoria/R03 commitment point or changes the macro route
   order before that phase;
2. carries packet resources for F39, H/Rockaway, M branch, G/7/Astoria/Flushing,
   Franklin, Harlem, and final tail;
3. prices exact event-to-event handoffs during generation rather than after a
   loose master has selected rows;
4. can emit complete route candidates or compatible staged columns with shared
   resource states;
5. reconstructs every promising output into route JSON and validates it with
   `solver.py`.

Avoid another full session of blind LNS, local shortcut scans, or single-suffix
column tuning unless it tests a concrete solver-generated hypothesis.

## Validation Contract

Every meaningful candidate must be validated on the full graph:

```bash
venv/bin/python -m subway_challenge.solver validate solutions/candidate.json --record
```

Always surface the complete `RESULT` line. `--record` may update
`solutions/best.json`, but only when the candidate is valid and faster than the
current recorded best.

Solver objectives, compact graph scores, lower bounds, and route plausibility
are diagnostics. They are not evidence of success until the full validator
reports a valid 472-station route.

## Working Loop

For each optimization session:

1. Inspect repo state and validate the current best.
2. State the optimization formulation or solver hypothesis.
3. Build the smallest useful model/export/adapter.
4. Run a smoke test on a tiny instance.
5. Run a bounded experiment and write artifacts under `reports/optimization_runs`
   or another ignored `reports/` directory.
6. Reconstruct promising solver outputs into route JSON.
7. Validate routes with `solver.py`; use `--record` only for real improvements.
8. Report the best `RESULT`, solver evidence, and what the experiment ruled in
   or out.

Useful commands:

```bash
venv/bin/python -m subway_challenge.solver best
venv/bin/python -m subway_challenge.search start-grid --terminal-runs
venv/bin/python -m subway_challenge.search lns --seed-from solutions/best.json --terminal-runs --run-radius 2500
venv/bin/python -m subway_challenge.solver validate solutions/candidate.json --record
```

## Completion

Keep going until the repository contains a candidate proven by the full validator
to satisfy:

```text
valid=true
stations=472/472
elapsed_s < 80050
elapsed < 22:14:10
```

Anything else is a baseline improvement, an experiment, a lower bound, or a
useful clue. It is not the finish line.
