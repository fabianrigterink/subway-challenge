# Subway Challenge - Codex Optimization Charter

This is the operating brief for Codex in this repository. Treat the checkout as
a clean starting point for solving the NYC Subway Challenge with **optimization
methods** over the real MTA schedule data.

## Mission

Use mathematical optimization, solver-backed search, and parallel computation to
find and validate a route that visits all **472 official NYC subway stations**
in strictly less than **22:14:10**.

The current repository best is **24:24:30**. Improving that baseline is useful,
but the mission is not complete unless the full validator proves:

```text
valid=true
stations=472/472
elapsed_s < 80050
elapsed < 22:14:10
```

A tie does not count. A heuristic estimate does not count. A solver objective
on a simplified model does not count. Only the full validator decides.

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

## Optimization Goal

Cast the problem as solver-friendly models while preserving a path back to the
real time-expanded graph. The useful abstraction is:

> Minimum elapsed-time walk in a cyclic, time-expanded, multimodal graph,
> visiting every official station label at least once, with incidental station
> coverage along paths.

This is related to time-dependent GTSP/TSP, graphical TSP, prize/coverage path
problems, set-covering path models, and column generation. Static TSP/VRP
solutions are only useful when they can be realized and validated on the actual
schedule.

## Current Frontier

The long optimization run has not found a record-beating route. The current
state of knowledge is:

```text
Best validated full route:      24:24:30, 472/472
World-record target:            22:14:10, 472/472
Best near-record partial route: 22:16:30, 454/472
```

The best exact partial route is:

```text
reports/optimization_runs/exact_cover_late_tail_branch_rlocal_04174_exact_top64_route.json
```

Its corresponding selected-column master is:

```text
reports/optimization_runs/exact_cover_late_tail_branch_rlocal_04174_exact_top64_require1_relaxed_start21601_end90000.json
```

It validates/replays at `454/472` stations and `22:16:30`. Its missing station
set is:

```text
108, 18, 193, 194, 195, 200, 201, 202, 203,
206, 207, 208, 209, 436, 437, 446, 7, 8
```

Treat this set as the current optimization frontier. It couples Middle Village,
both Rockaway branches, Manhattan holes around Canal/59 St, Harlem/145, and
Morris Park. The difficulty is not finding a fast partial route; it is covering
these stations without destroying coverage that the current late-middle route
already carries.

## Ruled-Out Families

Do not spend another full session repeating these unless there is a genuinely
new formulation or new data:

* Randomized incumbent-family LNS and start-grid variants plateau around
  `24:24:30`.
* Local splices, shortcut scans, exact endpoint-window replacement, and
  first-visit reinsertion did not break the incumbent plateau.
* Phase variants and phase-window variants can be locally faster, but unforced
  exact masters mostly ignore them; forced variants are worse.
* Fixed `J31S@87630 -> 227S@101970` final-tail continuation plus an open bridge
  into J31 is provably stuck at `454/472`.
* Rockaway split fragments before fixed J31 cover only already-covered trunk
  stations in the useful exact solution.
* Later J/M handoffs such as `J31N@89970`, `M12S@90210`, and `J31S@90030`
  close local residual clusters but lose too much final-tail coverage.
* Big early all-residual suffixes, such as `A47S@68220 -> 227S@96720` and
  `S01N@67680 -> 227S@96720`, cover the miss set locally but collapse the
  master/replay to roughly `394-395/472` because they delete too much unique
  mid-late coverage.

Important supporting diagnostics are documented in:

```text
OPTIMIZATION_STRATEGY.md
reports/optimization_runs/README.md
```

## Next Continuation Goal

The next useful attack is a protected multi-corridor branch-and-price model:

1. Start from the `454/472` exact partial route as the incumbent/hint.
2. Protect multiple coverage corridors at once, not only a final suffix.
3. Price replacement columns for the coverage lost before any earlier suffix
   handoff, especially the 7/F/G/Brooklyn and other middle-late blocks that the
   all-residual suffixes delete.
4. Generate alternative final-tail continuations and handoff states together,
   using exact connector costs and time buckets.
5. Use stabilized set-partitioning/column-generation ideas: dual smoothing,
   resource-constrained shortest-path pricing, and dynamic time-state
   refinement.

In practical repo terms: build small exact neighborhoods that can prove or
disprove multi-column replacements around the late middle of the route. A
single huge suffix or a single Rockaway patch is no longer the right unit of
search.

## Available External Capabilities

Use `.env` for local secrets and machine-specific settings. It is gitignored.
Never print secret values.

Smoke-tested capabilities:

* `NEOS_EMAIL` works for NEOS XML-RPC submissions. A tiny `milp:Cbc:AMPL` job
  was submitted and solved successfully.
* `NVIDIA_API_KEY` works against NVIDIA endpoints, and the account can see the
  active `ai-nvidia-cuopt` function.
* `OPT_ARTIFACT_DIR` is writable and defaults to `reports/optimization_runs`.
* OR-Tools is installed locally.

Current caveats:

* Python's default SSL certificate discovery may fail on this machine. Use
  `certifi` explicitly for HTTPS/XML-RPC clients.
* Local AMPL, `amplpy`, Gurobi, CPLEX, and MOSEK are not currently installed or
  configured unless `.env` is updated.

## Preferred Architecture

Build optimization infrastructure in layers:

1. **Transition oracle.** Given a start event and target station/block, return
   earliest-arrival path, elapsed time, modes, and incidentally covered stations.
   Cache aggressively.
2. **Compact/macro model.** Work on terminal blocks, branch blocks, decision
   nodes, and time buckets before trying huge models.
3. **Solver model.** Export ATSP/GTSP/MIP/CP-SAT/column-generation pilots to
   OR-Tools, NEOS, cuOpt, SCIP/HiGHS/Gurobi, or AMPL as available.
4. **Route reconstruction.** Convert solver macro orders or selected columns
   back into concrete path JSON.
5. **Validation and promotion.** Validate with `solver.py`; promote only useful
   candidate milestones into `solutions/`.

## Search Priorities

Prefer optimization work that can change the macro route shape:

* Protected multi-corridor set-partitioning/column-generation models where
  columns are schedule-realized path segments that cover station sets.
* Resource-constrained shortest-path pricing for replacement fragments around
  the late-middle route, especially fragments that recover coverage lost by
  earlier all-residual suffixes.
* Dynamic discretization: solve on sparse time buckets, then add event times
  where the incumbent route, exact connector costs, or duals show pressure.
* Time-bucketed ATSP/GTSP over station or branch-block orders.
* cuOpt/OR-Tools/Gurobi-style matrix routing as fast macro-order generators.
* NEOS/AMPL/Pyomo experiments for compact MILP/CP formulations and lower bounds.
* Parallel LNS portfolios seeded by solver-generated macro orders.

Avoid spending a full session only on blind local polishing around the incumbent
unless it tests a concrete optimizer-generated hypothesis.

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
