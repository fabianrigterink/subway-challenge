# Optimization Strategy

This repo's current best validated route is **24:24:30**. The target is a fully
validated route under **22:14:10**. A solver objective, simplified route score,
or lower bound is not enough; every serious candidate must round-trip through
`subway_challenge.solver`.

## Current Frontier

The latest optimization campaign did not find a record-beating 472-station
route. It did improve the world-record-capped partial frontier:

| Route family | Status |
|---|---|
| Best full route | valid `472/472`, `24:24:30` |
| Best record-capped partial | exact replay `456/472`, `22:12:30` |
| Prior over-cap partial | exact replay `456/472`, `22:14:30` |
| Human world record target | `472/472`, `22:14:10` |

The strongest record-capped partial misses:

```text
18, 108, 139, 141, 142, 148, 193, 194, 195,
200, 201, 202, 203, 253, 436, 437
```

The important result is not just the station count. The experiments repeatedly
showed that these misses are coupled to coverage that the partial route already
needs to preserve: G/7/Astoria/Flushing, H/Rockaway, M branch, Franklin,
Harlem/north Manhattan, Pelham/Dyre, and the final tail.

## Model Lessons

### Static TSP/VRP Is Only A Macro-Order Generator

The first OR-Tools station-level ATSP pilot produced complete legal routes, but
static cost did not predict schedule-realized elapsed time:

| Experiment | Static objective | Validated elapsed | Takeaway |
|---|---:|---:|---|
| Seeded from `solutions/best.json` | 17:00:07 | 25:57:30 | Useful macro order, weak time model. |
| Unseeded | 16:53:33 | 29:00:00 | Better static score can realize worse. |
| LNS from seeded ATSP route | n/a | 25:37:00 | LNS can repair, but not enough. |

Use ATSP/VRP engines to suggest route skeletons, not as the final model.

### Exact Timed Columns Are The Right Primitive

Column extraction and exact-cover path masters are calibrated well enough to
reproduce the 24:24:30 incumbent and to create near-record partials. The useful
column primitive is a schedule-realized path slice with:

- exact start/end stop events;
- elapsed time and mode mix;
- canonical station coverage;
- optional full path for exact route reconstruction.

Loose set-cover or proxy connector objectives are too optimistic. Exact
connector pricing, replay cuts, and final validation are mandatory.

### Hard Elapsed Caps Changed The Frontier

Adding `--max-total-elapsed` to exact-cover made record-capped diagnostics
meaningful. Under the hard cap, the best partial is 456/472 at 22:12:30. The
prior 456/472 at 22:14:30 is useful because it is only 20 seconds over the
record cap, but local prefix/end repairs failed to close that gap without losing
coverage.

### The Remaining Misses Are Packet Conflicts

The campaign repeatedly found the same structural tradeoffs:

- recovering `253/F39` tends to lose G/7/Flushing/Astoria or other mid-route
  packet coverage;
- recovering Franklin and inner H/Rockaway tends to miss exact handoffs into the
  incumbent tail;
- preserving broad `193-209` H/Rockaway coverage with M branch and final-tail
  coverage is infeasible in the tested active pools;
- B14/B23/F/G/R suffix resource chains cannot preserve both M branch and
  A/Rockaway-family coverage under the strict cap;
- R08 is early enough to generate M+A/Rockaway rows, but those rows destroy too
  much other coverage unless split, and split fragments still do not reach the
  frontier;
- G24/G22/R03 stage-chain chronologies can make large local F39/tail/H/M rows,
  but cannot append Flushing/Astoria before the cap and are infeasible in the
  global relaxed master;
- A47 H-first split rows produced clean Astoria and F39 fragments, but forcing a
  representative pair was infeasible under the record cap.

This points away from one-off local repairs and toward a model that carries
multiple branch packets as resources before committing to the R/Astoria/R03
phase.

## Durable Tooling

`subway_challenge.columns` is the main optimization workbench.

Useful master-model flags:

- `--max-total-elapsed`: hard route elapsed cap.
- `--uncovered-penalty-s`: allow relaxed coverage with a station miss penalty.
- `--min-covered-count`: turn relaxed coverage into a frontier feasibility test.
- `--uncovered-penalty-groups`: add extra penalties for critical station groups.
- `--protect-stations` and `--protect-stations-file`: require specific stations.
- `--protect-station-groups`: require at least `MIN` hits from station groups.
- `--require-column-id`, `--require-column-id-prefix`, and
  `--require-column-id-prefix-group`: force exact rows or generated families.
- `--exclude-column-id` and `--exclude-column-id-prefix`: remove brittle rows or
  families from a diagnostic pool.
- `--required-arc-top-k` and `--protected-arc-top-k`: exact-price additional
  connector arcs around important rows without widening the global graph.
- `--stop-after-first-solution`: stop large relaxed diagnostics once a feasible
  row order is found.

Useful generators and screens:

- `active-pool`: merge large JSONL libraries into a coverage-complete active
  pool.
- `exact-cover`: select and sequence timed columns with proxy or exact
  connectors.
- `block-replace`, `pair-replace`, `chain-replace`: bounded replacement screens
  around an incumbent selected-column order.
- `price-stage-chains`: generate exact columns through ordered resource groups.
- `price-resource-chains`: generate exact columns with unordered resource-state
  beam search.
- `split-columns`: split long priced rows into smaller exact subcolumns at
  gateway stations.
- `price-cuts`: generate bridge columns from failed proxy connector cuts.

Generated JSON/JSONL artifacts should stay under `reports/optimization_runs/`.
That directory is ignored; promote only small, documented, reproducible files.

## Closed Branches

Avoid rerunning these as standalone efforts unless the formulation changes.

### Incumbent Local Search

Randomized LNS, start-grid variants, shortcut scans, exact endpoint-window
replacement, shared-event splicing, and first-visit reinsertion did not improve
24:24:30. Local endpoint-fixed repairs are mostly score-neutral because elapsed
time is determined by start/end events.

### Static Or Sequencing-Free Masters

Pure set cover, static sequencing, and proxy-only connector masters produce
optimistic objectives that collapse on exact replay. They are useful for
diagnostics but not enough to claim progress.

### Single Huge Late Suffixes

Large late suffix rows can cover many current misses locally, but they delete
unique mid-late and final-tail coverage. R32/J31/H15/F39 whole-tail beams all
showed this pattern.

### Single Station Or Single Packet Inserts

Forced single-station probes from the 452/472 and 456/472 basins fell well below
the frontier. Recovering one missed station usually swaps out a different
critical packet rather than creating net coverage.

### Current Suffix Resource Chain

Starting from B14, B23, F36, F26, G35, G24, G22, or R03 within the current
over-cap skeleton does not generate a strict-cap suffix preserving both M branch
and A/Rockaway-family resources. R08 can generate such rows, but they are
globally destructive as opaque suffixes.

### Current Packet Stage Chronologies

G24/G22/R03 stage chains in both tail-first and F39-first orders fail in
complementary ways: tail-first cannot add F39; F39-first cannot add
Flushing/Astoria; neither can be forced into the global relaxed master.

## Recommended Next Plan

The next serious attempt should be a protected packet-state branch-and-price or
complete macro-route generator.

1. Move the decision point earlier than `R08/R03/G24`, or deliberately change
   the macro order before that phase.
2. Represent station packets as resources, not as optional late rewards:
   F39, H/Rockaway, M branch, Franklin, G/7/Astoria/Flushing, Harlem, Pelham,
   Dyre, and final tail.
3. Price exact event-to-event handoffs during generation. Do not rely on a loose
   master to discover connector infeasibility after selection.
4. Emit either complete macro-route candidates or staged columns that share a
   resource-state interface.
5. Use exact-cover only as a compatibility/regression check once rows are
   already event-compatible.
6. Reconstruct and validate every promising candidate with:

```bash
venv/bin/python -m subway_challenge.solver validate solutions/candidate.json --record
```

## Reproducibility Notes

- `reports/optimization_runs/` is intentionally ignored and may contain large
  scratch pools. Do not assume those artifacts exist on another checkout.
- `.env` is gitignored for NEOS/NVIDIA/solver credentials.
- `REPRODUCIBILITY.md` is the source of truth for GTFS/cache hashes.
- The exact validator, station mapping, and run model must remain unchanged
  unless the project intentionally changes the problem definition.

## Commit Hygiene

For a clean research commit:

```bash
git status --short
git diff --check
venv/bin/python -m py_compile subway_challenge/*.py
venv/bin/python -m subway_challenge.solver best
git ls-files -o --exclude-standard
```

Expected current validation baseline:

```text
RESULT valid=true stations=472/472 elapsed_s=87870 elapsed=24:24:30
```
