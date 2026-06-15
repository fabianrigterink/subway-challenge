# Optimization Strategy

This repo's current best validated route is **24:24:30**. The target is a
fully validated route under **22:14:10**. The solver objective, a simplified
route score, or a lower bound is not enough; every serious candidate must
round-trip through `subway_challenge.solver`.

## What The First Solver Pilot Showed

The first OR-Tools station-level ATSP pilot can produce complete, valid routes,
but the abstraction is too static:

| Experiment | Static objective | Validated elapsed | Takeaway |
|---|---:|---:|---|
| Seeded from `solutions/best.json` | 17:00:07 | 25:57:30 | Useful macro order, weak time model. |
| Unseeded | 16:53:33 | 29:00:00 | Better static score can realize worse. |
| LNS from seeded ATSP route | n/a | 25:37:00 | LNS can repair, but not enough. |

Conclusion: use TSP/VRP engines as **macro-order generators**, not as the final
model. The cost model has to know departure time, wait phase, and incidental
coverage.

## What The First Column Master Showed

The first seed-column extractor produced **539** schedule-realized route
fragments from the curated `solutions/*.json` portfolio. Their union covers all
472 stations.

The sequencing-free CP-SAT cover master solved optimally on this seed pool:

| Cover variant | Selected columns | Raw covered-fragment time |
|---|---:|---:|
| no column penalty | 39 | 18:10:30 |
| 10 min/column penalty | 32 | 18:41:30 |
| 30 min/column penalty | 31 | 18:52:30 |

Those numbers are optimistic because they do not pay real ordering, connector,
or timetable phase costs. Sequencing the selected columns as an open ATSP and
realizing the resulting station-anchor order gave:

| Cover variant | Connector objective | Validated elapsed | Result |
|---|---:|---:|---|
| 30 min/column penalty, static connectors | 22:28:38 | 26:17:00 | valid 472/472 |
| 30 min/column penalty, dynamic connectors | 24:19:00 | 26:17:00 | same realized route |
| 10 min/column penalty, static connectors | 23:44:02 | 28:55:00 | valid 472/472 |
| no penalty, static connectors | 23:28:35 | 29:03:30 | valid 472/472 |
| LNS repair from best column sequence | n/a | 26:30:00 | repair did not help |

Conclusion: the seed pool contains strong local fragments, but a pure
set-cover-then-sequence decomposition is too loose. The next column-generation
step should make sequencing part of the master, or add connector/phase cuts
while selecting columns.

The first combined path-cover master selected columns and connector arcs in the
same CP-SAT model, using sparse top-k connector neighborhoods:

| Path-cover variant | Selected columns | Master objective | Validated elapsed |
|---|---:|---:|---:|
| top-12 connectors, 10 min/column penalty | 36 | 29:11:34 | 26:29:00 |
| top-24 connectors, 10 min/column penalty | 36 | 28:59:34 | 25:51:30 |
| top-48 connectors, 10 min/column penalty | 34 | 28:39:36 | 26:44:00 |
| LNS repair from top-24 route | n/a | n/a | 25:51:30 |

This is better than set-cover-then-sequence, but still behind the 24:24:30
incumbent. Wider static connector neighborhoods improved the master objective
without reliably improving realized schedule time. The next step should attach
time buckets to columns/connectors so the master cannot choose locally cheap
fragments that line up poorly in the actual timetable.

The exact-time path-cover master then chained actual timed column slices instead
of discarding them back to station anchors:

| Exact-cover variant | Columns in pool | Selected columns | Master / validated elapsed |
|---|---:|---:|---:|
| coarse pool, top-16 phase arcs | 539 | 43 | 24:24:30 |
| coarse pool, top-32 phase arcs | 539 | 39 | 24:24:30 |
| coarse pool, top-64 phase arcs | 539 | 38 | 24:24:30 |
| fine-12 pool, top-64 phase arcs | 863 | 63 | 24:24:30 |
| fine-8 pool, top-48 phase arcs plus replay cuts | 1250 | 85 | 26:30:30 |
| fine-12 plus 2149 priced terminal columns, top-16 | 3012 | 60 | 24:24:30 |
| fine-12 plus 2149 priced terminal columns, top-32 plus replay cut | 3012 | 60 | 24:24:30 |
| incumbent all-target pricing, top-8 | 13734 | 61 | 24:24:30 |
| 101S-basin all-target pricing, top-8 plus replay cut | 13872 | 62 | 24:24:30 |
| expensive-window pricing, top-8 plus 8 replay cuts | 7070 | 64 | 24:24:30 |
| expensive-window pricing, slack-filtered top-8 | 7070 | 65 | 25:15:30 |
| expensive-window pricing, lazy exact-refined top-8 | 7070 | 66 | 24:24:30 |
| expensive-window pricing, lazy exact-refined top-12 | 7070 | 75 | no valid replay |
| failed-connector bridge pricing, top-12 | 7237 | 62 | 24:24:30 |
| failed-connector bridge pricing, top-16 | 7237 | 66 | no valid replay |
| LP dual target pricing, top-8 | 7636 | 78 | 31:28:00 |
| LP dual target pricing, top-12 | 7636 | 60 | no valid replay |
| LP dual beam pricing, top-12 | 7332 | 63 | no valid replay |
| 900-column active pool, exact top-12 | 900 | 62 | 24:24:30 |
| 900-column active pool, exact top-24 | 900 | 62 | 24:24:30 |
| 1500-column active pool, exact top-12 | 1500 | 69 | 24:24:30 |
| 1500-column active pool, exact top-24 | 1500 | 66 | 24:24:30 |
| connector-priced 1500 pool, exact top-12 | 1500 | 71 | 24:24:30 |
| connector-priced 1500 pool, exact top-24 | 1500 | 66 | 24:24:30 |
| connector-priced 2000 pool, exact top-24 forced | 2000 | 66 | 24:24:30 |
| window-event 2000 pool, exact top-24 forced | 2000 | 70 | 24:24:30 |
| LNS from window-event exact route | n/a | n/a | 25:10:00 |
| alternate-basin morning start-grid | n/a | n/a | 25:10:00 |
| incumbent high-freeze LNS seeds 8000-8002 | n/a | n/a | 25:16:00 |
| revisit-anchor LNS smoke | n/a | n/a | 27:05:00 |

The top-64 coarse run reconstructs the incumbent exactly. The exact-time model
is therefore calibrated: the master objective and validator agree. But widening
and refining columns did not find an improvement from the existing route
portfolio. A first terminal-target pricing pass added 2149 new graph-generated
columns, and all-station pricing from both the incumbent and the 101S alternate
basin generated more than 10k new columns per pass, but the exact master still
tied the incumbent. Expensive-window pricing produced below-incumbent proxy
solutions, including a 24:16:30 proxy objective, but those depended on
exact-infeasible platform connectors. Lazy exact refinement now catches selected
infeasible arcs before replay; it can validate another 24:24:30 tie, but the
top-12 below-incumbent proxy route still collapses after exact cuts. Pricing
167 feasible bridge columns from failed connector cuts makes the top-12 master
validate cleanly, but still only at the incumbent tie. The next
column-generation step must price **targeted negative-reduced-cost fragments
with exact connector feasibility**, not broad fragments from sampled route
points.

The first LP-guided column-generation pass is now implemented. A GLOP
set-cover relaxation over `seed_columns_windows_cut_priced.jsonl` gives a loose
14:34:08 lower proxy and 235 positive station duals. A direct high-dual target
pricer added 399 negative-reduced-cost columns, and a bounded beam pricer over
the time-expanded graph added 95 more clustered reward columns with stronger
median reduced cost. These pricing subproblems are producing plausible columns,
but the exact-cover master still falls back to worse legal routes or no replay
route. The bottleneck has moved from "can we generate new columns?" to "can the
master reason about exact connector feasibility before it spends effort on
proxy-infeasible chains?"

That exact-connector step is now practical for active pools. The `active-pool`
command merges the base/window/cut/dual/beam libraries, forces in known
selected incumbent and near-incumbent columns, ranks the remainder by dual
reduced cost, and repairs station coverage. Exact connector pricing now prunes
node-to-node Dijkstra searches by the target event time, which turns
top-k exact arc pricing from an open-ended wait into a repeatable run. Active
pools of 900 and 1500 columns with top-12/top-24 exact connector neighborhoods
all solved to optimality and validated, but all tied **24:24:30**. This is a
stronger negative than the proxy experiments: within the current priced column
library, exact recombination is not enough. The next improvement probably needs
new columns with different endpoints/departure phases, or duals from a
sequenced master rather than the loose station-cover LP.

Connector-aware pricing now targets the expensive handoffs in exact active-pool
solutions directly. It mines adjacent selected-column connectors, starts from
the source column's end or start event, and generates downstream-window columns
scored by dual reward plus connector-time credit. The first broad pass generated
1613 connector-priced columns; exact-cover selected 6-10 of them in later
solutions, so they are useful to the master. However, the best validated result
still tied **24:24:30**. Larger active pools also exposed a sparse-arc issue:
with many endpoint variants, top-k connector neighborhoods can omit the known
incumbent backbone. `exact-cover --force-solution-arcs` now exact-prices and
adds adjacent arcs from known selected-column solutions before solving, which
restores feasibility without returning to proxy connectors.

Multi-column event-window pricing was also tested. It mines windows of several
selected columns, prices exact paths from the source window start/end event to
the target window end/next-start event, and scores them by dual reward plus
window-time credit. The first high-score pass generated 199 exact-event
replacement columns with median coverage of 22 stations, but the 2000-column
active exact master selected none of them and again tied **24:24:30**. A
10-minute LNS repair seeded from that exact route produced only **25:10:00**.
The plateau is therefore not just one bad local connector: the current pricing
families are not generating an alternate station order/phase basin that can
beat the incumbent.

The latest alternate-basin route-search passes also failed to escape. A
morning start-grid over `A02S`, `726N`, `101S`, `R45N`, and `A09N` evaluated
200 rotated-seed configurations before interruption and topped out at
**25:10:00**. A 10-minute high-freeze incumbent LNS exploitation pass over
splits 0.65-0.9 and seeds 8000-8002 topped out at **25:16:00**. These are not
worth adding to the main solution portfolio. Future route-search effort should
either allocate substantially longer time to already-known promising basins, or
change the neighborhoods/anchor representation rather than reusing the same
rotated incumbent anchor order.

A first attempt at richer route-shape anchors was also negative. `search.py`
now supports `--anchor-mode first|all|revisit` and
`--no-skip-visited-anchors`, but a smoke LNS run using 573 revisit anchors from
the incumbent validated at only **27:05:00**. Plain repeated station anchors
overconstrain the realization without preserving the original platform/time
phase. If this branch continues, the richer representation should be
platform/time-aware segment anchors or branch-window moves, not duplicate
station IDs.

Local exact repair around the incumbent also failed to move the score.
`subway_challenge.repair` now has reproducible probes for exact route-window
replacement, shared-event route splicing, and exact-prefix anchor reinsertion:

| Probe | Scope | Best valid result | Lesson |
|---|---:|---:|---|
| Exact endpoint windows | 401 prioritized node windows, all-run 5000 m | 24:24:30 | fixed event endpoints have fixed elapsed time |
| Shared-event splicing | 16605 opportunities across 13 solution routes | 24:24:30 | route diversity did not produce a below-baseline prefix/suffix splice |
| Anchor reinsertion from best | 1314 feasible first-visit reinsertions | 24:24:30 | small phase-changing local windows remain on the incumbent plateau |
| Anchor reinsertion from tied 670-node route | 822 feasible first-visit reinsertions | 24:24:30 | alternate tied geometry did not expose an easy local move |

The important modeling correction is that a time-expanded path's score is
determined by its start and end event timestamps. Any local repair that fixes
both endpoints is almost guaranteed to be score-neutral. Useful repair/pricing
must alter endpoint phase, start phase, or the sequence of branch coverage
events, not only find a shorter path between exact events.

The first phase-altering column generator is now implemented as
`subway_challenge.columns price-phase-variants`. It replays useful column
anchor patterns from alternate departure events and tags the generated rows with
`pricing_phase_variant`. A calibration constraint was also added to exact-cover:
`--require-pricing-kind` plus `--min-required-pricing`.

| Probe | Scope | Best valid result | Lesson |
|---|---:|---:|---|
| Phase variants from top 300 active columns | 1600 variants, 828 locally faster | 24:24:30 | exact-cover selected 0 variants |
| Phase variants from 70 selected incumbent columns | 779 variants, 192 locally faster | 24:24:30 | repaired proxy route tied but selected 0 variants |
| Require 1 phase variant, proxy arcs | 1900-column pool | no valid route | best proxy solution was 24:57:30 before UNKNOWN |
| Require 1 phase variant, exact arcs top-4 | 1200-column pool, 3614 exact arcs | 25:55:00 | forced phase substitutions are connectable but expensive |
| Phase-window variants from selected incumbent windows | 2200 variants, 347 locally faster | 24:24:30 | exact-cover tied with 639 steps but selected 0 phase windows |
| Require phase-window variant, proxy arcs | 2600-column pool | 25:33:30 | moving blocks together is better than isolated variants but still dominated |
| Large phase-corridor variants | 1704 variants, median 95 stations covered | 24:24:30 | unforced exact-cover tied and selected 0 corridors |
| Require one large phase corridor | 2400-column pool | 24:47:00 | best forced phase-shifted result so far |
| All-run LNS from forced corridor | 360 s, seed at A02S@20430 | 25:51:30 | first-visit anchors lose corridor structure |
| Start bucket 05:30-06:00 | 18 eligible starts | 24:47:00 | naturally selects the large early corridor |
| Start bucket 06:00-06:30 | 142 eligible starts | 24:24:30 | incumbent phase, selects 0 corridors |
| Start buckets 05:00-05:30 / 06:30-07:00 | 22 / 30 eligible starts | infeasible | active pool lacks coherent chains |
| Incumbent start + record-ish end bucket | 142 starts, 141 ends | infeasible | no record-compatible chain in current pool |
| Record-end targeted pricing | 597 new columns, 341 eligible record-band ends | infeasible | final bucket alone is not enough |
| Middle-to-record bridge pricing | 401 generated bridge columns, 742 record-band ends in full pool | infeasible | even top-64 connector graph cannot cover all stations in record bucket |
| Relaxed record-end master | 1,000,000 s uncovered penalty | 447/472 covered | misses compact late cluster around J/M and A/H/Rockaway |
| Missing-station connector pricing | 1659 repair columns from relaxed skeleton | unknown hard solve | repairs part of J/M miss set but not A/H/Rockaway sequencing |
| Cluster-corridor station pricing | 100 focused corridors, 24/25 miss-station union | infeasible | unforced relaxed master still selects old final phase |
| Forced anchor corridors v2 | exact ordered chains ending at 441 by 99120 | 430/472 relaxed | improves forced case but sacrifices final-phase stations |
| Forced anchor corridors v3 | richer final-tail chains, best forced row ends 102720 | 446/472 relaxed | nearly matches old relaxed skeleton but misses M/R and A/H local gaps |
| Pre-final repair into old final start | sampled exact chains into M04N@85050 | no useful candidate | cannot simply repair A/H/J/M before old final window |
| Late-tail beam pricer | exact prize-collecting beam from F39S@77580 to 441 | 440/472 relaxed unforced | generated strong A/H columns, but master ignored them unless forced |
| Split plus M-branch continuation | A/H/J prefix to J31S@87630, then M-branch-aware continuation | 452/472 forced relaxed, optimal | first relaxed improvement over the old 447/472 record-bucket diagnostic |
| Branch/R-local continuation | fixed J31S@87630 handoff plus wider branch/R-local final continuation | 455/472 forced relaxed, optimal | best proxy soft master so far, but replay cuts collapse because selected connectors are proxy-infeasible |
| Exact branch/R-local master | same 04174 pool with exact top-12 connector arcs | 453/472 replayed route, 22:16:30 | near-record partial route; remaining task is 19-station coverage without losing exact feasibility |

This confirms that endpoint phase needs to be modeled jointly with sequencing.
Adding phase-shifted fragments to an otherwise event-fixed path-cover model
does not help: the connector penalty cancels local savings unless surrounding
columns move with the phase too. Moving adjacent columns as one phase-window
block helps directionally, but the unconstrained master still ignores those
blocks. Large 10-20-column corridors are better again, reaching **24:47:00**
when forced, but are still dominated unless the master is required to select
one. Phase therefore needs to become a state dimension of the master rather than
optional replacement-column metadata.

The first explicit bucket-state hook now exists: exact-cover accepts
`--start-time-window` and `--end-time-window`. Bucket sweeps show the current
active pool supports only two coherent start phases: the incumbent 06:00-06:30
bucket at **24:24:30** and an early 05:30-06:00 corridor bucket at **24:47:00**.
The 05:00-05:30 and 06:30-07:00 buckets are infeasible, and the incumbent start
bucket cannot pair with a record-compatible end bucket. The next generation step
should therefore be bucket-conditioned pricing: generate columns specifically
for missing/infeasible bucket corridors and record-compatible end buckets.
Filtering generation to the record-compatible end bucket increased eligible end
columns from 141 to 341 but the incumbent-start/record-end master remained
infeasible, so the next target is the middle corridor, not just the final end
state.

The first middle-to-record bridge pass generated 401 timed fragments that start
mid-route and end in the record-compatible bucket. A full 3001-column, top-64
connector master was still infeasible for the incumbent start bucket plus
90000-102000 end bucket. A soft-coverage diagnostic covered **447/472** and
repeatedly missed the same late cluster:

```text
93-96, 108, 185-209
```

Artificial dual pricing against those missing stations produced 1659 connector
repair columns, but the repaired hard master timed out and relaxed diagnostics
continued to miss variants of the same cluster. This points to a sequencing
conflict, not a lack of single-station target columns. The next pricer should
therefore generate exact replacement corridors spanning the A/H/Rockaway
cluster and J/M tail together, with exact handoff costs into the final
record-compatible phase.

The first version of that pricer now exists as
`price-cluster-corridors`. Station-destination mode generated focused
late-cluster corridors whose union covered **24/25** repeated miss stations.
Manual ordered anchor corridors can physically run from the late F segment
through A/H, J/M, and final-tail gateways to `441` by **99120** week-seconds.
Forcing the best v2 anchor corridor improved the forced relaxed model from
**421/472** to **430/472**, but the unforced relaxed model still picks the old
final phase window at **447/472**. The remaining problem is now sharper: a
single replacement corridor must preserve more of the final phase coverage, or
the master needs compatible late-corridor decomposition so it can combine A/H
repair with final-tail coverage instead of choosing one.

Richer v3 anchor corridors explicitly reinserted the missing final-tail
branches. The best strict record-bucket v3 row ends at **100320** and reduces
old-final misses from 31 to 23; the best near-record row ends at **102720** and
reduces old-final misses to 15. When forced with a looser `90000-103000` end
bucket, v3 reaches **446/472**, almost matching the old relaxed skeleton while
covering the A/H/J/M repair cluster. But hard full-cover remains infeasible,
and exact sampled attempts to repair A/H/J/M before the old final window start
`M04N@85050` found no useful chain. The late subproblem is no longer "find any
A/H repair"; it is an ordering/resource-constrained path problem over the
remaining late stations:

```text
102-114, 7-19, 139/141/142, 193-195, 200-203
```

This should be solved as a dedicated late-tail RCSP/prize-collecting path from
around `F39S@77580` to the final terminals with a hard end-time cap, or by a
master that can combine multiple compatible late replacement corridors.

The first dedicated late-tail beam pricer now exists. It successfully generated
exact columns under a `102000` end cap, including A/H-first rows with **39-41**
rewarded target hits and **22-25** A/H hits. But the master behavior is
instructive: forced A/H-first beam rows cover only **411/472**, while the
unforced relaxed model ignores beam rows and reaches **440/472** via another
late phase window. So the beam subproblem is viable, but the reward structure is
not yet aligned with the master. The next pricing iteration needs either higher
explicit rewards for final-tail stations that beam rows drop
(`360-367`, `416-431`, `433-477`, `59-70`) or shorter late-tail subcolumns with
compatibility states so exact-cover can combine A/H repair and final-tail
coverage.

The shorter-subcolumn direction produced the first real relaxed-master
improvement. `split-columns` now cuts exact priced paths at gateway stations.
Splitting A/H-first and balanced beam rows alone reached **444/472**, but the
useful pattern emerged after pricing a continuation from the exact split
endpoint `J31S@87630`. A final-tail continuation preserved tail coverage but
skipped the Middle Village M branch; a J/M-first continuation still missed the
branch; an M-branch-first continuation generated rows that hit `108-114` and
gave an **OPTIMAL 452/472** relaxed diagnostic under the `21601-23400` start
bucket and `90000-102000` end bucket. The selected late repair chain is:
`A52S@81360 -> H06S@83460 -> H11S@84030 -> J31S@87630 -> 227S@101970`.
The remaining uncovered set is:

```text
108, 139, 141, 142, 18, 193-195, 200-203, 436-437, 445-446, 475-477, 8
```

This is still not a valid full route, but it changes the next pricing problem
from "repair the whole late tail" to "preserve the M-branch continuation while
adding a small residual repair for Franklin shuttle, Rockaway gaps, and final
tail fragments."

Follow-up residual probes confirmed that these misses are coupled, not
independent. Forcing `108` in the M-branch continuation drops the relaxed cover
to **446/472**. Forcing the only generated continuation that contains both the
cheap M branch (`109-114`) and Franklin (`139/141/142`) drops to **406/472**.
Forcing a Rockaway prefix that covers `200-203` before the same J31 handoff
drops to **438/472**. The fixed `J31S@87630` handoff is therefore now the
limiting state. The next promising move is to generate alternative nearby
handoff states around J/M (`J31`, `M12`, `M20`) and Rockaway (`A52`, `A63`,
`H12`, `H11`), then price compatibility into the M-branch continuation instead
of holding J31 fixed.

That moved-handoff hypothesis was tested next. A prefix pricer from
`A52S@81360` generated 160 alternate handoffs to `J31`, `M12`, `M18`, and
`M20`, including:

```text
A52S@81360 -> J31N@89970 covering 200-203 and 109-114
A52S@81360 -> M12S@90210 covering 193-195, 200-203, and 113-114
A52S@81360 -> J31S@90030 covering 200-203 and 139/141/142
```

Compatible continuations from those later handoff events were then priced and
fed back to the relaxed master. The unrestricted merged master selected none of
the moved-handoff rows and covered **447/472**. Focused forced-pair tests were
also negative: `J31N@89970` reached **404/472**, `M12S@90210` reached
**438/472**, and `J31S@90030` reached **421/472**. The M12 pair was rerun with
top-64 connector neighborhoods and a 900-second solve budget; it proved
**OPTIMAL at 438/472**, confirming that this moved handoff is genuinely
dominated. The conclusion is sharper: moving the handoff later closes local
residual clusters but consumes too much of the final-tail time budget. The next
state space should support earlier or multi-handoff decompositions: Rockaway
repair, J/M or M-branch continuation, and final-tail preservation as separate
compatible event-fixed stages.

A fixed-handoff continuation variant then improved the soft master again. A
wider branch/R-local continuation from `J31S@87630` to `227S@101970`, anchored
on `109, 98, 475, 360, 445, 416, final:441`, raises the forced relaxed
record-bucket diagnostic to **455/472**. The remaining miss set is:

```text
108, 139, 141, 142, 18, 193-195, 200-203, 436-437, 446, 7, 8
```

Ranking all compatible `J31S@87630 -> 227S` continuation rows showed this row
is the best already in the pool. A deeper focused continuation pass generated
235 more rows, but the best fixed-skeleton coverage was only **446/472**:
adding `436/437/446` or `7/8/18` still costs too much M-branch/final-tail
coverage.

The important caveat is connector feasibility. Direct replay of the 455
selected order fails on a proxy connector
`priced_connector:08586 A48N@42750 -> priced_dual:07340 A50N@43110`. A replay
cut loop quickly drops the saved soft master from **455/472** to **432/472**
while exposing additional infeasible proxy arcs. Therefore the current best
late-tail row is useful, but the next model improvement should be exact
connector pricing or a smaller exact-arc active pool around this skeleton, not
more proxy-only continuation pricing.

That exact-arc active pool is now working. With exact top-12 connector arcs, the
same `04174` pool solves to **OPTIMAL 453/472** and reconstructs to a concrete
route:

```text
valid=false stations=453/472 elapsed=22:16:30
```

The miss set is:

```text
108, 18, 193-195, 200-203, 206-210, 436-437, 446, 7, 8
```

This is not challenge-valid, but it is the first route-shaped optimization
artifact close to the world-record duration. The immediate priority is now to
recover those 19 stations under exact connector arcs: widen exact neighborhoods,
price targeted residual columns from the top-12 route's late endpoints, and
penalize uncovered stations in a way that preserves the `04174` final-tail
skeleton instead of replacing it.

Follow-up exact runs sharpened that priority. Exact top-24 connector
neighborhoods proved **OPTIMAL at the same 453/472 and 22:16:30**, so the miss
set is not caused by a too-narrow exact arc graph. Residual corridors through
the 19 misses were generated from the top-24 route endpoints; unforced, the
master selected none, and forced it dropped to **450/472**. Forcing the proxy
455 H06S/H11/A55/J31 chain under exact arcs dropped further to **442/472**.
Splitting the J31 final continuations into 2919 subcolumns and adding 86
internal detours, mostly for Morris Park and the 59 St pair, still returned the
same **453/472, 22:16:30** route and selected no useful new continuation
detours.

Exact connector pricing is now cached at
`reports/optimization_runs/exact_connector_cache.jsonl`. The cache is keyed by
source/target event plus run policy, run radius, pop cap, and max-cost cap. A
repeat top-12 exact run went from 6258 new Dijkstra proofs to **0 misses** while
reproducing the same 453/472 master, so wider exact neighborhoods are now
practical.

Using that cache, the `04174` pool was widened further:

| Exact neighborhood | Status | Replayed coverage | Elapsed | Takeaway |
|---:|---|---:|---:|---|
| top-36 | OPTIMAL | 453/472 | 22:16:30 | no improvement over top-24 |
| top-48 | OPTIMAL | 454/472 | 22:16:30 | recovers Norwood-205 St (`210`) |
| top-64 | OPTIMAL | 454/472 | 22:16:30 | no further improvement |

The best exact partial route is now **454/472 at 22:16:30**. The remaining miss
set is:

```text
108, 18, 193-195, 200-203, 206-209, 436-437, 446, 7, 8
```

The conclusion is now stronger: the exact 453/472 route is a stable basin, and
even wider exact connector neighborhoods only recover one additional station.
Small late-tail detours cannot recover the missing branches. The next
optimization pass should change an earlier macro phase or generate a different
late macro-order that covers Rockaway/Far Rockaway and Bronx branch terminals
before entering the J31-to-final skeleton. It should also try later start
buckets or earlier final continuations, because the current exact partial route
is already 140 seconds slower than the `22:14:10` record before adding the 19
missing stations.

A Rockaway-prioritized whole-tail beam tested that macro-order hypothesis from
the exact `F39S@77580` late source. It generated 764 final-reaching columns,
including rows that hit up to eight residual stations. But the exact top-48
master selected none of them and remained **454/472 at 22:16:30**. Forcing the
best high-preservation row, which covers `18`, `193-195`, `437`, `446`, `7`,
and `8`, collapsed the master to **427/472** and replayed at **428/472**. This
rules out monolithic `F39S -> final` Rockaway-first replacements in the current
pool: they cover missing stations by throwing away too much J/M/M-branch and
final-tail coverage.

The next structural model should decompose the late problem into compatible
stages instead of one whole-tail column:

```text
F39 / A-H source
  -> Rockaway or Far Rockaway repair
  -> J/M/M-branch preservation
  -> Bronx/Pelham/Dyre/final-tail preservation
  -> 227S final
```

That suggests a small late-tail branch-and-price/RCSP subproblem with handoff
states, not more single-column whole-tail beams.

This matches the standard decomposition literature for huge routing problems:
set-partitioning or set-covering masters are typically paired with shortest
path pricing subproblems with resource constraints, and large time-expanded
networks often need bounded dynamic programming, stabilization, and targeted
state-space restriction rather than monolithic MIP. Relevant references include
Irnich and Desaulniers on SPPRC pricing, Desrosiers/Soumis-style VRPTW column
generation, Engineer/Nemhauser/Savelsbergh on large resource-constrained
time-expanded pricing, and Pessoa/Sadykov/Uchoa/Vanderbeck on stabilization.

## Problem Framing

The Subway Challenge is closest to a time-dependent generalized graph tour:

* physical station labels must be covered at least once;
* paths through the timetable cover incidental intermediate stations;
* arc cost depends on departure time, not just endpoints;
* runs/transfers are legal mode changes with physical penalties;
* the full week graph is large enough that naive exact time expansion is too
  expensive for monolithic MIP.

Good formulations should therefore keep a sparse connection to the full graph:
solve compact, reconstruct concrete events, validate, then refine where the
compact model lied.

## Recommended Architecture

1. **Transition oracle**

   Cache earliest-arrival labels of the form:

   ```text
   start node/time bucket + target station/block
       -> end node, elapsed, path slice, covered stations, mode counts
   ```

   This is the primitive for time-bucketed matrices, pricing subproblems, and
   reconstruction. Cache keys should include run policy and graph hash.

2. **Time-bucketed macro routing**

   Build several station/block matrices by departure bucket rather than one
   static matrix. Feed those matrices to OR-Tools/cuOpt as fast order generators,
   realize each order on the timetable, then add new buckets around the realized
   departure times that cause the largest model error.

3. **Column generation**

   Treat schedule-realized route segments as columns. A column has:

   ```text
   start gateway, end gateway, departure bucket, arrival bucket,
   elapsed, covered station set, path slice
   ```

   The restricted master starts with columns extracted from existing validated
   routes plus generated terminal/branch segments. The pricing subproblem is a
   resource-constrained shortest path over the timetable with dual rewards for
   newly covered stations.

4. **Dynamic discretization**

   Start with coarse time buckets. Whenever a compact solution realizes poorly,
   add the actual departure/arrival times from the failed reconstruction and
   rerun. This is a better fit than building a full time-expanded MIP upfront.

5. **Async solver portfolio**

   Use NEOS for compact AMPL/LP/MIP experiments and lower-bound probes. Use
   cuOpt or OR-Tools for high-throughput matrix routing. Feed every promising
   macro order back into local LNS and the full validator.

## External Solver Roles

| Tool | Best use | Caution |
|---|---|---|
| OR-Tools Routing | fast local TSP/VRP/LNS macro-order generation | static matrices need timetable repair |
| cuOpt Routing | GPU/batch macro-order portfolio over many matrices | service solves matrix VRP/TSP, not the real GTFS graph |
| NEOS | async compact MIP/AMPL/LP experiments and bounds | submit only compact models; full graph is too large |
| HiGHS/SCIP/Gurobi/CPLEX | restricted master, set cover, compact arc models | install/configure as available |

## Immediate Implementation Backlog

1. Time-bucketed transition cache: start with coarse hourly buckets at terminals
   and decision nodes, then refine around incumbent event times. Exact-cover now
   has start/end bucket constraints; the next cache should generate columns
   conditioned on those bucket states, especially buckets that were infeasible
   or missed the record-compatible end range. Final-end-only and broad
   middle-to-record targeted passes did not restore feasibility. First
   cluster-corridor passes can repair A/H/J/M but lose too much final-tail
   coverage. v3 corridors nearly close that gap when forced, so prioritize a
   dedicated late-tail resource-constrained pricing subproblem or a
   compatibility-state master for multiple late corridors.
2. Column pricing: generate new exact timed fragments with dual rewards for
   stations and penalties for awkward endpoints/time buckets. Broad terminal and
   all-target pricing tied the incumbent; expensive-window pricing found
   below-incumbent proxy artifacts; LP target and beam pricing now generate
   negative-reduced-cost fragments but do not yet improve exact sequencing. The
   phase-variant pricer generates locally useful alternate departure phases, but
   exact-cover ignores them unless forced, and forced exact routes are worse.
   The next pricer should include endpoint/time-bucket penalties, exact connector
   cost to/from active gateway events, and multi-column replacement windows that
   move adjacent phases together, not isolated phase-shifted fragments.
3. Time-bucketed sequenced master: promote priced fragments into active exact
   pools, exact-price connector arcs up front, and derive dual-like penalties
   from the sequenced master. Replay cuts are still useful, but active exact
   pricing should be the default for serious candidate pools. The phase-window
   experiments suggest the next master should duplicate high-value gateway
   states by time bucket and connect only bucket-compatible fragments, instead
   of relying on the path-cover solver to discover phase coherence from a flat
   pool of event-fixed columns. `exact-cover --uncovered-penalty-s` is now
   available as a diagnostic for identifying the station clusters that make a
   bucketed master infeasible. `exact-cover --hint-solution` should be used for
   enlarged exact pools; the staged Rockaway split pool was large enough that
   CP-SAT found weaker 451/472 incumbents until seeded with the known 454/472
   route. A fixed-prefix/open-bridge neighborhood with the final
   `J31S@87630 -> 227S@101970` continuation held fixed proved OPTIMAL at
   454/472, so the next late-tail model must move the J31 handoff or generate
   multiple compatible final-tail continuations instead of only inserting
   Rockaway repair fragments before J31. Early all-residual suffix rows were
   also tested exactly: forcing `A47S@68220 -> 227S@96720` or
   `S01N@67680 -> 227S@96720` covers the current miss set locally but collapses
   the master/replay to roughly 394-395/472 because the suffix deletes too much
   unique mid-late coverage. The next useful neighborhood should therefore
   protect multiple corridors at once, or price replacement columns for the
   coverage lost before any earlier suffix handoff.
4. cuOpt/OR-Tools portfolio: solve many bucketed matrices, reconstruct, then
   launch bounded LNS repairs.
5. NEOS jobs: export compact sequenced masters for async MILP experiments and
   lower-bound probes.
6. Route-search neighborhoods: if continuing LNS, add richer destroy/repair
   moves over branch windows and start-phase blocks; the existing rotated-anchor
   start-grid and high-freeze incumbent sweeps repeatedly fall back to
   25-hour routes. Revisit station anchors alone are too coarse; preserve
   platform/time segment structure instead. `search.py` now supports
   `--run-mode all`, but all-run LNS from the forced phase-corridor route still
   degraded to 25:51:30, so more run options do not fix the anchor abstraction.

## Sources

* NEOS XML-RPC supports job submission, polling, nonblocking final-result fetch,
  and solver templates through `submitJob`, `getJobStatus`,
  `getFinalResultsNonBlocking`, and `getSolverTemplate`.
* NVIDIA cuOpt exposes routing APIs around distance/time matrices, fleets, tasks,
  time windows, and solver time limits; it is appropriate as a routing engine for
  compact matrix instances.
* Dynamic discretization discovery is designed for time-dependent TSP variants
  where full time-expanded models become prohibitively large.
* Column generation and branch-and-price are standard approaches for huge VRP
  formulations, typically with set-covering masters and resource-constrained
  shortest-path pricing.
