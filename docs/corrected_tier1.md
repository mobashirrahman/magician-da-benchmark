# Corrected Tier 1 protocol: tier1-v3

Protocol revision: 2026-10-06. This revision precedes interpretation of any full
benchmark result. It is committed and tagged `preregistration-tier1-v3` before the
full run starts; the bundle records that revision in `revision.json`.

Validation completed: **86 regression tests passed**, and all **264 R method runs
across 12 datasets succeeded**. Reporting and conditional analysis also completed.
Null cases contain no positive design labels; all validation recommendations are
marked `smoke_only`. Evidence is retained in
`results_paper_tier1_v3_validation/provenance/validation.json`.

## Amendment, 2026-10-07: radEmu excluded from Tier 1

Tagged `preregistration-tier1-v3.1`. The first tier1-v3 launch was stopped after
six hours, with 43 of 4,060 datasets evaluated. radEmu took a median of 56 minutes
per dataset (maximum 2.4 hours) against under 4 minutes for the other 18 methods
together: about 3,800 core-hours, more than three weeks on the one shared 16-core
machine available. Tier 1 therefore runs **18 method configurations**, and radEmu
is reported as "excluded from Tier 1: compute cost", not as a failure. It remains
active in the registry. The decision used run times and completion statuses
only; no score from the stopped launch was inspected. The input bundle
(`cache/paper_tier1_v3`, generated at `preregistration-tier1-v3`) is unchanged;
the partial results were deleted and the workflow restarted from the bundle.
The validation counts below include radEmu.

## Why tier1-v2 was withdrawn

The `tier1-v2` run was stopped at about 1 % and deleted. Compact records are in
`archive/tier1_v2_20261006`; the original run's are in `archive/tier1_original_20261005`.

- Truth labelled every non-zero shift as differential. All 300 features were
  positive on the CLR lane in every spiked case, and on the relative-abundance and
  read-fraction lanes under one-directional change. Sixteen of 19 methods could not
  make a false discovery there, and recall rewarded calling everything.
- The G2 donor table was a Dirichlet CAMISIM design around one mean profile with
  30 implanted effects, not measured profiles.
- Simulated counts had no zeros, so filters and the prevalence-dependent zero cell
  did nothing.
- Null cases were regenerated for every effect cell although null data do not
  depend on effect factors: 7,000 of 9,800 datasets carried no signal.

## Corrections in tier1-v3

- **Design truth is primary.** Discoveries are scored against the implanted set,
  which is the same for every method (`fdr`, `recall`, `source_recall`, average
  precision, eligibility and ranking). Closure shifts of unselected features are
  recorded (`is_relative_shifted`, `true_log2fc`) but are not design truth.
- **Endpoint lanes are secondary.** Each method is also scored on the lane of its
  declared endpoint (`endpoint_fdr`, `endpoint_recall`). `endpoint_fdr` is
  unavailable when that scale has no null feature, instead of being zero.
  `spillover_share` is the share of a method's discoveries that are real shifts on
  its own scale but lie outside the design truth; it maps compositional bias.
- **G2 resamples measured profiles.** Donors are species profiles of 999 healthy
  adults from one study and body site (curatedMetagenomicData, AsnicarF_2021,
  MetaPhlAn 3; `tools/fetch_real_profiles.py`, pinned checksums, provenance in
  `resources/real_profiles/provenance.json`). Each case draws one distinct donor
  per sample, keeps measured zeros, and uses the most prevalent species. Synthetic
  feature expansion is disabled; the table has 636 species.
- **G1 and G3 are calibrated on the same table.** G1 is lognormal with the donor
  table's spread of mean abundances and its median within-feature log SD; its
  zeros come from sampling only. G3 is a hurdle model with each donor feature's
  prevalence and nonzero lognormal moments and prescribed block correlations. It
  is **synthetic**, not a fitted SparseDOSSA model.
- **Sparsity is realistic.** At 300 features about 69 % of G2 and G3 counts are
  zero and about a third of features fail the 10 % prevalence filter; G1 has
  4–10 % zeros. `diagnostics/realism.json` compares every null condition with
  donor profiles counted at the same library sizes.
- **Null cases are generated once per null condition** (generator, samples,
  features, depth, zero structure) and shared by every spiked cell with that
  condition. Failed null runs count against each spiked cell that depends on them.
- **Null criterion.** Eligibility needs a raw null p-value rejection rate of at
  most 0.07. A rate below 0.03 is valid and flagged `null_conservative`; the
  earlier 0.03–0.07 band rejected conservative methods.
- **Core effect size is log2FC 2.** A pilot of 30 datasets with 11 fast methods
  (not retained, not part of the run) showed that at log2FC 1 and 20 samples per
  group almost no method made a true discovery under donor-level variation. The
  pilot was used for this choice only.
- **edgeR** uses the classic `estimateDisp` + quasi-likelihood fit
  (`legacy = TRUE`). The edgeR 4.4.0 default fit aborts R on sparse counts when
  dispersions are supplied.
- DA jobs reserve 1 GB each for scheduling, so 16 run at once; measured use is
  below 1 GB. The earlier 2 GB reservation left half the cores idle.

Unchanged from tier1-v2: shared null/spiked random streams, balanced implants that
conserve the changed subset's mass, equal feature lengths, exact multinomial
library totals, exact binomial null intervals, paired bootstrap within conditions,
and the exclusion of batch confounding.

## Scope and evidence

The full run has **4,060 datasets**:

| Part | Cells | Replicates | Datasets |
|---|---|---|---|
| Core null conditions: 3 generators × 2 depths | 6 | 200 | 1,200 |
| Core spiked: × 2 fractions (10 %, 20 %) × 2 directions (balanced, all-up) | 24 | 100 | 2,400 |
| Survey null conditions not in the core | 8 | 20 | 160 |
| Survey spiked, one factor at a time around the G1 base cell | 15 | 20 | 300 |

The base cell is G1, 20 samples per group, 300 features, 2M read pairs, 10 %
changed, log2FC 2, balanced. The survey cannot identify factor interactions.
Its cells have 20 spiked replicates and stay exploratory; survey cells that share
a core null condition use that condition's 200 null runs.

Conditional eligibility: upper FDR confidence bound ≤0.10 against design truth,
null rejection rate ≤0.07, failure rate <5 %, positive recall, and no significant
unevaluable features. Recommendations require at least 100 null and 50 spiked
runs per condition.

Limits that remain. G2 inherits MetaPhlAn's detection floor, so its zeros do not
change with depth. Scoring one-directional cells against design truth counts real
closure shifts as false; `endpoint_fdr` and `spillover_share` show how much of a
method's error that is. Variance shares are descriptive and depend on factor
order; no mixed-effects model is fitted. Contrasts against the observed best
method are exploratory. The consensus rule is recorded but not evaluated. Tier 2
and real-data nulls are required for MAG and real-world reliability claims.

## Execution and retained evidence

```bash
python tools/fetch_real_profiles.py        # only to rebuild the committed donor table
python tools/paper_tier1.py --smoke --null-replicates 1 --spiked-replicates 1 --cores 16
python tools/paper_tier1.py --full --cores 16
```

Destinations are versioned and must be empty. The full command generates
`cache/paper_tier1_v3`, runs the workflow into `results_paper_tier1_v3`, then
writes the conditional analysis to `results_paper_tier1_v3/analysis`. Its log is
`cache/paper_tier1_v3.log` (workflow restarts: `cache/paper_tier1_v3.workflow.log`).
To resume an interrupted workflow:

```bash
python run_magician.py --configfile config/paper_tier1_v3.yaml --cores 16 --scheduler greedy
python tools/paper_analysis.py --tier1 results_paper_tier1_v3/benchmark --out results_paper_tier1_v3/analysis
```

Each bundle retains case-level generator records, full condition metadata, the
donor table's checksum, source hashes, a source-code snapshot and the git revision.
