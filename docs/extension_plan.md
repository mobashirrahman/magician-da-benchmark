# Extending the MAG differential abundance benchmark

Research checked on **2026-10-02**. This is an implementation and experiment plan;
the new adapters below are not enabled or validated yet. The current seven methods
remain the baseline. The larger artificial run uses [pilot20.yaml](../config/pilot20.yaml)
and completed all 195 jobs; see the
[pilot record](validation.md#twenty-genome-artificial-pilot).
User-selected scope: run that artificial pilot (completed 2026-10-02, see
[validation record](validation.md#twenty-genome-artificial-pilot)) and plan the real-genome study.

## Recommended additions

Start with MaAsLin 3, ALDEx3 and a scale-aware ALDEx2 variant. Add ADAPT and radEmu
next, alongside LinDA and ZicoSeq as established baselines. The priority reflects
scientific coverage, available software and integration cost; it is not a claim
that a particular method will win on MAGs.

| Priority | Method and evidence | Proposed input/endpoint | Main integration issue |
|---|---|---|---|
| 1 | **MaAsLin 3**, Nature Methods, January 2026; released package 1.4.0 | Length-corrected relative abundance; abundance first, prevalence separately | Select abundance rows explicitly; declare median correction and correction family |
| 1 | **ALDEx3**, CRAN 1.3.1, August 2026; related comparative paper September 2026 | Integer counts; CLR/scale-model contrasts | Streaming Monte Carlo, explicit scale assumptions, compiled C++17 dependencies |
| 1 | **ALDEx2 with scale uncertainty**, Genome Biology, May 2025 | Integer counts; sensitivity to normalization assumptions | Treat gamma settings as variants of one family; validate installed API |
| 2 | **ADAPT**, Bioinformatics 2024; release page reports 1.4.0 | Integer counts; adaptive-reference abundance contrast | Tobit censoring, reference selection, internal sample/feature filtering |
| 2 | **radEmu**, Biometrika 2026; CRAN 2.3.2.0, April 2026 | Counts or supported abundance measurements; contrast relative to typical taxon change | Explicit identification constraint and robust score tests; joint-model runtime |
| 2 | **LinDA**, Genome Biology 2022 | Counts or supported compositions; bias-corrected CLR contrast | Minority-change assumption and zero handling; preserve reference interpretation |
| 2 | **ZicoSeq**, Microbiome 2022; GUniFrac manual reports 1.8 | Counts initially; reference-normalized permutation tests | Permutation resolution, internal filtering and transformed effect units |
| 3 | **LOCOM**, PNAS 2022 | Counts; compositional logistic contrast | Permutation cost; retain taxon tests separately from community-wide test |
| 3 | **metaGEENOME**, BMC Bioinformatics, July 2025 | Counts; GEE with CTF normalization/CLR | Extract DA component from a broad reporting workflow; verify small-sample inference |

Release numbers are the versions observed in the linked primary distribution
pages, not promises of successful installation. Freeze sources and dependencies
only after the adapters pass their compatibility tests.

**MaAsLin 3** is the most direct extension of the existing MaAsLin2 comparison.
It models nonzero abundance and presence separately and supports more complex
designs. Its composition correction changes the hypothesis, so the relative
abundance lane will use `median_comparison_abundance=FALSE`. A separate corrected
lane will test change relative to the typical feature, with matching truth.
Evidence: [2026 paper](https://www.nature.com/articles/s41592-025-02923-9),
[release page](https://bioconductor.org/packages/release/bioc/html/maaslin3.html).

**ALDEx3** adds scale-model regression with linear and mixed-effects engines.
The September comparative paper supports investigating its computational benefits,
but used an earlier research version and does not establish a MAG-specific winner.
Keep ALDEx2 as a baseline when adding ALDEx3. Evidence:
[CRAN release](https://cran.r-project.org/web/packages/ALDEx3/index.html),
[September 2026 comparison](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014728).

**Scale-aware ALDEx2** offers a relatively small integration change using the
existing package family. Scale models express uncertainty in normalization
assumptions; gamma must be specified before inspecting benchmark discoveries.
Evidence: [2025 scale-model paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12100815/).

**ADAPT** treats zeros as censored observations and learns a reference set.
Test it especially under sampling dropout and sparse counts, without assuming
that its reference set is correct in every simulation. Evidence:
[original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11959182/),
[release package](https://bioconductor.org/packages/release/bioc/html/ADAPT.html).

**radEmu** estimates contrasts relative to typical changes and handles zeros
without a pseudocount. It does not identify absolute abundance in individual
samples from sequencing alone. Pin the CRAN release rather than the development
site's newer version, and record its constraint explicitly. Evidence:
[authors' documentation and 2026 citation](https://statdivlab.github.io/radEmu/),
[CRAN distribution](https://www.stat.ethz.ch/CRAN/web/packages/radEmu/index.html).

**LinDA, ZicoSeq and LOCOM** are useful established comparisons rather than new
2026 discoveries. They expand the model families already represented. Evidence:
[LinDA](https://link.springer.com/article/10.1186/s13059-022-02655-5),
[ZicoSeq](https://link.springer.com/article/10.1186/s40168-022-01320-0),
[LOCOM](https://pmc.ncbi.nlm.nih.gov/articles/PMC9335309/).

**metaGEENOME** is a promising optional GEE comparison. Its published package
combines DA with extensive downstream reporting. Its README also has differing
repository names in its installation example; resolve and pin the paper-linked
source. Disable unnecessary plotting/diversity work and verify inference with
small numbers of independent subjects before ranking it. Evidence:
[2025 paper](https://link.springer.com/article/10.1186/s12859-025-06217-x),
[paper-linked source](https://github.com/M-Mysara/metaGEENOME).

Do not add classifiers, community-wide association tests, consensus wrappers or
prevalence-only preprints as if they were interchangeable taxon abundance tests.
An ensemble would need its own predefined rule and independent evaluation seeds.
Longitudinal models require a genuine repeated-subject generator, not duplicated
independent samples. Optional fastANCOM can be considered later as a speed baseline.

## Adapter and workflow changes

1. Add a method registry under `config/` with accepted input types, hypothesis,
   effect units, package/environment, parameters, seed, thread cap and citations.
   Use separate `method_family` and `variant` fields so multiple gamma choices do
   not masquerade as independent methods. Extend the JSON schema and validate
   compatibility before constructing jobs.
2. Split `workflow/scripts/run_da.R` into a shared input/output contract and small
   adapters under `workflow/scripts/da/`. Track the selected adapter and registry
   settings as Snakemake dependencies. Each job continues to have its own log,
   benchmark, environment, threads, memory, scratch and explicit status.
3. Separate durable raw quantification from filtering/normalization. Currently
   these share `build_matrices`; changing filtering can require deleted upstream
   counts. Retain raw counts, lengths and library totals, then derive filtered
   inputs in a cheap rule. Add a matrix-input execution mode so new methods and
   matrix simulations never require FASTQ, assembly or binning rules.
4. Extend result rows with `endpoint`, `contrast`, `effect_scale`, `reference`,
   native p/q-values, `adjustment_family`, per-feature fit status and exclusion
   reason. Preserve every original feature ID and every untested feature.
   Join all input sample IDs strictly before calling packages; reject silent
   intersection, dropped groups, reversed contrasts or altered names.
5. Keep native adjusted values; never apply BH again to package q-values. Raw
   p-values remain the common average-precision score when available. If a future
   adapter exposes only q-values, label a separate q-ranked metric or report AP
   unavailable; do not reverse-engineer p-values.
6. Add method settings and their measured resources to the report, stratified
   by endpoint and input. Empty catalogues, convergence failures, timeouts and
   unavailable packages stay visible and make the affected method ineligible.

| Adapter | Implementation contract to validate |
|---|---|
| MaAsLin 3 | `maaslin3(..., formula='~ group', normalization='NONE', transform='LOG')` on genomic relative abundance; Control baseline; extract group abundance rows, `pval_individual`/`qval_individual` and declared correction family. Use native joint/prevalence output only in the matching lane. Disable summary/association plots, model and plot RDS storage; no implicit mirai workers. |
| ALDEx3 | `aldex(Y, ~group, metadata, method='lm', scale=clr.sm, gamma=..., n.cores=1, test='t.HC3')`; request `p.val` and `p.val.adj` from the fitted object, which are documented even though the printed summary emphasizes q-values. Fix seed and Monte Carlo count; set `streamsize` to a measured chunk limit, initially 256 MB, and omit persisted draws. |
| ALDEx2 scale | Validate `aldex.clr` scale/gamma support in the installed release; compare predefined gamma 0 and 0.5, with optional 0.1 sensitivity. Use 128 draws for screening and 512 or more for confirmation after Monte Carlo stability checks. |
| ADAPT | Build a minimal phyloseq object; `adapt(..., cond.var='group', base.cond='Control', prev.filter=0, depth.filter=0, alpha=0.05)` after common filtering. Preserve zeros, censor setting and learned reference IDs. Extract all results, not just discoveries. Native `log10foldchange` converts to log2 by multiplying by `log2(10)`, while retaining the reference endpoint. |
| radEmu | Samples by features; `emuFit(formula=~group, data=metadata, Y=..., test_kj=...)` tests only the treatment coefficient for retained features. Use robust score inference, one core initially, and record convergence and reference/typical-change constraint. Confirm log base from the pinned manual before unit conversion. |
| LinDA | `MicrobiomeStat::linda` with explicit input type/formula; disable additional prevalence/library filtering after common filtering where supported. Record adaptive zero handling and bias correction. Verify effect coordinates and sign on a deterministic fixture. |
| ZicoSeq | `GUniFrac::ZicoSeq` count input; explicit square-root link, posterior sampling and reference settings; `perm.no>=999`, `return.feature.dat=FALSE`. Retain `p.raw` and native `p.adj.fdr`. Its regression coefficient on transformed ratios is not a genomic log2 fold change. |
| LOCOM | `locom` with samples by features, binary group, explicit seed, one core and permutation budget. Retain native taxon p/q-values and reference-scaled effect; report `p.global` separately. |
| metaGEENOME | Call only the GEE DA component if public API permits; otherwise isolate the official wrapper in temporary scratch with optional analyses disabled. Pin a commit and test returned tables; defer if the broad wrapper cannot run without taxonomy/irrelevant outputs. |

API details above are grounded in the official
[MaAsLin 3 manual](https://bioconductor.org/packages/release/bioc/vignettes/maaslin3/inst/doc/maaslin3_manual.html),
[ALDEx3 manual](https://cran.r-project.org/web/packages/ALDEx3/ALDEx3.pdf),
[ADAPT tutorial](https://bioconductor.org/packages/release/bioc/vignettes/ADAPT/inst/doc/ADAPT-manual.html),
[radEmu examples](https://statdivlab.github.io/radEmu/),
[ZicoSeq manual](https://search.r-project.org/CRAN/refmans/GUniFrac/html/ZicoSeq.html),
and [LOCOM source](https://github.com/yijuanhu/LOCOM). Validate all calls against
the actual frozen package, rather than treating documentation as an executed test.

## Truth and fair comparison

The current generator defines **relative genomic abundance**, allocates reads in
proportion to abundance times genome length, and conserves the changed subset's
relative mass. The current fixture has equal genome lengths. Real genomes have
different lengths: count proportions are DNA/read fractions, whereas the retained
relative matrix is length corrected. MAG length, incomplete recovery and changing
catalogue composition add further differences. Define truth on the measurement
scale being tested instead of labelling a correct read-fraction association false.

Add four explicit truth lanes: relative genomic abundance; expected read fraction;
change relative to a declared reference/typical feature; and absolute abundance
when simulated or measured load information exists. Compute each truth from the
latent design, including genome lengths and the specified denominator. Separate
presence probability and abundance conditional on presence for hurdle models.
An unrecovered MAG is a recovery failure, not evidence of biological absence.

Preprocessing also needs a declared unit-aware policy. The current 0.5 pseudocount
changes CLR results differently for counts, TPM and proportions, as the pilot
shows. Test a count-based Dirichlet/multiplicative zero approach and predefined
measurement-scale pseudocount sensitivity, rather than treating a rescaled matrix
as an independent biological dataset. Compute expected truth using the latent
generator, accounting for sample-level closure/noise and finite-sample variation.

Apply common truth-blind filtering, retain raw totals, and log any unavoidable
method-specific exclusions. For Dirichlet methods, evaluate a reserved `other`
category containing filtered counts to preserve totals; it is not a tested
biological discovery. Adding that category changes the reference geometry, so
make this a declared variant and compute truth in its matching coordinates.
Unmapped reads and unbinned contigs require their own recorded library accounting;
do not infer them from source truth inside a method.

Report FDR, precision, recall, raw-p AP, direction accuracy, fitting success,
runtime and RSS. Compute fold-change error only for matching effect units and
truth coordinates. Distinguish recall among recovered MAGs from source-level
recall, where unrecovered true changes count as misses; deduplicate multiple bins
from the same source. Report unassigned significant bins separately.

For repeated experiments, report mean spiked FDR and its uncertainty, and null
probability of any false discovery (FWER under the global null). Bootstrap whole
independent seeds, keeping method comparisons paired. Add exact/binomial intervals
for the null rate. Five seeds satisfy the current minimum software label but
cannot substantiate 5% error control. Use at least 100 inexpensive null matrix
replicates per primary design for calibration, and show their intervals. Select
settings on development seeds and compare frozen methods on held-out seeds.
Allow ties and no eligible method; publish a tradeoff frontier rather than force
one universal winner.

## Larger experiments

| Stage | Design | Purpose and retained evidence |
|---|---|---|
| Artificial pilot (completed) | 20 genomes of 350 kb, 50 contigs each; 10 samples/group; null + spiked; seed 1001; 200,000 nominal pairs/sample; biological CV 0.30; four changed sources | Exercise all current seven methods and 52 compatible fits, MAG reconstruction and cleanup. One artificial seed remains smoke evidence. |
| Matrix screening after adapters | 50, 200 and 500 features; 10/20/40 samples/group; effects 2x/4x; changed fraction 0.1/0.2/0.4; low-count/dropout and load-shift designs | Predefine a small balanced subset of combinations, not the full product. Use 100 null and 30 spiked independent replicates per selected design; retain matrices, truth, scores and resources. |
| Real-genome feasibility | 20 diverse real genomes, then 40–50; 10 samples/group; one paired null/spiked seed per pilot | Verify depth, mapping ambiguity, usable MAG recovery and measured disk/RAM before repeated execution. |
| Real-genome confirmation | Two primary conditions chosen before held-out evaluation; 10 independent paired null/spiked seeds each, increased to 20 if intervals remain broad | All validated methods share each catalogue/matrix. Each pair uses shared baseline/library noise; different seeds are independent. Retain per-seed recovery and paired performance differences. |
| Later correlated designs | Repeated subjects with explicit subject effects and balanced batches | Compare mixed/GEE/permutation-strata methods only after corresponding truth and metadata exist. |

Prepare a versioned real-genome manifest with accession/version, local FASTA path,
checksum, length and provenance. Choose public benign reference communities with
different genome lengths/GC and include a predefined related-genome subset to
test ambiguous mapping. Reuse suitable archived local references by path after
checking provenance and completeness; do not copy whole historical collections
or download a large database. Preserve an immutable manifest and selection rule.
Real-genome downloads and execution are future steps; this iteration runs the
artificial pilot and records the real-genome design.

For 40 genomes averaging 3 Mb and 20 samples, 500,000 nominal pairs/sample gives
roughly 25x mean aggregate sequence coverage before library variation, assuming
uniform abundance. Rare genomes receive less coverage. Start lower if the storage
estimate rejects the design; increase depth only after measuring recovery.
Doubling sample count at the same depth doubles reads and can exceed 100 GB.
Record failures to recover low-abundance sources as part of the benchmark rather
than silently restricting truth to successful bins.

## Resource and reproducibility plan

Keep the existing heavy-case sequencing and `temp()` cleanup. Share one assembly,
one binning catalogue and the resulting matrices among all methods for each case.
Use one core per DA job initially, including BLAS/OpenMP/TBB and package workers.
Benchmark new jobs at 2 GB declared RAM, increasing only when measured use requires
it. Set explicit timeouts and record timeout status; never silently retry with a
different statistical method.

Reuse the existing Python and biology environments. Resolve a shared modern R
environment for compatible additions; use a second environment only for actual
incompatibilities. Bioconductor releases and their R versions must match: the
MaAsLin 3 1.4.0 release page refers to Bioconductor 3.23, so do not assume it drops
into the existing R 4.5 environment. Keep frozen baseline environments for comparison.
Freeze Conda builds/checksums and CRAN tarballs or Git commit SHAs; installation
happens before benchmark execution, with no network work inside DA rules.
Check the environment reserve against measured installed sizes before full runs.

ALDEx3's documented streaming threshold defaults to about 8,000 MB, larger than
the current DA reservation. Set it explicitly and retain aggregate coefficients
and p/q-values, not Monte Carlo arrays. Suppress per-feature plots and fitted-model
objects for every package; normalized input copies and temporary package reports
belong in scratch. Keep package messages, version/session records and one unified
report for audit.

The existing conservative preflight estimates roughly `8 * raw FASTQ bytes +
20 * reference bases + retained tables + 17 GB` for environment reserve/headroom.
At 150 bp reads, raw paired FASTQ is estimated at 900 bytes/pair. Twenty samples
with 500,000 pairs each therefore need about 72 GB work allowance plus references,
retained outputs and the reserve; actual lognormal allocations can push the total
near 100 GB. Preflight must evaluate those actual allocations. Use a 200 GB budget
for deeper selected cases, or reduce depth/sample count. More methods/seeds reuse
matrices and serial scratch; they primarily add compact tables and CPU time.
The artificial pilot's preflight estimate is about **45.7 GB**, within 100 GB.
The completed pilot retained about 20 MB and its largest job used 982 MB RSS; see
[validation record](validation.md#twenty-genome-artificial-pilot).

The pilot also exposes a scaling cost in exact simulation: a separate wgsim call
per source contig can require 1,000 subprocesses/sample here. Measured: 40 simulation
jobs took 12,377 of 16,082 job seconds (77%), about five minutes per sample. A future optimization
should batch simulation while preserving exact allocations and deterministic
record counts, and compare its output distribution against the current reference
implementation. Do not sacrifice allocation correctness for speed.

Archive moves preserve historical data without duplication; they do not free
space. Historical archives are outside the active run's budget, but filesystem
free-space checks still include their footprint. See [storage policy](storage.md)
and [archive manifest](../archive/2026-10-02/manifest.json).

## Implementation sequence and acceptance gates

1. **Registry, result schema, matrix-only entry point and truth lanes.** Demonstrate
   that changing method/filter parameters schedules analysis only, with zero
   simulation/assembly jobs after heavy temporary files are removed.
2. **MaAsLin 3, ALDEx3 and ALDEx2-scale.** Pin environments and test actual adapters
   on common fixtures: positive/negative effects, zeros, zero libraries, one/empty
   MAG, unsafe feature names, permuted sample order and missing metadata. Validate
   each hypothesis against its own truth and explicitly test adjustment families.
3. **ADAPT, radEmu, LinDA and ZicoSeq.** Add reference/zero-heavy fixtures and
   performance tests on 500 features/80 samples. Prespecify permutation counts:
   at least 999 for screening, then enough for the smallest relevant multiple-test
   threshold (often 9,999 or more), checking resolution rather than choosing only
   on speed. Validate failure and timeout reporting.
4. **Matrix calibration, then real-genome feasibility.** Freeze primary settings
   and held-out seeds. Permit full MAG confirmation only after adequate recovery
   and a storage estimate below budget with headroom.
5. **Repeated real-genome confirmation and report.** Require all job logs/resources,
   source/MAG truth accounting, uncertainty intervals, explicit exclusions and
   cleanup evidence. Archive obsolete runs using manifests after completion.
   Add LOCOM/metaGEENOME only if the earlier gates establish useful coverage
   within compute limits.

Deliverables: validated adapters and locks; versioned method registry and genome
manifest; scenario/truth tables; compact matrices; complete per-method predictions
and statuses; FDR/recall/resource comparisons with uncertainty; a reproducible
report stating when evidence supports a candidate, a tie or no recommendation.

## Implementation status

Research checked on **2026-10-04** for this update. Stages 1 to 4 of the sequence
below are implemented; stage 5 is designed and its tooling is in place, but the
repeated real-genome execution has not been run. See
[validation record](validation.md#extension-plan-implementation) for executed evidence.

| Item | State |
|---|---|
| Method registry with `method_family`/`variant`, declared inputs, units, references, adjustment families, environments, seeds, thread caps and citations | implemented in [config/methods.yaml](../config/methods.yaml) and `src/magician/registry.py` |
| Shared R contract and one small adapter per method | implemented in `workflow/scripts/da/` |
| Result rows with endpoint, contrast, effect scale, reference, native p/q, adjustment family, per-feature fit status and exclusion reason | implemented |
| Durable raw quantification separated from filtering and normalization | implemented (`quantify_features`, `derive_matrices`) |
| Matrix-only execution mode | implemented (`input_mode: matrix`, `src/magician/matrix.py`) |
| Four truth lanes plus the CLR geometry per declared zero policy | implemented in `src/magician/truth.py` |
| Unit-aware zero policies and the reserved `other` category | implemented in [config/zero_handling.yaml](../config/zero_handling.yaml) |
| Method settings and measured resources stratified by endpoint and input | implemented (`benchmark/method_settings.tsv`, `benchmark/method_resources.tsv`) |
| Original seven methods validated against real packages | done, 0 failures |
| MaAsLin 3, ALDEx3, ALDEx2-scale, ADAPT, radEmu, LinDA, ZicoSeq, LOCOM adapters | written; environments and locks being frozen, per-method validation pending |
| Matrix screening designs and replicate generator | implemented (`src/magician/screening.py`, `tools/matrix_screening.py`); 12 predeclared designs |
| Calibration replicates | 40 replicates (20 null, 20 spiked) across 4 designs ran all 16 active methods end to end (964/964 jobs); the 100-null-per-design scale uses the same command with larger replicate counts |
| Twenty-genome artificial pilot under the new schema | re-run end to end (205/205 jobs); endpoint-stratified candidates recorded |
| Versioned real-genome manifest, selection rule and feasibility record | implemented (`src/magician/genomes.py`, `tools/real_genomes.py`) |
| Repeated real-genome confirmation | not run; requires a frozen manifest and the stage 4 storage check |

One scientific decision changed during implementation and is recorded here because it
affects interpretation. The count-based negative-binomial models (DESeq2, edgeR,
limma-voom) fit normalized counts, so their estimand is the read or DNA fraction, not
the length-corrected relative abundance. Their registry endpoint is therefore
`read_fraction` and they are scored against the expected read fraction lane. With
equal genome lengths the two lanes coincide, so the artificial pilot numbers are
unchanged; with real genomes of different lengths they differ, and the benchmark now
reports which lane each method was tested on instead of labelling a correct
read-fraction association a false positive.
