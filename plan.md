# Experiment plan: a publishable benchmark of differential abundance on MAGs

Protocol correction, 2026-10-06: Tier 1 runs as `tier1-v3`. The original run and
`tier1-v2` were stopped and deleted after audits of the implants, truth labels,
the G2 donor table and sparsity. See [docs/corrected_tier1.md](docs/corrected_tier1.md)
for the operative design: 24 spiked core cells on 6 shared null conditions plus a
one-factor-at-a-time survey; batch confounding excluded. That document overrides
Sections 2 and 4 below where they differ: discoveries are scored against the
implanted design truth, G2 resamples measured profiles from one healthy cohort,
and the null criterion is a rejection rate of at most 0.07. The original sections
are retained as the broader research plan, not a record that all analyses are
implemented.

Status: tooling implemented 2026-10-04 (reduced scope Tier 1+2+3a first);
full compute not yet run. Resolved: reduced scope first, Phase 0 auto-picks
IBD/CRC cohorts, keep strain axis + second binner.
Compute cap: **60 cores** (the server has 128; leave the rest for other users).

## 1. Why this is publishable

Existing DA benchmarks (e.g. Calgaro 2020, Nearing 2022, Yang & Chen 2022, Hawinkel
2019) test taxon tables from 16S or marker-based profilers. None of them carries
errors from *genome reconstruction* into the DA step. MAG-based studies are now
routine, and their feature tables differ in ways nobody has quantified:

- bins are fragmentary, chimeric or merged across strains;
- low-abundance genomes are missed in a way that depends on coverage, so the zeros
  are not random;
- counts must be length-corrected, and the catalogue itself changes with the data.

Our repo already has the pieces most benchmarks lack: end-to-end simulation to
MAGs, a source-genome oracle, endpoint-matched truth lanes, a method registry with
families versus variants, exact null intervals, and refusal to let empty-handed
methods win.

**Claims we want to be able to make** (each is a hypothesis with a pre-declared test):

1. **MAG cost.** MAG-based DA loses measurable power and gains false discoveries
   relative to the same data analysed on source genomes, and the loss grows
   as coverage falls and community complexity rises.
2. **Truth definition matters.** Several apparent "false positives" (the current
   DESeq2 spiked FDR of 0.29, edgeR's null rate) may be mismatches between the
   hypothesis tested and the truth scored. Scoring each method on its own endpoint
   lane changes the ranking. We must verify this before claiming it; the
   calibration run could equally show real anti-conservatism.
3. **Calibration is not power.** Methods that control null error on parametric
   simulation fail on real-data nulls, and the failures concentrate in sparse,
   low-coverage features.
4. **Compositional bias is design-dependent.** Under one-directional change or
   total-load change, CLR/reference-based methods and count-based methods
   diverge. We map where, rather than declaring a winner.
5. **A decision guide.** The output is a pre-registered, conditional
   recommendation (by sample size, depth, sparsity, expected effect structure),
   not a single winner.

Target venues: *Genome Biology*, *Microbiome*, *Bioinformatics* (benchmark
article). Post a preprint on bioRxiv when Tier 1 and Tier 2 finish.

## 2. Design: four tiers, cheapest first

| Tier | What | Truth | Scale | Answers |
|---|---|---|---|---|
| 1 | Matrix-only simulation | Known by construction | ~6,000 cases | Calibration, power and robustness over a wide design grid |
| 2 | Full read-level pipeline to MAGs | Latent design + source genomes | ~60 experiments | MAG cost (claim 1), catalogue artefacts |
| 3 | Real-metagenome implants and real nulls | Implanted reads / label permutation | ~6 cohorts-sets | Realism of claims 2-4 |
| 4 | Real disease cohorts, replication | None; cross-cohort reproducibility | 2 diseases x 2-3 cohorts | External validity, not accuracy |

Rule: nothing in Tier 3-4 is used to tune methods or thresholds. Thresholds and
analysis code are frozen and committed (Section 8) before Tier 2 starts.

### Tier 1: matrix-only grid

Replace the current Dirichlet-multinomial-only generator with **three generators**
so conclusions do not depend on one simulator:

- **G1 parametric**: the existing Dirichlet-multinomial screening designs
  (baseline, overdispersed, low_count, dropout).
- **G2 empirical-resampled**: per-feature abundance and prevalence distributions
  taken from real profiles (the `data/from_real_v1..v3` tables are recoverable
  from git history, commit 6a1328e onwards, since they are deleted in the working
  tree), with sample-wise correlation structure kept by resampling whole real
  samples, then implanting effects.
- **G3 copula/SparseDOSSA-style**: fitted marginals plus a correlation structure,
  as an independent third generator (use the published package if it installs
  cleanly in a frozen env; otherwise a documented Gaussian copula).

Factorial axes (levels are proposals; trim by a fractional design if the full
cross exceeds budget):

| Axis | Levels |
|---|---|
| Samples per group | 5, 10, 20, 40 |
| Features | 100, 300, 1000 |
| Depth (read pairs/sample, scales counts) | 0.5M, 2M, 10M |
| Fraction DA | 0, 5, 10, 20, 40 % |
| Effect size (log2FC) | 0.5, 1, 2, 4 |
| Direction | balanced, all-up, all-down |
| Total-load change | none, group B load x2 (absolute-abundance lane) |
| Zero structure | none, sampling, dropout (structural), prevalence-dependent |
| Confounder | none, batch correlated with group (r = 0.3, 0.6) with adjusted model |

Replicates: **200 independent nulls and 100 spiked per cell** on a core subset of
~30 cells, plus 20 per cell on the remaining fractional grid. Null calibration is
judged at the **p-value level** (QQ plots, false-positive rate at alpha 0.05, KS
on uniformity), which yields thousands of null tests per cell, and by the
any-false-discovery rate only on the core cells. Use cluster (per-experiment)
bootstrap because features within an experiment are not independent.

### Tier 2: full pipeline to MAGs

Genomes come from a frozen manifest of real isolate/RefSeq genomes (`tools/real_genomes.py
freeze`), not synthetic ones. Abundances come from G2-style real profiles, so the
community structure is not artificial.

| Axis | Levels |
|---|---|
| Community size | 50, 150, 300 genomes |
| Strain heterogeneity | species-unique, 20 % species with 2-3 related strains |
| Depth | 2M, 10M read pairs per sample (40 samples per experiment) |
| Scenario | null, spiked (20 % DA, balanced) |
| Seeds | 5 independent per cell |

That is 3 x 2 x 2 x 2 x 5 = 120 experiments; run **60** (drop the strain axis
for the 300-genome cell) if the Phase 0 timing says the budget does not fit.

For every experiment, score **three feature tables** side by side from the same reads:

1. **Oracle**: reads mapped to the true source genomes (no binning).
2. **MAG**: one co-assembly per experiment, MetaBAT2 catalogue (existing path).
3. **Reference-based**: mapping to a deliberately incomplete reference, with 70 %
   of the genomes present, to emulate the usual database-gap situation.

Ablations on a subset (one design cell, 5 seeds): a second binner (SemiBin2 or
COMEBin), a completeness/contamination filter (none vs >=50/<10 vs >=90/<5), and
per-sample versus co-assembly if time allows.

Outputs for claim 1: DA metrics on each table, and MAG-quality covariates (recovery
by true abundance, completeness, contamination, share of DA genomes not recovered,
bins mixing DA and non-DA sources). The key figure is power and FDR as a function of
*true abundance and coverage*, which shows where MAG construction creates
missing-not-at-random zeros.

### Tier 3: real metagenomes with known truth

Two designs that avoid synthetic reads entirely:

- **3a. Real nulls.** Take healthy-cohort metagenomes (target >= 100 samples from
  one study, one body site, similar depth). Subsample reads to a common depth,
  assemble and bin, then draw **200 random label splits** (balanced, and confounded
  with sequencing batch where metadata allows). Any discovery is false. This
  measures type I error with real structure at no assumption cost.
- **3b. Read-level implants.** In the group-B samples, up-sample or down-sample the
  reads aligned to chosen MAGs by a known factor (extending Nearing-style
  implants to MAG space), conserving library size. Truth is exact by
  construction, and everything else (assembly, bins, noise, batch) is real. Use
  4 effect sizes x 3 DA fractions x 10 seeds.
- **3c. Mock communities** (ZymoBIOMICS-type with log-distributed ratios, public CAMI
  data) as a sanity anchor that quantification recovers known ratios. Small and
  cheap; not used for method ranking.

Dataset accessions are **not** pinned here. Phase 0 selects and records them with
URL, licence, checksum and read counts after checking they are public and
redistributable. I have not verified any accession.

### Tier 4: real disease cohorts

One IBD set (e.g. IBDMDB/HMP2 plus one independent IBD cohort) and one CRC set (>= 3
public cohorts). Without ground truth, report only (i) **cross-cohort replication**
of each method's discoveries (rank agreement and replicated-discovery rate
between independent cohorts, per method), and (ii) agreement between methods.
Replication is a proxy for reliability, not accuracy, and the paper must say so.
Also report the effect of the raw-data null from Tier 3a on how many "discoveries"
a method makes with permuted labels in the same cohort.

## 3. Methods and fairness rules

- Methods and variants come from `config/methods.yaml` (16 active). **Freeze the
  active list and parameters before Tier 2.** Variants (e.g. ALDEx2/3 gamma 0 vs
  0.5) are declared up front and reported as variants of a family.
- Defaults are the package defaults unless the registry already records a
  justified setting; no per-tier tuning.
- LOCOM is blocked (fails on its own example data). Report as "excluded: package
  failure", do not substitute.
- Add a **predefined consensus rule** (e.g. discovery if >= k of the eligible
  methods agree) only if written down before Tier 2 and scored with the same
  independent seeds.
- Include trivial baselines that any method must beat: a plain t-test/Wilcoxon on
  relative abundance, and a random-score ranker. Without a naive baseline a
  ranking is uninterpretable.
- A method that fails, times out or returns empty is a counted outcome and cannot
  win. Report the failure rate by condition.

## 4. Metrics (frozen before running)

Primary:

1. **Null calibration**: false-positive rate at alpha 0.05 on null p-values;
   any-false-discovery rate at nominal 0.05 with exact binomial interval.
2. **FDR control**: realised FDP vs nominal 0.05 on spiked experiments, cluster
   bootstrap CI.
3. **Power at realised FDR**: sensitivity at nominal 0.05, and AP from raw
   p-values (not q-values, which would not be comparable across packages).

Secondary: effect-size error (MAE of log2FC on the matching endpoint), rank
stability across seeds, recall on source genomes versus MAGs, runtime and memory,
and failure rate.

**Pre-registered eligibility**, applied per method x endpoint x condition: the
upper bound of the FDR confidence interval must not exceed 0.10 (twice nominal),
null FP rate must be within 0.03-0.07, and the failure rate must be below 5 %. A method
that passes only on part of the design space is reported as conditionally
eligible, with the region where it holds, not as a global winner.

Statistics: a mixed-effects model (method x factor interactions, random effect
for experiment) on FDP and sensitivity, with multiplicity-adjusted pairwise
contrasts. Report which axis explains the most variance (variance decomposition
on method, depth, DA fraction, direction, generator, tier). This is the headline
quantification of "what matters". Rankings are shown with rank-uncertainty
intervals from the bootstrap, never as bare ranks.

## 5. Compute plan on 60 cores

Always run with `--cores 60` and cap memory below 400 GB (the machine has 503 GB;
other users share it). Pin per-job thread limits so total demand never exceeds 60.

Estimates below rest on the existing runs and must be re-measured in Phase 0.
Calibration run: 964 jobs, 49,241 core-seconds for 40 cases = **~20 core-min per case**
(16 methods, sources plus MAGs). The 20-genome pilot: simulation ~5 min per sample and
assembly ~5.5 min for 20 samples of 200k read pairs.

| Tier | Work | Core-hours | Wall time at 60 cores |
|---|---|---|---|
| 1 | ~6,000 cases x 20 core-min (sources-only cases cost less) | ~2,000 | **~35 h** |
| 2 | 60 experiments, ~70 core-h each (simulation 40 samples, one co-assembly with 16 threads, mapping, binning, DA on three tables) | ~4,200 | **~3 days**, 3-4 experiments in parallel |
| 3a | 1 assembly + 200 splits x DA | ~800 | ~14 h |
| 3b | Implants, 120 DA sets on one assembly | ~600 | ~10 h |
| 4 | Real cohorts, assembly of subsampled reads, DA | ~1,500 | ~1 day |
| Re-runs, ablations, slack (30 %) | | ~2,700 | ~2 days |
| **Total** | | **~11,800** | **~8-9 days of continuous compute** |

Constraints to plan around:

- **Disk is the binding resource**: the home filesystem is 97 % full with 2.2 TB
  free. A 40-sample Tier 2 experiment at 10M read pairs writes roughly 100+ GB of
  FASTQ. Run at most 3 experiments at a time, delete reads after counting
  (already the repo's policy), and keep only raw counts, lengths, library totals
  and compact audit evidence. Ask for a scratch volume before Tier 2/4.
- Tier 1 is embarrassingly parallel and should start first, while the Tier 2
  manifests are prepared.
- Assembly is the long pole; schedule it before the cheap DA jobs so cores don't
  idle behind it.
- Keep a hard `storage.budget_gb` in each config and let `src/magician/storage.py`
  reject over-budget runs, as it already does.

## 6. Phases and timeline (about 10 weeks, compute is ~9 days of it)

| Week | Phase | Deliverable |
|---|---|---|
| 1 | **0. Freeze and measure** | Timing runs for every tier; pick datasets and record accessions; freeze `methods.yaml`, metrics, eligibility rules, analysis code. Commit a tagged `preregistration` revision. |
| 1-2 | **1. Generators** | G2 and G3 implemented and tested, including a realism check (feature prevalence, mean-variance, sparsity vs real data, with figures). |
| 2-3 | **2. Tier 1 run** | Full grid with 60 cores; calibration plots; first look at claim 2 (truth lanes) and claim 4. |
| 3-5 | **3. Tier 2 run** | Oracle/MAG/reference tables for all experiments; claim 1 analysis. |
| 5-7 | **4. Tier 3** | Real nulls, implants, mock. |
| 7-8 | **5. Tier 4** | Cohort replication analysis. |
| 8-9 | **6. Analysis and figures** | Variance decomposition, decision guide, ablations. |
| 9-10 | **7. Writing and archive** | Manuscript, Zenodo release (code, frozen env locks, matrices, results), container image. |

## 7. Risks and how to handle them

| Risk | Handling |
|---|---|
| Simulated reads (wgsim) lack realistic error and coverage bias | Compare against InSilicoSeq or CAMISIM error models on one cell; Tier 3 uses real reads. |
| Simulation circularity (generator favours methods sharing its model) | Three generators plus real-data tiers; report generator as a factor in the variance decomposition. |
| Findings turn out to be a truth-lane artefact | Run every method on all applicable lanes and publish the full method x lane table, so the sensitivity is visible. |
| MAG-to-source assignment errors confound "MAG cost" | Report assignment ambiguity separately; use the oracle table as the reference rather than the assignment. |
| Real datasets are unsuitable (depth, access) | Phase 0 decides this before any run; fallback is fewer cohorts, with the limitation stated. |
| Disk exhaustion | Section 5 constraints; hard storage budget. |
| Researcher degrees of freedom | Pre-registration tag; any post-hoc analysis is labelled exploratory. |
| Method installation failure | Frozen environments from `tools/bootstrap_envs.py`; excluded methods listed with reasons. |

## 8. Reproducibility package (required for acceptance)

- Frozen conda locks and checksums for every method environment (already supported).
- Tagged `preregistration` and `final` git revisions; run manifests with seeds.
- A public archive with matrices, truth, scores and logs (compact; reads excluded
  because they are regenerable from seeds), plus a container image.
- A single command per tier, e.g.
  `python run_magician.py --configfile config/paper_tier1.yaml --cores 60`.
- Figures regenerated from the archive by one script.

## 9. Moving to the experiment server

The runs happen on another server, so the repository must be self-sufficient:

1. Clone, then `conda env create -n magician --file environment.yaml`.
2. `python tools/bootstrap_envs.py --all` builds the frozen method environments
   (needs network once; later runs touch no network).
3. `pytest -q` and `python run_magician.py --dry-run` must pass before any tier starts.
4. Record in the Phase 0 report: core count, RAM, filesystem free space, Snakemake,
   Python and R versions. Re-measure the timings in Section 5 there and rescale.
5. Set `CONDA_PKGS_DIRS` and `XDG_CACHE_HOME` to a large local volume.
6. Prior results (`results_*`, `archive/`) are not in git; they stay on the old
   server and are only needed for the Phase 0 comparison.

## 10. Decisions I need from you

1. **Datasets for Tiers 3-4**: do you have preferred public cohorts, or should
   Phase 0 choose them (IBD and CRC are my default suggestion)?
2. **Scratch space**: can we get a volume with >= 5 TB? Without it, Tier 2 must
   shrink to ~3 experiments in flight.
3. **Scope**: full plan (~9 days compute), or a reduced paper using Tier 1 + Tier 2
   + Tier 3a only (~5 days)? I recommend the reduced version first, and adding
   Tier 3b/4 only if claim 1 holds.
4. **Strain axis and second binner**: keep or drop to save about 30 % of Tier 2.
