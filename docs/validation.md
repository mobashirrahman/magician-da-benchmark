# Validation record

The four-genome execution below is historical. Its results now reside in
[archive/2026-10-02/results/smoke4](../archive/2026-10-02/results/smoke4), and its
console logs in archive/2026-10-02/validation/smoke4. Original paths in the commands
describe the run before archiving. Current validation logs are in docs/validation.
The larger pilot is configured in [config/pilot20.yaml](../config/pilot20.yaml); its
outcome is recorded under [Twenty-genome artificial pilot](#twenty-genome-artificial-pilot).

Validation was performed on Linux x86_64 on 2026-10-02 with Snakemake 9.11.2.
The reproducible default fixture uses four artificial, fragmented genomes, six
independent samples, a null experiment and a spiked experiment. Its rankings are
software smoke evidence, not a scientific recommendation for real communities.

## Regression tests and workflow checks

```bash
pytest -q
snakemake -s workflow/Snakefile --workflow-profile profiles/testing --lint
python run_magician.py --dry-run
python run_magician.py --configfile config/production.yaml --dry-run
```

The 13 regression tests cover reproducible independent noise, conserving non-DA
abundance, exact integer allocations, fragment counting, abundance normalization,
empty MAG matrices, ambiguous/multiple genome matches, tied average precision,
unrecovered DA sources, failed-method eligibility, the literal null scenario,
zero-recall exclusion, configuration validation and storage rejection. The
Snakemake lint check passed; both default and production DAGs were constructible.
CI repeats the regression, lint and default DAG checks without downloading genomes.

## Actual method adapters

```bash
python tests/integration/validate_adapters.py --conda-prefix cache/conda
python tests/integration/validate_adapters.py --rscript r_core=/path/to/r_core/bin/Rscript
```

The harness resolves one Rscript per registry environment, preferring an
environment whose R actually carries the method's package. Every fixture is
designed to break something specific: positive and negative effects, structural
zeros and a zero library, a reserved `other` category, a single feature, an empty
catalogue, permuted sample order, reversed group labels, incomplete metadata, an
undeclared input metric and an explicit timeout.

The original seven adapters passed on the same reproducible 80-feature, 12-sample
count fixture, with no failures. Assertions check successful fitting, preservation
of original feature IDs including hyphens, declared endpoints and units, and the
direction of known positive effects. An empty MAG input produced the valid,
explicit `empty` result. Count-only methods refused a TPM input with a recorded
`failed` status and a non-zero exit. Reversed group labels and incomplete sample
metadata were rejected rather than silently accepted. A job whose declared budget
could not fit its work was reported as `timeout`, never replaced by another
method. This exercises the actual packages, not substitute implementations or
mocked predictions.

| Adapter | Validated package version |
|---|---|
| Wilcoxon CLR | Base R 4.5.3 |
| DESeq2 | 1.46.0 |
| edgeR | 4.8.2 |
| limma-voom | 3.66.0 |
| ALDEx2 (gamma 0 and 0.5) | 1.42.0 |
| ANCOM-BC2 (ANCOMBC) | 2.14.0 |
| MaAsLin2 | 1.22.0 |
| ALDEx3 (gamma 0 and 0.5) | 1.3.1 |
| radEmu | 2.3.2.0 |
| ZicoSeq (GUniFrac) | 1.9 |
| LinDA (MicrobiomeStat) | 1.4 |
| ADAPT | 1.4.0 |
| MaAsLin 3 (abundance and corrected lanes) | 1.4.0 |

The adapter code for MaAsLin 3, ALDEx3, ADAPT, radEmu, LinDA, ZicoSeq and LOCOM is
written against the frozen distributions named in the method registry, and each
adapter validates the installed release's own signature instead of assuming
argument names. The API mismatches found while reading those distributions are
recorded here so they are not reintroduced: ALDEx2 1.42 calls its count matrix
`reads`; ALDEx3 1.3.1 takes `aldex(Y, X, ...)` with `nsample`, `streamsize` and
gamma inside the scale function; ZicoSeq 1.9 takes `feature.dat` as features by
samples and returns per-feature F statistics with `p.raw` and `p.adj.fdr`;
MicrobiomeStat 1.4 exposes `linda(feature.dat, meta.dat, formula, ...)` with
`adaptive` zero handling and a `log2FoldChange`/`pvalue`/`padj` output table; ADAPT
1.4.0 returns an S4 `DAresult` whose full table is the `details` slot; radEmu
2.3.2.0 returns a `coef` table keyed by taxon index and covariate and requires an
explicit `test_kj`; LOCOM takes samples-by-features counts with a binary outcome
and returns `p.otu`/`q.otu` alongside a separate `p.global`; MaAsLin 3 1.4.0 uses
`normalization`, `transform`, `correction`, `median_comparison_abundance` and
`evaluate_only`, with natural-log coefficients under `transform = "LOG"`.
Per-method validation results are recorded below as each pinned environment
completes. metaGEENOME remains deferred: its published wrapper exposes no isolated
differential abundance entry point. LOCOM is blocked: the pinned checkout fails with
`subscript out of bounds` inside its own permutation machinery on its own throat
example data in the frozen R 4.4 environment, so it cannot be scheduled until the
package works in a frozen build. Its adapter uses the verified `locom(otu.table,
Y, ...)` API and stays `planned`.

## Extension-plan implementation

The method registry, extended result schema, split R adapters, durable raw
quantification, matrix-only execution mode, four truth lanes, declared zero
policies and the reserved `other` category, stratified reporting, matrix screening
designs and the real-genome manifest tooling are implemented and covered by 46
regression tests. The 8-replicate matrix screening bundle runs the full matrix-only
pipeline end to end with no simulation, assembly or binning rule: import, derive,
truth lanes, methods, scoring, aggregation and report. Gate-by-gate status lives in
the [extension plan](extension_plan.md#implementation-status).

The scientific Python environment contains Python 3.11.13, NumPy 1.26.4,
pandas 2.2.3 and matplotlib 3.9.4. Genome processing used wgsim 1.0, MEGAHIT 1.2.9,
minimap2 2.28, samtools 1.21 and MetaBAT2 2.17. Complete resolved package builds
and checksums are in [workflow/envs/locks](../workflow/envs/locks/README.md).

## End-to-end execution

```bash
python run_magician.py --cores 4
```

The run exercises controlled reads, shared assembly/binning, separate source and
MAG quantification, truth matching, all compatible method/metric combinations,
evaluation, ranking and the HTML/figure/resource report. Historical per-rule logs
and benchmarks are retained in archive/2026-10-02/results/smoke4/logs and benchmarks.
Each method retains a status JSON and R session information. Small catalogues can
make a statistical fit impossible; such outcomes remain failures in the report
and cannot improve a method's ranking by omitting the difficult experiment.

The full pinned-environment execution completed **111/111 jobs**. All **52 DA jobs**
succeeded: 13 compatible method/metric combinations for each of two catalogues
(sources and MAGs) across both experiments. The null experiment recovered four
MAGs/four sources. The spiked experiment recovered three MAGs/three sources,
including both truly changed sources. Both temporary work trees were empty after
normal completion; retained results, logs and benchmarks total approximately
**1.9 MB**. The managed cache is approximately **8.2 decimal GB**, including all
four environments and shared packages; obsolete validation environments were removed.
The largest sampled job RSS was **423.35 MB** for this tiny fixture, not a production
memory estimate or a measurement of concurrent total RAM.

## Resume and analysis-only changes

After completion, the default dry run reported nothing to do. An overlay retaining
only wilcoxon_clr and edger scheduled five jobs: two evaluations, aggregation, the
report and its final target. Changing alpha from 0.05 to 0.10 scheduled the 52 DA
jobs and five summary/target jobs. Neither change scheduled reference preparation,
simulation, assembly, mapping, binning or matrix construction. Preflight uses
retained lengths and designs, with analysis settings excluded from its dependency
fingerprint. Each new case waits for the previous pair of abundance matrices.
Report formatting is independent of preflight's implementation fingerprint.

The historical commands' console logs are retained locally under
archive/2026-10-02/validation/smoke4; generated logs and scientific results are
ignored by Git. Exact historical outcomes and resources are in
archive/2026-10-02/results/smoke4/benchmark/report.html, scores.tsv,
rule_resources.tsv and storage.json.

## Twenty-genome artificial pilot

```bash
python run_magician.py --configfile config/pilot20.yaml --cores 8
```

The pilot was re-run end to end under the registry schema (205/205 jobs, no
failures). Endpoints are now reported separately, so each catalogue and input names
one candidate per endpoint: ALDEx2 for the counts CLR lane, ANCOM-BC2 for the
counts relative-abundance lane, limma-voom for the source read-fraction lane, and
Wilcoxon CLR or MaAsLin2 for the normalized CLR lanes. No method made a null false
discovery on sources or MAGs. DESeq2 and edgeR keep their characteristic
high-recall, higher-FDR profile on counts. With one seed these picks remain
`smoke_only`; the stratification by endpoint is the point, not a new ranking.

Run on 2026-10-02 with 20 artificial 350 kb genomes (50 contigs each), 10 samples
per group, 200,000 nominal read pairs per sample, one null and one spiked case
(seed 1001, four changed sources, nominal 4x multiplier). The console log is
docs/validation/pilot20_run.log; results are in results_pilot20 (ignored by Git).

Redistribution conserves the changed subset's relative mass, so actual effects
are log2 fold changes +0.511332 for synthetic_005 and synthetic_012, and
−3.488668 for synthetic_014 and synthetic_020. synthetic_012 is the unrecovered
changed source. Each case contains 3,961,605 actual read pairs, for 7,923,210 total.

The run completed **195/195 jobs** in about 85 minutes wall time with no failed or
retried job. All **52 DA jobs** reported success. A subsequent dry run reported
nothing to do, and the 13 regression tests still pass.

| Case | Source genomes | MAGs | Matched MAGs | Recovered sources | Recovered changed sources |
|---|---|---|---|---|---|
| null_s1001 | 20 | 18 | 18 | 18 | 0 of 0 |
| spiked_s1001 | 20 | 17 | 17 | 17 | 3 of 4 |

No MAG was ambiguous, unassigned or matched to several sources. One truly changed
source was not recovered as a MAG in the spiked case, so MAG-level source recall
cannot exceed 0.75 there.

No method made a false discovery in the null case, on sources or MAGs. In the
spiked case, false discovery proportions at alpha 0.05 were:

| Input | Method | Sources: FDP / recall | MAGs: FDP / recall among recovered |
|---|---|---|---|
| counts | limma-voom | 0.00 / 0.75 | 0.25 / 1.00 |
| counts | ANCOM-BC2 | 0.00 / 0.50 | 0.00 / 0.67 |
| counts | DESeq2 | 0.20 / 1.00 | 0.25 / 1.00 |
| counts | edgeR | 0.25 / 0.75 | 0.25 / 1.00 |
| counts | ALDEx2 | 0.33 / 1.00 | 0.50 / 1.00 |
| counts | Wilcoxon CLR | 0.43 / 1.00 | 0.57 / 1.00 |
| counts | MaAsLin2 | 0.40 / 0.75 | 0.71 / 0.67 |
| relative | Wilcoxon CLR | 0.00 / 0.50 | 0.00 / 0.67 |
| relative | MaAsLin2 | 0.00 / 0.25 | 0.00 / 0.33 |

A method is eligible only if every fit succeeded and both its mean spiked FDP and
its null probability of any false discovery are at most alpha (0.05). Eligible
methods are then ordered by source recall. recommendations.json therefore names
ANCOM-BC2 for MAG counts, limma-voom for source counts and Wilcoxon CLR for
relative abundance on both catalogues. TPM and RPKM have no eligible method,
because Wilcoxon CLR and MaAsLin2 both exceeded 0.05 there; their rows are in
results_pilot20/benchmark/ranking.tsv.

Every significant true-positive call had the correct direction. With 13 to 20 tested features
and at most four true changes, the first false call alone raises the proportion
from 0 to 0.20 or more, so eligibility here rests on a single call in a single
seed. The picks are labelled `smoke_only` and do not establish a general method ranking. MaAsLin2
tested fewer features than the other methods (13 to 18 against 17 to 20) because
of its internal filtering.

Measured resources:

| Rule | Jobs | Total seconds | Longest job (s) | Largest RSS (MB) |
|---|---|---|---|---|
| simulate_reads | 40 | 12,377 | 367 | 175 |
| differential_abundance | 52 | 2,377 | 63 | 417 |
| assemble | 2 | 653 | 331 | 982 |
| map_contigs | 40 | 356 | 14 | 621 |
| map_sources | 40 | 221 | 10 | 431 |
| all other rules | 19 | 98 | 52 | 124 |

Read simulation took 77% of the 16,082 measured job seconds, because exact
allocation runs one wgsim call per source contig (1,000 per sample). The largest
job RSS was **981.81 MB** (assembly), which is a per-job figure, not concurrent RAM.
Both work trees were empty after completion: no FASTQ, BAM, PAF or assembly files
remain. Final retained outputs total **20,412,772 bytes (20.4 MB)** against a
preflight allowance of 40 MB. Logs account for 19.4 MB; matrices, predictions,
evaluations and reports occupy under 1 MB. The reusable managed environment/cache
snapshot is **8.169 decimal GB**. The conservative preflight total estimate is
45.7 GB including reserves. Peak scratch occupancy was not recorded, so the 28.7 GB peak-work estimate
is untested; the storage guard, which stops jobs above the 100 GB budget, did not
trigger.

The current CLR adapters use a fixed 0.5 pseudocount across input units. Its effect
differs between proportions and count/TPM values. The
[extension plan](extension_plan.md) includes unit-aware preprocessing sensitivity
and matching truth coordinates before real-genome rankings are interpreted.
New methods described there remain planned; this pilot tested the existing seven.

## Extension plan implementation

Executed on 2026-10-04 on Linux x86_64 with Snakemake 9.11.2, Python 3.11.13 and
R 4.5.3/4.4.3 environments. This section records what was run for the
[extension plan](extension_plan.md); the historical sections above describe the runs
they supersede.

### Regression tests and workflow checks

```bash
pytest -q
snakemake -s workflow/Snakefile --lint
python run_magician.py --dry-run
python run_magician.py --configfile config/pilot20.yaml --dry-run
python run_magician.py --configfile cache/screening/config.yaml --dry-run
```

**46 regression tests pass.** They cover reproducible independent designs, exact
integer allocations, fragment counting, library accounting, length-corrected metrics,
empty MAG matrices, the raw/derived split, reserved `other` category totals, ambiguous
and multiple genome matches, tied average precision, endpoint-matched truth lanes,
unmatched effect units, q-only adapters, exact binomial null intervals, method
families versus variants, registry and zero-policy validation, matrix bundle import
and rejection, screening replicate independence, and genome manifest selection.

The lint check passes. The simulation DAG builds **121 jobs**, of which **52 are DA
jobs** — the same 13 compatible method/metric combinations per catalogue and case as
before, because ALDEx2's package default gamma is the variant that was retained.

### Real adapters

```bash
python tests/integration/validate_adapters.py --conda-prefix cache/conda
```

All previously validated adapters pass against the real packages, with **0 failures**:

| Adapter | Validated package version | Endpoint | Tests |
|---|---|---|---|
| Wilcoxon CLR | base R 4.5.3 | `clr` | 80 |
| DESeq2 | 1.46.0 | `read_fraction` | 80 |
| edgeR | 4.8.2 | `read_fraction` | 80 |
| limma-voom | 3.66.0 | `read_fraction` | 80 |
| ALDEx2 (gamma 0.5) | 1.42.0 | `clr` | 80 |
| ALDEx2 (gamma 0) | 1.42.0 | `clr` | 80 |
| ANCOM-BC2 (ANCOMBC) | 2.14.0 | `relative_abundance` | 80 |

Each adapter was additionally exercised on fixtures designed to break something
specific, all passing:

| Fixture | Outcome |
|---|---|
| Structural zeros and a zero library | success, 80 tests, no fabricated values |
| Reserved `other` category | success, reserved rows present and never tested |
| Single feature | explicit status, no invented features |
| Empty catalogue | `empty`, valid empty result |
| Permuted sample order | explicit reordering, no samples dropped |
| Group labels disagreeing with sample order | **rejected**: reversed contrast |
| Metadata missing a sample | **rejected**: no silent intersection |
| Undeclared input metric | **rejected** with the accepted metrics listed |
| Explicit timeout on 4000 features | **timeout** recorded, not replaced |

Two failures found this way were real defects and are fixed. ALDEx2 1.42 names its
count argument `reads`, not `x`, so the adapter now reads the installed release's own
signature. The runner enforces its deadline with a supervised fork and records a
`timeout` status, because `setTimeLimit` does not interrupt a busy adapter.

### New adapters

Each adapter below was written against its frozen distribution and validated
against the real package on the same fixture family, with **0 failures** in the
final full run. Environments are frozen by `tools/bootstrap_envs.py`; the exact
resolved builds, checksums and pinned sources are in `workflow/envs/locks/`.

| Adapter | Frozen package | Endpoint | Notes |
|---|---|---|---|
| ALDEx3 (gamma 0 and 0.5) | 1.3.1 (CRAN) | `clr` | `aldex(Y, X, ...)` with `nsample`, `streamsize`, gamma in the scale function; effect is the mean of the Treatment coefficient draws |
| MaAsLin 3 (abundance) | 1.4.0 (Bioconductor 3.23) | `relative_abundance` | `evaluate_only="abundance"`, natural-log coefficients converted to log2 |
| MaAsLin 3 (corrected) | 1.4.0 (Bioconductor 3.23) | `reference_relative` | median comparison against the typical feature |
| ADAPT | 1.4.0 (Bioconductor) | `reference_relative` | S4 `DAresult`, full table from the `details` slot, learned reference set reported |
| radEmu | 2.3.2.0 (CRAN) | `reference_relative` | explicit `test_kj` for the Treatment coefficient only; natural-log estimates converted to log2 |
| LinDA (MicrobiomeStat) | 1.4 (CRAN) | `clr` | `linda(feature.dat, meta.dat, formula, ...)`, per-variable output table, adaptive zero handling |
| ZicoSeq (GUniFrac) | 1.9 (CRAN) | `clr` | features-by-samples matrix, per-feature F statistics with `p.raw`/`p.adj.fdr` |

Confirmed directions on the validation fixture: MaAsLin 3 abundance reports about
+2.3/−3.4 log2 for 4x/0.25x effects; ADAPT reports +2.0/−2.0 against the typical
change; radEmu's converted estimates match the simulated log2 effects. ZicoSeq's
statistic is non-negative, so direction is meaningless for it and is not asserted.

LOCOM is blocked: the pinned checkout fails with `subscript out of bounds` inside
its own permutation machinery on its own throat example data in the frozen R 4.4
environment, so it cannot be scheduled until the package works in a frozen build.
Its adapter uses the verified `locom(otu.table, Y, ...)` API and stays `planned`.

### Matrix-only mode

```bash
python tools/matrix_screening.py --designs baseline low_count \
    --null-replicates 2 --spiked-replicates 2
python run_magician.py --configfile cache/screening/config.yaml --cores 4
```

The DAG contains **100 jobs and no simulation, assembly, mapping or binning rule**:
8 bundle imports, 8 derivations, 8 truth lanes, 8 merges, 56 DA fits, 8 evaluations,
the aggregation and the report. All 100 jobs completed. Screening designs are
predeclared in `src/magician/screening.py`; replicates are independent Dirichlet
multinomial draws with disjoint seeds per design and scenario.

Observed on 8 replicates at alpha 0.05, with the exact binomial interval on the null
rate:

| Input | Endpoint | Method | Spiked FDR | Source recall | Null any-FD | Eligible |
|---|---|---|---|---|---|---|
| counts | clr | wilcoxon_clr | 0.027 | 0.850 | 0.00 [0.00, 1.00] | yes |
| counts | read_fraction | edger | 0.047 | 0.975 | 0.00 [0.00, 1.00] | yes |
| counts | read_fraction | limma_voom | 0.024 | 0.969 | 0.00 [0.00, 1.00] | yes |
| counts | read_fraction | deseq2 | 0.040 | 0.969 | 0.25 [0.01, 1.00] | no |
| counts | clr | aldex2 | 0.000 | 0.000 | 0.00 [0.00, 1.00] | no |
| counts | relative_abundance | ancombc2 | 0.000 | 0.781 | 0.00 [0.00, 1.00] | yes |
| relative | clr | wilcoxon_clr | 0.025 | 0.888 | 0.00 [0.00, 1.00] | yes |

These are four replicates per scenario: far too few to substantiate 5% error
control, which is exactly why the null rate carries an exact binomial interval
reaching 1.00. ALDEx2's zero recall with a high average precision is a genuine
result on these designs, not a scoring artefact: it declared no discoveries, so it
fails the positive-source-recall requirement and cannot win by discovering nothing.

### Calibration across screening designs

```bash
python tools/matrix_screening.py --designs baseline overdispersed low_count dropout \
    --null-replicates 5 --spiked-replicates 5
python run_magician.py --configfile cache/calibration/config.yaml --cores 8
```

40 independent replicates (20 null, 20 spiked) across four designs ran all 16
active methods end to end: **964/964 jobs, no failures**. Every method fit every
case. Findings at alpha 0.05, with paired-bootstrap FDR intervals and exact
binomial null intervals:

| Endpoint | Eligible | Not eligible, and why |
|---|---|---|
| `clr` | wilcoxon_clr, maaslin2, aldex2_g0, aldex3_g0, aldex3, aldex2, zicoseq | linda: null any-false-discovery rate 0.10 |
| `read_fraction` | limma_voom | deseq2: spiked FDR 0.29; edger: null rate 0.10 |
| `reference_relative` | rademu, adapt | maaslin3_corrected: null rate 0.15 |
| `relative_abundance` | ancombc2 | maaslin3_abundance: null rate 0.15 |

The gamma comparison is informative rather than decorative: ALDEx2 with gamma 0.5
declares almost nothing on these designs (source recall 0.004, FDR 0.000) while
gamma 0 recovers half the true changes at FDR 0.004, and ALDEx3 shows the same
pattern. DESeq2's high recall comes with a 0.29 spiked FDR here. With 20 null
replicates the exact null intervals remain wide; the planned 100-null-per-design
scale is what tightens them, and the machinery for it is the same command with
larger replicate counts.

### Real-genome manifests

`tools/real_genomes.py check` and `freeze` are implemented and tested against
fixtures: short references, unrecognised accessions, missing files, group sizes that
cannot be satisfied, and duplicate IDs are all rejected explicitly, and every frozen
manifest records accession, version, path, length, SHA256, provenance, the selection
rule and the manifest's own checksum. No genome has been downloaded; real-genome
execution is still ahead.
