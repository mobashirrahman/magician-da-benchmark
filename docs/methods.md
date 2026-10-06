# Experimental and statistical contract

The default target is **relative genomic abundance**, with independent Control and
Treatment samples. Each seed produces null and spiked experiments with the same
baseline, noise and library sizes. Allocations, expected group compositions, hashes
and seeds are retained. Sample and library noise are separate lognormal draws. Read
allocation is proportional to genomic abundance times genome length, using largest
remainders. wgsim adds substitution errors and no source mutations. Its simplified
uniform error model does not reproduce every platform. Default artificial genomes
test software.

The original read-level generator redistributes abundance within the DA subset while conserving its total mass.
Non-DA genomes retain their expected relative abundance. Truth records the actual
final log2 fold change after adjustment. Finite-sample means differ from expectations.
Samples are independent biological draws, not duplicated technical replicates.

Corrected Tier 1 uses the versioned protocol in [corrected_tier1.md](corrected_tier1.md).
Balanced implants conserve the changed subset's mass. One-directional implants
instead close the whole community, so unselected features also change on the
relative-abundance scale. Null and spiked cases share a noise mechanism; Tier 1
uses equal feature lengths, multinomial library accounting, and excludes unsupported
batch confounding. G2 resamples measured species profiles; G1 and G3 are synthetic
models calibrated on the same donor table. Tier 1 discoveries are scored against the
implanted set, with each method's endpoint lane reported beside it.

## Endpoints and truth lanes

A method is scored against the truth lane for the endpoint it declares. A correct
association on one measurement scale is not a false positive when another scale is
being tested, and the reverse. Every lane is computed from the latent design,
including genome lengths and the declared denominator, in
[src/magician/truth.py](../src/magician/truth.py).

| Lane | Endpoint | Effect | Denominator |
|---|---|---|---|
| `relative_genomic_abundance` | `relative_abundance` | log2 relative abundance | length-corrected total over tested features |
| `expected_read_fraction` | `read_fraction` | log2 read/DNA fraction | total expected read or DNA mass |
| `clr_log_ratio` | `clr` | log2 CLR difference | geometric mean of log2 values per sample |
| `reference_relative_change` | `reference_relative` | log2 change against the reference | median feature, learned set, or a declared feature |
| `absolute_abundance` | `absolute_abundance` | log2 copies | expected genome copies per sample |

Labels (`is_da`) come from the design's **target** in the endpoint's coordinates, never from the
noisy per-sample draws around it: under the null the targets are identical and every
measured difference is finite-sample variation. Effect sizes come from the
**expectation of the measured statistic**, which accounts for sample-level closure,
library variation and the finite counts actually drawn.

CLR labels include changes caused by a shifted geometric reference; read-fraction
truth includes its total-DNA denominator. Tier 1 source recall uses the matching
endpoint's positive sources. Null p-value calibration uses raw p-values, while
the rate of any adjusted discovery is reported separately.

The absolute abundance lane is only available when load information exists
(`experiment.load_information`). Otherwise it is emitted and reported as unavailable
rather than silently dropped.

Consequences for the count models: DESeq2, edgeR and limma-voom fit normalized
counts, so their estimand is the **read/DNA fraction** and their registry endpoint is
`read_fraction`. With equal genome lengths this coincides with relative abundance;
with real genomes of different lengths the two lanes differ, and each method is
scored against its own. ANCOM-BC2 applies explicit sample-scaled compositional
correction, so it keeps the `relative_abundance` endpoint.

## Declared zero handling

Zero handling is a declared, tracked variant in
[config/zero_handling.yaml](../config/zero_handling.yaml), never an implicit choice.

| Policy | Zero model | Pseudocount |
|---|---|---|
| `pseudocount_measurement` | additive, in the metric's own units | 0.5 |
| `pseudocount_count` | additive, in expected-count units | 0.5 |
| `dirichlet` | posterior mean of a Dirichlet fitted to the same sample | none |

A fixed 0.5 offset means three different things on proportions, TPM and counts; that
is a unit problem, not a biological difference, so the count-scale variant exists to
show it. Each registry entry names exactly one policy, and the CLR truth lane is
computed in the matching coordinates for every declared policy.

`analysis.reserve_other` adds a single `other` category holding the filtered counts,
so library totals are preserved for Dirichlet-based methods. It changes the reference
geometry of every CLR method, so it is a declared variant, its truth lane is computed
with the same geometry, and the category is never a tested discovery.

## Quantification and recovery

MEGAHIT reads all libraries directly, without a pooled FASTQ copy. One assembly and
one MetaBAT2 run per experiment use all sample BAMs in a multi-sample depth table.
Source reads are mapped separately. MAG counts reuse contig alignments. Counts are
properly paired primary first-mate fragments passing the mapping-quality threshold.
Secondary, supplementary, unmapped and low-quality alignments are excluded.

Raw quantification is durable: `quantify_features` writes counts, lengths and library
totals once, including per-sample mapped, catalogue, unassigned and unmapped pairs.
`derive_matrices` recomputes filtering, the reserved category and normalization from
those counts, so changing a threshold never needs deleted upstream reads.

TPM and relative genomic abundance divide counts by feature length and normalize
within samples. RPKM uses total simulated pairs as the library denominator. MAG
fractions are conditional on the recovered catalogue. All methods/metrics use the
same truth-blind count/prevalence filter; filtered rows remain untested in outputs.
Normalized abundance is never converted to pseudo-counts.

minimap2 aligns MAG contigs to known sources for evaluation only. Assignments use
identity and unioned query intervals. Low coverage, tied source coverage and no
alignment are explicit statuses. Secondary hits preserve ambiguity, with a 50-hit
cap: highly redundant communities may need more exhaustive matching. Aligned source
breadth and MAG fraction are sequence-based measures, not CheckM estimates.
Multiple MAGs assigned to one source remain visible.

## Methods

Every method is declared once in [config/methods.yaml](../config/methods.yaml) with
its accepted inputs, hypothesis, endpoint, effect units, reference, correction
family, environment, seed policy, thread cap and citation. `method_family` groups
statistical models so that several settings of one model are never reported as
independent methods. Compatibility between configured metrics and each method is
validated before any job is built.

| Method | Family | Inputs | Endpoint | Effect interpretation |
|---|---|---|---|---|
| Wilcoxon CLR | `clr_wilcoxon` | counts/TPM/RPKM/relative | `clr` | Log2 CLR difference |
| MaAsLin2 | `linear_clr` | counts/TPM/RPKM/relative | `clr` | Log2 CLR coefficient |
| DESeq2 | `count_nb` | counts | `read_fraction` | Normalized log2 fold change |
| edgeR | `count_nb` | counts | `read_fraction` | Normalized log2 fold change |
| limma-voom | `count_nb` | counts | `read_fraction` | Normalized log2 fold change |
| ANCOM-BC2 | `count_compositional` | counts | `relative_abundance` | Natural-log effect converted to log2 |
| ALDEx2 (gamma 0, 0.5) | `dirichlet_clr` | counts | `clr` | Log2 CLR difference |
| ALDEx3 (gamma 0, 0.5) | `scale_dirichlet_clr` | counts | `clr` | Log2 CLR difference under a scale model |
| MaAsLin 3 (abundance) | `hurdle_linear` | relative | `relative_abundance` | Non-zero log relative abundance |
| MaAsLin 3 (corrected) | `hurdle_linear` | relative | `reference_relative` | Change against the typical feature |
| ADAPT | `censored_reference` | counts | `reference_relative` | Censored log2 ratio, learned reference set |
| radEmu | `typical_change` | counts | `reference_relative` | Log2 change against the typical taxon |
| LinDA | `bias_corrected_clr` | counts/relative | `clr` | Bias-corrected log ratio |
| ZicoSeq | `compositional_permutation` | counts | `clr` | Slope on transformed ratios, **not** a log2 fold change |
| LOCOM | `compositional_logistic` | counts | `reference_relative` | Log ratio to a declared reference feature |

The shared contract in
[workflow/scripts/da/contract.R](../workflow/scripts/da/contract.R) enforces, for every
method: exact sample and feature ID agreement with no silent intersection, an
explicit reordering when sample order differs, exactly two groups with at least two
samples each, a Treatment-against-Control factor that cannot be silently reversed,
integer counts for count-only methods, a declared zero policy and transform, finite
prepared input, and stable internal feature IDs so packages cannot rewrite them.

Native adjusted values are never re-adjusted. Where a package returns no adjusted
values the contract applies the declared correction over its own raw p-values; where a
package returns raw p-values, average precision over raw p-values is the common
metric. An adapter that can expose only q-values is labelled with a separate
q-ranked average precision, and raw p-values are never reverse engineered.

Each job has its own log, benchmark, environment, thread, memory, scratch and status.
Statuses are `success`, `empty`, `failed`, `unavailable` and `timeout`. A timeout is
enforced by the runner and reported as a timeout; it is never replaced by another
statistical method. Empty catalogues, convergence failures and unavailable packages
stay visible and make the affected method ineligible.

ALDEx2 and ALDEx3 keep gamma as a declared variant of one Dirichlet family, set before
any discovery is inspected, and are validated against the installed release's own
signature rather than assumed. ZicoSeq and LOCOM declare permutation counts as
variants and are checked for resolution rather than chosen on speed. radEmu's
identification constraint is declared: it cannot recover absolute abundance in
individual samples from sequencing alone. metaGEENOME is deferred until its
differential abundance component can be called in isolation.

## Scoring and selection

At alpha, report TP/FP/FN/TN, FDR, precision, recall, false-positive rate and
direction agreement. Threshold-grouped average precision handles tied scores without
ordering advantages, and is computed over raw p-values when available and over
adjusted values otherwise, with the source recorded. Null runs report any false
discovery. Fold-change error is computed only where a method's reported units match
its truth lane's coordinates; otherwise it is unavailable, not rescaled.

Feature scoring is conditional on uniquely assigned MAGs. Source recall includes
unrecovered true DA sources as misses. Source discoveries are deduplicated when
multiple bins match. Unassigned significant MAGs are counted separately. Empty
catalogues are valid, with empty matrices and no eligible MAG method.

Candidates require all runs successful/evaluable, positive source recall, mean spiked
FDR and null probability of any false discovery at most alpha, and no significant
unassigned MAGs. Ordering uses source recall, average precision over adjusted values,
then FDR; equally scoring methods are reported as tied candidates instead of selecting
a winner by name, and variants of one family are reported as a single candidate.
Bootstrap intervals resample whole independent seeds, keeping method comparisons
paired. The null rate carries an exact binomial interval. Artificial runs or fewer
than five seeds are labelled `smoke_only`. These are empirical criteria, not proofs of
control. Repeated real-genome simulations still describe their tested scenarios, not a
universally best method.

## Primary references

- [Snakemake rules/resources/temporary outputs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)
- [MEGAHIT](https://github.com/voutcn/megahit), [minimap2](https://github.com/lh3/minimap2), [wgsim](https://github.com/lh3/wgsim)
- [ANCOM-BC2 tutorial](https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html)
- [ALDEx2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html)
- [MaAsLin2](https://github.com/biobakery/Maaslin2), [MaAsLin 3 manual](https://bioconductor.org/packages/release/bioc/vignettes/maaslin3/inst/doc/maaslin3_manual.html)
- [ALDEx3](https://cran.r-project.org/package=ALDEx3)
- [ADAPT manual](https://bioconductor.org/packages/release/bioc/vignettes/ADAPT/inst/doc/ADAPT-manual.html)
- [radEmu](https://statdivlab.github.io/radEmu/)
- [LinDA](https://link.springer.com/article/10.1186/s13059-022-02655-5)
- [ZicoSeq](https://link.springer.com/article/10.1186/s40168-022-01320-0)
- [LOCOM](https://pmc.ncbi.nlm.nih.gov/articles/PMC9335309/)
