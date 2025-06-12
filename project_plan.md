# Comprehensive Experiment Design — Benchmarking Differential-Abundance (DA) Methods on MAG-Coverage Tables

*Designed for ANCOM-BC2 focus, but extensible to all major DA tools*

## 0. Objectives & Hypotheses

- **Primary**: Quantify how input scale (raw counts, length-normalised counts, TPM/RPKM) affects FDR and power of ANCOM-BC2 and competing DA tools on MAG data.
- **Secondary**: Evaluate whether the ANCOM-BC2 pseudocount-sensitivity filter improves result reliability.
- **Exploratory**: Identify which DA method(s) remain robust under:
  - Extreme library-size variation
  - Genome-length bias
  - High sparsity/zero-inflation
  - Longitudinal (paired) sampling

## 1. Data Sets

| Track | Samples | Ground truth | Why include |
|-------|---------|--------------|-------------|
| Synthetic (core) | 20 × 2 groups per scenario | Spiked log2-fold-changes in 15% of MAGs | Allows full control; replicate 100× per scenario |
| Mock spike-in | 6 (Zymo mock + stool, 1% vs 5%) | Known eight-strain fold-changes | Noise from wet-lab & assembly captured |
| Public cohort | ≥30 (e.g. PRJNA389280 IBD WGS) | Unknown | Stress-test methods on real complexity |

## 2. Synthetic Community Generation (CAMISIM)
*Source: [microbiomejournal.biomedcentral.com](https://microbiomejournal.biomedcentral.com)*

- **Genome library**: 150 bacterial genomes (1–8 Mb) from GTDB; stratify by length quartile.
- **Community design**: Dirichlet-multinomial abundances (α = 0.1) to mimic over-dispersion.
- **Differential design**: Random 15% MAGs receive fold-changes of 1.5×, 2×, 4× (factor).
- **Depth profiles**: Log-normal read counts; coefficient of variation (CV) = 0.2, 0.5, 1.0 (factor).
- **Paired scenario**: Use STEMSIM to introduce within-strain evolution over time for 10 hosts, generating pre/post samples with correlated counts. *Source: [academic.oup.com](https://academic.oup.com)*

*Replicate every (depth × effect × sparsity) cell 100 times.*

## 3. Assembly & Binning

- **Assembler**: MEGAHIT with meta presets; retain contigs ≥ 1 kb.
- **Binning**: MetaBAT2 + VAMB; refined with DAS Tool; dereplicate at 95% ANI using dRep.
- **Taxonomy**: GTDB-Tk.
- Keep an assembly log to identify contigs/MAGs missing from some samples → structural zeros.

## 4. Coverage Quantification (CoverM)
*Source: [academic.oup.com](https://academic.oup.com)*

```bash
coverm genome --bam-files {sample}.bam \
              --genome-fasta-directory MAGs/ \
              --methods count,tpm,rpkm \
              --min-covered-fraction 0.9 \
              --output {sample}.tsv
```

Outputs three matrices per dataset:
- Counts
- Counts_per_kbp
- TPM

## 5. Feature Matrix Curation

- Filter MAGs present in ≥10% of samples or ≥10 mapped reads (per raw counts).
- Store genome length, GC%, completeness, and contamination for covariate analyses.

## 6. DA Methods & Configurations

| Tool | Input(s) | Special settings |
|------|----------|------------------|
| [ANCOM-BC2](https://nature.com) | counts, counts/kbp, TPM | p_adj_method="BH", offset="library_size"; run with and without sensitivity=TRUE filter |
| DESeq2 | counts | sfType="poscounts" for zero handling |
| edgeR | counts | robust dispersion (estimateGLMRobustDisp) |
| ALDEx2 | counts | 128 MCMC instances |
| metagenomeSeq | counts | CSS normalisation |
| corncob | counts | test="wald", zero_model=TRUE |
| ZicoSeq | counts | defaults |
| LinDA / MaAsLin2 (paired scenario) | counts | random-effect for subject |

All methods orchestrated via benchdamic to ensure identical preprocessing. *Source: [academic.oup.com](https://academic.oup.com)*

## 7. Experimental Factors (full-factorial grid)

| Factor | Levels |
|--------|--------|
| Sequencing-depth CV | 0.2 / 0.5 / 1.0 |
| Genome-length spread | narrow (1–3 Mb) / wide (1–8 Mb) |
| Sparsity (proportion zeros) | 30% / 60% / 80% |
| Effect size | 1.5× / 2× / 4× |
| Sample size | 10 / 20 / 40 per group |
| Correlation structure | independent / paired (pre–post) |
| Normalisation | counts / counts/kbp / TPM |
| Pseudocount filter (ANCOM-BC2) | off / on |

Total synthetic scenarios = 3 × 2 × 3 × 3 × 3 × 2 × 3 × 2.
Each scenario → 100 replicates → robust power/FDR estimates.

## 8. Workflow Automation (Snakemake sketch)

```yaml
rule all:
  input: expand("eval/{scenario}/{rep}/metrics.tsv", ...)

rule simulate_reads:     # CAMISIM
rule assemble:
rule bin_mags:
rule map_reads:          # minimap2
rule coverm_quant:
rule build_matrices:     # counts / kbp / TPM
rule run_benchdamic:     # wrapper calls each DA method
rule evaluate:
rule aggregate_results:
```

Docker/Conda image stores all versions; GitHub repo with CI tests.

## 9. Evaluation Metrics

| Category | Measure |
|----------|---------|
| Error control | FDR at q ≤ 0.05 |
| Sensitivity | TPR (recall), auPRC |
| Effect accuracy | Pearson r between estimated and true logFC |
| Robustness | ΔFDR/ΔTPR before vs. after ANCOM-BC2 sensitivity filter |
| Concordance | Jaccard & Kendall-W between method result sets |
| Paired fit | Type-I error under null in correlated scenario |

**Statistical comparisons**:
- Wilcoxon signed-rank (counts vs TPM within method)
- McNemar (pseudocount filter)
- Linear mixed models to quantify impact of each factor on FDR and TPR

## 10. Visual & Report Outputs

- Box-plots of FDR/TPR per method across scenarios
- Heat-map of method agreement
- Effect-size scatter (true vs. estimated)
- Upset plots for overlap of significant MAGs across methods
- R Markdown → benchmark_report.html + slide deck for stakeholders

## 11. Validation Steps

- **Mock spike-in**: Expect detection of eight Zymo strains at correct fold-change; check method ranking stability.
- **Public cohort**: Compare findings across tools; highlight disagreements explained by depth/length effects.

## 12. Resources & Timeline

- **Compute**: 64 CPU-hours per 100-rep synthetic block (parallelisable); 128 GB scratch.
- **Week 1–2**: Build pipeline; pilot 10 rep runs; QC assemblies.
- **Week 3–4**: Full factorial runs; evaluate; visualise.
- **Week 5**: Mock & public dataset analysis; draft report.
- **Week 6**: Review, optimise slides, release Zenodo snapshot.

## Key Design Choices & Rationale

- Ground-truth control via CAMISIM assures unbiased power/FDR measurement. *Source: [microbiomejournal.biomedcentral.com](https://microbiomejournal.biomedcentral.com)*
- benchdamic standardises method invocation, avoiding implementation drift. *Source: [academic.oup.com](https://academic.oup.com)*
- CoverM produces consistent count, length-normalised, and TPM outputs for direct comparison. *Source: [academic.oup.com](https://academic.oup.com)*
- ANCOM-BC2 sensitivity filter explicitly tested because its developers warn about depth-normalised inputs. *Source: [nature.com](https://nature.com)*
- Factorial design ensures we can model how each bias (depth, sparsity, length) alters each tool's behaviour.

*With this step-by-step plan you can generate reproducible evidence on which DA method—and which preprocessing—delivers the most trustworthy results for MAG-coverage studies, while directly addressing Alejandra's concerns about ANCOM-BC2 power and its pseudocount robustness.*










Sources

