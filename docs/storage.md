# Storage and compute policy

Default budget: **100 decimal GB**; set storage.budget_gb to 200 for 200 GB.
The new output directory and managed environment/package cache count towards it.
Existing unrelated data and external source libraries are preserved. Source lengths
are charged in the peak work estimate. The historical results symlink is preserved
in archive/2026-10-02/results_legacy_link; use a new output directory.

Archive moves do not duplicate or compress data and therefore do not reduce its
disk footprint. The dated archive is preserved historical material outside the
active run's managed budget. Filesystem free-space checks still reflect its space.

Preflight checks actual generated allocations, source lengths, conservative scratch
overhead, retained tables, and a 15 GB environment reserve plus 2 GB headroom. Oversized designs fail
before simulation. Heavy processing of cases runs sequentially: each case's reads
wait for both previous abundance matrices, after the previous large intermediates
have been consumed. DA fits may overlap the next case. Increasing seeds does not
multiply peak FASTQ/BAM storage. Reference copies, FASTQs, BAMs, assemblies, depth tables, MAG FASTAs
and PAF are temp() outputs, deleted after consumers finish. One genome scratch pair
is reused during simulation; SAM is streamed; no pooled FASTQ copy is made.

The local profile caps concurrent declared memory, disk reservations and IO jobs.
Filesystem locks stay enabled. Avoid --nolock in production and --notemp unless
deliberately debugging. Heavy external commands check managed usage/free space every
five seconds and terminate their process groups on excess; transient overshoot is
possible. Jobs share a cache-footprint snapshot, refreshed on environment/package
directory changes or every five minutes, avoiding repeated metadata scans on network
filesystems. Environments are installed before jobs; avoid installing additional
environments concurrently with a run.
Conda installation precedes jobs and relies on the reserve, not a quota.
A strict physical cap requires a filesystem/project quota. Resource reservations
also are not OS RAM limits; cluster limits belong in the executor.

For a 200 GB budget, increase the scheduler reservation with
--resources disk_mb=180000 if a job's conservative reservation exceeds the default
80000 MB. Keep space for environments and retained outputs. A larger storage budget
does not increase available RAM; reduce per-case depth/sample count when needed.

The launcher puts packages in cache/packages, environments in cache/conda, and small
application caches in cache/xdg. Direct Snakemake users should set CONDA_PKGS_DIRS
and XDG_CACHE_HOME accordingly. No shared cache is cleaned
automatically. Assembler/binning scratch is cleaned on success or Python exceptions.

Retained: truth, allocations, metadata, hashes/configuration; source/MAG matrices and
filter flags; catalogue assignments/matching evidence; per-method discoveries,
statuses/session versions; feature evaluations, scores, rankings, intervals; HTML,
SVG/PNG figures; logs and per-rule runtime/RSS benchmarks. These compact tables allow
DA/evaluation reruns without assembly. Resource summaries distinguish output/cache
sizes. Job max_rss is not concurrent total RAM. The final snapshot precedes the report
itself and its unfinished benchmark.
