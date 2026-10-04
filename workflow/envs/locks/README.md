# Exact resolved environments

linux-64/*.txt records the installed package builds with SHA256 checksums.
The workflow YAMLs pin the scientific tools; these explicit files archive the
complete resolved dependency graph from validation, including transitive packages.

Recreate an archived environment with:

```bash
conda create --prefix cache/replay_bio --file workflow/envs/locks/linux-64/bio.txt
```

For strict replay through Snakemake, copy each explicit file beside its YAML as
NAME.linux-64.pin.txt before environment creation (for example bio.linux-64.pin.txt).
Snakemake then installs those exact builds. Other platforms need their own resolution.

## Method environments

`tools/bootstrap_envs.py` builds one frozen environment per method family and writes
`ENV.lock.json` next to these files: the resolved conda builds, the pinned CRAN
tarball, Bioconductor release tarball or Git commit SHA, its SHA256, the installed R
version and the installed R package versions. A pinned source that failed to install
is recorded as `NOT INSTALLED` and the bootstrap fails, so an environment is never
reported complete because the script reached the end.

```bash
python tools/bootstrap_envs.py --all --force
```

Method environments are separate from the frozen baseline environments
(`r_core.txt`, `r_extended.txt`) because Bioconductor releases and their R versions
must match: MaAsLin 3 1.4.0 is a Bioconductor 3.23 release package and therefore does
not drop into the R 4.5 environment used by the older methods. Environments that
build R packages from source carry the conda-forge toolchain explicitly.
