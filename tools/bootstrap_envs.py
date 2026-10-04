#!/usr/bin/env python3
"""Build pinned method environments and record exactly what went into them.

Conda provides R and the compiled dependencies. Packages that are not distributed
through conda are installed from pinned CRAN tarballs or pinned Git commit SHAs into
the finished environment, so the frozen set is reproducible and the network is never
touched inside a DA rule.

  python tools/bootstrap_envs.py --env r_aldex3
  python tools/bootstrap_envs.py --all --dry-run
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tarfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
ENV_DIR = ROOT / "workflow/envs"
LOCK_DIR = ENV_DIR / "locks"
# Pinned sources, per environment. A CRAN entry is (package, version); a git entry is
# (package, repository, ref); a Bioconductor entry is (package, release, version).
# Versions, release branches and commits are recorded in the lock file. An environment
# receives only the sources it declares, so unrelated packages cannot leak in.
CRAN = {
    "ALDEx3": "1.3.1",
    "radEmu": "2.3.2.0",
    "GUniFrac": "1.9",
    "MicrobiomeStat": "1.4",
}
GIT = {
    # Pinned by commit SHA, not by branch: a SHA cannot move under us.
    "LOCOM": ("https://github.com/yijuanhu/LOCOM.git",
              "76f8c8d6f4d1920c020a4a183d59ef162a45cdb3"),
}
# Bioconductor release packages. Bioconductor releases and their R versions must
# match, so maaslin3 lives in its own environment.
BIOCONDUCTOR = {
    "maaslin3": ("3.23", "1.4.0"),
}
# Packages whose conda builds disagree on the Matrix ABI are rebuilt from the
# exact version conda resolved, so the version stays pinned and the ABI matches.
SOURCE_REINSTALL = {
    "r_linda": ["lme4"],
}
ENV_SOURCES = {
    "r_aldex3": ["ALDEx3"],
    "r_rademu": ["radEmu"],
    "r_zicoseq": ["GUniFrac"],
    "r_linda": ["MicrobiomeStat"],
    "r_locom": ["LOCOM"],
    "r_maaslin3": ["maaslin3"],
    "r_adapt": [],
}
# Packages conda already provides are recorded but not reinstalled from source.
CONDA_PROVIDED = {"ADAPT", "phyloseq"}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def environment(target):
    """A conda environment is not on PATH unless it is activated, and compiling an
    R package from source needs its toolchain and its shared libraries."""
    target = Path(target)
    env = dict(os.environ)
    env["PATH"] = f"{target / 'bin'}{os.pathsep}{env.get('PATH', '')}"
    env["R_HOME"] = str(target / "lib/R")
    for variable in ("LDFLAGS", "LIBRARY_PATH"):
        env.pop(variable, None)
    return env


def run(command, target=None, **kwargs):
    print("+", " ".join(str(part) for part in command), flush=True)
    if target is not None:
        kwargs.setdefault("env", environment(target))
    return subprocess.run([str(part) for part in command], check=True, **kwargs)


def bioconductor_url(package, release, version):
    return (f"https://bioconductor.org/packages/{release}/bioc/src/contrib/"
            f"{package}_{version}.tar.gz")


def cran_urls(package, version):
    """The current CRAN directory first, then the archive for superseded releases.

    R prints a dash in versions as a dot, but CRAN archive file names keep the
    dash, so both spellings are tried.
    """
    dotted = str(version)
    # R prints the release dash as a dot, but CRAN archive file names keep the dash
    # after the first two components (1.1-35.5 on disk, 1.1.35.5 in R).
    parts = dotted.split(".")
    dashed = f"{parts[0]}.{parts[1]}-{'.'.join(parts[2:])}" if len(parts) > 2 else dotted
    return [f"https://cran.r-project.org/src/contrib/{package}_{dotted}.tar.gz",
            f"https://cran.r-project.org/src/contrib/Archive/{package}/{package}_{dotted}.tar.gz",
            f"https://cran.r-project.org/src/contrib/{package}_{dashed}.tar.gz",
            f"https://cran.r-project.org/src/contrib/Archive/{package}/{package}_{dashed}.tar.gz"]


def current_source(package):
    """The current CRAN version, for recording when a release is not archived."""
    return f"https://cran.r-project.org/src/contrib/{package}_{{version}}.tar.gz"


def download(urls, target):
    if Path(target).exists():
        return Path(target)
    for url in ([urls] if isinstance(urls, str) else urls):
        try:
            print("download", url, flush=True)
            # Bioconductor and some mirrors reject the default urllib agent.
            request = urllib.request.Request(url, headers={"User-Agent": "magician-benchmark/1.0"})
            with urllib.request.urlopen(request) as response, Path(target).open("wb") as handle:
                handle.write(response.read())
            return Path(target)
        except urllib.error.HTTPError as exc:
            print(f"  {url} -> {exc.code}", flush=True)
    raise SystemExit(f"Could not download any of {urls}")


def env_name(stem):
    return f"magician_{stem}"


def create_conda_env(stem, prefix, dry_run, force=False):
    """Conda provides R and the compiled dependencies. Reuse an existing environment."""
    spec = ENV_DIR / f"{stem}.yaml"
    if not spec.exists():
        raise SystemExit(f"Missing environment specification: {spec}")
    target = Path(prefix) / env_name(stem)
    if (target / "bin/Rscript").exists() and not force:
        print(f"{target} already exists; reusing it", flush=True)
        return target
    if dry_run:
        print(f"would create {target} from {spec.name}", flush=True)
        return target
    run(["conda", "env", "create", "-y", "-p", str(target), "-f", str(spec)])
    return target


def install_sources(target, stem, dry_run):
    """Install only this environment's pinned sources into a finished environment."""
    target = Path(target)
    work = target / "vendor"
    work.mkdir(parents=True, exist_ok=True)
    wanted = ENV_SOURCES.get(stem)
    if wanted is None:
        raise SystemExit(f"No declared sources for environment {stem}; add it to ENV_SOURCES")
    installed, script = {}, ['options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 1L)']
    for package in wanted:
        if package in CRAN:
            version = CRAN[package]
            archive = download(cran_urls(package, version), work / f"{package}_{version}.tar.gz")
            installed[package] = dict(version=version, source=archive.name, sha256=sha256(archive))
            script.append(f'install.packages({str(archive)!r}, repos = NULL, type = "source")')
        elif package in BIOCONDUCTOR:
            release, version = BIOCONDUCTOR[package]
            archive = download(bioconductor_url(package, release, version),
                               work / f"{package}_{version}.tar.gz")
            installed[package] = dict(version=version, release=release, source=archive.name,
                                      sha256=sha256(archive))
            script.append(f'install.packages({str(archive)!r}, repos = NULL, type = "source")')
        elif package in GIT:
            repository, ref = GIT[package]
            checkout = work / package
            if not checkout.exists():
                if dry_run:
                    print(f"would clone {repository} at {ref}", flush=True)
                    continue
                run(["git", "clone", "--quiet", repository, str(checkout)])
            # Always re-pin: a reused checkout must not drift with its branch.
            if not dry_run:
                run(["git", "-C", str(checkout), "fetch", "--quiet", "origin", ref])
                run(["git", "-C", str(checkout), "checkout", "--quiet", ref])
            # A checkout can ship stale objects from another machine; a source build
            # must always start from clean C++.
            for stale in sorted(Path(checkout, "src").glob("*.o")) + sorted(Path(checkout, "src").glob("*.so")):
                stale.unlink()
            commit = subprocess.run(["git", "-C", str(checkout), "rev-parse", "HEAD"],
                                    capture_output=True, text=True, check=True).stdout.strip()
            installed[package] = dict(commit=commit, source=repository, sha256=None)
            # A Git checkout is installed from its local path, never from the URL.
            script.append(f'install.packages({str(checkout)!r}, repos = NULL, type = "source")')
        else:
            raise SystemExit(f"{stem}: unknown source {package}")
    reinstall = SOURCE_REINSTALL.get(stem, [])
    for package in reinstall:
        reported = subprocess.run(
            [target / "bin/Rscript", "--vanilla", "-e",
             f'cat(as.character(packageVersion("{package}")))'],
            capture_output=True, text=True, env=environment(target))
        version = reported.stdout.strip()
        if not version:
            raise SystemExit(f"{stem}: cannot read installed {package} version for reinstall")
        archive = download(cran_urls(package, version), work / f"{package}_{version}.tar.gz")
        print(f"pinned reinstall source for {package} {version}", flush=True)
    if reinstall:
        script.append(f'reinstall <- c({", ".join(repr(p) for p in reinstall)})')
        script.append('versions <- vapply(reinstall, function(p) as.character(packageVersion(p)), character(1))')
        script.append('for (i in seq_along(reinstall)) {')
        script.append(f'  archive <- file.path({str(work)!r}, paste0(reinstall[i], "_", versions[i], ".tar.gz"))')
        script.append('  if (!file.exists(archive)) stop("missing pinned reinstall tarball: ", archive)')
        script.append('  install.packages(archive, repos = NULL, type = "source")')
        script.append('}')
        script.append(f'writeLines(paste(reinstall, versions, sep = "="), {str(work / "reinstalled.txt")!r})')
    (work / "install.R").write_text("\n".join(script) + "\n")
    if dry_run:
        print("\n".join(script), flush=True)
        return installed
    if not wanted:
        print("no external sources declared; conda supplies every package", flush=True)
        return installed
    # Only pinned local sources are installed, so nothing else is fetched here.
    run([target / "bin/Rscript", "--vanilla", str(work / "install.R")], target=target)
    return installed


def record(target, stem, installed):
    """Resolve installed versions and the exact resolved conda builds for the lock.

    A package that failed to install is recorded as such: an environment is never
    reported complete because the script reached the end.
    """
    target = Path(target)
    packages = subprocess.run(
        [target / "bin/Rscript", "--vanilla", "-e",
         'cat(paste(sort(intersect(c("ALDEx3", "radEmu", "GUniFrac", "MicrobiomeStat", '
         '"maaslin3", "ADAPT", "LOCOM"), rownames(installed.packages()))), collapse="\n"))'],
        capture_output=True, text=True)
    installed_here = [name for name in packages.stdout.split() if name]
    versions = {}
    for name in installed_here:
        reported = subprocess.run(
            [target / "bin/Rscript", "--vanilla", "-e",
             f'cat(as.character(packageVersion("{name}")))'],
            capture_output=True, text=True)
        versions[name] = reported.stdout.strip() or "unknown"
    missing = sorted(set(installed) - set(installed_here))
    for name in missing:
        versions[name] = "NOT INSTALLED"
    resolved = subprocess.run(["conda", "list", "-p", str(target), "--json"],
                              capture_output=True, text=True)
    builds = []
    if resolved.returncode == 0:
        builds = sorted({f"{item['name']}={item['version']}={item['build_string']}"
                         for item in json.loads(resolved.stdout)})
    lock = dict(environment=stem, path=str(target),
                r_version=subprocess.run([target / "bin/Rscript", "--vanilla", "-e",
                                          'cat(as.character(getRversion()))'],
                                         capture_output=True, text=True).stdout.strip(),
                pinned_sources=installed, r_packages=versions,
                complete=not missing, missing=missing, conda_builds=builds)
    LOCK_DIR.mkdir(parents=True, exist_ok=True)
    (LOCK_DIR / f"{stem}.lock.json").write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n")
    print(f"wrote {LOCK_DIR / (stem + '.lock.json')}", flush=True)
    if missing:
        raise SystemExit(f"{stem}: pinned sources did not install: {', '.join(missing)}")
    return lock


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env", action="append", default=[], help="Environment stem, e.g. r_aldex3")
    parser.add_argument("--all", action="store_true", help="Every declared environment")
    parser.add_argument("--prefix", default=str(ROOT / "cache/envs"))
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true", help="Recreate environments from scratch")
    args = parser.parse_args()
    stems = args.env or ([p.stem for p in sorted(ENV_DIR.glob("r_*.yaml"))] if args.all else [])
    if not stems:
        parser.error("pass --env or --all")
    for stem in stems:
        print(f"=== {stem} ===", flush=True)
        target = create_conda_env(stem, args.prefix, args.dry_run, args.force)
        if args.dry_run:
            continue
        installed = install_sources(target, stem, args.dry_run)
        record(target, stem, installed)


if __name__ == "__main__":
    sys.exit(main())