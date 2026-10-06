#!/usr/bin/env python3
"""Fetch the empirical donor table for generator G2 and for G1/G3 calibration.

  python tools/fetch_real_profiles.py

Source: curatedMetagenomicData (MetaPhlAn 3 profiles), study AsnicarF_2021,
ExperimentHub EH5458. One study, one body site, one sample per subject:
healthy adult stool controls from Great Britain without current antibiotic use.
Species-level relative abundances are kept as measured, including their zeros.
Downloads are verified against pinned SHA256 digests; the committed table and
its provenance record are what the generators read, so runs need no network.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import subprocess
import sys
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from magician.io import write_json

STUDY = "AsnicarF_2021"
PROFILE_URL = ("https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/curatedMetagenomicData/"
               "2021-03-31/AsnicarF_2021/2021-03-31.AsnicarF_2021.relative_abundance.rda")
PROFILE_SHA256 = "89c635061b357583351ca33f520a72d0efced4a2963e73182e529228c8395c54"
METADATA_COMMIT = "f88ac2ee84df7094d04cbf225a9d4679fc4b8986"
METADATA_URL = ("https://github.com/waldronlab/curatedMetagenomicData/raw/"
                f"{METADATA_COMMIT}/data/sampleMetadata.rda")
METADATA_SHA256 = "168c8650d2893dcb49930b9b85a4187637775c9aae0a8023efb4cca52dc73444"
FILTER = ('study_name == "AsnicarF_2021" & body_site == "stool" & study_condition == "control" & '
          'age_category == "adult" & country == "GBR" & '
          '(is.na(antibiotics_current_use) | antibiotics_current_use == "no")')
EXTRACT = r'''
args <- commandArgs(trailingOnly = TRUE)
x <- get(load(args[1]))
load(args[2])
keep <- subset(sampleMetadata, %s)
if (anyDuplicated(keep$subject_id)) stop("more than one sample per subject")
species <- grepl("\\|s__", rownames(x)) & !grepl("\\|t__", rownames(x))
x <- x[species, colnames(x) %%in%% keep$sample_id, drop = FALSE]
x <- x[rowSums(x > 0) > 0, order(colnames(x)), drop = FALSE]
x <- x[order(rownames(x)), , drop = FALSE]
out <- data.frame(feature = sub(".*\\|s__", "", rownames(x)), x, check.names = FALSE)
write.table(out, args[3], sep = "\t", quote = FALSE, row.names = FALSE)
''' % FILTER


def download(url, digest, target):
    urllib.request.urlretrieve(url, target)
    found = hashlib.sha256(Path(target).read_bytes()).hexdigest()
    if found != digest:
        raise SystemExit(f"{url}: SHA256 {found} differs from the pinned {digest}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="resources/real_profiles")
    p.add_argument("--rscript", default="cache/envs/magician_r_core/bin/Rscript")
    args = p.parse_args()
    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        download(PROFILE_URL, PROFILE_SHA256, tmp / "profiles.rda")
        download(METADATA_URL, METADATA_SHA256, tmp / "metadata.rda")
        (tmp / "extract.R").write_text(EXTRACT)
        subprocess.run([str(ROOT / args.rscript), "--vanilla", str(tmp / "extract.R"),
                        str(tmp / "profiles.rda"), str(tmp / "metadata.rda"), str(tmp / "table.tsv")],
                       check=True)
        text = (tmp / "table.tsv").read_bytes()
    table = out / "asnicar_2021_species.tsv.gz"
    # A fixed mtime keeps the compressed file byte-identical between fetches.
    with open(table, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as handle:
        handle.write(text)
    lines = text.decode().splitlines()
    write_json(out / "provenance.json", dict(
        study=STUDY, resource="ExperimentHub EH5458", profiler="MetaPhlAn 3 (curatedMetagenomicData 2021-03-31)",
        profile_url=PROFILE_URL, profile_sha256=PROFILE_SHA256,
        metadata_url=METADATA_URL, metadata_sha256=METADATA_SHA256, sample_filter=FILTER,
        level="species", unit="relative abundance in percent of the profiled community; zeros as measured",
        n_samples=len(lines[0].split("\t")) - 1, n_features=len(lines) - 1,
        table=table.name, table_sha256=hashlib.sha256(text).hexdigest(),
        licence="curatedMetagenomicData is distributed under Artistic-2.0; cite Asnicar et al. 2021 "
                "and Pasolli et al. 2017 when using these profiles"))
    print(f"wrote {table}: {len(lines) - 1} species x {len(lines[0].split(chr(9))) - 1} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
