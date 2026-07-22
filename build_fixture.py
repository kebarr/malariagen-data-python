"""
Build a local file:// fixture for AnophelesSnpData, using real VCF/zarr data
already downloaded for sample_set 1276-AD-BD-ALAM-VMF00156 (AdirusWRAIR2 assembly).

Restricted to the 2 chromosome-scale contigs: KB672490, KB672491.
Genome *sequence* is a stub (placeholder bases, real contig lengths) since the
real AdirusWRAIR2 FASTA needs a login I don't have.
Site list (POS/REF/ALT) comes from one representative VCF (same for every
sample, since GATK was run with GENOTYPE_GIVEN_ALLELES/EMIT_ALL_SITES).
Genotypes are stacked directly from the 28 already-downloaded per-sample
.zarr.zip GT arrays (no need to re-parse any VCF for GT).
"""
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

REPO = Path("/root/malariagen-data-python")
DOWNLOADED = REPO / "downoaded_data"
FIXTURE_ROOT = REPO / "local_fixture"
BUCKET = FIXTURE_ROOT / "vo_adir_release"

CONTIGS = ["KB672490", "KB672491"]
CONTIG_LENGTHS = {"KB672490": 22947322, "KB672491": 6774089}
SAMPLE_SET = "1276-AD-BD-ALAM-VMF00156"
RELEASE = "1.0"
RELEASE_PATH = "v1"
REPRESENTATIVE_VCF = DOWNLOADED / "VBS46299-6321STDY9453299.vcf.gz"


def get_sample_ids():
    return sorted(p.name.split(".zarr.zip")[0] for p in DOWNLOADED.glob("*.zarr.zip"))


def write_config():
    config = {
        "PUBLIC_RELEASES": [RELEASE],
        "GENOME_FASTA_PATH": "reference/genome/adirwrair2/genome.fa",
        "GENOME_FAI_PATH": "reference/genome/adirwrair2/genome.fa.fai",
        "GENOME_ZARR_PATH": "reference/genome/adirwrair2/genome.zarr",
        "GENOME_REF_ID": "AdirusWRAIR2",
        "GENOME_REF_NAME": "Anopheles dirus (WRAIR2)",
        "CONTIGS": CONTIGS,
        "SITE_MASK_IDS": [],
    }
    BUCKET.mkdir(parents=True, exist_ok=True)
    with open(BUCKET / "config.json", "w") as f:
        json.dump(config, f, indent=4)


def write_manifest(sample_ids):
    release_dir = BUCKET / RELEASE_PATH
    release_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.DataFrame(
        {
            "sample_set": [SAMPLE_SET],
            "sample_count": [len(sample_ids)],
            "study_id": ["1276-AD-BD-ALAM"],
            "study_url": [
                "https://www.malariagen.net/network/where-we-work/1276-AD-BD-ALAM"
            ],
            "terms_of_use_expiry_date": ["2026-12-31"],
            "terms_of_use_url": [
                "https://malariagen.github.io/vector-data/vobs/vobs.html#terms-of-use"
            ],
        }
    )
    manifest.to_csv(release_dir / "manifest.tsv", sep="\t", index=False)


def write_sample_metadata(sample_ids):
    src = REPO / "metadata/general" / SAMPLE_SET / "samples.meta.csv"
    df = pd.read_csv(src)
    df_ds = df[df["sample_id"].isin(sample_ids)].reset_index(drop=True)
    assert len(df_ds) == len(sample_ids), "sample metadata subset mismatch"
    dst = BUCKET / RELEASE_PATH / "metadata/general" / SAMPLE_SET / "samples.meta.csv"
    dst.parent.mkdir(parents=True, exist_ok=True)
    df_ds.to_csv(dst, index=False)
    return df_ds["sample_id"].tolist()


def write_stub_genome():
    path = BUCKET / "reference/genome/adirwrair2/genome.zarr"
    path.parent.mkdir(parents=True, exist_ok=True)
    root = zarr.open(path, mode="w")
    rng = np.random.default_rng(42)
    bases = np.array([b"A", b"C", b"G", b"T"])
    for contig in CONTIGS:
        size = CONTIG_LENGTHS[contig]
        seq = rng.choice(bases, size=size, replace=True)
        root.create_dataset(name=contig, data=seq, chunks=(2_000_000,))
    zarr.consolidate_metadata(path)


def extract_sites(vcf_path, contigs):
    """Stream CHROM/POS/REF/ALT for `contigs` out of a whole-genome VCF,
    stopping as soon as we pass the target contigs (VCF body is sorted in
    contig-header order, and these are the first two, largest, contigs)."""
    awk_prog = (
        '$1=="%s" || $1=="%s" {print $1"\\t"$2"\\t"$4"\\t"$5; seen=1; next} '
        "seen==1 {exit}" % (contigs[0], contigs[1])
    )
    cmd = f"zcat {vcf_path} | awk -F'\\t' '{awk_prog}'"
    proc = subprocess.Popen(
        ["bash", "-c", cmd], stdout=subprocess.PIPE, text=True, bufsize=1 << 20
    )

    data = {c: {"pos": [], "ref": [], "alt": []} for c in contigs}
    for line in proc.stdout:
        chrom, pos, ref, alt = line.rstrip("\n").split("\t")
        d = data[chrom]
        d["pos"].append(int(pos))
        d["ref"].append(ref)
        d["alt"].append(alt.split(","))
    proc.stdout.close()
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"awk/zcat pipeline failed with code {ret}")

    out = {}
    for c in contigs:
        pos = np.array(data[c]["pos"], dtype="i4")
        ref = np.array(data[c]["ref"], dtype="S1")
        alt = np.array(data[c]["alt"], dtype="S1")  # (n, 3)
        out[c] = (pos, ref, alt)
        print(f"{c}: {len(pos)} sites extracted")
    return out


def write_sites_zarr(sites):
    path = BUCKET / RELEASE_PATH / "snp_genotypes/all/sites"
    root = zarr.open(path, mode="w")
    for contig, (pos, ref, alt) in sites.items():
        variants = root.require_group(contig).require_group("variants")
        variants.create_dataset(name="POS", data=pos)
        variants.create_dataset(name="REF", data=ref)
        variants.create_dataset(name="ALT", data=alt)
    zarr.consolidate_metadata(path)


def write_genotypes_zarr(sample_ids):
    path = BUCKET / RELEASE_PATH / "snp_genotypes/all" / SAMPLE_SET
    root = zarr.open(path, mode="w")
    samples = np.array(sample_ids, dtype="S")
    root.create_dataset(name="samples", data=samples)

    n_samples = len(sample_ids)
    for contig in CONTIGS:
        n_sites = CONTIG_LENGTHS_ACTUAL[contig]
        gt = np.empty(shape=(n_sites, n_samples, 2), dtype="i1")
        for i, sample_id in enumerate(sample_ids):
            src_path = DOWNLOADED / f"{sample_id}.zarr.zip"
            src = zarr.open_group(src_path, mode="r")
            gt[:, i, :] = src[sample_id][contig]["calldata"]["GT"][:, 0, :]
            print(f"  {contig}: merged {sample_id} ({i + 1}/{n_samples})")
        contig_grp = root.require_group(contig)
        calldata = contig_grp.require_group("calldata")
        gt_chunks = (n_sites // 20, n_samples, None)
        calldata.create_dataset(name="GT", data=gt, chunks=gt_chunks)
    zarr.consolidate_metadata(path)


if __name__ == "__main__":
    sample_ids = get_sample_ids()
    print(f"{len(sample_ids)} samples with GT data available")

    print("Writing config...")
    write_config()

    print("Writing manifest...")
    write_manifest(sample_ids)

    print("Writing sample metadata...")
    sample_ids = write_sample_metadata(sample_ids)  # re-order to metadata row order

    print("Writing stub genome zarr...")
    write_stub_genome()

    print("Extracting site list from representative VCF (streaming)...")
    sites = extract_sites(REPRESENTATIVE_VCF, CONTIGS)
    CONTIG_LENGTHS_ACTUAL = {c: len(sites[c][0]) for c in CONTIGS}

    print("Writing sites zarr...")
    write_sites_zarr(sites)

    print("Writing combined genotypes zarr (stacking per-sample GT arrays)...")
    write_genotypes_zarr(sample_ids)

    print("Done. Fixture at:", BUCKET)
