"""Draft: package amplicon genotype calls into an AIM-style Zarr store.

Scope
-----
This does NOT do read alignment or variant calling. It assumes those steps
have already happened upstream (e.g. via an amplicon-calling pipeline such
as iVar/GATK/a panel vendor's own caller) and that, per sample, you have:

  * genotype calls at each targeted site in the panel
  * per-site call quality (GQ) and allele depth (AD)
  * per-amplicon read coverage (for dropout / evenness QC)

`simulate_upstream_calls()` below fabricates data in that shape purely so
this script is runnable end-to-end without real inputs. In a real version,
that function would be replaced by a parser for whatever the upstream
caller actually emits (VCF, a caller-specific TSV, etc.) - everything
downstream of it is agnostic to where the arrays came from.

Design choices, and how they map onto conventions in this repo
----------------------------------------------------------------
- Genotype calls follow the same dimension/coordinate naming as
  `malariagen_data.anoph.aim_data`: dims `variants`/`samples`/`ploidy`/
  `alleles`, coords `sample_id`/`variant_contig`/`variant_position`, data
  variable `call_genotype`. This is the closest existing analogue to a
  small, curated marker panel (see `AnophelesAimData.aim_calls`).
- Per-amplicon coverage doesn't fit the per-variant calldata shape used for
  WGS SNPs (`GQ`/`AD`/`MQ` in `AnophelesSnpData`), because amplicon dropout
  is a property of the *amplicon*, not any single site within it. So it
  gets its own dimension (`amplicons`) and its own dense (amplicons,
  samples) arrays: `amplicon_depth` and the derived `amplicon_dropout`.
  A `variant_amplicon` coordinate on the `variants` dimension links each
  site back to the amplicon that produced it, the same way `variant_contig`
  links a site back to its chromosome.
- Kit/library-prep/platform metadata is NOT stored inside the Zarr arrays.
  It is per-sample provenance, so it belongs in a CSV table keyed by
  `sample_id`, exactly like `general_metadata()` /
  `samples.meta.csv` for WGS sample sets. Mixing provenance into the
  genotype array would make it awkward to update independently of calls.
- The panel definition (which sites belong to which amplicon, and their
  ref/alt alleles) is written as its own small Zarr store, independent of
  any one sample set - mirroring `aim_defs_{analysis}/{aims}.zarr`, which
  is separate from `aim_calls_{analysis}/{sample_set}/{aims}.zarr`.

What's deliberately left as a stub for a production version
-------------------------------------------------------------
- Storage: this writes to a local directory via a plain path. The real
  repo opens stores via `_init_zarr_store(fs, path)` over an `fsspec`
  filesystem (GCS in production; see `malariagen_data/anoph/base.py` and
  `malariagen_data/util.py`), so swap `zarr.DirectoryStore` /
  `xr.Dataset.to_zarr(path, ...)` for a `zarr.storage.FSStore` built from
  the resource's `self._fs`.
- Config: a real integration would add panel/analysis identifiers to the
  release's config JSON (e.g. `AMPLICON_PANEL_IDS`,
  `DEFAULT_AMPLICON_ANALYSIS`), the same way `AIM_IDS` /
  `DEFAULT_AIM_ANALYSIS` are declared in `v3-config.json`.
- Reader API: `read_amplicon_calls()` at the bottom sketches the read side
  (open Zarr, join sample metadata, optional query filter) in the same
  shape as `AnophelesAimData.aim_calls()`, so that adding this as a mixin
  method later is a small step, not a rewrite.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# Mirrors the dimension names in malariagen_data/util.py
# (DIM_VARIANT, DIM_SAMPLE, DIM_PLOIDY, DIM_ALLELE).
DIM_VARIANT = "variants"
DIM_SAMPLE = "samples"
DIM_PLOIDY = "ploidy"
DIM_ALLELE = "alleles"
DIM_AMPLICON = "amplicons"

PLOIDY = 2
ALLELE_CODES = np.array(["A", "C", "G", "T"])


def make_synthetic_panel(
    *,
    n_amplicons: int,
    sites_per_amplicon: int,
    contig: str,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Fabricate a panel definition: which sites belong to which amplicon.

    In production this table would come from the panel/kit design (primer
    BED file, manifest, etc.), not be randomly generated.
    """
    rows = []
    pos = 1000
    for amplicon_ix in range(n_amplicons):
        amplicon_id = f"AMP{amplicon_ix:03d}"
        pos += rng.integers(200, 500)  # amplicons don't overlap
        for _ in range(sites_per_amplicon):
            pos += rng.integers(5, 50)
            ref, alt = rng.choice(ALLELE_CODES, size=2, replace=False)
            rows.append(
                dict(
                    amplicon_id=amplicon_id,
                    contig=contig,
                    variant_position=pos,
                    ref_allele=ref,
                    alt_allele=alt,
                )
            )
    return pd.DataFrame(rows)


def make_synthetic_sample_sheet(
    *,
    n_samples: int,
    panel_id: str,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Fabricate per-sample provenance metadata.

    This is the amplicon-data analogue of `general_metadata()` - kept as a
    flat table keyed by `sample_id`, separate from the genotype arrays.
    """
    kits = ["AmpliconKitA_v2", "AmpliconKitB_v1"]
    platforms = ["Illumina MiSeq", "Oxford Nanopore MinION"]
    return pd.DataFrame(
        {
            "sample_id": [f"SAMP{i:04d}" for i in range(n_samples)],
            "panel_id": panel_id,
            "library_prep_kit": rng.choice(kits, size=n_samples),
            "sequencing_platform": rng.choice(platforms, size=n_samples),
        }
    )


def simulate_upstream_calls(
    *,
    panel_df: pd.DataFrame,
    sample_sheet: pd.DataFrame,
    dropout_rate: float,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """Stand in for parsing an upstream caller's output.

    Returns arrays already shaped as this script's Zarr layout expects:
    call_genotype (variants, samples, ploidy), call_gq (variants, samples),
    call_ad (variants, samples, alleles=2 i.e. ref/alt depth), and
    amplicon_depth (amplicons, samples).

    Replace this function with a real parser (e.g. reading per-sample VCFs
    and a coverage table) when wiring up real data; everything downstream
    only depends on this shape, not on how it was produced.
    """
    n_variants = len(panel_df)
    n_samples = len(sample_sheet)
    n_amplicons = panel_df["amplicon_id"].nunique()

    # Simulate amplicon-level coverage first, since dropout should make the
    # genotype calls at that amplicon's sites come out missing.
    amplicon_depth = rng.poisson(lam=200, size=(n_amplicons, n_samples)).astype(
        "float64"
    )
    dropout_mask_amplicon = rng.random(size=(n_amplicons, n_samples)) < dropout_rate
    amplicon_depth[dropout_mask_amplicon] = rng.uniform(
        0, 5, size=int(dropout_mask_amplicon.sum())
    )

    # Map each variant to its amplicon's dropout mask.
    amplicon_ids = panel_df["amplicon_id"].astype("category")
    variant_amplicon_ix = amplicon_ids.cat.codes.to_numpy()
    dropout_mask_variant = dropout_mask_amplicon[variant_amplicon_ix, :]

    call_genotype = rng.integers(
        0, 2, size=(n_variants, n_samples, PLOIDY), dtype="int8"
    )
    call_genotype[dropout_mask_variant] = -1  # missing call where amplicon dropped out

    call_gq = rng.integers(10, 60, size=(n_variants, n_samples)).astype("int16")
    call_gq[dropout_mask_variant] = 0

    ad_ref = rng.integers(5, 100, size=(n_variants, n_samples))
    ad_alt = rng.integers(0, 50, size=(n_variants, n_samples))
    call_ad = np.stack([ad_ref, ad_alt], axis=-1).astype("int32")
    call_ad[dropout_mask_variant] = 0

    return dict(
        call_genotype=call_genotype,
        call_gq=call_gq,
        call_ad=call_ad,
        amplicon_depth=amplicon_depth,
        variant_amplicon_ix=variant_amplicon_ix,
    )


def build_amplicon_calls_dataset(
    *,
    panel_df: pd.DataFrame,
    sample_sheet: pd.DataFrame,
    calls: dict[str, np.ndarray],
    dropout_depth_threshold: float,
) -> xr.Dataset:
    """Assemble the calls + coverage into one xarray Dataset.

    Mirrors the coordinate/data-variable naming used by
    `AnophelesAimData.aim_calls()`, extended with an `amplicons` dimension
    for coverage/dropout, which has no equivalent in the AIM or SNP data
    because those don't have a notion of "amplicon".
    """
    amplicon_ids = panel_df["amplicon_id"].astype("category").cat.categories.to_numpy()
    amplicon_depth = calls["amplicon_depth"]
    amplicon_dropout = amplicon_depth < dropout_depth_threshold

    ds = xr.Dataset(
        data_vars={
            "variant_allele": (
                (DIM_VARIANT, DIM_ALLELE),
                panel_df[["ref_allele", "alt_allele"]].to_numpy(),
            ),
            "call_genotype": (
                (DIM_VARIANT, DIM_SAMPLE, DIM_PLOIDY),
                calls["call_genotype"],
            ),
            "call_GQ": ((DIM_VARIANT, DIM_SAMPLE), calls["call_gq"]),
            "call_AD": ((DIM_VARIANT, DIM_SAMPLE, DIM_ALLELE), calls["call_ad"]),
            "amplicon_depth": ((DIM_AMPLICON, DIM_SAMPLE), amplicon_depth),
            "amplicon_dropout": ((DIM_AMPLICON, DIM_SAMPLE), amplicon_dropout),
        },
        coords={
            "sample_id": (DIM_SAMPLE, sample_sheet["sample_id"].to_numpy()),
            "variant_position": (
                DIM_VARIANT,
                panel_df["variant_position"].to_numpy(),
            ),
            "variant_amplicon": (
                DIM_VARIANT,
                amplicon_ids[calls["variant_amplicon_ix"]],
            ),
            "amplicon_id": (DIM_AMPLICON, amplicon_ids),
        },
        attrs={
            "contigs": sorted(panel_df["contig"].unique().tolist()),
            "panel_id": sample_sheet["panel_id"].iloc[0],
            "dropout_depth_threshold": dropout_depth_threshold,
        },
    )
    return ds


def write_amplicon_store(ds: xr.Dataset, path: Path) -> None:
    """Write the calls dataset to a local Zarr store with consolidated metadata.

    In production, `path` would be replaced by an `fsspec`-backed store
    (see `_init_zarr_store` in `malariagen_data/util.py`) pointing at
    `{release_path}/amplicon_calls_{analysis}/{sample_set}/{panel_id}.zarr`.
    """
    ds.to_zarr(path, mode="w", consolidated=True)


def write_panel_definition(panel_df: pd.DataFrame, path: Path) -> None:
    """Write the panel definition as its own small Zarr store.

    Kept separate from any one sample set's calls, mirroring
    `aim_defs_{analysis}/{aims}.zarr` vs.
    `aim_calls_{analysis}/{sample_set}/{aims}.zarr`.
    """
    amplicon_ids = panel_df["amplicon_id"].astype("category")
    ds = xr.Dataset(
        data_vars={
            "variant_allele": (
                (DIM_VARIANT, DIM_ALLELE),
                panel_df[["ref_allele", "alt_allele"]].to_numpy(),
            ),
        },
        coords={
            "variant_position": (DIM_VARIANT, panel_df["variant_position"].to_numpy()),
            "variant_amplicon": (DIM_VARIANT, amplicon_ids.to_numpy()),
        },
        attrs={"contigs": sorted(panel_df["contig"].unique().tolist())},
    )
    ds.to_zarr(path, mode="w", consolidated=True)


def write_sample_metadata(sample_sheet: pd.DataFrame, path: Path) -> None:
    """Write per-sample provenance as a flat CSV, analogous to samples.meta.csv."""
    sample_sheet.to_csv(path, index=False)


def read_amplicon_calls(
    *,
    calls_path: Path,
    metadata_path: Path,
    sample_query: str | None = None,
) -> xr.Dataset:
    """Sketch of the read side, shaped like `AnophelesAimData.aim_calls()`.

    Opens the calls store directly with `xr.open_zarr` (no manual Dask
    wiring needed, since it was written in xarray's native Zarr encoding),
    then optionally filters samples using a pandas query against the
    metadata table, the same pattern `aim_calls()` uses via
    `_filter_sample_dataset()`.
    """
    ds = xr.open_zarr(calls_path, consolidated=True)
    if sample_query is not None:
        df_samples = pd.read_csv(metadata_path)
        loc_samples = df_samples.eval(sample_query).to_numpy()
        ds = ds.isel({DIM_SAMPLE: loc_samples})
    return ds


def main() -> None:
    rng = np.random.default_rng(42)
    out_dir = Path(__file__).parent / "amplicon_demo_output"
    out_dir.mkdir(exist_ok=True)

    panel_df = make_synthetic_panel(
        n_amplicons=20, sites_per_amplicon=3, contig="2L", rng=rng
    )
    sample_sheet = make_synthetic_sample_sheet(
        n_samples=50, panel_id="demo_panel_v1", rng=rng
    )
    calls = simulate_upstream_calls(
        panel_df=panel_df, sample_sheet=sample_sheet, dropout_rate=0.08, rng=rng
    )
    ds = build_amplicon_calls_dataset(
        panel_df=panel_df,
        sample_sheet=sample_sheet,
        calls=calls,
        dropout_depth_threshold=10.0,
    )

    calls_path = out_dir / "amplicon_calls_demo" / "demo_panel_v1.zarr"
    panel_path = out_dir / "amplicon_defs_demo" / "demo_panel_v1.zarr"
    metadata_path = out_dir / "samples.meta.csv"

    write_amplicon_store(ds, calls_path)
    write_panel_definition(panel_df, panel_path)
    write_sample_metadata(sample_sheet, metadata_path)

    print(f"Wrote calls to: {calls_path}")
    print(f"Wrote panel definition to: {panel_path}")
    print(f"Wrote sample metadata to: {metadata_path}")
    print()
    print(ds)
    print()

    dropout_rate_per_sample = ds["amplicon_dropout"].mean(dim=DIM_AMPLICON)
    print("Dropout rate per sample (first 5):")
    print(
        pd.DataFrame(
            {
                "sample_id": ds["sample_id"].to_numpy()[:5],
                "dropout_rate": dropout_rate_per_sample.to_numpy()[:5],
            }
        )
    )

    print()
    print("Round-trip read check via read_amplicon_calls():")
    ds_read = read_amplicon_calls(
        calls_path=calls_path,
        metadata_path=metadata_path,
        sample_query="sequencing_platform == 'Illumina MiSeq'",
    )
    print(f"Filtered to {ds_read.sizes[DIM_SAMPLE]} of {ds.sizes[DIM_SAMPLE]} samples")


if __name__ == "__main__":
    main()
