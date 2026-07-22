from pathlib import Path
from malariagen_data.anoph.ld import AnophelesLdAnalysis

BUCKET = Path("/root/malariagen-data-python/local_fixture/vo_adir_release")


class MinimalAnophelesData(AnophelesLdAnalysis):
    pass


ds_res = MinimalAnophelesData(
    url=BUCKET.as_uri(),
    public_url=BUCKET.as_uri(),
    config_path="config.json",
    pre=False,
    major_version_number=1,
    major_version_path="v1",
)

print("releases:", ds_res.releases)
print("sample sets:", ds_res.sample_sets())
print("contigs:", ds_res.contigs)

df_meta = ds_res.general_metadata()
print("n samples in metadata:", len(df_meta))
print(df_meta[["sample_id", "country", "location"]].head())

# N Snps based on references mentioning 506 ancestry informative markers
gt_dask = ds_res.biallelic_snp_calls_ld_pruned(
    region="KB672490:1-1,000,000", n_snps=500, site_mask=None
)
