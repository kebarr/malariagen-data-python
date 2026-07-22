from pathlib import Path
import allel
from malariagen_data.anoph.snp_data import AnophelesSnpData

BUCKET = Path("/root/malariagen-data-python/local_fixture/vo_adir_release")


class MinimalAnophelesData(AnophelesSnpData):
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

gt_dask = ds_res.snp_genotypes(region="KB672490:1-1,000,000", field="GT")
gt = gt_dask.compute()
print("genotype array shape (variants, samples, ploidy):", gt.shape)


gt_allel = allel.GenotypeArray(gt)
ac = gt_allel.count_alleles()
print("allele counts shape:", ac.shape)
print("n biallelic segregating sites in this 1Mb window:", ac.count_segregating())
