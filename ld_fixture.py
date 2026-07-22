from pathlib import Path
from malariagen_data.anoph.ld import AnophelesLdAnalysis
from malariagen_data.util import _dask_compress_dataset
import bed_reader
import numpy as np
import allel

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

## then think i want it similar to to_plink for admixture, then just call to complete
# or this: https://github.com/jmcastelo/mixtum/blob/main/gui/core.py

# both return xr.dataset so should be OK to steal from to_plink

gt_asc = gt_dask["call_genotype"].data
gn_ref = allel.GenotypeDaskArray(gt_asc).to_n_ref(fill=-127)
gn_ref = gn_ref.compute()

loc_var = np.any(gn_ref != gn_ref[:, 0, np.newaxis], axis=1)
ds_snps_final = _dask_compress_dataset(gt_dask, loc_var, dim="variants")
gn_ref_final = gn_ref[loc_var]
val = gn_ref_final.T
# with self._spinner("Prepare output data"):
alleles = ds_snps_final["variant_allele"].values
properties = {
    "iid": ds_snps_final["sample_id"].values,
    "chromosome": ds_snps_final["variant_contig"].values,
    "bp_position": ds_snps_final["variant_position"].values,
    "allele_1": alleles[:, 0],
    "allele_2": alleles[:, 1],
}

# local bed file
bed_file_path = "/root/malariagen-data-python/local_fixture/test.bed"
bed_reader.to_bed(
    filepath=bed_file_path,
    val=val,
    properties=properties,
    count_A1=True,
)

#  ../admixture_linux-1.4.0/admixture_linux-1.4.0/admixture local_fixture/test.bed 2
# now runs and produces bed, bim and fam
# but none of the options work- they do I was just being tired, should set threads i.e. jN
# consider projection approach

## to consider
# the number of markers needed to resolve populations in this kind of analysis is inversely proportional to the genetic distance (FST ) betweeen the populations.
# As a rule of thumb, we have found that 10,000 markers suffice to perform GWAS correction
# for continentally separated populations (for example, African, Asian, and European pop-
# ulations FST > .05) while more like 100,000 markers are necessary when the populations
# are within a continent (Europe, for instance, FST < 0.01).
