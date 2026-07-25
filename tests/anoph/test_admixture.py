import numpy as np
import pytest
from pytest_cases import parametrize_with_cases

from malariagen_data import af1 as _af1
from malariagen_data import ag3 as _ag3
from malariagen_data import adir1 as _adir1
from malariagen_data import as1 as _as1

from malariagen_data.anoph.admixture import Admixture
from malariagen_data.anoph import ld_params

import os
import bed_reader

from numpy.testing import assert_array_equal


@pytest.fixture
def ag3_sim_api(ag3_sim_fixture):
    return Admixture(
        url=ag3_sim_fixture.url,
        public_url=ag3_sim_fixture.url,
        config_path=_ag3.CONFIG_PATH,
        major_version_number=_ag3.MAJOR_VERSION_NUMBER,
        major_version_path=_ag3.MAJOR_VERSION_PATH,
        pre=True,
        aim_metadata_dtype={
            "aim_species_fraction_arab": "float64",
            "aim_species_fraction_colu": "float64",
            "aim_species_fraction_colu_no2l": "float64",
            "aim_species_gambcolu_arabiensis": object,
            "aim_species_gambiae_coluzzii": object,
            "aim_species": object,
        },
        gff_gene_type="gene",
        gff_gene_name_attribute="Name",
        gff_default_attributes=("ID", "Parent", "Name", "description"),
        default_site_mask="gamb_colu_arab",
        results_cache=ag3_sim_fixture.results_cache_path.as_posix(),
        taxon_colors=_ag3.TAXON_COLORS,
        virtual_contigs=_ag3.VIRTUAL_CONTIGS,
    )


@pytest.fixture
def af1_sim_api(af1_sim_fixture):
    return Admixture(
        url=af1_sim_fixture.url,
        public_url=af1_sim_fixture.url,
        config_path=_af1.CONFIG_PATH,
        major_version_number=_af1.MAJOR_VERSION_NUMBER,
        major_version_path=_af1.MAJOR_VERSION_PATH,
        pre=False,
        gff_gene_type="protein_coding_gene",
        gff_gene_name_attribute="Note",
        gff_default_attributes=("ID", "Parent", "Note", "description"),
        default_site_mask="funestus",
        results_cache=af1_sim_fixture.results_cache_path.as_posix(),
        taxon_colors=_af1.TAXON_COLORS,
    )


@pytest.fixture
def adir1_sim_api(adir1_sim_fixture):
    return Admixture(
        url=adir1_sim_fixture.url,
        public_url=adir1_sim_fixture.url,
        config_path=_adir1.CONFIG_PATH,
        major_version_number=_adir1.MAJOR_VERSION_NUMBER,
        major_version_path=_adir1.MAJOR_VERSION_PATH,
        pre=False,
        gff_gene_type="protein_coding_gene",
        gff_gene_name_attribute="Note",
        gff_default_attributes=("ID", "Parent", "Note", "description"),
        default_site_mask="dirus",
        results_cache=adir1_sim_fixture.results_cache_path.as_posix(),
        taxon_colors=_adir1.TAXON_COLORS,
    )


@pytest.fixture
def as1_sim_api(as1_sim_fixture):
    return Admixture(
        url=as1_sim_fixture.url,
        public_url=as1_sim_fixture.url,
        config_path=_as1.CONFIG_PATH,
        major_version_number=_as1.MAJOR_VERSION_NUMBER,
        major_version_path=_as1.MAJOR_VERSION_PATH,
        pre=False,
        gff_gene_type="protein_coding_gene",
        gff_gene_name_attribute="Note",
        gff_default_attributes=("ID", "Parent", "Note", "description"),
        default_site_mask="stephensi",
        results_cache=as1_sim_fixture.results_cache_path.as_posix(),
        taxon_colors=_as1.TAXON_COLORS,
    )


# N.B., here we use pytest_cases to parametrize tests. Each
# function whose name begins with "case_" defines a set of
# inputs to the test functions. See the documentation for
# pytest_cases for more information, e.g.:
#
# https://smarie.github.io/python-pytest-cases/#basic-usage
#
# We use this approach here because we want to use fixtures
# as test parameters, which is otherwise hard to do with
# pytest alone.


def case_ag3_sim(ag3_sim_fixture, ag3_sim_api):
    return ag3_sim_fixture, ag3_sim_api


def case_af1_sim(af1_sim_fixture, af1_sim_api):
    return af1_sim_fixture, af1_sim_api


def case_as1_sim(as1_sim_fixture, as1_sim_api):
    return as1_sim_fixture, as1_sim_api


@parametrize_with_cases("fixture,api", cases=".")
def test_biallelic_snps_to_admixture(fixture, api: Admixture, tmp_path):
    # Parameters for selecting input data and filtering, before LD pruning.
    all_sample_sets = api.sample_sets()["sample_set"].to_list()

    data_params = dict(
        region=str(np.random.choice(api.contigs)),
        sample_sets=np.random.choice(all_sample_sets, size=2, replace=False).tolist(),
        site_mask=np.random.choice(list(api.site_mask_ids) + [None]),
        min_minor_ac=1,
        max_missing_an=1,
        thin_offset=1,
        random_seed=int(np.random.randint(1, 2001)),
    )

    # Load a ds containing the randomly generated samples and regions to get
    # the number of available biallelic snps to subset from (before LD pruning).
    ds = api.biallelic_snp_calls(
        **data_params,
    )

    n_snps_available = ds.sizes["variants"]
    n_snps = int(np.random.randint(1, n_snps_available + 1))

    # LD pruning parameters, kept at their defaults to keep the test robust -
    # biallelic_snps_to_admixture applies LD pruning on top of data_params,
    # unlike biallelic_snps_to_plink.
    ld_params_dict = dict(
        ld_window_size=ld_params.ld_window_size_default,
        ld_window_step=ld_params.ld_window_step_default,
        ld_threshold=ld_params.ld_threshold_default,
    )

    # Define admixture params.
    admixture_params_dict = dict(
        output_dir=str(tmp_path), n_snps=n_snps, **ld_params_dict, **data_params
    )

    # Make the admixture (plink-format) files.
    api.biallelic_snps_to_admixture(**admixture_params_dict)

    # Test to see if bed, bim, fam output files exist.
    file_path = f"{str(tmp_path)}/{data_params['region']}.{n_snps}.{data_params['min_minor_ac']}.{data_params['max_missing_an']}.{data_params['thin_offset']}"

    assert os.path.exists(f"{file_path}.bed")
    assert os.path.exists(f"{file_path}.bim")
    assert os.path.exists(f"{file_path}.fam")

    # Read bed, bim, and fam files (bed_reader searches for the .bim and .fam files matching the prefix of the .bed file).
    bed = bed_reader.open_bed(f"{file_path}.bed")

    # Load a ds containing the same LD-pruned data exported to ADMIXTURE to test against.
    ds_test = api.biallelic_snp_calls_ld_pruned(
        n_snps=n_snps,
        **ld_params_dict,
        **data_params,
    )

    # Test to make sure that the rows and columns (no. variants and no. samples) of the .bed file match.
    assert bed.shape[1] == ds_test.variant_position.shape[0]
    assert bed.shape[0] == ds_test.samples.shape[0]

    # Test to see if sample_id is exported correctly (stored in the .fam file).
    assert_array_equal(bed.iid, ds_test.sample_id.values)

    # Test to see if variant position is exported to the .bim correctly.
    assert set(bed.bp_position) == set(ds_test.variant_position.values)

    # Test to make sure chromosome ID is exported to the .bim file correctly (coerce to str to match types).
    assert set(bed.chromosome) == set(ds_test.variant_contig.values.astype(str))

    # Test to make sure that the major and minor allele are exported to the .bim file as expected (coerce to str to match types).
    assert set(bed.allele_1) == set(ds_test.variant_allele.values[:, 0].astype(str))
    assert set(bed.allele_2) == set(ds_test.variant_allele.values[:, 1].astype(str))
