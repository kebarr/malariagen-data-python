from .ld import AnophelesLdAnalysis
from ..util import _dask_compress_dataset
from . import base_params, ld_params, pca_params, admixture_params
from typing import Optional
import os
import subprocess as sp
import shutil
from datetime import datetime
import numpy as np
from numpy import random
from pathlib import Path
from numpydoc_decorator import doc  # type: ignore
import allel
import bed_reader


class Admixture(
    AnophelesLdAnalysis,
):
    def __init__(
        self,
        **kwargs,
    ):
        # N.B., this class is designed to work cooperatively, and
        # so it's important that any remaining parameters are passed
        # to the superclass constructor.
        super().__init__(**kwargs)

    @doc(
        summary="""
            Convert Anopheles biallelic SNP data to the ADMIXTURE file format.
            Run ADMIXTURE on output bed files using code pulled from: https://github.com/dportik/admixture-wrapper/tree/master
            Description from repo:
                    admixture-wrapper - A tool for automating analyses with the program admixture. A directory of
            ped files should be specified using the -i flag. The minimum and maximum K values, the
            number of replicates per K, and the cross-validation procedure folds value are set by the user.
            The number of threads can also be specified. Outputs from each replicate are written to a unique
            directory created for each ped file. Two main output files are produced per ped file, one which
            contains the cross-validation scores for every replicate, and one which contains the average
            cross-validation score per K value. The second file can be used to plot the CV scores with the
            associated R script. A log file is also produced, which contains the analysis settings and the
            commands used to execute admixture for all K value replicates.


            DEPENDENCIES: admixture (in path).
        """,
        extended_summary="""
            This function writes biallelic SNPs to the admixture binary file format. It enables
            subsetting to specific regions (`region`), selecting specific sample sets, or lists of
            samples, randomly downsampling sites, and specifying filters based on missing data and
            minimum minor allele count (see the docs for `biallelic_snp_calls` for more information).
            The `overwrite` parameter, set to true, will enable overwrite of data with the same
            SNP selection parameter values.

            Cheat is to do this: https://github.com/sophiemoss/Anopheles_darlingi_genome_wide_analysis_Rondonia_Brazil/blob/main/Admixture.sh
        """,
        returns="""
        Base path to files containing binary admixture output files. Append .bed,
        .bim or .fam to obtain paths for the binary genotype table file, variant
        information file and sample information file respectively.
        """,
        notes="""
            This computation may take some time to run, depending on your computing
            environment. Unless the `overwrite` parameter is set to `True`, results will be returned
            from a previous computation, if available.
        """,
    )
    def biallelic_snps_to_admixture(
        self,  # those are input from  ld
        output_dir: admixture_params.input_file_dir,
        region: base_params.regions,
        n_snps: base_params.n_snps,
        overwrite: admixture_params.overwrite = False,
        thin_offset: base_params.thin_offset = 0,
        sample_sets: Optional[base_params.sample_sets] = None,
        sample_query: Optional[base_params.sample_query] = None,
        sample_query_options: Optional[base_params.sample_query_options] = None,
        sample_indices: Optional[base_params.sample_indices] = None,
        site_mask: Optional[base_params.site_mask] = base_params.DEFAULT,
        min_minor_ac: Optional[
            base_params.min_minor_ac
        ] = pca_params.min_minor_ac_default,
        max_missing_an: Optional[
            base_params.max_missing_an
        ] = pca_params.max_missing_an_default,
        random_seed: base_params.random_seed = 42,
        inline_array: base_params.inline_array = base_params.inline_array_default,
        chunks: base_params.chunks = base_params.native_chunks,
        out: Optional[admixture_params.out] = None,
        ld_window_size: ld_params.ld_window_size = ld_params.ld_window_size_default,
        ld_window_step: ld_params.ld_window_step = ld_params.ld_window_step_default,
        ld_threshold: ld_params.ld_threshold = ld_params.ld_threshold_default,
    ):
        # Check that either sample_query xor sample_indices are provided.
        base_params._validate_sample_selection_params(
            sample_query=sample_query, sample_indices=sample_indices
        )

        # Use user-provided prefix or fall back to auto-generated default
        if out is not None:
            admixture_file_path = f"{output_dir}/{out}"
        else:
            admixture_file_path = f"{output_dir}/{region}.{n_snps}.{min_minor_ac}.{max_missing_an}.{thin_offset}"

        bed_file_path = f"{admixture_file_path}.bed"

        # Check to see if file exists and if overwrite is set to false, return existing file
        if os.path.exists(bed_file_path):
            if not overwrite:
                return admixture_file_path
        base_params._validate_sample_selection_params(
            sample_query=sample_query, sample_indices=sample_indices
        )

        # Validate LD parameters.
        if ld_window_size <= 0:
            raise ValueError(f"ld_window_size must be > 0, got {ld_window_size}")
        if ld_window_step <= 0:
            raise ValueError(f"ld_window_step must be > 0, got {ld_window_step}")
        if not (0 < ld_threshold <= 1):
            raise ValueError(f"ld_threshold must be in (0, 1], got {ld_threshold}")

        gt_dask = self.biallelic_snp_calls_ld_pruned(
            region=region,
            n_snps=n_snps,
            ld_window_size=ld_window_size,
            ld_window_step=ld_window_step,
            ld_threshold=ld_threshold,
            thin_offset=thin_offset,
            sample_sets=sample_sets,
            sample_query=sample_query,
            sample_query_options=sample_query_options,
            sample_indices=sample_indices,
            site_mask=site_mask,
            min_minor_ac=min_minor_ac,
            max_missing_an=max_missing_an,
            random_seed=random_seed,
            inline_array=inline_array,
            chunks=chunks,
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
        with self._spinner("Prepare output data"):
            alleles = ds_snps_final["variant_allele"].values
            properties = {
                "iid": ds_snps_final["sample_id"].values,
                "chromosome": ds_snps_final["variant_contig"].values,
                "bp_position": ds_snps_final["variant_position"].values,
                "allele_1": alleles[:, 0],
                "allele_2": alleles[:, 1],
            }

        # local bed file
        bed_reader.to_bed(
            filepath=bed_file_path,
            val=val,
            properties=properties,
            count_A1=True,
        )

        return admixture_file_path

    def write_log(
        self,
        output_dir: admixture_params.output_dir,
        text: str,
        argd: Optional[admixture_params.argd] = None,
    ):
        log_path = Path(output_dir, "admixture_wrapper.log")
        if argd:
            with open(log_path, "a") as fh:
                fh.write(
                    "Run executed: {}\n\nadmixture_wrapper settings:\n"
                    "-i:\t\t{}\n--kmin:\t{}\n--kmax:\t{}\n--reps:\t{}\n--cv:\t{}\n"
                    "-t:\t\t{}\n--seed:\t{}\n--method:\t{}\n--acceleration:\t{}\n"
                    "-C:\t{}\n-c:\t{}\n-B:\t{}\n\n".format(
                        datetime.now(),
                        argd.get("indir"),
                        argd.get("kmin"),
                        argd.get("kmax"),
                        argd.get("reps"),
                        argd.get("cv"),
                        argd.get("threads"),
                        argd.get("seed"),
                        argd.get("method"),
                        argd.get("acceleration"),
                        argd.get("major_convergence"),
                        argd.get("minor_convergence"),
                        argd.get("bootstrap"),
                    )
                )
        else:
            with open(log_path, "a") as fh:
                fh.write("{}".format(text))

    def get_beds(self, input_dir: admixture_params.input_file_dir):
        os.chdir(input_dir)
        beds = [f for f in os.listdir(".") if f.endswith(".bed")]
        if not beds:
            raise ValueError(
                "\n\n\nERROR: No bed files (X.bed) were "
                "located in input directory.\n\n\n"
            )
        else:
            print("\n\nFound {} ped files to run:".format(len(beds)))
            for p in beds:
                print("\t{}".format(p))
            return beds

    def run_admixture(
        self,
        p: str,
        input_dir: admixture_params.input_file_dir,
        output_dir: admixture_params.output_dir,
        kmin: admixture_params.kmin,
        kmax: admixture_params.kmax,
        reps: admixture_params.reps = admixture_params.reps_default,
        cv: admixture_params.cv = admixture_params.cv_default,
        threads: admixture_params.threads = admixture_params.threads_default,
        seed: Optional[admixture_params.seed] = None,
        method: admixture_params.method = admixture_params.method_default,
        acceleration: Optional[admixture_params.acceleration] = None,
        major_convergence: Optional[admixture_params.major_convergence] = None,
        minor_convergence: Optional[admixture_params.minor_convergence] = None,
        bootstrap: Optional[admixture_params.bootstrap] = None,
    ):
        os.chdir(input_dir)

        # Record the settings this run was executed with.
        self.write_log(
            output_dir=output_dir,
            text="",
            argd={
                "indir": input_dir,
                "kmin": kmin,
                "kmax": kmax,
                "reps": reps,
                "cv": cv,
                "threads": threads,
                "seed": seed,
                "method": method,
                "acceleration": acceleration,
                "major_convergence": major_convergence,
                "minor_convergence": minor_convergence,
                "bootstrap": bootstrap,
            },
        )

        # Build the extra ADMIXTURE flags that are only added when specified,
        # since ADMIXTURE has its own internal defaults for each of these.
        extra_flags = f"--method={method}"
        if acceleration is not None:
            extra_flags += f" --acceleration={acceleration}"
        if major_convergence is not None:
            extra_flags += f" -C={major_convergence}"
        if minor_convergence is not None:
            extra_flags += f" -c={minor_convergence}"
        if bootstrap is not None:
            extra_flags += f" -B{bootstrap}"

        kreps = []
        for i in range(kmin, kmax + 1):
            kreps.append([[i, x] for x in list(range(1, reps + 1))])

        for i in kreps:
            for j in i:
                tb = datetime.now()
                print("\n\n{}".format("-" * 50))
                print("Running: K{0} replicate {1}".format(j[0], j[1]))
                print("{}\n".format("-" * 50))

                rep_seed = seed if seed is not None else random.randint(5000)
                call_str = "admixture {0} {1} -j{2} --cv={3} -s {4} {5} | tee {6}.{1}.out".format(
                    p, j[0], threads, cv, rep_seed, extra_flags, p.split(".ped")[0]
                )
                self.write_log(
                    output_dir=output_dir,
                    text="{0}: K{1} replicate {2}: {3}\n".format(
                        datetime.now(), j[0], j[1], call_str
                    ),
                )
                print("{}\n".format(call_str))
                retcode = sp.call(call_str, shell=True)
                if retcode != 0:
                    self.write_log(
                        output_dir=output_dir,
                        text=f"Process {call_str} exited with returncode {retcode}",
                    )
                outs = [f for f in os.listdir(".") if f.endswith((".P", ".Q", ".out"))]
                for o in outs:
                    pieces = o.split(".")
                    if len(pieces) == 3:
                        shutil.move(
                            o,
                            os.path.join(
                                output_dir,
                                "{}.{}.{}.{}".format(
                                    o.split(".")[0],
                                    o.split(".")[1],
                                    j[1],
                                    o.split(".")[-1],
                                ),
                            ),
                        )
                    elif len(pieces) > 3:
                        out_prefix = "_".join(pieces[:-2])
                        shutil.move(
                            o,
                            os.path.join(
                                output_dir,
                                "{}.{}.{}.{}".format(
                                    out_prefix, pieces[-2], j[1], pieces[-1]
                                ),
                            ),
                        )

                tf = datetime.now()
                print("\n{}".format("-" * 50))
                print("Finished: K{0} replicate {1}".format(j[0], j[1]))
                print("Elapsed time: {} (H:M:S)".format(tf - tb))
                self.write_log(
                    output_dir=output_dir,
                    text="{0}: K{1} replicate {2}: Finished. Elapsed time: {3}\n\n".format(
                        datetime.now(), j[0], j[1], tf - tb
                    ),
                )
                print("{}\n".format("-" * 50))

    def summarize_outputs(
        self,
        output_dir: admixture_params.output_dir,
        kmin: admixture_params.kmin,
        kmax: admixture_params.kmax,
        prefix: admixture_params.prefix,
    ):
        #  Q (the ancestry fractions), and P (the allele frequencies of the inferred ancestral populations).
        os.chdir(output_dir)
        outs = [f for f in os.listdir(".") if f.endswith(".out")]

        if outs:
            tb = datetime.now()
            print("\n\n{}".format("-" * 50))
            print("\nSummarizing output files...")

            outall = "Cross_Validation_All_Replicates.txt"
            with open(outall, "a") as fh:
                fh.write("{}\t{}\t{}\n".format("K", "Rep", "CV"))

            outavg = "Cross_Validation_Averages.txt"
            with open(outavg, "a") as fh:
                fh.write("{}\t{}\t{}\n".format("K", "CV_Avg", "CV_Stdev"))

            for i in range(kmin, kmax + 1):
                kouts = sorted([f for f in outs if int(f.split(".")[1]) == i])
                cv_vals = []
                for f in kouts:
                    with open(f, "r") as fh:
                        # for line in fh.readlines():
                        cv = [
                            float(line.strip().split()[-1])
                            for line in fh
                            if line.startswith("CV error")
                        ]
                        cv_vals.append(cv[0])
                        with open(outall, "a") as fh:
                            fh.write("{}\t{}\t{}\n".format(i, f.split(".")[2], cv[0]))

                with open(outavg, "a") as fh:
                    fh.write(
                        "{}\t{}\t{}\n".format(
                            i,
                            np.round(np.mean(cv_vals), 4),
                            np.round(np.std(cv_vals), 4),
                        )
                    )

            tf = datetime.now()
            print("\tFinished.\n\tElapsed time: {} (H:M:S)".format(tf - tb))
            print("{}\n".format("-" * 50))

            shutil.move(
                outall,
                os.path.join(output_dir, "{}.CV_All.txt".format(prefix)),
            )
            shutil.move(
                outavg,
                os.path.join(output_dir, "{}.CV_Avg.txt".format(prefix)),
            )

        else:
            raise ValueError(
                "\n\n\nERROR: No output log files found in directory: {}\n\n\n".format(
                    output_dir
                )
            )
