"""Parameters for ADMIXTURE converter functions."""

from typing import Literal, Mapping

from typing_extensions import Annotated, TypeAlias

overwrite: TypeAlias = Annotated[
    bool,
    """
    A boolean indicating whether a previously written file with the same name ought
    to be overwritten. Default is False.
    """,
]

output_dir: TypeAlias = Annotated[
    str,
    """
    A string indicating the desired output file location.
    """,
]

input_file_dir: TypeAlias = Annotated[
    str,
    """
    A string indicating the desired input file location for ADMIXTURE.
    """,
]

out: TypeAlias = Annotated[
    str,
    """
    A string specifying the output file path prefix. The ADMIXTURE input files
    will be written as ``{output_dir}/{out}.bed``, ``{output_dir}/{out}.bim``,
    and ``{output_dir}/{out}.fam``. If not provided, a default prefix is
    generated from the SNP selection parameters (region, n_snps,
    min_minor_ac, max_missing_an, thin_offset).
    """,
]

# --- K sweep (admixture-wrapper specific; ADMIXTURE itself is run once per K) ---

kmin: TypeAlias = Annotated[
    int,
    "Minimum number of ancestral populations (K) to fit, inclusive.",
]

kmax: TypeAlias = Annotated[
    int,
    "Maximum number of ancestral populations (K) to fit, inclusive.",
]

reps: TypeAlias = Annotated[
    int,
    """
    Number of replicate ADMIXTURE runs to perform for each value of K, each
    with an independently-chosen random seed unless `seed` is set.
    """,
]

reps_default: reps = 1

# --- General options ---

threads: TypeAlias = Annotated[
    int,
    "Number of threads to use for computation (ADMIXTURE's ``-j`` option).",
]

threads_default: threads = 1

seed: TypeAlias = Annotated[
    int,
    """
    Random seed for initialization (ADMIXTURE's ``--seed`` option). If not
    provided, a random seed is chosen for each replicate.
    """,
]

# --- Algorithm options ---

method: TypeAlias = Annotated[
    Literal["em", "block"],
    """
    Algorithm used to compute ADMIXTURE estimates (ADMIXTURE's ``-m``/
    ``--method`` option). ADMIXTURE's own default is 'block'.
    """,
]

method_default: method = "block"

acceleration: TypeAlias = Annotated[
    str,
    """
    Acceleration scheme used to speed up convergence (ADMIXTURE's ``-a``/
    ``--acceleration`` option). One of 'none', 'sqs<X>' (SqS3 acceleration
    with X secant conditions, e.g. 'sqs3'), or 'qn<X>' (quasi-Newton
    acceleration with X secant conditions, e.g. 'qn3'). If not provided,
    ADMIXTURE's own default acceleration scheme is used.
    """,
]

# --- Convergence criteria ---

major_convergence: TypeAlias = Annotated[
    float,
    """
    Major convergence criterion used for point estimation (ADMIXTURE's
    ``-C`` option), expressed as a log-likelihood delta.
    """,
]

minor_convergence: TypeAlias = Annotated[
    float,
    """
    Minor convergence criterion used for bootstrap and cross-validation
    re-estimates (ADMIXTURE's ``-c`` option), expressed as a log-likelihood delta.
    """,
]

# --- Cross-validation and bootstrapping ---

cv: TypeAlias = Annotated[
    int,
    """
    Number of folds to use in the cross-validation procedure (ADMIXTURE's
    ``--cv`` option), used to help choose the best value of K.
    """,
]

cv_default: cv = 5

bootstrap: TypeAlias = Annotated[
    int,
    """
    Number of bootstrap replicates to perform in order to compute standard
    errors (ADMIXTURE's ``-B`` option). If not provided, no bootstrapping
    is performed.
    """,
]

# --- Bundle of all ADMIXTURE run arguments, for callers that want to build
# and pass around one settings object instead of many separate parameters. ---

argd: TypeAlias = Annotated[
    Mapping[str, object],
    """
    A dict of ADMIXTURE run arguments, with keys matching the parameter names
    above (any subset of: kmin, kmax, reps, threads, seed, method,
    acceleration, major_convergence, minor_convergence, cv, bootstrap).
    Used by `write_log` to record the settings a run was executed with.
    """,
]

prefix: TypeAlias = Annotated[
    str,
    "Prefix used when naming summarized cross-validation output files.",
]
