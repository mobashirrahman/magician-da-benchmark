"""Configuration validation shared by the launcher and workflow."""
from copy import deepcopy
from pathlib import Path
import re
import yaml

ROOT = Path(__file__).resolve().parents[2]
METRICS = {"counts", "tpm", "rpkm", "relative"}
ENV_DIR = ROOT / "workflow/envs"
ADAPTER_DIR = ROOT / "workflow/scripts/da"
INPUT_MODES = {"simulation", "matrix"}
ZERO_HANDLING = {"pseudocount", "dirichlet"}
SCALES = {"measurement", "count"}


def merge(base, update):
    result = deepcopy(base)
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = deepcopy(value)
    return result


def resolve(path):
    p = Path(path).expanduser()
    return p if p.is_absolute() else ROOT / p


def load(path=None, overrides=None):
    with (ROOT / "config/config.yaml").open() as handle:
        cfg = yaml.safe_load(handle)
    if path:
        with resolve(path).open() as handle:
            cfg = merge(cfg, yaml.safe_load(handle) or {})
    cfg = merge(cfg, overrides or {})
    validate(cfg)
    return cfg


def validate(cfg):
    with (ROOT / "config/config.yaml").open() as handle:
        template = yaml.safe_load(handle)
    def unknown_keys(value, expected, prefix=""):
        if isinstance(value, dict) and isinstance(expected, dict):
            for key in value:
                if key not in expected:
                    raise ValueError(f"Unknown configuration key: {prefix}{key}")
                unknown_keys(value[key], expected[key], f"{prefix}{key}.")
    unknown_keys(cfg, template)
    for key in ("output_dir", "cache_dir"):
        p = resolve(cfg[key])
        if p.is_symlink() or p.resolve() == ROOT.resolve() or p.resolve() == Path("/"):
            raise ValueError(f"{key} must be a dedicated directory, not a symlink or repository root")
    out, cache = (resolve(cfg[k]).resolve() for k in ("output_dir", "cache_dir"))
    if out == cache or out in cache.parents or cache in out.parents:
        raise ValueError("output_dir and cache_dir must be separate, non-nested directories")
    if cfg["input_mode"] not in INPUT_MODES:
        raise ValueError(f"input_mode must be one of {sorted(INPUT_MODES)}")
    if cfg["input_mode"] == "matrix" and not cfg["matrix_input"]["bundle"]:
        raise ValueError("Matrix input mode requires matrix_input.bundle")
    if cfg["input_mode"] == "simulation" and cfg["matrix_input"]["bundle"]:
        raise ValueError("matrix_input.bundle is only valid when input_mode is matrix")
    if cfg["input_mode"] == "simulation":
        validate_simulation(cfg)
    validate_analysis(cfg)
    validate_storage(cfg)


def validate_simulation(cfg):
    exp = cfg["experiment"]
    if not exp["scenarios"] or len(set(exp["scenarios"])) != len(exp["scenarios"]):
        raise ValueError("scenarios must be nonempty and unique")
    if set(exp["scenarios"]) - {"null", "spiked"}:
        raise ValueError("scenarios must be null or spiked")
    if not exp["seeds"] or len(set(exp["seeds"])) != len(exp["seeds"]):
        raise ValueError("seeds must be nonempty and unique")
    if any(type(s) is not int or s <= 0 or s >= 2**31 for s in exp["seeds"]):
        raise ValueError("seeds must be positive 31-bit integers")
    positive = [("experiment", "samples_per_group"), ("experiment", "read_pairs"),
                ("simulation", "read_length"), ("simulation", "insert_size"),
                ("resources", "threads"), ("assembly", "min_contig_length"),
                ("binning", "min_bin_size"), ("binning", "min_contig_length")]
    for section, key in positive:
        if type(cfg[section][key]) is not int or cfg[section][key] < 1:
            raise ValueError(f"{section}.{key} must be a positive integer")
    if exp["samples_per_group"] < 3:
        raise ValueError("At least three independent samples per group are required")
    if not 0 < exp["differential_fraction"] <= 1 or exp["effect_size"] <= 1:
        raise ValueError("differential_fraction must be in (0,1]; effect_size must exceed 1")
    for key in ("library_cv", "biological_cv", "load_cv"):
        if not 0 <= exp[key] <= 2:
            raise ValueError(f"{key} must be in [0,2]")
    if not 0 <= cfg["simulation"]["error_rate"] < 0.25:
        raise ValueError("error_rate must be in [0,0.25)")
    if cfg["simulation"]["insert_size"] < 2 * cfg["simulation"]["read_length"]:
        raise ValueError("insert_size must be at least twice read_length")
    if cfg["simulation"]["insert_sd"] < 0:
        raise ValueError("insert_sd must be nonnegative")
    if not cfg["genomes"]["manifest"]:
        if cfg["genomes"]["synthetic_count"] < 4 or cfg["genomes"]["synthetic_length"] < 10000:
            raise ValueError("Artificial fixtures need at least four genomes of at least 10000 bp")
        if type(cfg["genomes"]["synthetic_contigs"]) is not int or cfg["genomes"]["synthetic_contigs"] < 1:
            raise ValueError("synthetic_contigs must be a positive integer")
    if not 0 <= cfg["mapping"]["min_mapq"] <= 60:
        raise ValueError("min_mapq must be in [0,60]")
    for key in ("min_identity", "min_aligned_fraction", "ambiguity_margin"):
        if not 0 <= cfg["matching"][key] <= 1:
            raise ValueError(f"matching.{key} must be in [0,1]")
    k = cfg["assembly"]["k_list"]
    if not k or k != sorted(set(k)) or any(type(v) is not int or v < 15 or v % 2 == 0 for v in k):
        raise ValueError("k_list must contain increasing, unique, odd integers >=15")
    if max(k) >= cfg["simulation"]["read_length"]:
        raise ValueError("Assembly k-mers must be shorter than reads")


def zero_policies(cfg):
    """Declared zero-handling policies, loaded from their own tracked file."""
    path = resolve(cfg["analysis"]["zero_handling_policies"])
    if not path.exists():
        raise ValueError(f"Missing zero handling policy file: {path}")
    with path.open() as handle:
        policies = (yaml.safe_load(handle) or {}).get("policies")
    if not isinstance(policies, list):
        raise ValueError("Zero handling policy file must define a 'policies' list")
    for policy in policies:
        if set(policy) - {"name", "zero_handling", "pseudocount", "pseudocount_scale", "dirichlet_alpha"}:
            raise ValueError(f"Unknown key in zero handling policy {policy.get('name', '?')}")
        if policy["zero_handling"] not in ZERO_HANDLING:
            raise ValueError(f"zero_handling must be one of {sorted(ZERO_HANDLING)}")
        if policy["zero_handling"] == "pseudocount":
            if float(policy.get("pseudocount", -1)) < 0:
                raise ValueError("pseudocount must be nonnegative")
            if policy.get("pseudocount_scale", "measurement") not in SCALES:
                raise ValueError(f"pseudocount_scale must be one of {sorted(SCALES)}")
        if policy["zero_handling"] == "dirichlet" and float(policy.get("dirichlet_alpha", 0)) <= 0:
            raise ValueError("dirichlet_alpha must be positive")
    return policies


def validate_analysis(cfg):
    a = cfg["analysis"]
    policies = zero_policies(cfg)
    if not policies or len({p["name"] for p in policies}) != len(policies):
        raise ValueError("zero_handling_policies must be nonempty with unique names")
    if not a["metrics"] or set(a["metrics"]) - METRICS or "counts" not in a["metrics"]:
        raise ValueError("metrics must include counts and use counts/tpm/rpkm/relative")
    if len(set(a["metrics"])) != len(a["metrics"]):
        raise ValueError("metrics must be unique")
    if not 0 < a["alpha"] < 1 or not 0 <= a["min_prevalence"] <= 1 or a["min_total_count"] < 0:
        raise ValueError("Invalid alpha, prevalence, or count threshold")
    if type(a["mc_samples"]) is not int or a["mc_samples"] < 16:
        raise ValueError("mc_samples must be an integer >=16")
    if type(a["da_timeout_minutes"]) is not int or a["da_timeout_minutes"] < 1:
        raise ValueError("da_timeout_minutes must be a positive integer")
    if not a["enable_status"] or set(a["enable_status"]) - {"active", "planned"}:
        raise ValueError("enable_status may only contain active and planned")
    if a["truth_reference"] not in {"median_feature", "learned_set", "declared_reference_feature"}:
        raise ValueError("truth_reference must be median_feature, learned_set or declared_reference_feature")
    if a["truth_reference"] == "declared_reference_feature" and not a["truth_reference_feature"]:
        raise ValueError("truth_reference_feature is required for a declared reference feature")
    if not isinstance(a["allow_method_failures"], bool) or not isinstance(a["reserve_other"], bool):
        raise ValueError("allow_method_failures and reserve_other must be booleans")
    from .registry import load_registry, resolve_methods
    registry = load_registry(a["method_registry"], ENV_DIR, ADAPTER_DIR, cfg)
    resolve_methods(registry, a["methods"], a["enable_status"])
    if not registry.select(a["methods"], a["metrics"]):
        raise ValueError("No selected method accepts any configured metric")


def validate_storage(cfg):
    if cfg["storage"]["budget_gb"] <= cfg["storage"]["environment_reserve_gb"] + cfg["storage"]["headroom_gb"]:
        raise ValueError("Storage budget must exceed environment reserve plus headroom")
    for key, value in cfg["resources"].items():
        if type(value) is not int or value < 1:
            raise ValueError(f"resources.{key} must be a positive integer")


def registry(cfg):
    from .registry import load_registry
    return load_registry(cfg["analysis"]["method_registry"], ENV_DIR, ADAPTER_DIR, cfg)


def cases(cfg):
    return [f"{scenario}_s{seed}" for seed in cfg["experiment"]["seeds"]
            for scenario in cfg["experiment"]["scenarios"]]


def case_info(case):
    match = re.fullmatch(r"(null|spiked)_s([0-9]+)", case)
    if not match:
        raise ValueError(f"Invalid experiment identifier: {case}")
    return match[1], int(match[2])


def samples(cfg):
    return [f"{group}{i:02d}" for group in ("C", "T")
            for i in range(1, cfg["experiment"]["samples_per_group"] + 1)]


def combinations(cfg):
    """Compatible (metric, method) pairs, from the registry's declared inputs."""
    reg = registry(cfg)
    return reg.select(cfg["analysis"]["methods"], cfg["analysis"]["metrics"])