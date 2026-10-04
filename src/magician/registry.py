"""Versioned DA method registry: declared inputs, hypothesis, units and environment.

The registry is the single place a method is declared. `method_family` groups
statistical models so that several settings of one model (for example ALDEx2
gamma choices) are never reported as independent methods. Compatibility between
the configured metrics and each selected method is validated here, before any
job is constructed.
"""
from pathlib import Path
from copy import deepcopy
import json
import yaml
from .config import resolve
from .io import sha256

REQUIRED = ("id", "method_family", "variant", "adapter", "environment", "package",
            "package_version", "accepted_inputs", "transform", "zero_policy",
            "endpoint", "hypothesis", "effect_scale", "reference", "contrast",
            "adjustment_family", "native_q", "parameters", "thread_cap", "seed",
            "citation", "status")
# Per-method resource overrides; the configured default applies when absent.
OPTIONAL = ("mem_mb", "timeout_minutes", "notes", "truth_lane_policy")
INPUTS = {"counts", "tpm", "rpkm", "relative"}
TRANSFORMS = {"none", "clr"}
STATUSES = {"active", "planned", "deferred"}
CONTRAST = "Treatment / Control"


class Registry:
    """Immutable view over the registry file plus its repository location."""

    def __init__(self, payload, path, env_dir, adapter_dir, policies=None):
        self.path = Path(path)
        self.payload = payload
        self.env_dir = Path(env_dir)
        self.adapter_dir = Path(adapter_dir)
        self.digest = sha256(self.path)
        self.methods = {entry["id"]: entry for entry in payload["methods"]}
        self.endpoints = payload["endpoints"]
        self.policies = {policy["name"]: policy for policy in (policies or [])}

    def __contains__(self, method_id):
        return method_id in self.methods

    def entry(self, method_id):
        try:
            return self.methods[method_id]
        except KeyError:
            raise ValueError(f"Unknown DA method: {method_id}") from None

    def environment(self, method_id):
        """The yaml specification for a method's environment."""
        return self.env_dir / f"{self.entry(method_id)['environment']}.yaml"

    def execution_env(self, method_id, prefix_root):
        """An already frozen prefix if the bootstrap completed it, else the yaml.

        A frozen prefix carries its lock, so DA jobs prefer it and rerun exactly when
        the frozen package set changes. Before the bootstrap runs, jobs fall back to
        the yaml and report the package as explicitly unavailable.
        """
        from .io import read_json
        stem = self.entry(method_id)["environment"]
        prefix = Path(prefix_root) / f"magician_{stem}"
        lock = self.env_dir / "locks" / f"{stem}.lock.json"
        try:
            frozen = read_json(lock)
        except (OSError, ValueError):
            frozen = {}
        if (frozen.get("complete", False) and frozen.get("path") == str(prefix)
                and (prefix / "bin" / "Rscript").exists()):
            return prefix
        return self.environment(method_id)

    def adapter_path(self, method_id):
        return self.adapter_dir / f"{self.entry(method_id)['adapter']}.R"

    def settings(self, method_id):
        """Every declared value that changes a fit, for the report and logs."""
        declared = deepcopy(self.entry(method_id))
        return dict(method_id=method_id, method_family=declared["method_family"],
                    variant=declared["variant"], endpoint=declared["endpoint"],
                    effect_scale=declared["effect_scale"], reference=declared["reference"],
                    hypothesis=declared["hypothesis"], adjustment_family=declared["adjustment_family"],
                    package=declared["package"] or "base R", package_version=declared["package_version"],
                    environment=declared["environment"], thread_cap=declared["thread_cap"],
                    seed=declared["seed"], transform=declared["transform"],
                    zero_policy=declared["zero_policy"],
                    truth_lane_policy=declared.get("truth_lane_policy", declared["zero_policy"]),
                    parameters=declared["parameters"],
                    citation=declared["citation"], status=declared["status"])

    def select(self, method_ids, metrics):
        """Compatible (metric, method_id) pairs in declaration order."""
        pairs = []
        for metric in metrics:
            for method_id in method_ids:
                if metric in self.entry(method_id)["accepted_inputs"]:
                    pairs.append((metric, method_id))
        return pairs

    def lookup(self, method_id):
        """Registry fields a workflow may read without knowing the registry schema."""
        entry = self.entry(method_id)
        return dict(adapter=entry["adapter"], environment=entry["environment"],
                    endpoint=entry["endpoint"], method_family=entry["method_family"],
                    variant=entry["variant"], transform=entry["transform"],
                    zero_policy=entry["zero_policy"],
                    truth_lane_policy=entry.get("truth_lane_policy", entry["zero_policy"]),
                    thread_cap=entry["thread_cap"],
                    mem_mb=entry.get("mem_mb"), timeout_minutes=entry.get("timeout_minutes"))

    def entry_json(self, method_id):
        """The single argument the R runner is given: everything the job may rely on."""
        entry = self.entry(method_id)
        policy = self.policies.get(entry["zero_policy"])
        payload = {key: entry[key] for key in
                   ("id", "method_family", "variant", "adapter", "package", "endpoint",
                    "effect_scale", "reference", "contrast", "adjustment_family",
                    "native_q", "transform", "accepted_inputs", "parameters")}
        payload["method_id"] = entry["id"]
        payload["zero_policy"] = policy
        payload["truth_lane_policy"] = entry.get("truth_lane_policy", entry["zero_policy"])
        payload["registry_digest"] = self.digest
        return json.dumps(payload, sort_keys=True)

    def settings_table(self, method_ids, metrics):
        rows = []
        for metric, method_id in self.select(method_ids, metrics):
            row = self.settings(method_id)
            row["input_metric"] = metric
            rows.append(row)
        return rows


def load_registry(path, env_dir, adapter_dir, cfg=None):
    path = resolve(path)
    policies = []
    if cfg is not None:
        from .config import zero_policies
        policies = zero_policies(cfg)
    with path.open() as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict) or "methods" not in payload or "endpoints" not in payload:
        raise ValueError("Method registry must define 'methods' and 'endpoints'")
    if payload.get("schema_version") != 1:
        raise ValueError("Unsupported method registry schema_version")
    registry = Registry(payload, path, env_dir, adapter_dir, policies)
    validate(registry, cfg)
    return registry


def validate(registry, cfg=None):
    """Reject registry mistakes that would otherwise surface as silent run failures."""
    policies = set()
    if cfg is not None:
        policies = set(registry.policies)
    seen_ids, seen_labels = set(), set()
    for entry in registry.payload["methods"]:
        missing = [key for key in REQUIRED if key not in entry]
        if missing:
            raise ValueError(f"Registry entry {entry.get('id', '?')} is missing: {', '.join(missing)}")
        unknown = [key for key in entry if key not in REQUIRED and key not in OPTIONAL]
        if unknown:
            raise ValueError(f"Registry entry {entry['id']} has unknown keys: {', '.join(unknown)}")
        method_id = entry["id"]
        if method_id in seen_ids:
            raise ValueError(f"Duplicate registry id: {method_id}")
        seen_ids.add(method_id)
        if not entry["method_family"] or not entry["variant"]:
            raise ValueError(f"{method_id}: method_family and variant must be non-empty")
        label = (entry["method_family"], entry["variant"])
        if label in seen_labels:
            raise ValueError(f"Duplicate method_family/variant: {method_id}")
        seen_labels.add(label)
        if entry["status"] not in STATUSES:
            raise ValueError(f"{method_id}: status must be one of {sorted(STATUSES)}")
        if entry["contrast"] != CONTRAST:
            raise ValueError(f"{method_id}: only the Treatment / Control contrast is supported")
        if set(entry["accepted_inputs"]) - INPUTS or not entry["accepted_inputs"]:
            raise ValueError(f"{method_id}: accepted_inputs must be a nonempty subset of {sorted(INPUTS)}")
        if entry["endpoint"] not in registry.endpoints:
            raise ValueError(f"{method_id}: unknown endpoint {entry['endpoint']}")
        if entry["transform"] not in TRANSFORMS:
            raise ValueError(f"{method_id}: transform must be one of {sorted(TRANSFORMS)}")
        if entry["zero_policy"] != "none" and entry["zero_policy"] not in policies:
            raise ValueError(f"{method_id}: undeclared zero policy {entry['zero_policy']}")
        truth_policy = entry.get("truth_lane_policy", entry["zero_policy"])
        if truth_policy != "none" and truth_policy not in policies:
            raise ValueError(f"{method_id}: undeclared truth lane policy {truth_policy}")
        if entry["endpoint"] == "clr" and entry["transform"] == "none" and truth_policy == "none":
            # Acceptable: a package that centres raw counts itself, such as ALDEx with
            # gamma zero. The truth lane is then the CLR of the expected counts.
            pass
        if entry["transform"] == "clr":
            # A CLR endpoint must be tested on a measure that carries no feature
            # outside the tested composition; packages build their own geometry.
            if entry["endpoint"] != "clr" and entry["zero_policy"] == "none":
                raise ValueError(f"{method_id}: a CLR transform with no zero policy cannot fix its geometry")
        if type(entry["thread_cap"]) is not int or entry["thread_cap"] < 1:
            raise ValueError(f"{method_id}: thread_cap must be a positive integer")
        for key in ("mem_mb", "timeout_minutes"):
            if entry.get(key) is not None and (type(entry[key]) is not int or entry[key] < 1):
                raise ValueError(f"{method_id}: {key} must be a positive integer or null")
        if not isinstance(entry["parameters"], dict):
            raise ValueError(f"{method_id}: parameters must be a mapping")
        if not entry["citation"].strip():
            raise ValueError(f"{method_id}: a citation is required")
        if not registry.environment(method_id).exists():
            raise ValueError(f"{method_id}: missing environment {registry.environment(method_id)}")
        if not registry.adapter_path(method_id).exists():
            raise ValueError(f"{method_id}: missing adapter {registry.adapter_path(method_id)}")
    for name, endpoint in registry.endpoints.items():
        if not endpoint.get("lane") or not endpoint.get("description"):
            raise ValueError(f"Endpoint {name} needs a lane and a description")
    from .truth import LANES
    for name, endpoint in registry.endpoints.items():
        if endpoint["lane"] not in LANES:
            raise ValueError(f"Endpoint {name} refers to unknown truth lane {endpoint['lane']}")
    unused = sorted({name for entry in registry.payload["methods"]
                     for name in (entry["zero_policy"], entry.get("truth_lane_policy", entry["zero_policy"]))
                     if name != "none"} - policies)
    if unused:
        raise ValueError(f"Registry references undeclared zero policies: {', '.join(unused)}")


def resolve_methods(registry, requested, enabled):
    """Validate a requested method selection against the registry.

    Unknown ids, disabled entries and methods whose declared inputs exclude every
    configured metric are hard errors: a silently dropped method would look like a
    failed fit in the report.
    """
    if not requested or len(set(requested)) != len(requested):
        raise ValueError("methods must be nonempty and unique")
    unknown = [method_id for method_id in requested if method_id not in registry]
    if unknown:
        raise ValueError(f"Unknown DA method: {', '.join(unknown)}")
    blocked = [method_id for method_id in requested
               if registry.entry(method_id)["status"] not in enabled]
    if blocked:
        detail = ", ".join(f"{m} ({registry.entry(m)['status']})" for m in blocked)
        raise ValueError(f"Method status does not allow scheduling: {detail}")
    return list(requested)


def families(registry, method_ids):
    """Map each selected method to its family, for reporting shared models."""
    return {method_id: registry.entry(method_id)["method_family"] for method_id in method_ids}


def distinct_families(registry, method_ids):
    return sorted({registry.entry(m)["method_family"] for m in method_ids})