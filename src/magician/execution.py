"""External commands with bounded scratch, version records and disk checks."""
from pathlib import Path
import os
import shlex
import shutil
import signal
import subprocess
import time
from .config import resolve
from .io import tree_bytes, read_json, write_json
_last_check = (0.0, None, 0)


def cache_bytes(root):
    """Share the immutable environment footprint across jobs on slow filesystems."""
    snapshot = root / "usage" / "cache.json"
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    signature = []
    for path in (root, root / "conda", root / "packages"):
        signature.append(path.stat().st_mtime_ns if path.exists() else 0)
    try:
        previous = read_json(snapshot)
        if previous["signature"] == signature and 0 <= time.time() - previous["measured_at"] < 300:
            return previous["bytes"]
    except (OSError, ValueError, KeyError):
        pass
    size = tree_bytes(root) + 4096  # Include the small snapshot itself conservatively.
    write_json(snapshot, dict(signature=signature, bytes=size, measured_at=time.time()))
    return size


def check_storage(cfg, required=0, throttle=False):
    global _last_check
    roots = [resolve(cfg[k]) for k in ("output_dir", "cache_dir")]
    key = (tuple(roots), cfg["storage"]["budget_gb"])
    if throttle and not required and _last_check[1] == key and time.monotonic() - _last_check[0] < 5:
        return _last_check[2]
    limit = cfg["storage"]["budget_gb"] * 10**9
    # Installation precedes jobs. Refresh the shared cache footprint on directory
    # changes or every five minutes instead of rescanning all packages per job.
    used = tree_bytes(roots[0]) + cache_bytes(roots[1])
    if used + required > limit:
        raise RuntimeError(f"Managed storage budget exceeded: {used + required:,} > {limit:,} bytes")
    for root in roots:
        existing = root
        while not existing.exists():
            existing = existing.parent
        free = shutil.disk_usage(existing).free
        headroom = cfg["storage"]["headroom_gb"] * 10**9
        if free < required + headroom:
            raise RuntimeError(f"Insufficient free space on {existing}: {free:,} bytes")
    _last_check = (time.monotonic(), key, used)
    return used


def terminate(proc):
    if proc.poll() is None:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
            proc.wait(timeout=5)
        except (ProcessLookupError, subprocess.TimeoutExpired):
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()


def run(commands, cfg, stdout=None):
    """Run a command or pipeline without a shell; kill all stages on any failure."""
    if isinstance(commands[0], str):
        commands = [commands]
    check_storage(cfg, throttle=True)
    procs = []
    previous = None
    print("COMMAND:", " | ".join(shlex.join([str(v) for v in c]) for c in commands), flush=True)
    try:
        for i, command in enumerate(commands):
            proc = subprocess.Popen([str(v) for v in command], stdin=previous,
                                    stdout=stdout if i == len(commands)-1 else subprocess.PIPE,
                                    start_new_session=True)
            if previous is not None:
                previous.close()
            previous = proc.stdout
            procs.append(proc)
        last_check = time.monotonic()
        while any(p.poll() is None for p in procs):
            failed = [p.returncode for p in procs if p.poll() not in (None, 0)]
            if failed:
                raise subprocess.CalledProcessError(failed[0], commands)
            if time.monotonic() - last_check >= 5:
                check_storage(cfg)
                last_check = time.monotonic()
            active = [proc for proc in procs if proc.poll() is None]
            if active:
                try:
                    active[-1].wait(timeout=0.1)
                except subprocess.TimeoutExpired:
                    pass
        for proc in procs:
            if proc.returncode:
                raise subprocess.CalledProcessError(proc.returncode, commands)
    finally:
        for proc in procs:
            terminate(proc)
    check_storage(cfg, throttle=True)


def version(command):
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, timeout=15)
        return result.stdout.strip()[:2000]
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"unavailable: {exc}"
