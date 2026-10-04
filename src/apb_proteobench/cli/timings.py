"""Write optional, tool-owned timing files for an integrated benchmark run."""

from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Mapping, Sequence
from math import isfinite
from pathlib import Path


def write_tool_timings(
    target: Path,
    /,
    *,
    tool: str,
    operation: str,
    phases: Sequence[tuple[str, float]],
    levels: Sequence[Mapping[str, str | float]] = (),
) -> Path:
    """Atomically publish a version-1 timing file without replacing an existing file."""
    if target.exists():
        raise ValueError(f"timing output already exists: {target}")
    if not tool or not operation:
        raise ValueError("timing tool and operation must be nonempty")
    names = [name for name, _seconds in phases]
    if len(set(names)) != len(names) or any(not name for name in names):
        raise ValueError("timing phase names must be nonempty and unique")
    if any(not isfinite(seconds) or seconds < 0 for _name, seconds in phases):
        raise ValueError("timing phase durations must be finite and nonnegative")
    document: dict[str, object] = {
        "format": "apb-tool-timings",
        "format_version": 1,
        "tool": tool,
        "operation": operation,
        "phases": [{"name": name, "seconds": seconds} for name, seconds in phases],
    }
    if levels:
        document["levels"] = [dict(level) for level in levels]
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=target.parent,
        prefix=f".{target.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(document, handle, ensure_ascii=False, allow_nan=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if target.exists():
            raise ValueError(f"timing output already exists: {target}")
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)
    return target
