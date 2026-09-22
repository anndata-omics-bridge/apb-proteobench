"""Versioned timing artifacts remain separate from APB2 scientific results."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from apb_proteobench.io.tool_timings import write_tool_timings


def test_tool_timing_file_is_atomic_and_refuses_overwrite(tmp_path: Path) -> None:
    target = tmp_path / "timings" / "apb-fasta.verify-peptides.timings.json"

    written = write_tool_timings(
        target,
        tool="apb-fasta",
        operation="verify-peptides",
        phases=(("load_database", 0.25), ("verify_peptides", 1.5)),
    )

    assert written == target
    assert json.loads(target.read_text(encoding="utf-8"))["phases"] == [
        {"name": "load_database", "seconds": 0.25},
        {"name": "verify_peptides", "seconds": 1.5},
    ]
    assert not list(target.parent.glob("*.tmp"))
    with pytest.raises(ValueError, match="already exists"):
        write_tool_timings(
            target,
            tool="apb-fasta",
            operation="verify-peptides",
            phases=(("verify_peptides", 2.0),),
        )
    assert json.loads(target.read_text(encoding="utf-8"))["phases"][1]["seconds"] == 1.5
