"""Write the pMultiQC intermediate and matching ProteoBot datapoint."""

from __future__ import annotations

import json
import os
import re
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import pandas as pd

RESULT_PERFORMANCE_FILENAME = "result_performance.csv"
_SHA1 = re.compile(r"^[0-9a-f]{40}$")
_REQUIRED_COLUMNS = frozenset(
    {
        "precursor ion",
        "species",
        "log_Intensity_mean_A",
        "log_Intensity_mean_B",
        "log_Intensity_std_A",
        "log_Intensity_std_B",
        "CV_A",
        "CV_B",
        "log2_A_vs_B",
        "epsilon",
    }
)
_REQUIRED_PROTEOBOT_FIELDS = frozenset(
    {
        "id",
        "software_name",
        "software_version",
        "search_engine",
        "search_engine_version",
        "ident_fdr_psm",
        "ident_fdr_peptide",
        "ident_fdr_protein",
        "enable_match_between_runs",
        "precursor_mass_tolerance",
        "fragment_mass_tolerance",
        "enzyme",
        "allowed_miscleavages",
        "min_peptide_length",
        "max_peptide_length",
        "is_temporary",
        "intermediate_hash",
        "results",
        "median_abs_epsilon_global",
        "mean_abs_epsilon_global",
        "median_abs_epsilon_eq_species",
        "mean_abs_epsilon_eq_species",
        "median_abs_epsilon_precision_global",
        "mean_abs_epsilon_precision_global",
        "median_abs_epsilon_precision_eq_species",
        "mean_abs_epsilon_precision_eq_species",
        "nr_feature",
        "comments",
        "proteobench_version",
        "old_new",
        "semi_enzymatic",
        "fixed_mods",
        "variable_mods",
        "max_mods",
        "min_precursor_charge",
        "max_precursor_charge",
        "quantification_method",
        "protein_inference",
        "abundance_normalization_ions",
        "postprocessing_performed",
        "postprocessing_description",
        "submission_comments",
    }
)


@dataclass(frozen=True, slots=True)
class ResultPerformanceFiles:
    """The two files published by one result-performance export."""

    csv: Path
    proteobot_json: Path


class _Syncable(Protocol):
    """Minimal writable-file capability needed for durable staging."""

    def flush(self) -> None:
        """Flush buffered content."""
        ...

    def fileno(self) -> int:
        """Return the operating-system file descriptor."""
        ...


def write_result_performance(frame: pd.DataFrame, target: Path, /) -> Path:
    """Write one ion-level ProteoBench intermediate without overwriting.

    Args:
        frame: ProteoBench-compatible ion-level intermediate table.
        target: Exact output path, named ``result_performance.csv``.

    Returns:
        The written output path.

    Raises:
        ValueError: The target name or input columns are incompatible, or the target exists.
        OSError: The output directory or file cannot be written.
    """
    _validate_frame(frame, target)
    if target.exists():
        raise ValueError(f"result-performance output already exists: {target}")

    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = _stage_csv(frame, target)
    try:
        os.link(temporary, target)
    except FileExistsError as error:
        raise ValueError(f"result-performance output already exists: {target}") from error
    finally:
        temporary.unlink(missing_ok=True)
    return target


def write_result_performance_bundle(
    frame: pd.DataFrame,
    datapoint: Mapping[str, object],
    target: Path,
    /,
) -> ResultPerformanceFiles:
    """Stage and publish a pMultiQC CSV and its ProteoBot datapoint JSON.

    The JSON is named from the ProteoBench intermediate hash, matching the
    ``Proteobench/Results_quant_ion_DDA`` repository convention. If either final
    path already exists, neither file is written.

    Args:
        frame: ProteoBench-compatible ion-level intermediate table.
        datapoint: Complete JSON-compatible ProteoBot datapoint.
        target: Exact CSV output path, named ``result_performance.csv``.

    Returns:
        Both published paths.

    Raises:
        ValueError: The inputs are incompatible or either output already exists.
        OSError: The output directory or files cannot be written.
    """
    _validate_frame(frame, target)
    json_target = _proteobot_target(datapoint, target)
    existing = next((path for path in (target, json_target) if path.exists()), None)
    if existing is not None:
        raise ValueError(f"result-performance output already exists: {existing}")

    target.parent.mkdir(parents=True, exist_ok=True)
    csv_temporary: Path | None = None
    json_temporary: Path | None = None
    published: list[Path] = []
    try:
        csv_temporary = _stage_csv(frame, target)
        json_temporary = _stage_json(datapoint, json_target)
        for temporary, final in (
            (json_temporary, json_target),
            (csv_temporary, target),
        ):
            os.link(temporary, final)
            published.append(final)
    except FileExistsError as error:
        for path in published:
            path.unlink(missing_ok=True)
        raise ValueError(f"result-performance output already exists: {error.filename}") from error
    except Exception:
        for path in published:
            path.unlink(missing_ok=True)
        raise
    finally:
        if csv_temporary is not None:
            csv_temporary.unlink(missing_ok=True)
        if json_temporary is not None:
            json_temporary.unlink(missing_ok=True)
    return ResultPerformanceFiles(csv=target, proteobot_json=json_target)


def _validate_frame(frame: pd.DataFrame, target: Path) -> None:
    if target.name != RESULT_PERFORMANCE_FILENAME:
        raise ValueError(
            f"result-performance output must be named {RESULT_PERFORMANCE_FILENAME!r}; got {target}"
        )
    missing = sorted(_REQUIRED_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(
            "result-performance export requires an ion-level ProteoBench intermediate; "
            f"missing columns: {missing}"
        )


def _proteobot_target(datapoint: Mapping[str, object], target: Path) -> Path:
    missing = sorted(_REQUIRED_PROTEOBOT_FIELDS - set(datapoint))
    if missing:
        raise ValueError(f"ProteoBot datapoint is missing fields: {missing}")
    intermediate_hash = datapoint.get("intermediate_hash")
    if not isinstance(intermediate_hash, str) or _SHA1.fullmatch(intermediate_hash) is None:
        raise ValueError("ProteoBot intermediate_hash must be a lowercase SHA-1 hexadecimal value")
    return target.with_name(f"{intermediate_hash}.json")


def _stage_csv(frame: pd.DataFrame, target: Path) -> Path:
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            prefix=f".{target.name}.",
            suffix=".tmp",
            dir=target.parent,
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            frame.to_csv(handle, index=False, lineterminator="\n")
            _flush(handle)
        return temporary
    except Exception:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise


def _stage_json(datapoint: Mapping[str, object], target: Path) -> Path:
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            prefix=f".{target.name}.",
            suffix=".tmp",
            dir=target.parent,
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(dict(datapoint), handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
            _flush(handle)
        return temporary
    except Exception:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise


def _flush(handle: _Syncable) -> None:
    handle.flush()
    os.fsync(handle.fileno())
