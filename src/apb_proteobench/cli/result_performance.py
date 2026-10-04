"""Write the pMultiQC intermediate and matching ProteoBot datapoint."""

from __future__ import annotations

import hashlib
import json
import os
import re
import struct
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import numpy as np
import pandas as pd

from apb_proteobench.api import SubmissionContent

RESULT_PERFORMANCE_FILENAME = "result_performance.csv"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_HASH_VERSION = b"apb-proteobench-content-v2\0"
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
    *,
    content: SubmissionContent,
) -> ResultPerformanceFiles:
    """Stage and publish a pMultiQC CSV and its ProteoBot datapoint JSON.

    The JSON is named from a versioned hash of the selected layer and its scoring
    inputs. If either final path already exists, neither file is written.

    Args:
        frame: ProteoBench-compatible ion-level intermediate table.
        datapoint: JSON-compatible ProteoBot fields excluding writer-owned identity.
        target: Exact CSV output path, named ``result_performance.csv``.
        content: Selected abundance, feature, protein, and sample values.

    Returns:
        Both published paths.

    Raises:
        ValueError: The inputs are incompatible or either output already exists.
        OSError: The output directory or files cannot be written.
    """
    _validate_frame(frame, target)
    if target.exists():
        raise ValueError(f"result-performance output already exists: {target}")
    document = _complete_datapoint(datapoint, content)
    json_target = _proteobot_target(document, target)
    existing = next((path for path in (target, json_target) if path.exists()), None)
    if existing is not None:
        raise ValueError(f"result-performance output already exists: {existing}")

    target.parent.mkdir(parents=True, exist_ok=True)
    csv_temporary: Path | None = None
    json_temporary: Path | None = None
    published: list[Path] = []
    try:
        csv_temporary = _stage_csv(frame, target)
        json_temporary = _stage_json(document, json_target)
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


def _complete_datapoint(
    datapoint: Mapping[str, object], content: SubmissionContent
) -> dict[str, object]:
    if {"id", "intermediate_hash"} & datapoint.keys():
        raise ValueError("ProteoBot identity fields belong to the result-performance writer")
    software_name = datapoint.get("software_name")
    if not isinstance(software_name, str):
        raise ValueError("ProteoBot software_name must be a string")
    comments = datapoint.get("submission_comments", "")
    if not isinstance(comments, str):
        raise ValueError("ProteoBot submission_comments must be a string")
    digest = content_hash(content)
    return {
        **datapoint,
        "id": f"{software_name.replace(' ', '_')}_{digest[:12]}",
        "intermediate_hash": digest,
        "submission_comments": (
            f"{comments}\n\nDataset URL: https://proteobench.cubimed.rub.de/datasets/{digest}/"
        ),
    }


def content_hash(content: SubmissionContent) -> str:
    """Hash one canonical, order-independent scientific submission input."""
    matrix = content.matrix
    rows, columns = matrix.shape
    if rows != len(content.raw_files) or rows != len(content.conditions):
        raise ValueError("submission content has misaligned observations")
    if columns != len(content.feature_ids) or columns != len(content.reported_proteins):
        raise ValueError("submission content has misaligned features")
    if len(set(content.raw_files)) != rows or not content.feature_ids.is_unique:
        raise ValueError("submission content requires unique raw files and feature IDs")

    digest = hashlib.sha256(_HASH_VERSION)

    def put_text(value: object) -> None:
        if value is None or value is pd.NA or (isinstance(value, float) and np.isnan(value)):
            digest.update(b"\x00")
            return
        encoded = str(value).encode("utf-8")
        digest.update(b"\x01")
        digest.update(struct.pack("<Q", len(encoded)))
        digest.update(encoded)

    settings = json.dumps(content.settings, sort_keys=True, separators=(",", ":"), allow_nan=False)
    put_text(settings)
    observation_order = sorted(range(rows), key=content.raw_files.__getitem__)
    feature_order = sorted(range(columns), key=content.feature_ids.__getitem__)
    digest.update(struct.pack("<QQ", rows, columns))
    for row in observation_order:
        put_text(content.raw_files[row])
        put_text(content.conditions[row])
    for column in feature_order:
        put_text(content.feature_ids[column])
        put_text(content.reported_proteins.iloc[column])

    values = np.array(matrix[np.ix_(observation_order, feature_order)], dtype="<f8", order="C")
    values.view("<u8")[np.isnan(values)] = 0x7FF8000000000000
    values[values == 0] = 0.0
    digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


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
    if not isinstance(intermediate_hash, str) or _SHA256.fullmatch(intermediate_hash) is None:
        raise ValueError(
            "ProteoBot intermediate_hash must be a lowercase SHA-256 hexadecimal value"
        )
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
