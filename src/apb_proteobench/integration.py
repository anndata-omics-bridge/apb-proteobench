"""Translate between APB storage-neutral results and ProteoBench calculation values."""

from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import dataclass, replace
from typing import cast

import numpy as np
import pandas as pd
import polars as pl
from apb2.result_facade import (
    JsonValue,
    ParsedLevel,
    ParsedLevelName,
    ParsedLevels,
    quantitative_layer_values,
)

from apb_proteobench.calculation.contracts import QuantitativeLevelInput
from apb_proteobench.calculation.metrics import (
    PROTEOBENCH_COMPATIBILITY_VERSION,
    PROTEOBENCH_SOURCE_REVISION,
)
from apb_proteobench.configuration.schema import ModuleSettings
from apb_proteobench.workflow import ProteoBenchResult

_STORAGE_KEY = "proteobench"


@dataclass(frozen=True, slots=True)
class ResolvedRoles:
    """Logical APB columns and layer selected for one calculation."""

    proteins: str
    feature: str
    intensity: str

    def as_json(self) -> dict[str, JsonValue]:
        """Return portable role locations."""
        return {
            "Proteins": f"var:{self.proteins}",
            "feature": f"var:{self.feature}",
            "Intensity": f"layer:{self.intensity}",
            "Sample name": "obs:sample_name",
            "Condition": "obs:condition",
        }


@dataclass(frozen=True, slots=True)
class ExtractedProteoBenchLevel:
    """Calculation input paired with its persisted APB role resolution."""

    name: ParsedLevelName
    calculation: QuantitativeLevelInput
    roles: ResolvedRoles


def embedded_configuration(parsed: ParsedLevels, /) -> ModuleSettings:
    """Read the normalized module configuration contributed during annotation."""
    annotation = _object(parsed.metadata.get("annotation"), "APB annotation section")
    proteobench = _object(annotation.get(_STORAGE_KEY), "ProteoBench annotation section")
    details = _object(proteobench.get("metadata"), "ProteoBench annotation metadata")
    configuration = _object(details.get("configuration"), "ProteoBench configuration")
    document = dict(configuration)
    sample_columns = _object(document.get("samples"), "ProteoBench sample configuration")
    raw_files = _string_list(sample_columns.get("raw_file"), "sample raw_file")
    sample_names = _string_list(sample_columns.get("sample_name"), "sample sample_name")
    conditions = _string_list(sample_columns.get("condition"), "sample condition")
    if not (len(raw_files) == len(sample_names) == len(conditions)):
        raise ValueError("embedded ProteoBench sample columns have different lengths")
    document["samples"] = [
        {"raw_file": raw_file, "sample_name": sample_name, "condition": condition}
        for raw_file, sample_name, condition in zip(
            raw_files,
            sample_names,
            conditions,
            strict=True,
        )
    ]
    return ModuleSettings.model_validate(document)


def extract_level(
    parsed: ParsedLevels,
    configuration: ModuleSettings,
    /,
) -> ExtractedProteoBenchLevel:
    """Extract exactly the configured level's typed scientific inputs."""
    name = configuration.general.level
    try:
        level = parsed.levels[name]
    except KeyError as error:
        raise ValueError(
            f"embedded ProteoBench configuration selects unavailable level {name!r}"
        ) from error
    _require_available_storage(level)
    feature = _single_feature_key(level)
    proteins = _protein_role(level)
    matrix = (
        quantitative_layer_values(level, level.primary_layer_name)
        .to_numpy()
        .astype(np.float64, copy=False)
        .T
    )
    return ExtractedProteoBenchLevel(
        name=name,
        calculation=QuantitativeLevelInput(
            observations=level.obs.frame.to_pandas(),
            matrix=matrix,
            feature_ids=pd.Index(level.var.frame.get_column(feature).cast(pl.String).to_list()),
            reported_proteins=level.var.frame.get_column(proteins).to_pandas(),
            level=name,
        ),
        roles=ResolvedRoles(
            proteins=proteins,
            feature=feature,
            intensity=level.primary_layer_name,
        ),
    )


def persist_result(
    parsed: ParsedLevels,
    extracted: ExtractedProteoBenchLevel,
    result: ProteoBenchResult,
    /,
) -> ParsedLevels:
    """Return a copy with feature diagnostics and APB metadata attached."""
    level = parsed.levels[extracted.name]
    _require_available_storage(level)
    diagnostics = pl.from_pandas(
        result.diagnostics.varm.reset_index(drop=True),
        include_index=False,
    )
    if diagnostics.height != level.var.frame.height:
        raise ValueError("ProteoBench diagnostics do not align to the APB variable axis")
    record = _result_record(extracted, result)
    varm = dict(level.varm)
    varm[_STORAGE_KEY] = diagnostics
    level_metadata = deepcopy(level.metadata)
    level_metadata[_STORAGE_KEY] = record
    levels = dict(parsed.levels)
    levels[extracted.name] = replace(level, varm=varm, metadata=level_metadata)
    return ParsedLevels(
        levels=levels,
        uns=deepcopy(parsed.uns),
        metadata=deepcopy(parsed.metadata),
    )


def _result_record(
    extracted: ExtractedProteoBenchLevel,
    result: ProteoBenchResult,
) -> dict[str, JsonValue]:
    score_document = json.loads(result.scores.model_dump_json())
    mapping_document = result.diagnostics.protein_mapping.model_dump(mode="json")
    if not isinstance(score_document, dict):
        raise TypeError("ProteoBench result serialization did not produce JSON objects")
    return {
        "schema_version": "1",
        "compatibility_version": PROTEOBENCH_COMPATIBILITY_VERSION,
        "source_revision": PROTEOBENCH_SOURCE_REVISION,
        "diagnostic_method": cast(dict[str, JsonValue], result.diagnostic_method),
        "scoring_method": cast(dict[str, JsonValue], result.scoring_method),
        "column_roles": extracted.roles.as_json(),
        "protein_mapping": cast(dict[str, JsonValue], mapping_document),
        "scores": cast(dict[str, JsonValue], score_document),
    }


def _single_feature_key(level: ParsedLevel) -> str:
    if len(level.var.key_columns) != 1:
        raise ValueError(
            "ProteoBench scoring currently requires one APB variable key; "
            f"got {list(level.var.key_columns)}"
        )
    return level.var.key_columns[0]


def _protein_role(level: ParsedLevel) -> str:
    roles = _object(level.uns.get("column_roles"), "APB column_roles")
    value = roles.get("protein_assignment")
    if not isinstance(value, str) or value not in level.var.frame.columns:
        raise ValueError("ProteoBench scoring requires a stored protein_assignment column role")
    return value


def _require_available_storage(level: ParsedLevel) -> None:
    if _STORAGE_KEY in level.varm:
        raise ValueError("varm['proteobench'] already exists; refusing to overwrite diagnostics")
    if _STORAGE_KEY in level.metadata:
        raise ValueError("uns['apb']['proteobench'] already exists; refusing to overwrite scores")


def _object(value: JsonValue | None, role: str) -> dict[str, JsonValue]:
    if not isinstance(value, dict):
        raise ValueError(f"{role} is absent or is not an object; run annotation first")
    return value


def _string_list(value: JsonValue | None, role: str) -> list[str]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError(f"{role} is absent or is not a string list")
    return cast(list[str], value)
