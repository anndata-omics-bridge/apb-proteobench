"""Translate between APB storage-neutral results and ProteoBench calculation values."""

from __future__ import annotations

import json
from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass, replace
from typing import Literal, Protocol, cast
from urllib.parse import quote

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
_DIAGNOSTICS_PREFIX = f"{_STORAGE_KEY}:"


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
class ResolvedLayerSelection:
    """Persistable resolution of one requested layer-selection policy."""

    mode: Literal["primary", "named_abundance", "all_abundance"]
    layer_names: tuple[str, ...]
    requested_layer: str | None = None
    fallback: Literal["primary_missing_abundance_roles"] | None = None

    def as_json(self) -> dict[str, JsonValue]:
        """Return common provenance; selected quantities are the persisted scoring keys."""
        document: dict[str, JsonValue] = {"selection_mode": self.mode}
        if self.fallback is not None:
            document["selection_fallback"] = self.fallback
        return document


class LayerSelection(Protocol):
    """Select quantitative layers from one configured APB level."""

    def resolve(self, level: ParsedLevel, /) -> ResolvedLayerSelection:
        """Validate and return the ordered selected layer names."""
        ...


@dataclass(frozen=True, slots=True)
class PrimaryLayer:
    """Select the APB primary layer represented by AnnData ``X``."""

    def resolve(self, level: ParsedLevel, /) -> ResolvedLayerSelection:
        """Return the validated primary layer."""
        _require_layer(level, level.primary_layer_name)
        return ResolvedLayerSelection(mode="primary", layer_names=(level.primary_layer_name,))


@dataclass(frozen=True, slots=True)
class NamedAbundanceLayer:
    """Select one named layer carrying the semantic abundance role."""

    name: str

    def resolve(self, level: ParsedLevel, /) -> ResolvedLayerSelection:
        """Return the named layer after validating its semantic role."""
        _require_layer(level, self.name)
        abundance = _abundance_layers(level)
        if abundance is None or self.name not in abundance:
            raise ValueError(f"layer {self.name!r} does not carry the stored abundance role")
        return ResolvedLayerSelection(
            mode="named_abundance",
            layer_names=(self.name,),
            requested_layer=self.name,
        )


@dataclass(frozen=True, slots=True)
class AllAbundanceLayers:
    """Select every stored abundance layer, falling back to primary for old results."""

    def resolve(self, level: ParsedLevel, /) -> ResolvedLayerSelection:
        """Return abundance layers in their authored order."""
        abundance = _abundance_layers(level)
        if not abundance:
            _require_layer(level, level.primary_layer_name)
            return ResolvedLayerSelection(
                mode="all_abundance",
                layer_names=(level.primary_layer_name,),
                fallback="primary_missing_abundance_roles",
            )
        return ResolvedLayerSelection(mode="all_abundance", layer_names=abundance)


PRIMARY_LAYER = PrimaryLayer()
ALL_ABUNDANCE_LAYERS = AllAbundanceLayers()


@dataclass(frozen=True, slots=True)
class ExtractedProteoBenchLayer:
    """One layer's calculation input paired with its persisted APB role resolution."""

    level_name: ParsedLevelName
    layer_name: str
    calculation: QuantitativeLevelInput
    roles: ResolvedRoles


@dataclass(frozen=True, slots=True)
class ScoredLayerResult:
    """One selected layer's roles, persisted slot, and completed analysis."""

    level_name: ParsedLevelName
    layer_name: str
    diagnostics_slot: str
    roles: ResolvedRoles
    analysis: ProteoBenchResult


def embedded_configuration(parsed: ParsedLevels, /) -> ModuleSettings:
    """Read the normalized module configuration contributed during annotation."""
    proteobench = _object(parsed.metadata.get(_STORAGE_KEY), "ProteoBench metadata")
    provenance = _object(proteobench.get("provenance"), "ProteoBench provenance")
    details = _object(provenance.get("annotation"), "ProteoBench annotation provenance")
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


def resolve_layer_selection(
    parsed: ParsedLevels,
    configuration: ModuleSettings,
    selection: LayerSelection = ALL_ABUNDANCE_LAYERS,
    /,
) -> ResolvedLayerSelection:
    """Resolve one selection policy against the configured APB level."""
    level_name = configuration.general.level
    try:
        level = parsed.levels[level_name]
    except KeyError as error:
        raise ValueError(
            f"embedded ProteoBench configuration selects unavailable level {level_name!r}"
        ) from error
    _require_available_storage(level)
    resolved = selection.resolve(level)
    if not resolved.layer_names:
        raise ValueError("layer selection resolved to no layers")
    return resolved


def extract_layer(
    parsed: ParsedLevels,
    configuration: ModuleSettings,
    layer_name: str,
    /,
) -> ExtractedProteoBenchLayer:
    """Extract one layer from the configured level as typed scientific input."""
    level_name = configuration.general.level
    try:
        level = parsed.levels[level_name]
    except KeyError as error:
        raise ValueError(
            f"embedded ProteoBench configuration selects unavailable level {level_name!r}"
        ) from error
    _require_available_storage(level)
    _require_layer(level, layer_name)
    feature = _single_feature_key(level)
    proteins = _protein_role(level)
    matrix = (
        quantitative_layer_values(level, layer_name).to_numpy().astype(np.float64, copy=False).T
    )
    return ExtractedProteoBenchLayer(
        level_name=level_name,
        layer_name=layer_name,
        calculation=QuantitativeLevelInput(
            observations=level.obs.frame.to_pandas(),
            matrix=matrix,
            feature_ids=pd.Index(level.var.frame.get_column(feature).cast(pl.String).to_list()),
            reported_proteins=level.var.frame.get_column(proteins).to_pandas(),
            level=level_name,
        ),
        roles=ResolvedRoles(
            proteins=proteins,
            feature=feature,
            intensity=layer_name,
        ),
    )


def persist_results(
    parsed: ParsedLevels,
    selection: ResolvedLayerSelection,
    results: Mapping[str, ScoredLayerResult],
    /,
) -> ParsedLevels:
    """Return a copy with every selected layer's diagnostics and scores attached."""
    if not results:
        raise ValueError("cannot persist an empty ProteoBench layer result")
    first = next(iter(results.values()))
    level = parsed.levels[first.level_name]
    _require_available_storage(level)
    if tuple(results) != selection.layer_names:
        raise ValueError("ProteoBench analyses do not match the selected layers")
    varm = dict(level.varm)
    layer_records: dict[str, JsonValue] = {}
    for selected in results.values():
        result = selected.analysis
        diagnostics = pl.from_pandas(
            result.diagnostics.varm.reset_index(drop=True),
            include_index=False,
        )
        if diagnostics.height != level.var.frame.height:
            raise ValueError(
                f"ProteoBench diagnostics for layer {selected.layer_name!r} "
                "do not align to the APB variable axis"
            )
        varm[selected.diagnostics_slot] = diagnostics
        layer_records[_metadata_layer_key(selected.layer_name)] = _layer_result_record(
            selected,
        )
    level_metadata = deepcopy(level.metadata)
    tool = _tool_section(level_metadata)
    tool["scoring"] = layer_records
    root_metadata = deepcopy(parsed.metadata)
    root_tool = _tool_section(root_metadata)
    provenance = root_tool.setdefault("provenance", {})
    if not isinstance(provenance, dict):
        raise ValueError("ProteoBench provenance must be an object")
    if "scoring" in provenance:
        raise ValueError("ProteoBench scoring provenance already exists; refusing to overwrite")
    provenance["scoring"] = _result_record(selection, results)
    levels = dict(parsed.levels)
    levels[first.level_name] = replace(level, varm=varm, metadata=level_metadata)
    return replace(
        parsed,
        levels=levels,
        uns=deepcopy(parsed.uns),
        metadata=root_metadata,
        annotation_tables=deepcopy(parsed.annotation_tables),
        feature_relations=deepcopy(parsed.feature_relations),
    )


def diagnostics_slot(layer_name: str, /) -> str:
    """Return the logical APB ``varm`` slot for one layer's diagnostics."""
    return f"{_DIAGNOSTICS_PREFIX}{layer_name}"


def _result_record(
    selection: ResolvedLayerSelection,
    results: Mapping[str, ScoredLayerResult],
) -> dict[str, JsonValue]:
    first_result = next(iter(results.values())).analysis
    return {
        "schema_version": "3",
        "compatibility_version": PROTEOBENCH_COMPATIBILITY_VERSION,
        "source_revision": PROTEOBENCH_SOURCE_REVISION,
        **selection.as_json(),
        "diagnostic_method": cast(dict[str, JsonValue], first_result.diagnostic_method),
        "scoring_method": cast(dict[str, JsonValue], first_result.scoring_method),
    }


def _layer_result_record(
    result: ScoredLayerResult,
) -> dict[str, JsonValue]:
    score_document = json.loads(result.analysis.scores.model_dump_json())
    mapping_document = result.analysis.diagnostics.protein_mapping.model_dump(mode="json")
    if not isinstance(score_document, dict):
        raise TypeError("ProteoBench result serialization did not produce JSON objects")
    return {
        "layer_name": result.layer_name,
        "diagnostics": f"varm:{result.diagnostics_slot}",
        "column_roles": result.roles.as_json(),
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
    if any(name == _STORAGE_KEY or name.startswith(_DIAGNOSTICS_PREFIX) for name in level.varm):
        raise ValueError("ProteoBench diagnostics already exist; refusing to overwrite them")
    tool = level.metadata.get(_STORAGE_KEY, {})
    if not isinstance(tool, dict):
        raise ValueError("ProteoBench metadata must be an object")
    if "scoring" in tool:
        raise ValueError("ProteoBench scoring already exists; refusing to overwrite scores")


def _tool_section(metadata: dict[str, JsonValue]) -> dict[str, JsonValue]:
    tool = metadata.setdefault(_STORAGE_KEY, {})
    if not isinstance(tool, dict):
        raise ValueError("ProteoBench metadata must be an object")
    return tool


def _metadata_layer_key(layer_name: str, /) -> str:
    """Encode one logical layer name as an HDF5-safe metadata mapping key."""
    return quote(layer_name, safe="")


def _require_layer(level: ParsedLevel, name: str, /) -> None:
    if name not in level.layers:
        raise ValueError(f"configured APB level has no layer {name!r}")


def _abundance_layers(level: ParsedLevel, /) -> tuple[str, ...] | None:
    raw_roles = level.uns.get("layer_roles")
    if raw_roles is None:
        return None
    if not isinstance(raw_roles, dict):
        raise ValueError("APB layer_roles is not an object")
    raw_abundance = raw_roles.get("abundance")
    if raw_abundance is None:
        return None
    if not isinstance(raw_abundance, list) or not all(
        isinstance(item, str) for item in raw_abundance
    ):
        raise ValueError("APB abundance layer role is not a string list")
    abundance = tuple(cast(list[str], raw_abundance))
    if len(abundance) != len(set(abundance)):
        raise ValueError("APB abundance layer role contains duplicate layer names")
    missing = [name for name in abundance if name not in level.layers]
    if missing:
        raise ValueError(f"APB abundance layer role references missing layers: {missing}")
    return abundance


def _object(value: JsonValue | None, role: str) -> dict[str, JsonValue]:
    if not isinstance(value, dict):
        raise ValueError(f"{role} is absent or is not an object; run annotation first")
    return value


def _string_list(value: JsonValue | None, role: str) -> list[str]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError(f"{role} is absent or is not a string list")
    return cast(list[str], value)
