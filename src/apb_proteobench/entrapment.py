"""Score ProteoBench entrapment on canonical APB2 results."""

from __future__ import annotations

import json
import math
from copy import deepcopy
from dataclasses import dataclass, replace
from typing import cast

import polars as pl
from apb2.api import JsonValue, ParsedLevels
from apb_catalog.api import Catalog, attach_snapshot
from apb_fasta.api import FastaAnnotator

from apb_proteobench.annotation import PROTEOBENCH_SCHEMA_VERSION
from apb_proteobench.calculation.entrapment import (
    ENTRAPMENT_SOURCE_REVISION,
    TIE_RULE,
    EntrapmentError,
    EntrapmentScores,
    fasta_pairs,
    label_precursors,
    score_entrapment,
)
from apb_proteobench.calculation.mapping import render_proteobench_features
from apb_proteobench.configuration.entrapment import EntrapmentModuleSettings
from apb_proteobench.configuration.schema import QuantificationLevel

CATALOGUE = "proteobench_entrapment"
"""The apb-catalog set naming each precursor q-value kind a result offers."""
_STORAGE_KEY = "proteobench"
DIAGNOSTICS_SLOT = f"{_STORAGE_KEY}:entrapment"
"""The ``varm`` slot holding each precursor's label, pair and best q-value per kind."""
_PEPTIDE = "ProForma_peptide"
_PEPTIDOFORM = "ProForma_peptidoform"


@dataclass(frozen=True, slots=True)
class EntrapmentAnalysisResult:
    """A complete entrapment analysis and its scored APB2 result."""

    parsed: ParsedLevels
    configuration: EntrapmentModuleSettings
    precursors: pl.DataFrame
    """Labelled precursors: variable keys, ``peptide``, ``sequence``, best q-value per kind,
    ``label`` and ``pair``."""
    scores: dict[str, EntrapmentScores]
    """Scores per q-value kind, such as ``q_value`` or ``library_q_value``."""


class EntrapmentAnalyzer:
    """Bind one entrapment module and the target/entrapment pairs of its FASTA."""

    __slots__ = ("_configuration", "_pairs")

    def __init__(self, module: EntrapmentModuleSettings, fasta: FastaAnnotator) -> None:
        """Derive the target/entrapment pairs from the entrapment FASTA.

        Args:
            module: A packaged entrapment module.
            fasta: ProteoBench's entrapment FASTA, read by ``FastaAnnotator.read``.

        Raises:
            EntrapmentError: The frame lacks ``is_entrapment``, or two adjacent entries are not a
                target followed by its same-length entrapment.
        """
        self._configuration = module
        self._pairs = fasta_pairs(fasta.proteins)

    def analyze(self, parsed: ParsedLevels) -> EntrapmentAnalysisResult:
        """Label precursors and score every q-value kind the result offers, without I/O."""
        settings = self._configuration
        catalog = Catalog(parsed, CATALOGUE)
        precursors, kinds = _precursors(parsed, settings.level, catalog)
        labelled = label_precursors(
            precursors, self._pairs, max_missing_fraction=settings.max_missing_fraction
        )
        scores = {
            kind: score_entrapment(
                labelled.with_columns(pl.col(kind).alias("q_value")),
                intervals=settings.curve_intervals,
            )
            for kind in kinds
        }
        stored = _persist(parsed, settings, labelled, scores)
        return EntrapmentAnalysisResult(
            parsed=attach_snapshot(stored, catalog.snapshot()),
            configuration=settings,
            precursors=labelled,
            scores=scores,
        )


def _precursors(
    parsed: ParsedLevels, level_name: QuantificationLevel, catalog: Catalog
) -> tuple[pl.DataFrame, tuple[str, ...]]:
    """Return one row per precursor with its best q-value across runs for each kind."""
    kinds = catalog.layer(level_name, concept="confidence")
    if kinds is None:
        raise EntrapmentError(
            f"the {level_name} level offers no precursor q-value catalogued in {CATALOGUE!r}"
        )
    level = parsed.levels[level_name]
    keys = list(level.var.key_columns)
    rendered = render_proteobench_features(level.var.frame.get_column(_PEPTIDOFORM).to_pandas())
    frame = level.var.frame.select(
        *keys,
        pl.col(_PEPTIDE).alias("peptide"),
        pl.Series("sequence", rendered.to_list(), dtype=pl.String),
    )
    for kind in kinds:
        table = catalog.layer(level_name, concept="confidence", kind=kind)
        if table is None:
            raise EntrapmentError(f"the catalogue offers kind {kind!r} but resolves no layer")
        best = level.var.frame.select(keys).hstack(
            table.values.select(pl.min_horizontal(pl.all()).cast(pl.Float64).alias(kind))
        )
        frame = frame.join(best, on=keys, how="left", maintain_order="left")
    return frame, kinds


def _persist(
    parsed: ParsedLevels,
    settings: EntrapmentModuleSettings,
    labelled: pl.DataFrame,
    scores: dict[str, EntrapmentScores],
) -> ParsedLevels:
    """Return a copy with labels in ``varm`` and scores and provenance in metadata."""
    level = parsed.levels[settings.level]
    level_metadata = deepcopy(level.metadata)
    tool = _section(level_metadata, _STORAGE_KEY)
    tool_result = _section(tool, "result")
    if DIAGNOSTICS_SLOT in level.varm or "entrapment" in tool_result:
        raise ValueError("ProteoBench entrapment results already exist; refusing to overwrite")
    keys = list(level.var.key_columns)
    kinds = list(scores)
    diagnostics = (
        level.var.frame.select(keys)
        .join(
            labelled.select(*keys, "label", "pair", *kinds),
            on=keys,
            how="left",
            maintain_order="left",
        )
        .drop(keys)
    )
    # Through JSON text, so an undefined FDP is stored as null rather than NaN.
    tool_result["entrapment"] = {
        kind: cast(JsonValue, json.loads(score.model_dump_json())) for kind, score in scores.items()
    }
    _list(tool, "summary").extend(
        entry for kind, score in scores.items() for entry in _summary(kind, score)
    )
    _list(tool, "details").append({"slot": "varm", "name": DIAGNOSTICS_SLOT})
    root_metadata = deepcopy(parsed.metadata)
    root_tool = _section(root_metadata, _STORAGE_KEY)
    root_tool["schema_version"] = PROTEOBENCH_SCHEMA_VERSION
    provenance = _section(root_tool, "provenance")
    provenance["entrapment"] = {
        "module": settings.model_dump(mode="json"),
        "source_revision": ENTRAPMENT_SOURCE_REVISION,
        "ties": TIE_RULE,
        "kinds": list(kinds),
    }
    levels = dict(parsed.levels)
    levels[settings.level] = replace(
        level, varm={**level.varm, DIAGNOSTICS_SLOT: diagnostics}, metadata=level_metadata
    )
    return replace(parsed, levels=levels, metadata=root_metadata)


def _summary(kind: str, score: EntrapmentScores) -> list[JsonValue]:
    """One q-value kind's identified features and FDP estimates, keyed by the kind.

    An FDP needs attention unless its entrapment category is ``valid``.
    """
    entries: list[JsonValue] = [
        {
            "name": "identified_features",
            "label": "Identified precursors",
            "value": score.nr_id_features,
            "unit": "features",
            "status": "ok",
            "layer": kind,
        }
    ]
    for name, label, value, category in (
        ("combined_fdp", "Combined FDP", score.combined_FDP, score.category_combined),
        ("paired_fdp", "Paired FDP", score.paired_FDP, score.category_paired),
    ):
        entries.append(
            {
                "name": name,
                "label": label,
                "value": value if math.isfinite(value) else None,
                "unit": "fraction",
                "status": "ok" if category == "valid" else "attention",
                "layer": kind,
            }
        )
    return entries


def _list(document: dict[str, JsonValue], key: str) -> list[JsonValue]:
    value = document.setdefault(key, [])
    if not isinstance(value, list):
        raise ValueError(f"metadata {key!r} must be a list")
    return value


def _section(document: dict[str, JsonValue], key: str) -> dict[str, JsonValue]:
    section = document.setdefault(key, {})
    if not isinstance(section, dict):
        raise ValueError(f"metadata {key!r} must be an object")
    return section
