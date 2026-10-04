"""Score ProteoBench entrapment on canonical APB2 results."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, replace

import polars as pl
from apb2.api import JsonValue, ParsedLevels
from apb_catalog.api import Catalog, attach_snapshot

from apb_proteobench.calculation.entrapment import (
    ENTRAPMENT_SOURCE_REVISION,
    TIE_RULE,
    EntrapmentError,
    EntrapmentScores,
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
    """Bind one entrapment module and the pair table derived from its FASTA."""

    __slots__ = ("_configuration", "_pairs")

    def __init__(self, module: EntrapmentModuleSettings, /, *, pairs: pl.DataFrame) -> None:
        """Create an analyzer; ``pairs`` comes from :func:`fasta_pairs`."""
        self._configuration = module
        self._pairs = pairs

    def analyze(self, parsed: ParsedLevels, /) -> EntrapmentAnalysisResult:
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
    if DIAGNOSTICS_SLOT in level.varm or "entrapment" in tool:
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
    tool["entrapment"] = {kind: score.model_dump(mode="json") for kind, score in scores.items()}
    root_metadata = deepcopy(parsed.metadata)
    provenance = _section(_section(root_metadata, _STORAGE_KEY), "provenance")
    provenance["entrapment"] = {
        "module": settings.model_dump(mode="json"),
        "source_revision": ENTRAPMENT_SOURCE_REVISION,
        "ties": TIE_RULE,
        "kinds": list(kinds),
        "diagnostics": f"varm:{DIAGNOSTICS_SLOT}",
    }
    levels = dict(parsed.levels)
    levels[settings.level] = replace(
        level, varm={**level.varm, DIAGNOSTICS_SLOT: diagnostics}, metadata=level_metadata
    )
    return replace(parsed, levels=levels, metadata=root_metadata)


def _section(document: dict[str, JsonValue], key: str) -> dict[str, JsonValue]:
    section = document.setdefault(key, {})
    if not isinstance(section, dict):
        raise ValueError(f"metadata {key!r} must be an object")
    return section
