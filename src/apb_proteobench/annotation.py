"""Bind one ProteoBench module to the configured APB quantification level."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

import polars as pl
from apb2.annotation_extension import (
    AnnotationError,
    AnnotationFileOrigin,
    AnnotationMatches,
    AnnotationResult,
    RequireCompleteAnnotation,
    annotation_matching_for,
    make_annotation_table,
    match_annotation,
    record_annotation_provenance,
)
from apb2.result_facade import ParsedLevels

from apb_proteobench.configuration.load import LoadedModule, load_module


@dataclass(frozen=True, slots=True)
class ProteoBenchAnnotation:
    """A complete ProteoBench sample design bound to one dataset level."""

    parsed: ParsedLevels
    matches: AnnotationMatches
    module: LoadedModule

    def annotate(self) -> AnnotationResult:
        """Attach the validated samples and store normalized module evidence."""
        level_name = self.module.settings.general.level
        selected = ParsedLevels(
            levels={level_name: self.parsed.levels[level_name]},
            uns=deepcopy(self.parsed.uns),
            metadata=deepcopy(self.parsed.metadata),
        )
        applied = RequireCompleteAnnotation().apply(selected, self.matches)
        recorded = record_annotation_provenance(
            applied,
            "proteobench",
            AnnotationFileOrigin(Path(self.module.source.name)),
            metadata={"schema_version": "1", **self.module.metadata()},
        )
        levels = dict(self.parsed.levels)
        levels[level_name] = recorded.parsed.levels[level_name]
        return AnnotationResult(
            parsed=ParsedLevels(
                levels=levels,
                uns=deepcopy(self.parsed.uns),
                metadata=recorded.parsed.metadata,
            ),
            reports=recorded.reports,
        )


@dataclass(frozen=True, slots=True)
class ProteoBenchAnnotationParser:
    """A validated module source ready to bind to one APB result."""

    module: LoadedModule

    @classmethod
    def from_path(cls, path: Path, /) -> ProteoBenchAnnotationParser:
        """Decode and validate the complete module document once."""
        return cls(load_module(path))

    def parse(self, parsed: ParsedLevels, /) -> ProteoBenchAnnotation:
        """Validate complete one-to-one dataset coverage before constructing annotation."""
        level_name = self.module.settings.general.level
        if level_name not in parsed.levels:
            raise AnnotationError(
                f"ProteoBench module selects unavailable level {level_name!r}; "
                f"available={list(parsed.levels)}"
            )
        level = parsed.levels[level_name]
        table = make_annotation_table(
            _sample_frame(self.module),
            ("__match_raw_file",),
            ("__match_raw_file_alias",),
            AnnotationFileOrigin(Path(self.module.source.name)),
        )
        selected = ParsedLevels(
            levels={level_name: level},
            uns=parsed.uns,
            metadata=parsed.metadata,
        )
        matches = match_annotation(
            table,
            selected,
            {level_name: annotation_matching_for(level)},
        )
        RequireCompleteAnnotation().validate(matches)
        coverage = matches.levels[level_name].coverage
        if coverage.annotation_only_count:
            raise AnnotationError(
                "ProteoBench module contains samples absent from quantification; "
                f"count={coverage.annotation_only_count}, "
                f"examples={list(coverage.annotation_only_examples)}"
            )
        return ProteoBenchAnnotation(parsed=parsed, matches=matches, module=self.module)


def _sample_frame(module: LoadedModule) -> pl.DataFrame:
    samples = module.settings.samples
    values: dict[str, list[str | None]] = {
        "__match_raw_file": [sample.raw_file for sample in samples],
        "raw_file": [sample.raw_file for sample in samples],
        "sample_name": [sample.sample_name for sample in samples],
        "condition": [sample.condition for sample in samples],
    }
    values["__match_raw_file_alias"] = [sample.raw_file_alias for sample in samples]
    return pl.DataFrame(values)
