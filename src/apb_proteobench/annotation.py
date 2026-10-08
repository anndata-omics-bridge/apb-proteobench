"""Bind one ProteoBench module to the configured APB quantification level."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, replace
from pathlib import Path

import polars as pl
from apb2.api import AnnotationCompiler, AnnotationError, ParsedLevels

from apb_proteobench.configuration.load import LoadedModule, load_module

_STORAGE_KEY = "proteobench"
PROTEOBENCH_SCHEMA_VERSION = "4"
"""The version of the ``proteobench`` record every ProteoBench step writes."""


@dataclass(frozen=True, slots=True)
class ProteoBenchAnnotation:
    """A dataset annotated with every module sample, awaiting its module provenance."""

    annotated: ParsedLevels
    module: LoadedModule

    def annotate(self) -> ParsedLevels:
        """Store the normalized module beside APB2's own annotation provenance."""
        metadata = deepcopy(self.annotated.metadata)
        section = metadata.setdefault(_STORAGE_KEY, {})
        if not isinstance(section, dict):
            raise AnnotationError("ProteoBench metadata must be an object")
        provenance = section.setdefault("provenance", {})
        if not isinstance(provenance, dict):
            raise AnnotationError("ProteoBench provenance must be an object")
        provenance["annotation"] = dict(self.module.metadata())
        section["schema_version"] = PROTEOBENCH_SCHEMA_VERSION
        return replace(self.annotated, metadata=metadata)


@dataclass(frozen=True, slots=True)
class ProteoBenchAnnotationParser:
    """A validated module source ready to bind to one APB result."""

    module: LoadedModule

    @classmethod
    def from_path(cls, path: Path, /) -> ProteoBenchAnnotationParser:
        """Decode and validate the complete module document once."""
        return cls(load_module(path))

    def parse(self, parsed: ParsedLevels, /) -> ProteoBenchAnnotation:
        """Match every run to exactly one module sample and every sample to a run."""
        level_name = self.module.settings.general.level
        if level_name not in parsed.levels:
            raise AnnotationError(
                f"ProteoBench module selects unavailable level {level_name!r}; "
                f"available={list(parsed.levels)}"
            )
        annotation = AnnotationCompiler("error").compile(_sample_frame(self.module)).parse(parsed)
        coverage = annotation.matches.levels[level_name].coverage
        if coverage.annotation_only_count:
            raise AnnotationError(
                "ProteoBench module contains samples absent from quantification; "
                f"count={coverage.annotation_only_count}, "
                f"examples={list(coverage.annotation_only_examples)}"
            )
        return ProteoBenchAnnotation(annotated=annotation.annotate().parsed, module=self.module)


def _sample_frame(module: LoadedModule) -> pl.DataFrame:
    """The module samples as a prolfquapp table keyed by raw file, with its aliases."""
    samples = module.settings.samples
    return pl.DataFrame(
        {
            "raw_file": [sample.raw_file for sample in samples],
            "raw_file_aliases": [
                list(
                    dict.fromkeys(
                        identifier
                        for identifier in (*sample.raw_file_aliases, sample.sample_name)
                        if identifier != sample.raw_file
                    )
                )
                for sample in samples
            ],
            "sample_name": [sample.sample_name for sample in samples],
            "condition": [sample.condition for sample in samples],
        },
        schema={
            "raw_file": pl.String,
            "raw_file_aliases": pl.List(pl.String),
            "sample_name": pl.String,
            "condition": pl.String,
        },
    )
