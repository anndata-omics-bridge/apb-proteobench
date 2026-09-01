"""Decode one ProteoBench module document and retain portable source evidence."""

from __future__ import annotations

import hashlib
import json
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import cast

from apb2.result_facade import JsonValue

from apb_proteobench.configuration.schema import ModuleSettings


@dataclass(frozen=True, slots=True)
class ModuleSource:
    """Portable identity of the authored ProteoBench module document."""

    name: str
    sha256: str
    format: str = "proteobench-module-toml"

    def as_json(self) -> dict[str, JsonValue]:
        """Return the JSON-compatible provenance record."""
        return {"name": self.name, "sha256": self.sha256, "format": self.format}


@dataclass(frozen=True, slots=True)
class LoadedModule:
    """Validated module settings paired with their source identity."""

    settings: ModuleSettings
    source: ModuleSource

    def metadata(self) -> dict[str, JsonValue]:
        """Return normalized configuration and portable source provenance."""
        configuration = json.loads(self.settings.model_dump_json(by_alias=True))
        if not isinstance(configuration, dict):
            raise TypeError("ProteoBench module serialization did not produce an object")
        configuration["samples"] = {
            "raw_file": [sample.raw_file for sample in self.settings.samples],
            "sample_name": [sample.sample_name for sample in self.settings.samples],
            "condition": [sample.condition for sample in self.settings.samples],
        }
        return {
            "source": self.source.as_json(),
            "configuration": cast(dict[str, JsonValue], configuration),
        }


def load_module(path: Path, /) -> LoadedModule:
    """Load and validate one complete ProteoBench module TOML."""
    source = path.expanduser().resolve()
    payload = source.read_bytes()
    settings = ModuleSettings.model_validate(tomllib.loads(payload.decode("utf-8")))
    return LoadedModule(
        settings=settings,
        source=ModuleSource(name=source.name, sha256=hashlib.sha256(payload).hexdigest()),
    )
