"""Decode one ProteoBench module document and retain portable source evidence."""

from __future__ import annotations

import hashlib
import json
import tomllib
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import cast

from apb2.result_facade import JsonValue

from apb_proteobench.configuration.schema import ModuleSettings

SUPPORTED_MODULE_NAMES = (
    "dda_astral",
    "dda_peptidoform",
    "dda_qexactive",
    "dia_aif",
    "dia_astral",
    "dia_diapasef",
    "dia_singlecell",
    "dia_zenotof",
)
"""Stable names of modules supported by the quantitative scorer."""

PACKAGED_MODULE_NAMES = (
    *SUPPORTED_MODULE_NAMES,
    "dia_plasma",
    "denovo_dda_hcd",
    "entrapment_dia_astral",
)
"""Stable names of all upstream ProteoBench module documents in the package."""


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
    return _load_module(payload, source.name)


def available_modules() -> tuple[str, ...]:
    """Return the packaged modules supported by the quantitative scorer."""
    return SUPPORTED_MODULE_NAMES


def packaged_module_names() -> tuple[str, ...]:
    """Return the stable names of all packaged ProteoBench module documents."""
    return PACKAGED_MODULE_NAMES


def load_packaged_module(name: str, /) -> LoadedModule:
    """Load one packaged quantitative benchmark module by its stable name.

    Args:
        name: A value returned by :func:`available_modules`.

    Returns:
        The validated settings and source identity of the packaged module.

    Raises:
        ValueError: The name does not identify a supported packaged module.
    """
    if name in PACKAGED_MODULE_NAMES and name not in SUPPORTED_MODULE_NAMES:
        raise ValueError(
            f"packaged ProteoBench module {name!r} is not supported by the quantitative scorer"
        )
    if name not in SUPPORTED_MODULE_NAMES:
        raise ValueError(
            f"unknown packaged ProteoBench module {name!r}; available: "
            f"{list(SUPPORTED_MODULE_NAMES)}"
        )
    resource = files("apb_proteobench.data.modules").joinpath(f"{name}.toml")
    return _load_module(resource.read_bytes(), resource.name)


def _load_module(payload: bytes, name: str) -> LoadedModule:
    settings = ModuleSettings.model_validate(tomllib.loads(payload.decode("utf-8")))
    return LoadedModule(
        settings=settings,
        source=ModuleSource(name=name, sha256=hashlib.sha256(payload).hexdigest()),
    )
