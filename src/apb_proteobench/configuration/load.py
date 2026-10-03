"""Decode one ProteoBench module TOML and SDRF pair and retain portable source evidence."""

from __future__ import annotations

import hashlib
import json
import tomllib
from dataclasses import dataclass
from importlib.resources import as_file, files
from pathlib import Path
from typing import cast

from apb2.api import JsonValue, SdrfSource

from apb_proteobench.configuration.design import compose_module_settings
from apb_proteobench.configuration.schema import ModuleDocument, ModuleSettings

SUPPORTED_MODULE_NAMES = (
    "dda_astral",
    "dda_peptidoform",
    "dda_qexactive",
    "dia_aif",
    "dia_astral",
    "dia_diapasef",
    "dia_plasma",
    "dia_singlecell",
    "dia_zenotof",
)
"""Stable names of modules supported by the quantitative scorer."""

PACKAGED_MODULE_NAMES = (
    *SUPPORTED_MODULE_NAMES,
    "denovo_dda_hcd",
    "entrapment_dia_astral",
)
"""Stable names of all upstream ProteoBench module documents in the package."""


@dataclass(frozen=True, slots=True)
class SourceFile:
    """Portable identity of one authored file."""

    name: str
    sha256: str

    def as_json(self) -> dict[str, JsonValue]:
        """Return the JSON-compatible provenance record."""
        return {"name": self.name, "sha256": self.sha256}


@dataclass(frozen=True, slots=True)
class ModuleSource:
    """Portable identity of the authored module TOML and its SDRF."""

    name: str
    sha256: str
    sdrf: SourceFile
    format: str = "proteobench-module-toml-sdrf"

    def as_json(self) -> dict[str, JsonValue]:
        """Return the JSON-compatible provenance record."""
        return {
            "name": self.name,
            "sha256": self.sha256,
            "format": self.format,
            "sdrf": self.sdrf.as_json(),
        }


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
    """Load and validate one module TOML and the SDRF it names, relative to its directory."""
    source = path.expanduser().resolve()
    payload = source.read_bytes()
    document = _document(payload)
    return _load_module(payload, source.name, document, source.parent / document.sdrf)


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
    resources = files("apb_proteobench.data.modules")
    resource = resources.joinpath(f"{name}.toml")
    payload = resource.read_bytes()
    document = _document(payload)
    with as_file(resources.joinpath(document.sdrf)) as sdrf_path:
        return _load_module(payload, resource.name, document, sdrf_path)


def _document(payload: bytes) -> ModuleDocument:
    return ModuleDocument.model_validate(tomllib.loads(payload.decode("utf-8")))


def _load_module(
    payload: bytes,
    name: str,
    document: ModuleDocument,
    sdrf_path: Path,
) -> LoadedModule:
    sdrf = SdrfSource.read(sdrf_path)
    return LoadedModule(
        settings=compose_module_settings(document, sdrf),
        source=ModuleSource(
            name=name,
            sha256=hashlib.sha256(payload).hexdigest(),
            sdrf=SourceFile(
                name=sdrf_path.name,
                sha256=hashlib.sha256(sdrf_path.read_bytes()).hexdigest(),
            ),
        ),
    )
