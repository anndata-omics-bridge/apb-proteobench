"""Entrapment module settings, read from ProteoBench's packaged module document."""

from __future__ import annotations

import hashlib
import tomllib
from importlib.resources import as_file, files

from apb2.api import SdrfSource
from pydantic import BaseModel, ConfigDict, Field

from apb_proteobench.configuration.schema import QuantificationLevel

ENTRAPMENT_MODULE_NAMES = ("entrapment_dia_astral",)
"""Stable names of modules supported by the entrapment scorer."""


class EntrapmentModuleSettings(BaseModel):
    """What one entrapment module scores, and ProteoBench's fixed scoring constants."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    name: str
    level: QuantificationLevel
    mapping_file: str = Field(description="Where ProteoBench publishes the pair file")
    sdrf: str = Field(description="Packaged module SDRF describing the runs")
    sdrf_sha256: str
    max_missing_fraction: float = Field(
        default=0.01, description="ProteoBench's tolerated fraction of unmapped peptides"
    )
    curve_intervals: int = Field(
        default=10, description="ProteoBench's evenly spaced FDP-curve thresholds"
    )


def load_packaged_entrapment_module(name: str, /) -> EntrapmentModuleSettings:
    """Load one packaged entrapment module.

    Raises:
        ValueError: The name is not a packaged entrapment module.
    """
    if name not in ENTRAPMENT_MODULE_NAMES:
        raise ValueError(
            f"unknown packaged entrapment module {name!r}; packaged: {list(ENTRAPMENT_MODULE_NAMES)}"
        )
    resources = files("apb_proteobench.data.modules")
    document = tomllib.loads(resources.joinpath(f"{name}.toml").read_text("utf-8"))
    sdrf = resources.joinpath(document["sdrf"])
    with as_file(sdrf) as path:
        SdrfSource.read(path)
    general = document["general"]
    return EntrapmentModuleSettings(
        name=name,
        level=general["level"],
        mapping_file=general["mapping_file"],
        sdrf=document["sdrf"],
        sdrf_sha256=hashlib.sha256(sdrf.read_bytes()).hexdigest(),
    )
