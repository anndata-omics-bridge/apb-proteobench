"""Compose resolved module settings from the module TOML and its SDRF."""

from __future__ import annotations

import re
from dataclasses import dataclass

import polars as pl
from apb2.api import SdrfSource

from apb_proteobench.configuration.schema import (
    ExpectedRatio,
    ModuleDocument,
    ModuleSettings,
    SampleSettings,
)

_DATA_FILE = "comment[data file]"
_ASSAY_NAME = "assay name"
_SPIKED_COMPOUND = "characteristics[spiked compound]"
_FACTOR_PREFIX = "factor value["
_QUANTITY = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*(\S+)\s*$")


@dataclass(frozen=True, slots=True)
class _Quantity:
    value: float
    unit: str


def compose_module_settings(document: ModuleDocument, sdrf: SdrfSource, /) -> ModuleSettings:
    """Project SDRF rows and species quantities onto the scoring configuration.

    Each SDRF row becomes one sample: its data-file basename is the run identifier, its assay
    name the sample name, and its single factor value the condition. Each species' expected
    A/B ratio is the quotient of its ``characteristics[spiked compound]`` quantities.
    """
    data_files = _strings(sdrf, sdrf.column(_DATA_FILE))
    unknown = sorted(set(document.run_aliases) - set(data_files))
    if unknown:
        raise ValueError(f"[run_aliases] names data files absent from the SDRF: {unknown}")
    basenames = sdrf.data_file_basenames().to_list()
    conditions = _strings(sdrf, _factor_column(sdrf))
    samples = [
        SampleSettings(
            raw_file=str(basename),
            raw_file_aliases=document.run_aliases.get(data_file, []),
            sample_name=sample_name,
            condition=condition,
        )
        for data_file, basename, sample_name, condition in zip(
            data_files,
            basenames,
            _strings(sdrf, sdrf.column(_ASSAY_NAME)),
            conditions,
            strict=True,
        )
    ]
    quantities = _species_quantities(sdrf, conditions)
    scored = {species.organism.lower() for species in document.species.values()}
    unscored = sorted(set(quantities) - scored)
    if unscored:
        raise ValueError(f"SDRF spiked organisms have no [species] entry: {unscored}")
    return ModuleSettings(
        species_expected_ratio={
            name: ExpectedRatio(
                A_vs_B=_ratio(quantities, species.organism),
                color=species.color,
            )
            for name, species in document.species.items()
        },
        species_mapper={species.suffix: name for name, species in document.species.items()},
        general=document.general,
        samples=samples,
    )


def _factor_column(sdrf: SdrfSource) -> str:
    headers = list(
        dict.fromkeys(
            header
            for header in sdrf.source.headers
            if header.strip().lower().startswith(_FACTOR_PREFIX)
        )
    )
    if len(headers) != 1:
        raise ValueError(
            f"ProteoBench module SDRF needs exactly one factor value column; found {headers}"
        )
    return sdrf.column(headers[0])


def _species_quantities(
    sdrf: SdrfSource,
    conditions: list[str],
) -> dict[str, dict[str, _Quantity]]:
    quantities: dict[str, dict[str, _Quantity]] = {}
    rows = zip(
        conditions,
        *(_strings(sdrf, column) for column in sdrf.columns(_SPIKED_COMPOUND)),
        strict=True,
    )
    for condition, *cells in rows:
        organisms: set[str] = set()
        for cell in cells:
            organism, quantity = _spiked(cell)
            if organism in organisms:
                raise ValueError(f"SDRF row repeats spiked organism {organism!r}")
            organisms.add(organism)
            known = quantities.setdefault(organism, {}).setdefault(condition, quantity)
            if known != quantity:
                raise ValueError(
                    f"SDRF rows of condition {condition!r} disagree on {organism!r}: "
                    f"{known.value:g} {known.unit} and {quantity.value:g} {quantity.unit}"
                )
    return quantities


def _spiked(cell: str) -> tuple[str, _Quantity]:
    fields = {
        key.strip().upper(): value.strip()
        for key, separator, value in (part.partition("=") for part in cell.split(";"))
        if separator
    }
    organism = fields.get("SP")
    match = _QUANTITY.match(fields.get("QY", ""))
    if not organism or match is None:
        raise ValueError(f"spiked compound needs SP=<organism> and QY=<amount unit>: {cell!r}")
    return organism.lower(), _Quantity(value=float(match.group(1)), unit=match.group(2))


def _ratio(quantities: dict[str, dict[str, _Quantity]], organism: str) -> float:
    by_condition = quantities.get(organism.lower(), {})
    missing = [condition for condition in ("A", "B") if condition not in by_condition]
    if missing:
        raise ValueError(f"SDRF has no {organism!r} quantity for condition(s) {missing}")
    a, b = by_condition["A"], by_condition["B"]
    if a.unit != b.unit:
        raise ValueError(f"{organism!r} quantities use different units: {a.unit!r}, {b.unit!r}")
    if b.value == 0:
        raise ValueError(f"{organism!r} quantity in condition 'B' must be positive")
    return a.value / b.value


def _strings(sdrf: SdrfSource, column: str) -> list[str]:
    values = sdrf.source.frame.get_column(column).cast(pl.String).to_list()
    if any(value is None for value in values):
        raise ValueError(f"SDRF column {column!r} must not contain empty cells")
    return [str(value) for value in values]
