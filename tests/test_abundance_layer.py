"""A module's abundance normalisation picks the scored layer through apb-catalog."""

from __future__ import annotations

from pathlib import Path
from typing import cast

import pytest
from apb2.api import (
    JsonValue,
    ParsedLevels,
    packaged_rule_declarations,
    read_parsed_levels,
    write_parsed_levels,
)

from apb_proteobench.api import ProteoBenchAnalysisResult, ProteoBenchAnalyzer
from apb_proteobench.cli.app import app
from apb_proteobench.configuration.load import load_module
from conftest import parsed_result, write_module


def _rule_json(rule: str) -> str:
    return packaged_rule_declarations()[(rule, "ion")][0]


def _diann() -> ParsedLevels:
    """DIA-NN 1.8's ion level: normalised quantity primary, un-normalised beside it."""
    return parsed_result(
        ("Precursor_Normalised", "Precursor_Quantity"), _rule_json("diann/v1_8/rules.json")
    )


def _analyze(
    folder: Path,
    parsed: ParsedLevels,
    /,
    *,
    normalisation: str | None,
    layers: tuple[str, ...] | None = None,
) -> ProteoBenchAnalysisResult:
    folder.mkdir(exist_ok=True)
    module = folder / "module.toml"
    write_module(module, normalisation=normalisation)
    return ProteoBenchAnalyzer(load_module(module), layers=layers).analyze(parsed)


def _entry(result: ProteoBenchAnalysisResult) -> dict[str, JsonValue] | None:
    """The level summary entry for the module's abundance normalisation, if written."""
    tool = result.parsed.levels["ion"].metadata["proteobench"]
    assert isinstance(tool, dict)
    summary = tool["summary"]
    assert isinstance(summary, list)
    entries = [
        entry
        for entry in summary
        if isinstance(entry, dict) and entry["name"] == "abundance_normalisation"
    ]
    return entries[0] if entries else None


def _general(result: ProteoBenchAnalysisResult, layer: str) -> dict[str, object]:
    return cast(dict[str, object], result.submission(layer).settings["general"])


def test_the_module_scores_the_layer_the_catalogue_names(tmp_path: Path) -> None:
    result = _analyze(tmp_path, _diann(), normalisation="none")

    assert list(result.layers) == ["Precursor_Quantity"]
    entry = _entry(result)
    assert entry is not None
    assert (entry["value"], entry["status"], entry["layer"]) == (True, "ok", "Precursor_Quantity")
    catalog = result.parsed.metadata["catalog"]
    assert isinstance(catalog, dict)
    assert "abundance" in catalog


@pytest.mark.parametrize("rule", ["peaks/rules.json", None], ids=["missing", "unknown"])
def test_without_a_catalogued_layer_the_primary_layer_needs_attention(
    tmp_path: Path, rule: str | None
) -> None:
    """PEAKS writes no un-normalised layer; a result bound to no rule leaves the catalogue unsure."""
    parsed = parsed_result(("Normalized_Area",), None if rule is None else _rule_json(rule))

    result = _analyze(tmp_path, parsed, normalisation="none")

    assert list(result.layers) == ["Normalized_Area"]
    entry = _entry(result)
    assert entry is not None
    assert (entry["value"], entry["status"], entry["layer"]) == (
        False,
        "attention",
        "Normalized_Area",
    )


def test_a_named_layer_overrides_the_module_normalisation(tmp_path: Path) -> None:
    result = _analyze(tmp_path, _diann(), normalisation="none", layers=("Precursor_Normalised",))

    assert list(result.layers) == ["Precursor_Normalised"]
    assert _entry(result) is None
    assert "catalog" not in result.parsed.metadata


def test_an_unset_normalisation_stays_out_of_the_hashed_settings(tmp_path: Path) -> None:
    """Modules without the setting keep the ProteoBot identity they had before it existed."""
    plain = _analyze(tmp_path / "plain", parsed_result(), normalisation=None)
    plasma = _analyze(tmp_path / "plasma", _diann(), normalisation="none")

    assert "abundance_normalisation" not in _general(plain, "Intensity")
    assert _general(plasma, "Precursor_Quantity")["abundance_normalisation"] == "none"


def test_cli_layer_x_is_the_modules_layer(tmp_path: Path) -> None:
    source = tmp_path / "source.parquet"
    module = tmp_path / "module.toml"
    target = tmp_path / "scored.parquet"
    write_module(module, normalisation="none")
    write_parsed_levels(_diann(), source)

    status = app(
        ["benchmark", str(source), str(module), str(target)],
        exit_on_error=False,
        result_action="return_value",
    )

    assert status == 0
    assert "proteobench:Precursor_Quantity" in read_parsed_levels(target).levels["ion"].varm
