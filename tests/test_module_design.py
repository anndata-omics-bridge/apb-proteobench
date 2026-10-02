"""Composition of resolved module settings from the module TOML and its SDRF."""

from __future__ import annotations

from pathlib import Path

import pytest

from apb_proteobench.configuration.load import load_module
from conftest import write_module


def _replace(path: Path, old: str, new: str, /, *, count: int = -1) -> None:
    text = path.read_text(encoding="utf-8")
    assert old in text
    path.write_text(text.replace(old, new, count), encoding="utf-8")


def _module(tmp_path: Path, /, *, alias: bool = False) -> Path:
    module = tmp_path / "module.toml"
    write_module(module, alias=alias)
    return module


def test_module_settings_combine_sdrf_rows_with_toml_aliases(tmp_path: Path) -> None:
    loaded = load_module(_module(tmp_path, alias=True))

    sample = loaded.settings.samples[0]
    assert (sample.raw_file, sample.raw_file_aliases, sample.sample_name, sample.condition) == (
        "run_A1",
        ["run_A1_alias"],
        "A1",
        "A",
    )
    assert {
        species: ratio.a_vs_b for species, ratio in loaded.settings.species_expected_ratio.items()
    } == {"YEAST": 2.0, "ECOLI": 0.25, "HUMAN": 1.0}
    assert loaded.source.sdrf.name == "module.sdrf.tsv"
    assert loaded.metadata()["source"] == loaded.source.as_json()


def test_quantities_must_share_a_unit(tmp_path: Path) -> None:
    module = _module(tmp_path)
    _replace(
        module.with_suffix(".sdrf.tsv"),
        "SP=Saccharomyces cerevisiae;QY=2 ng",
        "SP=Saccharomyces cerevisiae;QY=2 pg",
    )

    with pytest.raises(ValueError, match="different units"):
        load_module(module)


def test_rows_of_one_condition_must_agree(tmp_path: Path) -> None:
    module = _module(tmp_path)
    _replace(
        module.with_suffix(".sdrf.tsv"),
        "SP=Saccharomyces cerevisiae;QY=2 ng",
        "SP=Saccharomyces cerevisiae;QY=3 ng",
        count=1,
    )

    with pytest.raises(ValueError, match="disagree on 'saccharomyces cerevisiae'"):
        load_module(module)


def test_every_scored_species_needs_quantities(tmp_path: Path) -> None:
    module = _module(tmp_path)
    _replace(
        module,
        "[general]",
        '[species.MOUSE]\norganism = "mus musculus"\nsuffix = "_MOUSE"\n[general]',
    )

    with pytest.raises(ValueError, match="no 'mus musculus' quantity"):
        load_module(module)


def test_every_spiked_organism_needs_a_species_entry(tmp_path: Path) -> None:
    module = _module(tmp_path)
    _replace(module, '[species.ECOLI]\norganism = "escherichia coli"\nsuffix = "_ECOLI"\n', "")

    with pytest.raises(ValueError, match=r"no \[species\] entry: \['escherichia coli'\]"):
        load_module(module)


def test_spiked_compounds_need_organism_and_quantity(tmp_path: Path) -> None:
    module = _module(tmp_path)
    _replace(module.with_suffix(".sdrf.tsv"), "QY=1 ng", "QY=unknown", count=1)

    with pytest.raises(ValueError, match="needs SP=<organism> and QY=<amount unit>"):
        load_module(module)


def test_module_sdrf_needs_exactly_one_factor_value(tmp_path: Path) -> None:
    module = _module(tmp_path)
    sdrf = module.with_suffix(".sdrf.tsv")
    lines = sdrf.read_text(encoding="utf-8").splitlines()
    sdrf.write_text(
        "\n".join([f"{lines[0]}\tfactor value[batch]", *(f"{line}\t1" for line in lines[1:])])
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="exactly one factor value column"):
        load_module(module)


def test_run_aliases_must_name_sdrf_data_files(tmp_path: Path) -> None:
    module = _module(tmp_path, alias=True)
    _replace(module, '"run_A1.raw"', '"elsewhere.raw"')

    with pytest.raises(ValueError, match=r"absent from the SDRF: \['elsewhere.raw'\]"):
        load_module(module)
