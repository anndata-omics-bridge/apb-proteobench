"""Packaged ProteoBench module catalogue."""

from __future__ import annotations

import hashlib
import subprocess
import sys
from importlib.resources import as_file, files
from pathlib import Path

import pytest

from apb_proteobench.configuration.entrapment import load_packaged_entrapment_module
from apb_proteobench.configuration.load import (
    PACKAGED_MODULE_NAMES,
    SUPPORTED_MODULE_NAMES,
    load_packaged_module,
)

EXPECTED_SUPPORTED_MODULE_HASHES = {
    "dda_astral": "f556f0aa761eca058300eb0874cd764e076015cf81030f1194dbaf966195cd7d",
    "dda_peptidoform": "f27fcdb89d53a06722282720fd16c89f627d9967036bb551dccc1f9a2cacce22",
    "dda_qexactive": "80fde2bd093de08b16e11bf1cdf74d32626573e7c752184de544c1ea0dfa251e",
    "dia_aif": "984dfa21f80fa4631bd6cbe120c662636014a7570f51ff0b36d92be322de51ae",
    "dia_astral": "35bf718c54b964c40b3b4ff2e352bc82e9fd492fdd1c6d4da9d51cd588475568",
    "dia_diapasef": "6072039f558cc979d41d2ca29c8134cab89a1b2eba1cd0171387464cca19d103",
    "dia_plasma": "50901ea040970407721a417e8ebf8decadf21a8cd9d4ae83b7f6657fd035bef5",
    "dia_singlecell": "de1f4e132aa5bd6ef7ed9538225a70fe8e69aac9ee2d7e2031b3a6dfac661eff",
    "dia_zenotof": "b8daad3330a1b93b8f8ba24b4c081c17a77a9aaab4a3ca5b8ad2e17e8c6d525b",
}

EXPECTED_SDRF_HASHES = {
    "dda_astral": "d63efd1e9f0894ce025128c5e4b2e2e4f64ed626968ecf0da0d37db4ee223494",
    "dda_peptidoform": "ac1addb6e6ea42e58f89eb63d4d5e321830c8ca7e15c7a8228e17dbb70849cb5",
    "dda_qexactive": "ac1addb6e6ea42e58f89eb63d4d5e321830c8ca7e15c7a8228e17dbb70849cb5",
    "dia_aif": "fdf5165d3e3ec79442a915193c24931920b9809262a12c466cdb198250783a76",
    "dia_astral": "9979dd3ebbdc04f6820a5b6866a139fb9aa61da95c402820f0f63dd8606ccb0e",
    "dia_diapasef": "b5b9f831b143539cc42dbd13f6c78a6096e069f3cfdf6071d6a99f24c3321171",
    "dia_plasma": "43c9142f635086db339357fa395a46c031cd57a350df0e7a681a9d31f48cf027",
    "dia_singlecell": "36ae5f085ef283eef9798760432bcfdcc12b52a5c3b3bafd0a475ff3d19431e3",
    "dia_zenotof": "08a2a0a964a998b7347a3ec38dd80789437ac24d612ded4f39d15722839ac7e9",
}

EXPECTED_PACKAGED_MODULE_HASHES = {
    **EXPECTED_SUPPORTED_MODULE_HASHES,
    "denovo_dda_hcd": "d53cd86c02228c57cc35aba996e44c4739ce49e4de8fe57bf556a266d18460b7",
    "entrapment_dia_astral": ("c6740cd8f2954b9b137bc870437331c53eb92c09a0e53650a30884cd79692db1"),
}


def test_every_supported_module_is_packaged_validated_and_pinned() -> None:
    assert tuple(EXPECTED_SUPPORTED_MODULE_HASHES) == SUPPORTED_MODULE_NAMES

    loaded = {name: load_packaged_module(name) for name in SUPPORTED_MODULE_NAMES}

    assert {name: module.source.sha256 for name, module in loaded.items()} == (
        EXPECTED_SUPPORTED_MODULE_HASHES
    )
    assert {name: module.source.sdrf.sha256 for name, module in loaded.items()} == (
        EXPECTED_SDRF_HASHES
    )
    assert {module.settings.general.level for module in loaded.values()} == {
        "ion",
        "peptidoform",
    }
    assert all(
        len(module.settings.samples) == module.settings.general.max_nr_observed
        for module in loaded.values()
    )
    assert {module.source.format for module in loaded.values()} == {"proteobench-module-toml-sdrf"}


def test_every_upstream_module_document_is_packaged_and_pinned() -> None:
    assert tuple(EXPECTED_PACKAGED_MODULE_HASHES) == PACKAGED_MODULE_NAMES

    module_resources = files("apb_proteobench.data.modules")
    hashes = {
        name: hashlib.sha256(module_resources.joinpath(f"{name}.toml").read_bytes()).hexdigest()
        for name in PACKAGED_MODULE_NAMES
    }

    assert hashes == EXPECTED_PACKAGED_MODULE_HASHES


def test_packaged_ratios_preserve_plot_colors() -> None:
    ratios = load_packaged_module("dda_qexactive").settings.species_expected_ratio

    assert ratios["YEAST"].color == "#88ccef"


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("dda_qexactive", {"YEAST": 2.0, "ECOLI": 0.25, "HUMAN": 1.0}),
        ("dia_singlecell", {"YEAST": 0.2, "HUMAN": 1.2}),
        ("dia_plasma", {"YEAST": 1 / 3, "ECOLI": 2.0, "HUMAN": 1.0}),
    ],
)
def test_packaged_ratios_are_quotients_of_sdrf_species_quantities(
    name: str,
    expected: dict[str, float],
) -> None:
    settings = load_packaged_module(name).settings

    assert {
        species: ratio.a_vs_b for species, ratio in settings.species_expected_ratio.items()
    } == expected
    assert list(settings.species_mapper) == [f"_{species}" for species in expected]


def test_peptidoform_module_names_the_qexactive_raw_files() -> None:
    peptidoform = load_packaged_module("dda_peptidoform").settings.samples
    qexactive = load_packaged_module("dda_qexactive").settings.samples

    assert [sample.raw_file for sample in peptidoform] == [sample.raw_file for sample in qexactive]
    assert peptidoform[0].raw_file_aliases == ["abundance_A_1", "A_1"]


@pytest.mark.parametrize("name", [*EXPECTED_SDRF_HASHES, "entrapment_dia_astral"])
def test_packaged_sdrf_passes_the_sdrf_proteomics_validator(name: str) -> None:
    resource = files("apb_proteobench.data.modules").joinpath(f"{name}.sdrf.tsv")
    templates = ["-t", "ms-proteomics"]
    if "Data-independent acquisition" in resource.read_text(encoding="utf-8"):
        templates += ["-t", "dia-acquisition"]
    validator = Path(sys.executable).with_name("parse_sdrf")

    with as_file(resource) as path:
        completed = subprocess.run(
            [str(validator), "validate-sdrf", "-s", str(path), *templates, "--use_ols_cache_only"],
            capture_output=True,
            check=False,
            text=True,
        )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Everything seems to be fine" in completed.stdout


@pytest.mark.parametrize(
    "name",
    ["denovo_dda_hcd", "entrapment_dia_astral"],
)
def test_packaged_but_unsupported_module_is_explicit(name: str) -> None:
    with pytest.raises(ValueError, match="is not supported by the quantitative scorer"):
        load_packaged_module(name)


def test_unknown_packaged_module_reports_available_names() -> None:
    with pytest.raises(ValueError, match="unknown packaged ProteoBench module 'missing'"):
        load_packaged_module("missing")


def test_entrapment_module_records_its_sdrf() -> None:
    module = load_packaged_entrapment_module("entrapment_dia_astral")

    assert module.sdrf == "entrapment_dia_astral.sdrf.tsv"
    assert module.sdrf_sha256 == "5ff1a9af225505a5631126e4753aa7ca32a377d69afe422c2fe58535894ca72f"
