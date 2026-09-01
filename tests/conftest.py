"""Shared synthetic HYE values for ProteoBench tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from apb2.result_facade import (
    FinalLayerTable,
    ObsFinal,
    ParsedLevel,
    ParsedLevels,
    VarFinal,
)
from numpy.typing import NDArray

from apb_proteobench.calculation.contracts import QuantitativeLevelInput
from apb_proteobench.configuration.schema import (
    ExpectedRatio,
    ModuleGeneral,
    ModuleSettings,
    SampleSettings,
)


def matrix_values() -> NDArray[np.float64]:
    """Return the small legacy-compatible quantitative matrix."""
    return np.asarray(
        [
            [10, 20, 5, 10, 10, 10],
            [10, 20, 5, 10, 10, 0],
            [10, 10, 20, 10, 10, np.nan],
            [10, 10, 20, 10, 10, -1],
        ],
        dtype=np.float64,
    )


def module_settings(*, hy: bool = False) -> ModuleSettings:
    """Return HYE or HY settings for the same config-driven calculation."""
    ratios = {"YEAST": ExpectedRatio(A_vs_B=2.0)}
    mapper = {"_YEAST": "YEAST"}
    if not hy:
        ratios["ECOLI"] = ExpectedRatio(A_vs_B=0.25)
        mapper["_ECOLI"] = "ECOLI"
    ratios["HUMAN"] = ExpectedRatio(A_vs_B=1.0)
    mapper["_HUMAN"] = "HUMAN"
    return ModuleSettings(
        species_expected_ratio=ratios,
        species_mapper=mapper,
        general=ModuleGeneral(min_count_multispec=1, level="ion"),
        samples=[
            SampleSettings(raw_file="run_A1", sample_name="A1", condition="A"),
            SampleSettings(raw_file="run_A2", sample_name="A2", condition="A"),
            SampleSettings(raw_file="run_B1", sample_name="B1", condition="B"),
            SampleSettings(raw_file="run_B2", sample_name="B2", condition="B"),
        ],
    )


def quantitative_input(*, hy: bool = False) -> QuantitativeLevelInput:
    """Return calculation values without an APB storage container."""
    proteins = [
        "P1_HUMAN",
        "P2_YEAST",
        "P3_HUMAN" if hy else "P3_ECOLI",
        "Cont_P4_HUMAN",
        "P5_HUMAN_YEAST",
        "P6_UNKNOWN",
    ]
    return QuantitativeLevelInput(
        observations=pd.DataFrame(
            {
                "Run": ["run_A1", "run_A2", "run_B1", "run_B2"],
                "sample_name": ["A1", "A2", "B1", "B2"],
                "condition": ["A", "A", "B", "B"],
            }
        ),
        matrix=matrix_values(),
        feature_ids=pd.Index(["H/2", "Y/2", "E/2", "C/2", "M/2", "N/2"]),
        reported_proteins=pd.Series(proteins),
        level="ion",
    )


def parsed_result() -> ParsedLevels:
    """Return an unannotated one-level numeric APB result."""
    inputs = quantitative_input()
    features = inputs.feature_ids.astype(str).tolist()
    proteins = inputs.reported_proteins.astype(str).tolist()
    matrix = matrix_values()
    level = ParsedLevel(
        obs=ObsFinal(
            frame=pl.DataFrame({"Run": ["run_A1", "run_A2", "run_B1", "run_B2"]}),
            key_columns=("Run",),
        ),
        var=VarFinal(
            frame=pl.DataFrame({"feature": features, "Protein_Ids": proteins}),
            key_columns=("feature",),
        ),
        primary_layer_name="Intensity",
        uns={
            "quantification_level": "ion",
            "software_name": "Synthetic",
            "matrix_values_projected": True,
            "column_roles": {"protein_assignment": "Protein_Ids"},
        },
        layers={
            "Intensity": FinalLayerTable(
                layer_name="Intensity",
                var_key_columns=("feature",),
                values=pl.DataFrame(
                    {
                        "feature": features,
                        **{f"obs_{index}": matrix[index, :] for index in range(matrix.shape[0])},
                    }
                ),
            )
        },
        obsm={},
        varm={},
        obsp={},
        varp={},
    )
    return ParsedLevels(levels={"ion": level}, uns={"produced_by": "apb2"})


def write_module(path: Path, /, *, alias: bool = False, hy: bool = False) -> None:
    """Write the small module in the existing ProteoBench TOML shape."""
    settings = module_settings(hy=hy)
    lines = []
    for species, ratio in settings.species_expected_ratio.items():
        lines.extend([f"[species_expected_ratio.{species}]", f"A_vs_B = {ratio.a_vs_b}"])
    lines.append("[species_mapper]")
    lines.extend(f'"{flag}" = "{species}"' for flag, species in settings.species_mapper.items())
    lines.extend(["[general]", 'level = "ion"', "min_count_multispec = 1"])
    for sample in settings.samples:
        lines.extend(
            [
                "[[samples]]",
                f'raw_file = "{sample.raw_file}"',
                *([f'raw_file_alias = "{sample.raw_file}_alias"'] if alias else []),
                f'sample_name = "{sample.sample_name}"',
                f'condition = "{sample.condition}"',
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
