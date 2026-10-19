"""Validated storage schema for ProteoBench module settings."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

QuantificationLevel = Literal["ion", "peptidoform"]


class _SettingsModel(BaseModel):
    model_config = ConfigDict(extra="ignore")


class SampleSettings(_SettingsModel):
    """One run in a ProteoBench experiment design."""

    raw_file: str = Field(min_length=1)
    raw_file_aliases: list[str] = Field(default_factory=list[str])
    sample_name: str = Field(min_length=1)
    condition: str = Field(min_length=1)


class ExpectedRatio(_SettingsModel):
    """Expected abundance ratio for one species."""

    a_vs_b: float = Field(alias="A_vs_B", gt=0)
    color: str | None = None


class ModuleGeneral(_SettingsModel):
    """Scoring fields from the module ``[general]`` table."""

    min_count_multispec: int = Field(ge=1)
    level: QuantificationLevel
    default_cutoff_min_feature: int = Field(default=1, ge=1)
    max_nr_observed: int = Field(default=6, ge=1)

    @model_validator(mode="after")
    def _validate_cutoffs(self) -> ModuleGeneral:
        if self.default_cutoff_min_feature > self.max_nr_observed:
            raise ValueError("default_cutoff_min_feature must not exceed max_nr_observed")
        return self


class SpeciesSettings(_SettingsModel):
    """One scored species: its SDRF organism, FASTA header suffix, and plot colour."""

    organism: str = Field(min_length=1)
    suffix: str = Field(min_length=1)
    color: str | None = None


class ModuleDocument(_SettingsModel):
    """The authored module TOML: scoring settings the module SDRF cannot hold.

    ``sdrf`` names the module SDRF relative to the TOML's directory. ``run_aliases`` maps an
    SDRF ``comment[data file]`` value to further run names that vendor tables report.
    """

    sdrf: str = Field(min_length=1)
    species: dict[str, SpeciesSettings] = Field(min_length=1)
    run_aliases: dict[str, list[str]] = Field(default_factory=dict[str, list[str]])
    general: ModuleGeneral

    @model_validator(mode="after")
    def _validate_species(self) -> ModuleDocument:
        organisms = [species.organism.lower() for species in self.species.values()]
        if len(organisms) != len(set(organisms)):
            raise ValueError("[species] organisms must be unique")
        suffixes = [species.suffix for species in self.species.values()]
        if len(suffixes) != len(set(suffixes)):
            raise ValueError("[species] suffixes must be unique")
        return self


class ModuleSettings(_SettingsModel):
    """Resolved experiment design and HYE scoring configuration.

    Loading composes it from a ``ModuleDocument`` and the module SDRF; it is also the
    normalized configuration recorded in result provenance.
    """

    species_expected_ratio: dict[str, ExpectedRatio]
    species_mapper: dict[str, str]
    general: ModuleGeneral
    samples: list[SampleSettings] = Field(min_length=1)

    @model_validator(mode="after")
    def _validate_design(self) -> ModuleSettings:
        species = list(self.species_mapper.values())
        if len(species) != len(set(species)):
            raise ValueError("species_mapper values must be unique")
        if set(species) != set(self.species_expected_ratio):
            raise ValueError(
                "species_mapper values must equal species_expected_ratio keys; "
                f"mapper={species}, ratios={list(self.species_expected_ratio)}"
            )

        raw_files = [sample.raw_file for sample in self.samples]
        if len(raw_files) != len(set(raw_files)):
            raise ValueError("[[samples]].raw_file values must be unique")
        sample_names = [sample.sample_name for sample in self.samples]
        if len(sample_names) != len(set(sample_names)):
            raise ValueError("[[samples]].sample_name values must be unique")

        conditions = {sample.condition for sample in self.samples}
        if not {"A", "B"} <= conditions:
            raise ValueError(
                "ProteoBench HYE scoring requires conditions 'A' and 'B'; "
                f"found {sorted(conditions)}"
            )
        return self
