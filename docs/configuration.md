# Module configuration

The authored configuration is the existing ProteoBench module TOML:

```toml
[species_expected_ratio.YEAST]
A_vs_B = 2.0

[species_expected_ratio.HUMAN]
A_vs_B = 1.0

[species_mapper]
_YEAST = "YEAST"
_HUMAN = "HUMAN"

[general]
min_count_multispec = 1
level = "ion"
default_cutoff_min_feature = 1
max_nr_observed = 6

[[samples]]
raw_file = "run_A1"
raw_file_alias = "run_A1_uncalibrated"
sample_name = "A1"
condition = "A"
```

- `species_expected_ratio` maps each reported species to its expected A/B abundance ratio.
- `species_mapper` maps protein-identifier patterns to those same species names.
- `min_count_multispec` controls exclusion of features assigned to multiple species.
- `level` selects the single APB2 quantification level to annotate and score.
- `default_cutoff_min_feature` selects the aggregate projection shown at the top of the score.
- `max_nr_observed` controls the complete cutoff-indexed score range.
- Each `samples` record supplies a matching `raw_file`, optional alias, displayed `sample_name`,
  and condition. Conditions A and B are required.

HY omits ECOLI from both species sections; HYE includes it. There is no separate preset catalogue.
Unknown upstream module sections are tolerated because the source document is shared with
ProteoBench, but only the validated scoring and sample subset is embedded into APB metadata.
The embedded sample design is normalized as a column mapping (`raw_file`, `sample_name`, and
`condition` arrays), rather than TOML's list of tables, so the same JSON-compatible value also
round-trips through AnnData's HDF5 representation. Authored aliases remain traceable to the
checksum-identified source module but are unnecessary after the match has been applied.
