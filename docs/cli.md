# CLI

## Annotate

```bash
apb-proteobench annotate SOURCE MODULE TARGET
```

`SOURCE` and `TARGET` may be APB2 H5AD, H5MU, Parquet-directory, or DuckDB results. `MODULE` is
the existing ProteoBench module TOML. The configured `general.level` must exist, every quantified
observation must match one sample, and every module sample must be used. The command writes
`raw_file`, `sample_name`, and `condition` into that level's `obs` and embeds the resolved module.

## Score

```bash
apb-proteobench score ANNOTATED TARGET
apb-proteobench score ANNOTATED TARGET --verbose
```

Scoring reads the module embedded by `annotate`. It refuses to overwrite the input, an existing
target, or existing ProteoBench diagnostics/scores. `--verbose` reports paths, level, species,
ratios, sample coverage, feature inclusion, mapper counts, cutoff, principal scores, and storage
locations. Normal output is one completion line.
