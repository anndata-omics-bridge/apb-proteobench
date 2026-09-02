# CLI

## Convert

Convert one quantification level from a vendor table to AnnData:

```bash
apb-proteobench convert report.tsv ion \
    --params search-parameters.txt \
    --software spectronaut \
    --output results/ion
```

This writes `results/ion.h5ad`. Omit the level to parse every compatible level into MuData:

```bash
apb-proteobench convert report.tsv \
    --params search-parameters.txt \
    --software spectronaut \
    --output results/all-levels
```

This writes `results/all-levels.h5mu`.

| Argument or option | Meaning |
| --- | --- |
| `DATA` | Vendor result table |
| `LEVEL` | Optional APB2 quantification level |
| `--params PATH` | Required vendor search-parameter file |
| `--software NAME` | Select and verify the packaged vendor rules |
| `--params-software NAME` | Select the parameter parser independently |
| `--output BASENAME` | Output basename without `.h5ad` or `.h5mu` |
| `--strict` | Promote layer-contract warnings to errors |

`--software` may be omitted when APB2 can identify one vendor from the table columns.

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
