# End-to-end workflow

An APB ProteoBench run has three explicit stages: convert, annotate, and score. Keeping the stages
separate allows an existing APB2 result to enter at annotation and makes each intermediate result
inspectable and reusable.

## What you need

To start from vendor files, provide:

- a vendor quantification table supported by [APB2](https://anndata-omics-bridge.github.io/apb2/supported_software/);
- the corresponding vendor search-parameter file; and
- a ProteoBench `module_settings.toml` describing the experiment and samples.

To start from an existing APB2 result, only the result and module TOML are required. The result may
be h5ad, h5mu, an APB2 Parquet directory, or DuckDB.

## Command-line workflow

### 1. Convert a vendor table

Omit the level to convert every compatible quantification level to MuData:

```bash
apb-proteobench convert report.tsv \
    --params search-parameters.txt \
    --software spectronaut \
    --output results/all-levels
```

This writes `results/all-levels.h5mu`. To convert only the ion level to AnnData:

```bash
apb-proteobench convert report.tsv ion \
    --params search-parameters.txt \
    --software spectronaut \
    --output results/ion
```

This writes `results/ion.h5ad`. APB2 can infer the vendor from unambiguous table columns, so
`--software` is optional when detection has exactly one answer.

### 2. Annotate the experiment

```bash
apb-proteobench annotate \
    results/all-levels.h5mu \
    module_settings.toml \
    results/annotated.h5mu
```

The module's `general.level` selects one level in the result. Annotation requires exact, one-to-one
sample coverage: every quantified observation must match one module sample and every module sample
must be used. The command adds `raw_file`, `sample_name`, and `condition` to that level and embeds
the normalized module configuration and source checksum.

### 3. Calculate diagnostics and scores

```bash
apb-proteobench score \
    results/annotated.h5mu \
    results/scored.h5mu \
    --verbose
```

Scoring reads the configuration stored by annotation. `--verbose` reports the selected level,
species and expected ratios, sample coverage, included features, protein mapping, cutoff, and key
scores. Omit it for one completion line.

## Start from an existing APB2 result

Skip conversion when another process already produced the APB2 result:

```bash
apb-proteobench annotate \
    existing-result.parquet \
    module_settings.toml \
    results/annotated.parquet
apb-proteobench score \
    results/annotated.parquet \
    results/scored.parquet
```

The same workflow works with `.h5ad`, `.h5mu`, `.parquet`, and `.duckdb` paths. APB2 selects the
result reader and writer from the suffix.

## Python workflow

The file-to-file API follows the same stages:

```python
from pathlib import Path

from apb_proteobench.api import annotate_result, convert_vendor_result, score_result

conversion = convert_vendor_result(
    Path("report.tsv"),
    Path("search-parameters.txt"),
    Path("results/all-levels.h5mu"),
    software="spectronaut",
)
print(list(conversion.parsed.levels))

annotation = annotate_result(
    Path("results/all-levels.h5mu"),
    Path("module_settings.toml"),
    Path("results/annotated.h5mu"),
)
print(annotation.reports)

scored = score_result(
    Path("results/annotated.h5mu"),
    Path("results/scored.h5mu"),
)
print(scored.analysis.scores.nr_feature)
```

Pass `level="ion"` to `convert_vendor_result()` and use an `.h5ad` target for a single level.
Omitting `level` requires an `.h5mu` target and converts every compatible level.

For in-memory annotation, packaged module selection, custom calculation methods, and typed return
values, continue to the [Python API reference](api.md).

## Output safety

Annotation and scoring require a target different from the source and refuse an existing target.
Scoring also refuses to replace existing `varm["proteobench"]` diagnostics or
`uns["apb"]["proteobench"]` scores. Choose a new path for each stage.
