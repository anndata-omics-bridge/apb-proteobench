# APB ProteoBench

ProteoBench annotation, mixed-species diagnostics, and scoring for storage-neutral APB2
results. HYE and HY use the same configuration-driven calculation.

## Installation

APB ProteoBench requires Python 3.13 or later.

```bash
pip install apb-proteobench
```

## Usage

Run the complete workflow directly from vendor files:

```bash
apb-proteobench run quant report.tsv proteins.fasta \
    --params search-parameters.txt \
    --module module_settings.toml \
    --software spectronaut \
    --level ion \
    --output results/scored.h5ad
```

This one call converts the selected APB2 quantification level, verifies modification-stripped peptide sequences against the FASTA database, applies the ProteoBench sample design, calculates diagnostics and scores, and writes one final APB2 result. Omit `--level` to convert every compatible level. Single-level H5AD and multi-level H5MU targets are supported alongside Parquet and DuckDB. Intermediate values remain storage-neutral `ParsedLevels`, and APB2's `write_parsed_levels` selects persistence from the target suffix. It scores the level named in `module_settings.toml` as the vendor table reports it and performs no quantitative aggregation.

The equivalent reusable workflow calls the existing tools directly:

```bash
apb2 convert report.tsv --params search-parameters.txt --software spectronaut --output results/all
apb-fasta verify-peptides results/all.h5mu proteins.fasta --output results/fasta-checked.h5mu
apb-aggregate ion protein sum results/fasta-checked.h5mu results/aggregated.h5mu
apb-proteobench benchmark results/aggregated.h5mu module_settings.toml results/scored.h5mu
```

The `apb-aggregate` step is optional: include it only when the scored level must be derived from a lower one. Omit it and `benchmark` reads `results/fasta-checked.h5mu` directly.

APB ProteoBench declares only `apb2` and `apb-fasta`, and reaches aggregation solely as a subprocess.

The CLI intentionally exposes only `run` and `benchmark`; APB2 owns conversion and persistence, while the Python API exposes in-memory ProteoBench analysis over canonical `ParsedLevels`. `benchmark` combines annotation and scoring for an existing APB2 result. Scoring uses one abundance layer: `--layer NAME`, by default `X`, the APB primary layer projected to AnnData `X`. Per-layer diagnostics live in `varm["proteobench:<layer-name>"]`, while selection provenance and layer-keyed scores live in `metadata["proteobench"]`. See the [documentation](https://anndata-omics-bridge.github.io/apb-proteobench/).

To feed the existing pMultiQC ProteoBench module and retain the matching ProteoBot result, add `--result-performance reports/result_performance.csv` to `run` or `benchmark`. The option publishes two staged files in the same directory: `result_performance.csv` plus `<intermediate_hash>.json` in the schema and naming convention used by `Proteobench/Results_quant_ion_DDA`. Each file is published atomically without overwrite, and a failed publication rolls back any file added by the same call. This export writes the scored layer, supports only an ion level, and refuses either existing target. Its SHA-256 submission hash is computed only during this export from the selected layer and scoring inputs. The new identity does not match historical SHA-1 uploads.

For optional operation-level timing files, pass `--timings-dir DIR` to `run`. It writes separate JSON files for APB2 conversion, FASTA verification, and ProteoBench benchmarking without changing the scored result or Studio's process-level runtime measurement.

The package owns all 11 ProteoBench module TOMLs: the nine quantitative HYE/HY and plasma modules
and the newer de novo and entrapment documents. Entrapment is scored by `run entrapment`; de novo
is packaged for planned support. Neither is accepted by the quantitative scorer. Every quantitative module's scores
include ProteoBench's plasma metrics; see [results](https://anndata-omics-bridge.github.io/apb-proteobench/results/). Stable names,
support status, and the Python loading API are documented under
[module configuration](https://anndata-omics-bridge.github.io/apb-proteobench/configuration/).

## Development

```bash
uv sync --group dev
make check
.venv/bin/pre-commit install --hook-type pre-commit --hook-type pre-push
```

All Python commands run from the synchronized project `.venv`.
