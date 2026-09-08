# End-to-end workflow

APB ProteoBench supports a direct one-call route and a reusable staged route. Both perform the same ordered work: APB2 conversion, FASTA peptide verification, ProteoBench annotation, and ProteoBench scoring.

Quantitative aggregation is not part of either route. APB ProteoBench depends only on APB2 and APB FASTA; when the scored level must be derived from a lower one, insert the separate `apb-aggregate` command as a staged step.

## What you need

To start from vendor files, provide:

- a vendor quantification table supported by [APB2](https://anndata-omics-bridge.github.io/apb2/supported_software/);
- the corresponding vendor search-parameter file;
- one or more protein FASTA files; and
- a ProteoBench `module_settings.toml` describing the experiment and samples.

## Direct one-call workflow

```bash
apb-proteobench run report.tsv proteins.fasta \
    --params search-parameters.txt \
    --module module_settings.toml \
    --software spectronaut \
    --output results/scored.h5mu \
    --verbose
```

The direct command compiles every compatible APB2 level, verifies modification-stripped peptides against the FASTA database, applies the sample design, calculates diagnostics and scores, and writes only the final `.h5mu`. Intermediate APB2 values stay in memory. It scores the level named by `module_settings.toml` as the vendor table reports it; it derives no new level. Use the staged route when aggregation is required.

## Staged workflow

Use the independent tools when each boundary should be cached or inspected:

```bash
apb2 convert report.tsv \
    --params search-parameters.txt \
    --software spectronaut \
    --output results/all-levels
apb-fasta verify-peptides \
    results/all-levels.h5mu proteins.fasta \
    --output results/fasta-checked.h5mu
apb-aggregate ion protein sum \
    results/fasta-checked.h5mu results/aggregated.h5mu
apb-proteobench benchmark \
    results/aggregated.h5mu module_settings.toml results/scored.h5mu
```

Aggregation is reached only through the `apb-aggregate` CLI, never as a library call. That is what keeps APB ProteoBench free of a dependency on it.

Include the step only when `module_settings.toml` names a level the vendor table does not report. Run it once per source level, chaining outputs when both `ion` and `fragment` must reach `protein`. Omit it and `benchmark` reads `fasta-checked.h5mu` directly.

## Start from an existing APB2 result

Skip conversion when another process already produced the APB2 result:

```bash
apb-fasta verify-peptides \
    existing-result.parquet proteins.fasta \
    --output results/fasta-checked.parquet
apb-aggregate ion protein sum \
    results/fasta-checked.parquet results/aggregated.parquet
apb-proteobench benchmark \
    results/aggregated.parquet module_settings.toml results/scored.parquet
```

The staged workflow accepts `.h5ad`, `.h5mu`, `.parquet`, and `.duckdb` paths. APB2 selects the
result reader and writer from the suffix.

## Fine-grained ProteoBench stages

Annotation and scoring remain separately callable:

```bash
apb-proteobench annotate \
    existing-result.h5mu module_settings.toml results/annotated.h5mu
apb-proteobench score \
    results/annotated.h5mu results/scored.h5mu --verbose
```

Annotation requires exact, one-to-one sample coverage and embeds the normalized module
configuration and source checksum. `score` reads that embedded configuration.

## Python workflow

The direct Python API follows the same composition:

```python
from pathlib import Path

from apb_proteobench.api import run_vendor_benchmark

result = run_vendor_benchmark(
    Path("report.tsv"),
    Path("search-parameters.txt"),
    (Path("proteins.fasta"),),
    Path("module_settings.toml"),
    Path("results/scored.h5mu"),
    software="spectronaut",
)

print(list(result.parsed.levels))
print(result.fasta_reports.peptide_levels)
print(result.scored.analysis.scores.nr_feature)
```

For separate file-to-file stages, in-memory parsers, packaged module selection, custom calculation
methods, and typed return values, continue to the [Python API reference](api.md).

## Output safety

All result-to-result operations require a target different from the source and refuse an existing
target. The direct vendor workflow requires an exact `.h5mu` target and also refuses to overwrite
its input. Choose a new path for each persisted stage.
