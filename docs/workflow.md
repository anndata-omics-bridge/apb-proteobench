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
apb-proteobench run quant report.tsv proteins.fasta \
    --params search-parameters.txt \
    --module module_settings.toml \
    --software spectronaut \
    --level ion \
    --output results/scored.h5ad \
    --result-performance reports/result_performance.csv \
    --verbose
```

The direct command compiles the requested APB2 quantification level, verifies modification-stripped peptides against the FASTA database, applies the sample design, calculates diagnostics and scores, and writes the final APB2 result plus the optional pMultiQC/ProteoBot export pair. Omit `--level` to compile every compatible level. Intermediate APB2 values stay in memory as storage-neutral `ParsedLevels`; `write_parsed_levels` selects H5AD, H5MU, Parquet, or DuckDB from the target suffix. It scores the level named by `module_settings.toml` as the vendor table reports it; it derives no new level. Use the staged route when aggregation is required.

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

## Export for pMultiQC

Both `run` and `benchmark` can write the current ion-level ProteoBench intermediate under the exact filename consumed by the existing pMultiQC ProteoBench module. The same option also writes the sibling `<intermediate_hash>.json` datapoint used by the ProteoBot result repository:

```bash
apb-proteobench benchmark \
    results/fasta-checked.h5mu module_settings.toml results/scored.h5mu \
    --result-performance reports/result_performance.csv
multiqc --proteobench-plugin reports -o reports/multiqc
```

The export writes the one scored layer, `--layer NAME` or the default `X`; it rejects non-ion levels, wrong CSV filenames, and either existing output. Both files are staged before publication; each final path is created atomically, and a failed publication rolls back files added by the same call. The CLI reads the selected in-memory APB2 layer once more to compute the versioned SHA-256 submission hash; scoring itself does not compute a hash. The CSV uses `index=False`, while the JSON uses ProteoBot's hash filename and top-level datapoint fields. pMultiQC and MultiQC require no changes.

## Python workflow

The direct Python API keeps package ownership explicit:

```python
from pathlib import Path

from apb2.api import ParseRuleCompiler, write_parsed_levels
from apb_fasta.api import FastaAnnotator
from apb_proteobench.api import ProteoBenchAnalyzer, load_module

compiler = ParseRuleCompiler(
    Path("report.tsv"),
    Path("search-parameters.txt"),
    software="spectronaut",
)
parsed = compiler.compile().parse()
verified = FastaAnnotator.read((Path("proteins.fasta"),)).verify_peptides(parsed)
result = ProteoBenchAnalyzer(load_module(Path("module_settings.toml"))).analyze(verified.parsed)
write_parsed_levels(result.parsed, Path("results/scored.h5mu"))

print(list(result.parsed.levels))
print(verified.reports.peptide_levels)
for layer_name, layer in result.layers.items():
    print(layer_name, layer.analysis.scores.nr_feature)
```

For packaged module selection, custom calculation methods, and typed return values, continue to the [Python API reference](api.md).

## Output safety

All result-to-result operations require a target different from the source and refuse an existing
target. The direct vendor workflow requires an exact `.h5mu`, `.parquet`, or `.duckdb` target and also refuses to overwrite its input. Choose a new path for each persisted stage.
