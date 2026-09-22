# CLI reference

`apb-proteobench` exposes two complete workflows:

| Command | Input | Outputs |
| --- | --- | --- |
| `apb-proteobench run` | vendor table, parameters, FASTA, and module TOML | scored APB2 result and optional pMultiQC/ProteoBot pair |
| `apb-proteobench benchmark` | existing APB2 result and module TOML | scored APB2 result and optional pMultiQC/ProteoBot pair |

In-memory annotation and scoring remain available through `ProteoBenchAnalyzer`, but are not separate CLI commands. APB2 owns conversion and persistence. Use `apb2 convert`, `apb-fasta verify-peptides`, and `apb-aggregate` when a staged shell workflow is needed.

Use `apb-proteobench --help` or a command's `--help` for the installed version's generated Cyclopts reference. Direct conversion supports the software, versions, inputs, parameter parsers, and levels listed in the [APB2 support matrix](https://anndata-omics-bridge.github.io/apb2/supported_software/).

## Outputs

Both commands write the scored APB2 result. Its ProteoBench score and provenance record is stored as JSON-compatible metadata under `uns["apb"]["proteobench"]`.

Both commands also accept `--result-performance PATH`. This option writes two ion-level compatibility artifacts in the same directory: the `result_performance.csv` intermediate consumed by the existing pMultiQC ProteoBench module and the `<intermediate_hash>.json` datapoint stored by ProteoBot in `Proteobench/Results_quant_ion_DDA`.

The supplied path must be named exactly `result_performance.csv`. Neither it nor the derived hash-named JSON may already exist; both files are staged before their final paths are created without overwrite, and the CSV has no DataFrame index. The export requires exactly one selected ion-level layer, so callers must pass `--x` for the APB primary/X layer or `--layer NAME` for one named abundance layer.

## `apb-proteobench run`

```text
apb-proteobench run DATA FASTA... --params PATH --module PATH --output RESULT [OPTIONS]
```

`run` converts vendor inputs, verifies peptides against FASTA, applies the module experiment design, scores the selected layer, and writes one final APB2 result:

```bash
apb-proteobench run report.tsv human.fasta contaminants.fasta \
    --params search-parameters.txt \
    --module module_settings.toml \
    --software spectronaut \
    --level ion \
    --x \
    --output results/scored.h5ad \
    --result-performance reports/result_performance.csv \
    --verbose
```

The main options are `--level LEVEL`, `--software`, `--strict`, `--backend`, `--il-equivalent`, `--protein-group-separator`, `--x`, `--layer NAME`, and `--timings-dir DIR`. `--level` converts one quantification level; omitting it converts every compatible level. Scoring includes every declared abundance layer by default; `--x` restricts scoring to the APB primary/X layer. `--timings-dir` optionally writes separate version-1 timing JSON files for APB2 conversion (`compile`, `read`, `parse`), FASTA verification (`load_database`, `verify_peptides`), and ProteoBench benchmarking (`load_module`, `analyze`, `write`, and `export` when requested). The files are independent of the scientific result and refuse existing targets. At least one FASTA is required. A single-level output may end in `.h5ad`; `.h5mu`, `.parquet`, and `.duckdb` support multiple levels. The target must differ from the vendor table and must not already exist.

`run` performs no quantitative aggregation and scores the level named in `module_settings.toml` as the vendor table reports it. To derive the configured level from a lower one, use the staged route and insert `apb-aggregate` before `benchmark`.

## `apb-proteobench benchmark`

```text
apb-proteobench benchmark SOURCE MODULE TARGET [OPTIONS]
```

`benchmark` reads an APB2 result, applies and embeds the module experiment design, calculates diagnostics and scores, and writes a new APB2 result:

```bash
apb-proteobench benchmark \
    results/fasta-checked.h5mu \
    module_settings.toml \
    results/scored.h5mu \
    --x \
    --result-performance reports/result_performance.csv \
    --verbose
```

With neither layer option, `benchmark` scores every declared abundance layer. `--x` selects only the APB primary/X layer; `--layer NAME` selects one named abundance layer. The APB2 target must differ from the source and must not already exist.

## Run pMultiQC

After either command writes the intermediate:

```bash
multiqc --proteobench-plugin reports -o reports/multiqc
```

The CLI passes the completed in-memory diagnostics table directly to the CSV writer and neither reopens nor recalculates the APB2 result. pMultiQC and MultiQC require no changes.

## Exit behavior

- `0`: operation completed successfully
- `1`: expected input, detection, parsing, annotation, scoring, result-I/O, or writing failure
- `2`: invalid option combination, including `--layer` with `--x` or `--result-performance` without either single-layer selector

Cyclopts reports invalid command-line usage. Unexpected programming errors remain visible with their traceback.
