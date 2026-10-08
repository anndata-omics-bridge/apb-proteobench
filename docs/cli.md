# CLI reference

`apb-proteobench` exposes three complete workflows:

| Command | Input | Outputs |
| --- | --- | --- |
| `apb-proteobench run quant` | vendor table, parameters, FASTA, and packaged module or module TOML | scored APB2 result and optional pMultiQC/ProteoBot pair |
| `apb-proteobench run entrapment` | vendor table, parameters, and the entrapment FASTA or its protein database | APB2 result scored for every precursor q-value kind |
| `apb-proteobench benchmark` | existing APB2 result and packaged module or module TOML | scored APB2 result and optional pMultiQC/ProteoBot pair |

In-memory annotation and scoring remain available through `ProteoBenchAnalyzer` and `EntrapmentAnalyzer`, but are not separate CLI commands. APB2 owns conversion and persistence. Use `apb2 convert`, `apb-fasta verify-peptides`, and `apb-aggregate` when a staged shell workflow is needed.

Use `apb-proteobench --help` or a command's `--help` for the installed version's generated Cyclopts reference. Direct conversion supports the software, versions, inputs, parameter parsers, and levels listed in the [APB2 support matrix](https://anndata-omics-bridge.github.io/apb2/supported_software/).

## Outputs

Both commands write the scored APB2 result. Its ProteoBench score and provenance record is stored as JSON-compatible metadata under `uns["apb"]["proteobench"]`.

Both commands also accept `--result-performance PATH`. This option writes two ion-level compatibility artifacts in the same directory: the `result_performance.csv` intermediate consumed by the existing pMultiQC ProteoBench module and the `<intermediate_hash>.json` datapoint stored by ProteoBot in `Proteobench/Results_quant_ion_DDA`.

`--scores PATH` writes that datapoint alone, at any level: the search-parameter fields ProteoBench records plus the scores, including `results` per cutoff, in the layout of the upstream datapoints, so APB scores and ProteoBench's can be compared field by field. It requires exactly one scored layer.

The supplied path must be named exactly `result_performance.csv`. Neither it nor the derived hash-named JSON may already exist; both files are staged before their final paths are created without overwrite, and the CSV has no DataFrame index. The export writes the one scored layer, which must be at the ion level.

## `apb-proteobench run quant`

```text
apb-proteobench run quant DATA FASTA... --params PATH --module NAME|PATH --output RESULT [OPTIONS]
```

`run quant` converts vendor inputs, verifies peptides against FASTA, applies the module experiment design, scores the selected layer, and writes one final APB2 result:

```bash
apb-proteobench run quant report.tsv human.fasta contaminants.fasta \
    --params search-parameters.txt \
    --module module_settings.toml \
    --software spectronaut \
    --level ion \
    --output results/scored.h5ad \
    --result-performance reports/result_performance.csv \
    --verbose
```

The main options are `--level LEVEL`, `--software`, `--strict`, `--backend`, `--protein-group-separator`, `--layer NAME`, and `--timings-dir DIR`. `--level` converts one quantification level; omitting it converts every compatible level. `--layer NAME` scores one abundance layer; the default `X` is the APB primary layer stored in AnnData `X`. `--timings-dir` optionally writes separate version-1 timing JSON files for APB2 conversion (`compile`, `read`, `parse`), FASTA verification (`load_database`, `verify_peptides`), and ProteoBench benchmarking (`load_module`, `analyze`, `write`, and `export` when requested). The files are independent of the scientific result and refuse existing targets. At least one FASTA is required. A single-level output may end in `.h5ad`; `.h5mu`, `.parquet`, and `.duckdb` support multiple levels. The target must differ from the vendor table and must not already exist.

`run quant` performs no quantitative aggregation and scores the level named in `module_settings.toml` as the vendor table reports it. To derive the configured level from a lower one, use the staged route and insert `apb-aggregate` before `benchmark`.

## `apb-proteobench run entrapment`

```text
apb-proteobench run entrapment DATA FASTA --params PATH --output RESULT [OPTIONS]
```

`run entrapment` converts vendor inputs at the module's level, verifies peptides against FASTA, labels each precursor target or entrapment from the same FASTA, and scores ProteoBench's entrapment metrics: lower-bound, combined and paired FDP, their categories, and the FDP curve. apb-catalog's `proteobench_entrapment` set names the precursor q-values the result offers (`q_value`, `library_q_value`, `global_q_value`); each kind is scored from its best value across runs:

```bash
apb-proteobench run entrapment report.parquet ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment_pep.fasta \
    --params report.log.txt \
    --output results/entrapment.h5ad
```

The FASTA argument may instead be the Parquet file `protein-fasta database entrapment.parquet ProteoBenchFASTA_Entrapment_….fasta` writes once; it loads about 70 times faster than parsing the 2.84 M-entry FASTA. Labels and pairs come from that FASTA rather than ProteoBench's pair file, and reproduce the pair file exactly; see [differences from ProteoBench](compatibility.md).

`--scores PATH` writes one ProteoBench entrapment datapoint per q-value kind, keyed by kind: the search-parameter fields plus `lower_bound_FDP`, `combined_FDP`, `paired_FDP`, their categories and the FDP curve, in the layout of the upstream entrapment datapoints.

`--module` defaults to the packaged `entrapment_dia_astral`. Scores per kind sit in the level's `metadata["proteobench"]["result"]["entrapment"]`, with identified precursors and combined and paired FDP per kind in its `summary`; each precursor's label, pair and best q-values sit in `varm["proteobench:entrapment"]`, and the catalogue snapshot in the root `metadata["catalog"]["proteobench_entrapment"]["result"]`. Precursors with equal q-values share a rank, so a tie never counts as an entrapment out-scoring its target ([ProteoBench#1159](https://github.com/Proteobench/ProteoBench/issues/1159)). `--timings-dir` writes the same three timing files as `run quant`, without `export`.

## `apb-proteobench benchmark`

```text
apb-proteobench benchmark SOURCE MODULE TARGET [OPTIONS]
```

`benchmark` reads an APB2 result, applies and embeds the module experiment design, calculates diagnostics and scores, and writes a new APB2 result:

```bash
apb-proteobench benchmark \
    results/fasta-checked.h5mu \
    dda_qexactive \
    results/scored.h5mu \
    --result-performance reports/result_performance.csv \
    --verbose
```

The module argument of `benchmark` and `--module` of `run` name a packaged module, such as `dda_qexactive`, or a module TOML path ending in `.toml`, whose SDRF is resolved beside it. Both commands score one abundance layer: `--layer NAME`, by default `X`, the APB primary layer. The APB2 target must differ from the source and must not already exist.

## Run pMultiQC

After either command writes the intermediate:

```bash
multiqc --proteobench-plugin reports -o reports/multiqc
```

The CLI passes the completed in-memory diagnostics table directly to the CSV writer and neither reopens nor recalculates the APB2 result. pMultiQC and MultiQC require no changes.

## Exit behavior

- `0`: operation completed successfully
- `1`: expected input, detection, parsing, annotation, scoring, result-I/O, or writing failure

Cyclopts reports invalid command-line usage. Unexpected programming errors remain visible with their traceback.
