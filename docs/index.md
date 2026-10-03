# APB ProteoBench

`apb-proteobench` connects [APB2](https://anndata-omics-bridge.github.io/apb2/) results to the
[ProteoBench platform](https://proteobench.cubimed.rub.de/). It converts vendor tables through
APB2, applies a ProteoBench experiment design, calculates mixed-species diagnostics and scores,
and writes a portable APB2 result.

## Choose an interface

| Interface | Start here | Best suited to |
| --- | --- | --- |
| Command line | [CLI workflow](#command-line-interface) or [complete CLI reference](cli.md) | shell use, scripts, and workflow engines |
| Python | [Python workflow](#python-api) or [complete API reference](api.md) | libraries, notebooks, and custom pipelines |

The interfaces expose two primary routes:

```text
direct:  vendor table + parameters + FASTA + module -> checked + scored MuData

staged:  vendor table -> APB2 result -> FASTA check -> [aggregation] -> scored result
         apb2 convert   apb-fasta       apb-aggregate    apb-proteobench
```

Neither route aggregates. APB ProteoBench declares only APB2 and APB FASTA; the optional
aggregation step is reached through the separate `apb-aggregate` CLI, never as a library import.

Use the direct route for one final artifact. Use the staged route when intermediates need to be inspected, cached, or reused, or when the scored level must be derived from a lower one. Existing APB2 results can enter at either FASTA checking or ProteoBench benchmarking. The command line stays at these two complete workflows; the Python API provides in-memory annotation and scoring, while APB2 owns conversion and persistence.

## Command-line interface

Run everything from the raw vendor files:

```bash
apb-proteobench run quant report.tsv proteins.fasta \
    --params search-parameters.txt \
    --module module_settings.toml \
    --software spectronaut \
    --output results/scored.h5mu
```

Or retain each boundary as an artifact:

```bash
apb2 convert report.tsv --params search-parameters.txt --output results/all
apb-fasta verify-peptides results/all.h5mu proteins.fasta --output results/checked.h5mu
apb-aggregate ion protein sum results/checked.h5mu results/aggregated.h5mu
apb-proteobench benchmark results/aggregated.h5mu module_settings.toml results/scored.h5mu
```

The aggregation call is conditional and belongs to `apb-aggregate`, not to this package: include it only when the level named in `module_settings.toml` is not the level the vendor table reports, chaining one call per source level. Otherwise `benchmark` reads `results/checked.h5mu` directly.

The [end-to-end guide](workflow.md) explains both routes; the [CLI reference](cli.md) lists every
argument and option.

For the existing pMultiQC ProteoBench module, add `--result-performance reports/result_performance.csv` to `run` or `benchmark`. The compatibility export accepts one ion-level layer and publishes both the pMultiQC CSV and the sibling `<intermediate_hash>.json` ProteoBot datapoint while leaving pMultiQC and MultiQC unchanged.

## Python API

The Python API analyzes canonical APB2 results in memory; callers compose conversion and persistence through APB2:

```python
from pathlib import Path

from apb2.api import read_parsed_levels, write_parsed_levels
from apb_proteobench.api import ProteoBenchAnalyzer
from apb_proteobench.configuration.load import load_module

parsed = read_parsed_levels(Path("results/fasta-checked.h5mu"))
analyzer = ProteoBenchAnalyzer(load_module(Path("module_settings.toml")))
result = analyzer.analyze(parsed)
write_parsed_levels(result.parsed, Path("results/scored.h5mu"))
print(result.layers["Intensity"].analysis.scores.nr_feature)
```

The analyzer owns annotation, layer selection, diagnostics, scoring, and attachment of results. It owns no paths, vendor conversion, FASTA loading, or persistence. See the [complete Python API](api.md).

## Inputs and outputs

Direct vendor conversion uses the same packaged rules and parameter parsers as APB2. Consult the
[APB2 support matrix](https://anndata-omics-bridge.github.io/apb2/supported_software/) for supported
software, versions, file types, and quantification levels.

| Operation | Accepted input | Output |
| --- | --- | --- |
| `run` | vendor table, parameters, one or more FASTAs, and module TOML | FASTA-checked, scored H5MU, Parquet, or DuckDB result |
| `benchmark` | APB2 result plus module TOML | annotated and scored APB2 result |

Both commands can additionally emit the canonical `result_performance.csv` pMultiQC input and matching hash-named ProteoBot JSON when requested. Their APB2 result also retains the ProteoBench score and provenance record as JSON-compatible metadata under `uns["apb"]["proteobench"]`.

The annotation stage checks that the module describes every observation exactly once. It adds `raw_file`, `sample_name`, and `condition` to the configured level. HYE and HY are module configurations consumed by the same calculation rather than separate hard-coded modes.

Scoring uses one layer listed under the APB `abundance` role: `--layer NAME`, by default `X`, the APB primary layer represented by AnnData `X`. The Python API can still score every abundance layer. Per-layer feature diagnostics live in `varm["proteobench:<layer-name>"]`; selection provenance, aggregate scores, method identities, role resolution, compatibility versions, and species-mapping provenance live in `uns["apb"]["proteobench"]`.

The [module guide](configuration.md) lists all packaged module documents and their current support
status. The [result-layout guide](results.md) documents the stored annotation, diagnostics, scores,
and round-trip guarantees.
