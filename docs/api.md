# Python API

The Python API provides two levels of composition:

- file-to-file functions mirror complete CLI operations; and
- parser and workflow objects expose APB2's storage-neutral `ParsedLevels` between annotation,
  calculation, and persistence.

## Run the complete vendor workflow

`run_vendor_benchmark()` is the one-call API from raw vendor files to one final scored APB2 result:

```python
from pathlib import Path

from apb_proteobench.api import run_vendor_benchmark
from apb_proteobench.integration import ALL_ABUNDANCE_LAYERS

result = run_vendor_benchmark(
    Path("report.tsv"),
    Path("search-parameters.txt"),
    (Path("human.fasta"), Path("contaminants.fasta")),
    Path("module_settings.toml"),
    Path("results/scored.parquet"),
    software="spectronaut",
    selection=ALL_ABUNDANCE_LAYERS,
)

print(result.fasta_reports.peptide_levels)
for layer_name, layer in result.scored.layers.items():
    print(layer_name, layer.analysis.scores.nr_feature)
```

It constructs APB2's `ParseRuleCompiler`, compiles and parses canonical `ParsedLevels`, passes those levels directly to APB FASTA, applies the ProteoBench design, scores the configured level, and persists once at the end through APB2's `write_parsed_levels`. `ParseRuleCompiler.parameters` supplies ProteoBot parameter fields; parameters are never reconstructed from parsed provenance. Pass `level="ion"` to parse only that quantification level and enable a `.h5ad` target; all-level targets may be `.h5mu`, `.parquet`, or `.duckdb`.

It performs no quantitative aggregation, and APB ProteoBench imports no aggregation package. The scored level must already exist in the vendor result.

When the scored level must be derived from a lower one, aggregate as a separate staged step through the `apb-aggregate` CLI:

```bash
apb-aggregate ion protein sum \
    results/fasta-checked.h5mu results/aggregated.h5mu
```

Reaching aggregation only as a subprocess is deliberate: it keeps APB ProteoBench's declared dependencies to APB2 and APB FASTA.

## Benchmark an existing result

`benchmark_result()` combines annotation and scoring while preserving the APB2 result boundary:

```python
from pathlib import Path

from apb_proteobench.api import benchmark_result
from apb_proteobench.integration import NamedAbundanceLayer

scored = benchmark_result(
    Path("results/fasta-checked.h5mu"),
    Path("module_settings.toml"),
    Path("results/scored.h5mu"),
    selection=NamedAbundanceLayer("LFQ_Intensity"),
)
print(scored.layers["LFQ_Intensity"].analysis.scores.nr_feature)
```

The pMultiQC compatibility writer is deliberately a separate, storage-neutral boundary. Pass it the completed ion-level diagnostics table; it validates the canonical filename and columns, refuses overwrites, and writes atomically without a DataFrame index:

```python
from pathlib import Path

from apb_proteobench.io.result_performance import write_result_performance

layer = scored.layers["LFQ_Intensity"]
write_result_performance(
    layer.analysis.diagnostics.legacy,
    Path("reports/result_performance.csv"),
)
```

The writer performs no APB2 reads or calculations. Callers must select exactly one ion-level `ScoredLayerResult`; the CLI enforces that constraint for `run` and `benchmark`.

## Convert vendor results

### File-to-file facade

`convert_vendor_result()` parses a vendor table through APB2's packaged rules, writes through APB2's suffix-selected persistence boundary, and returns the same parsed result in memory:

```python
from pathlib import Path

from apb_proteobench.api import convert_vendor_result

conversion = convert_vendor_result(
    Path("report.tsv"),
    Path("search-parameters.txt"),
    Path("results/all-levels.parquet"),
    software="spectronaut",
)

print(conversion.software)
print(conversion.software_version)
print(list(conversion.parsed.levels))
```

Omitting `level` converts every compatible level and accepts `.h5mu`, `.parquet`, or `.duckdb`. Select one level explicitly; `.h5ad` is then also valid:

```python
conversion = convert_vendor_result(
    Path("report.tsv"),
    Path("search-parameters.txt"),
    Path("results/ion.h5ad"),
    level="ion",
    software="spectronaut",
    checks="strict",
)
```

`parameters_software` selects the parameter parser independently. When `software` is omitted, APB2
detects the vendor from the result-table columns and parameter evidence.

For direct access to APB2's compiler/parser boundary, rule documents, and parser output, see the
[APB2 Python API](https://anndata-omics-bridge.github.io/apb2/api/#convert-vendor-results).

## Annotate a result

### File-to-file facade

```python
from pathlib import Path

from apb_proteobench.api import annotate_result

annotation = annotate_result(
    Path("results/all-levels.h5mu"),
    Path("module_settings.toml"),
    Path("results/annotated.h5mu"),
)

for level, report in annotation.reports.items():
    print(level, report.coverage)
```

`annotate_result()` reads the source, validates complete module-sample coverage, writes a new result,
and returns APB2's typed `AnnotationResult`.

### Parser and in-memory result

Use `ProteoBenchAnnotationParser` to keep `ParsedLevels` in memory:

```python
from pathlib import Path

from apb2.result_facade import read_parsed_levels, write_parsed_levels
from apb_proteobench.annotation import ProteoBenchAnnotationParser

parsed = read_parsed_levels(Path("results/all-levels.h5mu"))
parser = ProteoBenchAnnotationParser.from_path(Path("module_settings.toml"))
annotation = parser.parse(parsed)

for level, match in annotation.matches.levels.items():
    print(level, match.coverage)

annotated = annotation.annotate().parsed
write_parsed_levels(annotated, Path("results/annotated.parquet"))
```

The parser validates and decodes the module once. `parse(parsed)` binds it to one dataset and raises
before constructing an annotation when the selected level or complete sample coverage is invalid.
`annotate()` applies the stored match without recomputing it.

## Use packaged modules

Eight quantitative HYE/HY modules are loadable without an external ProteoBench checkout:

```python
from pathlib import Path

from apb2.result_facade import read_parsed_levels
from apb_proteobench.annotation import ProteoBenchAnnotationParser
from apb_proteobench.configuration.load import (
    available_modules,
    load_packaged_module,
    packaged_module_names,
)

print(available_modules())
print(packaged_module_names())

module = load_packaged_module("dia_singlecell")
parsed = read_parsed_levels(Path("results/all-levels.h5mu"))
annotation = ProteoBenchAnnotationParser(module).parse(parsed)
annotated = annotation.annotate().parsed
```

`packaged_module_names()` inventories all 11 module TOMLs in the distribution. Plasma, de novo,
and entrapment are packaged for planned support but deliberately rejected by
`load_packaged_module()`. See [Module configuration](configuration.md).

## Score a result

### File-to-file facade

```python
from pathlib import Path

from apb_proteobench.api import score_result
from apb_proteobench.integration import ALL_ABUNDANCE_LAYERS

scored = score_result(
    Path("results/annotated.h5mu"),
    Path("results/scored.h5mu"),
    selection=ALL_ABUNDANCE_LAYERS,
)

for layer_name, layer in scored.layers.items():
    print(layer_name)
    print(layer.analysis.scores.nr_feature)
    print(layer.analysis.scores.median_abs_epsilon_global)
```

The default `ALL_ABUNDANCE_LAYERS` selection preserves the declared abundance-layer order and falls back to primary/X when older input has no abundance metadata, recording that fallback in the result. Pass `PRIMARY_LAYER` to score only `ParsedLevel.primary_layer_name`, which is projected to AnnData `X`; it does not require an `abundance` role. `NamedAbundanceLayer(name)` requires the named layer to exist and carry that role.

The returned `ScoredResult.layers` is an ordered layer-keyed dictionary. Each `ScoredLayerResult` retains resolved roles, diagnostics location, complete diagnostics, and aggregate scores, but not the large calculation matrix. See [Result layout](results.md) for persisted locations.

### Storage-neutral workflow

Use the lower-level workflow to inspect or transform values before writing:

```python
from pathlib import Path

from apb2.result_facade import read_parsed_levels, write_parsed_levels
from apb_proteobench.integration import (
    ALL_ABUNDANCE_LAYERS,
    ScoredLayerResult,
    diagnostics_slot,
    embedded_configuration,
    extract_layer,
    persist_results,
    resolve_layer_selection,
)
from apb_proteobench.workflow import (
    MixedSpeciesDiagnostics,
    ProteoBenchCompatibleScoring,
    analyze_level,
)

parsed = read_parsed_levels(Path("results/annotated.h5mu"))
configuration = embedded_configuration(parsed)
selection = resolve_layer_selection(parsed, configuration, ALL_ABUNDANCE_LAYERS)
layers = {}
for layer_name in selection.layer_names:
    extracted = extract_layer(parsed, configuration, layer_name)
    analysis = analyze_level(
        extracted.calculation,
        configuration,
        MixedSpeciesDiagnostics(),
        ProteoBenchCompatibleScoring(),
    )
    layers[layer_name] = ScoredLayerResult(
        level_name=extracted.level_name,
        layer_name=layer_name,
        diagnostics_slot=diagnostics_slot(layer_name),
        roles=extracted.roles,
        analysis=analysis,
    )

scored = persist_results(parsed, selection, layers)
write_parsed_levels(scored, Path("results/scored.duckdb"))
```

The calculation consumes `QuantitativeLevelInput`, not AnnData, MuData, or `ParsedLevels`.
`extract_layer()` materializes one selected layer at a time, and `persist_results()` attaches all completed layer outputs to a copy of the APB2 result.

## Substitute diagnostic and scoring methods

`score_result()` accepts implementations of the client-owned `DiagnosticMethod` and
`ScoringMethod` protocols:

```python
scored = score_result(
    Path("results/annotated.h5mu"),
    Path("results/custom-scored.h5mu"),
    diagnostic_method=my_diagnostics,
    scoring_method=my_scoring,
)
```

A diagnostic method receives `QuantitativeLevelInput` plus `ModuleSettings` and returns an
`IntermediateResult`. A scoring method reduces that intermediate result to `ProteoBenchScores`.
Both methods provide a persistable `identity()`.

## Errors and output safety

File-to-file annotation and scoring raise `ValueError` when source and target resolve to the same
path, the target exists, required annotation is absent, or ProteoBench output already exists.
Pydantic `ValidationError` reports invalid module documents. APB2 result readers and writers report
format and persistence errors through their own documented exception hierarchy.

## API objects

::: apb_proteobench.api

::: apb_proteobench.annotation

::: apb_proteobench.workflow

::: apb_proteobench.integration

::: apb_proteobench.configuration.schema

::: apb_proteobench.configuration.load
