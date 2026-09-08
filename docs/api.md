# Python API

The Python API provides two levels of composition:

- file-to-file functions mirror complete CLI operations; and
- parser and workflow objects expose APB2's storage-neutral `ParsedLevels` between annotation,
  calculation, and persistence.

## Run the complete vendor workflow

`run_vendor_benchmark()` is the one-call API from raw vendor files to final scored MuData:

```python
from pathlib import Path

from apb_proteobench.api import run_vendor_benchmark

result = run_vendor_benchmark(
    Path("report.tsv"),
    Path("search-parameters.txt"),
    (Path("human.fasta"), Path("contaminants.fasta")),
    Path("module_settings.toml"),
    Path("results/scored.h5mu"),
    software="spectronaut",
)

print(result.fasta_reports.peptide_levels)
print(result.analysis.scores.nr_feature)
```

It compiles and runs every compatible APB2 parser, verifies modification-stripped peptide sequences against the FASTA database, applies the ProteoBench design, scores the configured level, and persists once at the end. The target must be an exact `.h5mu` path.

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

scored = benchmark_result(
    Path("results/fasta-checked.h5mu"),
    Path("module_settings.toml"),
    Path("results/scored.h5mu"),
)
print(scored.analysis.scores.nr_feature)
```

## Convert vendor results

### File-to-file facade

`convert_vendor_result()` parses a vendor table through APB2's packaged rules, writes h5ad or h5mu,
and returns the same parsed result in memory:

```python
from pathlib import Path

from apb_proteobench.api import convert_vendor_result

conversion = convert_vendor_result(
    Path("report.tsv"),
    Path("search-parameters.txt"),
    Path("results/all-levels.h5mu"),
    software="spectronaut",
)

print(conversion.software)
print(conversion.software_version)
print(list(conversion.parsed.levels))
```

Omitting `level` converts every compatible level and requires an `.h5mu` target. Select one level
and an `.h5ad` target explicitly:

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

scored = score_result(
    Path("results/annotated.h5mu"),
    Path("results/scored.h5mu"),
)

print(scored.extracted.name)
print(scored.analysis.scores.nr_feature)
print(scored.analysis.scores.median_abs_epsilon_global)
```

The returned `ScoredResult` retains the validated configuration, extracted typed calculation input,
complete diagnostics, aggregate scores, selected methods, and input/output paths. See
[Result layout](results.md) for persisted locations.

### Storage-neutral workflow

Use the lower-level workflow to inspect or transform values before writing:

```python
from pathlib import Path

from apb2.result_facade import read_parsed_levels, write_parsed_levels
from apb_proteobench.integration import embedded_configuration, extract_level, persist_result
from apb_proteobench.workflow import (
    MixedSpeciesDiagnostics,
    ProteoBenchCompatibleScoring,
    analyze_level,
)

parsed = read_parsed_levels(Path("results/annotated.h5mu"))
configuration = embedded_configuration(parsed)
extracted = extract_level(parsed, configuration)

analysis = analyze_level(
    extracted.calculation,
    configuration,
    MixedSpeciesDiagnostics(),
    ProteoBenchCompatibleScoring(),
)

scored = persist_result(parsed, extracted, analysis)
write_parsed_levels(scored, Path("results/scored.duckdb"))
```

The calculation consumes `QuantitativeLevelInput`, not AnnData, MuData, or `ParsedLevels`.
`persist_result()` attaches the calculation output to a copy of the APB2 result.

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
