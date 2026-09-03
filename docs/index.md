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

Both interfaces expose the same three explicit operations:

```text
vendor table + search parameters
    -> convert
    -> APB2 result

APB2 result + ProteoBench module TOML
    -> annotate
    -> annotated APB2 result

annotated APB2 result
    -> score
    -> scored APB2 result
```

Conversion is optional when an APB2 result already exists. Annotation stores the validated module
with the result, so scoring needs neither a preset name nor the module file a second time.

## Command-line interface

Start directly from a vendor table:

```bash
apb-proteobench convert report.tsv \
    --params search-parameters.txt \
    --software spectronaut \
    --output results/converted
apb-proteobench annotate \
    results/converted.h5mu \
    module_settings.toml \
    results/annotated.h5mu
apb-proteobench score \
    results/annotated.h5mu \
    results/scored.h5mu \
    --verbose
```

With no level, `convert` writes every compatible level to `results/converted.h5mu`. Pass a level
after the vendor table—for example, `report.tsv ion`—to write one h5ad instead. The
[end-to-end guide](workflow.md) covers both starting points; the [CLI reference](cli.md) lists every
argument and option.

## Python API

File-to-file functions mirror the complete CLI operations and return typed results for inspection:

```python
from pathlib import Path

from apb_proteobench.api import annotate_result, convert_vendor_result, score_result

conversion = convert_vendor_result(
    Path("report.tsv"),
    Path("search-parameters.txt"),
    Path("results/converted.h5mu"),
    software="spectronaut",
)
print(conversion.software, list(conversion.parsed.levels))

annotate_result(
    Path("results/converted.h5mu"),
    Path("module_settings.toml"),
    Path("results/annotated.h5mu"),
)
scored = score_result(
    Path("results/annotated.h5mu"),
    Path("results/scored.h5mu"),
)
print(scored.analysis.scores.nr_feature)
```

The API also exposes the annotation parser, storage-neutral calculation workflow, and replaceable
diagnostic and scoring protocols. See the [complete Python API](api.md).

## Inputs and outputs

Direct vendor conversion uses the same packaged rules and parameter parsers as APB2. Consult the
[APB2 support matrix](https://anndata-omics-bridge.github.io/apb2/supported_software/) for supported
software, versions, file types, and quantification levels.

| Operation | Accepted input | Output |
| --- | --- | --- |
| `convert` | supported vendor table plus search-parameter file | one-level `.h5ad` or all-level `.h5mu` |
| `annotate` | APB2 h5ad, h5mu, Parquet, or DuckDB result plus module TOML | APB2 h5ad, h5mu, Parquet, or DuckDB result |
| `score` | annotated APB2 h5ad, h5mu, Parquet, or DuckDB result | scored APB2 h5ad, h5mu, Parquet, or DuckDB result |

The annotation command checks that the module describes every observation exactly once. It
adds `raw_file`, `sample_name`, and `condition` to the configured level. HYE and HY are module
configurations consumed by the same calculation rather than separate hard-coded modes.

The resulting feature table lives in `varm["proteobench"]`. Aggregate scores, method identities,
role resolution, compatibility versions, and mapper provenance live in
`uns["apb"]["proteobench"]`.

The [module guide](configuration.md) lists all packaged module documents and their current support
status. The [result-layout guide](results.md) documents the stored annotation, diagnostics, scores,
and round-trip guarantees.
