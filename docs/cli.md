# CLI reference

`apb-proteobench` exposes three commands:

| Command | Purpose | Guide |
| --- | --- | --- |
| `apb-proteobench convert` | parse a vendor table through APB2 | [End-to-end workflow](workflow.md#1-convert-a-vendor-table) |
| `apb-proteobench annotate` | bind a ProteoBench module to an APB2 result | [End-to-end workflow](workflow.md#2-annotate-the-experiment) |
| `apb-proteobench score` | calculate diagnostics and scores from embedded configuration | [End-to-end workflow](workflow.md#3-calculate-diagnostics-and-scores) |

Use `apb-proteobench --help` or a command's `--help` for the installed version's generated Cyclopts
reference. Direct conversion supports the software, versions, inputs, parameter parsers, and levels
listed in the [APB2 support matrix](https://anndata-omics-bridge.github.io/apb2/supported_software/).

## `apb-proteobench convert`

```text
apb-proteobench convert DATA [LEVEL] [OPTIONS]
```

| Argument or option | Meaning |
| --- | --- |
| `DATA` | vendor quantification table |
| `LEVEL` | optional APB2 level: `ion`, `peptidoform`, `peptide`, `protein`, or `fragment` |
| `--params PATH` | required vendor search-parameter file |
| `--software NAME` | select and verify the packaged vendor rule document |
| `--params-software NAME` | select the parameter parser independently of the result-table rules |
| `--output BASENAME` | output basename without `.h5ad` or `.h5mu` |
| `--strict` | promote APB2 layer-contract warnings to errors |

An explicit level writes `.h5ad`; an omitted level writes every compatible level to `.h5mu`. The
command appends the suffix, so `--output` must be a basename. If `--output` is omitted, the vendor
table's suffix is replaced with `.h5ad` or `.h5mu`.

Convert only the ion level:

```bash
apb-proteobench convert report.tsv ion \
    --params search-parameters.txt \
    --software spectronaut \
    --output results/ion
```

This writes `results/ion.h5ad`.

Convert every compatible level:

```bash
apb-proteobench convert report.tsv \
    --params search-parameters.txt \
    --software spectronaut \
    --output results/all-levels
```

This writes `results/all-levels.h5mu`.

`--software` may be omitted when APB2 can identify exactly one vendor from the table columns. Use
`--params-software` when the parameter file needs a different parser selection from the result
table's packaged rules.

## `apb-proteobench annotate`

```text
apb-proteobench annotate SOURCE MODULE TARGET
```

| Argument | Meaning |
| --- | --- |
| `SOURCE` | existing APB2 h5ad, h5mu, Parquet, or DuckDB result |
| `MODULE` | ProteoBench `module_settings.toml` |
| `TARGET` | new APB2 result; format selected from its suffix |

Example:

```bash
apb-proteobench annotate \
    results/all-levels.h5mu \
    module_settings.toml \
    results/annotated.h5mu
```

The module's `general.level` must exist in the source. Every quantified observation must match one
module sample, and every module sample must be used. Successful annotation adds `raw_file`,
`sample_name`, and `condition` and embeds the normalized module plus its source checksum.

The target must differ from the source and must not exist. Source and target may use any supported
APB2 result formats; using another suffix changes storage format while preserving the logical
result.

## `apb-proteobench score`

```text
apb-proteobench score SOURCE TARGET [--verbose]
```

| Argument or option | Meaning |
| --- | --- |
| `SOURCE` | APB2 result previously processed by `apb-proteobench annotate` |
| `TARGET` | new scored APB2 result; format selected from its suffix |
| `--verbose` | report configuration, coverage, feature, mapping, and score details |

Examples:

```bash
apb-proteobench score results/annotated.h5mu results/scored.h5mu
apb-proteobench score results/annotated.h5mu results/scored-verbose.h5mu --verbose
```

Scoring reads the module configuration embedded by annotation; there is no module argument or
benchmark-name option. It refuses to overwrite the source, an existing target, or existing
ProteoBench diagnostics and scores. See [Result layout](results.md) for the values written.

## Exit behavior

- `0`: operation completed successfully
- `1`: expected input, detection, parsing, annotation, scoring, result-I/O, or writing failure
- `2`: `convert --output` already includes the suffix that the command appends

Cyclopts reports invalid command-line usage. Unexpected programming errors are not broadly
swallowed; they remain visible with their traceback.
