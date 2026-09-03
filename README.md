# APB ProteoBench

ProteoBench annotation, mixed-species diagnostics, and scoring for storage-neutral APB2
results. HYE and HY use the same configuration-driven calculation.

Start directly from a vendor result table:

```bash
apb-proteobench convert report.tsv \
    --params search-parameters.txt \
    --software spectronaut \
    --output converted
```

This writes `converted.h5mu` with every compatible quantification level. Pass a level such as
`ion` after `report.tsv` to write `converted.h5ad` instead.

Continue from that result—or from an existing APB2 result—with annotation and scoring:

```bash
apb-proteobench annotate converted.h5mu module_settings.toml annotated.h5mu
apb-proteobench score annotated.h5mu scored.h5mu --verbose
```

Annotation validates complete one-to-one sample coverage and embeds the normalized module
configuration. Scoring reads that embedded configuration, writes feature diagnostics to
`varm["proteobench"]`, and writes scores and provenance to
`uns["apb"]["proteobench"]`. See the [documentation](docs/index.md).

The package owns all 11 ProteoBench module TOMLs: the eight quantitative HYE/HY modules used by
legacy APB and the newer plasma, de novo, and entrapment documents. The latter three are explicitly
packaged for planned support but are not accepted by the current quantitative scorer. Stable names,
support status, and the Python loading API are documented under
[module configuration](docs/configuration.md).

## Development

```bash
uv sync --group dev
make check
.venv/bin/pre-commit install --hook-type pre-commit --hook-type pre-push
```

All Python commands run from the synchronized project `.venv`.
