# APB ProteoBench

ProteoBench annotation, mixed-species diagnostics, and scoring for storage-neutral APB2
results. HYE and HY use the same configuration-driven calculation.

```bash
apb-proteobench annotate converted.h5mu module_settings.toml annotated.h5mu
apb-proteobench score annotated.h5mu scored.h5mu --verbose
```

Annotation validates complete one-to-one sample coverage and embeds the normalized module
configuration. Scoring reads that embedded configuration, writes feature diagnostics to
`varm["proteobench"]`, and writes scores and provenance to
`uns["apb"]["proteobench"]`. See the [documentation](docs/index.md).

## Development

```bash
uv sync --group dev
make check
.venv/bin/pre-commit install --hook-type pre-commit --hook-type pre-push
```

All Python commands run from the synchronized project `.venv`.
