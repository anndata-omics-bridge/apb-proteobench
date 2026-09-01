# Python API

The result-level convenience operations are:

```python
from pathlib import Path

from apb_proteobench.api import annotate_result, score_result

annotate_result(
    Path("converted.h5mu"),
    Path("module_settings.toml"),
    Path("annotated.h5mu"),
)
result = score_result(Path("annotated.h5mu"), Path("scored.h5mu"))
```

`score_result` accepts explicit `diagnostic_method` and `scoring_method` implementations for
programmatic composition. The storage-independent calculation entry point is
`apb_proteobench.workflow.analyze_level`.

::: apb_proteobench.api

::: apb_proteobench.workflow

::: apb_proteobench.configuration.schema
