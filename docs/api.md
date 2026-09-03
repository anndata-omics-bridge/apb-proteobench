# Python API

Vendor conversion uses APB2's compiler/parser boundary and returns the parsed in-memory result as
well as persisting it:

```python
from pathlib import Path

from apb_proteobench.api import convert_vendor_result

conversion = convert_vendor_result(
    Path("report.tsv"),
    Path("search-parameters.txt"),
    Path("results/all-levels.h5mu"),
    software="spectronaut",
)
parsed = conversion.parsed
```

Pass `level="ion"` and an `.h5ad` target for one quantification level. Omitting `level` compiles
and parses every compatible level and requires an `.h5mu` target.

The packaged quantitative module catalogue is available without an external ProteoBench checkout.
`available_modules()` reports the eight modules supported by the scorer; `packaged_module_names()`
inventories all 11 upstream documents, including the three retained for planned support:

```python
from apb_proteobench.configuration.load import (
    available_modules,
    load_packaged_module,
    packaged_module_names,
)

print(packaged_module_names())

for name in available_modules():
    module = load_packaged_module(name)
    print(name, module.settings.general.level)
```

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

::: apb_proteobench.configuration.load
