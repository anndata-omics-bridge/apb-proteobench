# Python API

The public API analyzes canonical APB2 values in memory. It does not open vendor files, load FASTA files, choose output formats, write results, log progress, or reproduce CLI workflows.

## Analyze an APB2 result

Construct one complete analyzer from a validated ProteoBench module, then pass it `ParsedLevels`:

```python
from pathlib import Path

from apb2.api import read_parsed_levels, write_parsed_levels
from apb_proteobench.api import ProteoBenchAnalyzer, load_module

parsed = read_parsed_levels(Path("results/fasta-checked.h5mu"))
module = load_module(Path("module_settings.toml"))
analyzer = ProteoBenchAnalyzer(module)

result = analyzer.analyze(parsed)

print(list(result.layers))
print(result.layers["Intensity"].analysis.scores.nr_feature)
write_parsed_levels(result.parsed, Path("results/scored.h5mu"))
```

`load_packaged_module(name)` loads a packaged module instead; see [Configuration](configuration.md).

`analyze()` validates complete sample coverage, attaches the normalized module annotation, resolves the requested quantitative layers, calculates diagnostics and scores, and returns a new storage-neutral APB2 result. It does not mutate the input or perform physical I/O.

`ProteoBenchAnalysisResult` contains only analysis-owned values:

- `parsed`: annotated and scored `ParsedLevels`
- `configuration`: validated module settings
- `layers`: ordered results by logical layer name

Input paths, output paths, vendor detection, search parameters, FASTA reports, and persistence outcomes belong to the calling application and are not part of this result.

## Select layers

Every declared abundance layer is scored by default; `layers` names the ones to score instead:

```python
from apb_proteobench.api import ProteoBenchAnalyzer

every = ProteoBenchAnalyzer(module).analyze(parsed)
primary = ProteoBenchAnalyzer(
    module, layers=(parsed.levels["ion"].primary_layer_name,)
).analyze(parsed)
named = ProteoBenchAnalyzer(module, layers=("LFQ_Intensity",)).analyze(parsed)
```

APB2's `ParsedLevel.abundance_layers()` resolves `layers`: `None` keeps declared order and requires at least one layer with the abundance role; a named layer must exist and carry that role. Scoring provenance records the scored layer names under `layers`.

## Substitute calculation methods

The analyzer accepts implementations of the client-owned `DiagnosticMethod` and `ScoringMethod` protocols at construction:

```python
analyzer = ProteoBenchAnalyzer(
    module,
    diagnostic_method=my_diagnostics,
    scoring_method=my_scoring,
)
result = analyzer.analyze(parsed)
```

The selected objects perform the calculation directly. The analyzer does not inspect method names or branch on implementation types.

## Compose vendor conversion explicitly

Applications starting from vendor files compose the owning packages at their boundary:

```python
from pathlib import Path

from apb2.api import ParseRuleCompiler, write_parsed_levels
from apb_fasta.api import FastaAnnotator
from apb_proteobench.api import ProteoBenchAnalyzer, load_module

compiler = ParseRuleCompiler(
    Path("report.tsv"),
    Path("search-parameters.txt"),
    requested_levels=("ion",),
    software="spectronaut",
)
parsed = compiler.compile().parse()

verified = FastaAnnotator.read((Path("proteins.fasta"),)).verify_peptides(parsed)

module = load_module(Path("module_settings.toml"))
result = ProteoBenchAnalyzer(module).analyze(verified.parsed)
write_parsed_levels(result.parsed, Path("results/scored.h5ad"))
```

This is the same ownership sequence used by the CLI: APB2 converts, APB FASTA verifies, APB ProteoBench analyzes, and APB2 persists. No APB ProteoBench wrapper duplicates APB2 conversion or storage.

## Score entrapment

`EntrapmentAnalyzer` scores ProteoBench's entrapment module. It takes the packaged module and the entrapment FASTA as apb-fasta reads it, from which it derives the target/entrapment pairs:

```python
from pathlib import Path

from apb_fasta.api import FastaAnnotator
from apb_proteobench.api import EntrapmentAnalyzer, load_packaged_entrapment_module

fasta = FastaAnnotator.read((Path("entrapment.fasta"),))
module = load_packaged_entrapment_module("entrapment_dia_astral")
result = EntrapmentAnalyzer(module, fasta).analyze(fasta.verify_peptides(parsed).parsed)
print(result.scores["q_value"].paired_FDP)
```

`analyze()` labels every precursor, scores each q-value kind the result offers, and returns a new APB2 result with the labels in `varm`. It does not perform physical I/O.
