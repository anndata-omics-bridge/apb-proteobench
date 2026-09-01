"""Loguru presentation for completed ProteoBench workflows."""

from __future__ import annotations

from collections import Counter

from loguru import logger

from apb_proteobench.api import ScoredResult


def report_score(result: ScoredResult, /, *, verbose: bool) -> None:
    """Log a concise completion or the stable detailed score summary."""
    if not verbose:
        logger.info(
            "scored level={} diagnostics=varm['proteobench'] scores=uns['apb']['proteobench'] "
            "output={}",
            result.extracted.name,
            result.output_path,
        )
        return
    configuration = result.configuration
    diagnostics = result.analysis.diagnostics
    frame = diagnostics.varm
    scores = result.analysis.scores
    conditions = Counter(sample.condition for sample in configuration.samples)
    accession = diagnostics.protein_mapping.accession_mapper
    included = int(frame["included"].sum())
    logger.info("ProteoBench score summary")
    logger.info("input={} output={}", result.input_path, result.output_path)
    logger.info("level={} layer={}", result.extracted.name, result.extracted.roles.intensity)
    logger.info(
        "species={} expected_A_vs_B={}",
        list(configuration.species_expected_ratio),
        {species: ratio.a_vs_b for species, ratio in configuration.species_expected_ratio.items()},
    )
    logger.info("samples={} conditions={}", len(configuration.samples), dict(conditions))
    logger.info(
        "features total={} included={} excluded={}",
        len(frame),
        included,
        len(frame) - included,
    )
    logger.info(
        "accessions mapper_entries={} matched={} unmatched={}",
        accession.entries,
        accession.matched_token_occurrences,
        accession.unmatched_token_occurrences,
    )
    logger.info(
        "cutoff={} nr_feature={} median_abs_error={} median_abs_precision={} "
        "cv_median={} roc_auc={}",
        configuration.general.default_cutoff_min_feature,
        scores.nr_feature,
        scores.median_abs_epsilon_global,
        scores.median_abs_epsilon_precision_global,
        scores.results[str(configuration.general.default_cutoff_min_feature)].root["CV_median"],
        scores.results[str(configuration.general.default_cutoff_min_feature)].root["roc_auc"],
    )
    logger.info("stored diagnostics=varm['proteobench'] scores=uns['apb']['proteobench']")
