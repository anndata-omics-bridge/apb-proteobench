"""Typed in-memory values shared by ProteoBench calculations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.sparse import csc_array, csc_matrix, csr_array, csr_matrix

from apb_proteobench.configuration.schema import QuantificationLevel

# Callers may hand over either precision, because a deserialized AnnData layer is often
# float32. Calculations always widen to float64: storage precision is not arithmetic
# precision, and narrowing costs half the mantissa for no benefit.
type FloatArray = NDArray[np.float32] | NDArray[np.float64]
type CompressedSparseMatrix[ScalarT: (np.float32, np.float64)] = (
    csr_matrix[ScalarT] | csc_matrix[ScalarT] | csr_array[ScalarT] | csc_array[ScalarT]
)
type SparseQuantMatrix = CompressedSparseMatrix[np.float32] | CompressedSparseMatrix[np.float64]
type QuantMatrix = FloatArray | SparseQuantMatrix


@dataclass(frozen=True, slots=True)
class QuantitativeLevelInput:
    """Exact level values consumed by the mixed-species calculation."""

    observations: pd.DataFrame
    matrix: QuantMatrix
    feature_ids: pd.Index
    reported_proteins: pd.Series
    level: QuantificationLevel
