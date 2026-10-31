from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from cc_perturbseq.manuscript_feature_umaps import (
    FIGURE2E_SPECS,
    FIGURE3B_SPECS,
    FIGURE4D_SPECS,
    SOURCE_CMAP_COLORS,
    SOURCE_EXECUTED_P99,
    SOURCE_POINT_SIZE,
    SOURCE_UPPER_PERCENTILE,
    source_program_vector,
    upper_clip_only,
)


class SourceFeatureUmapTests(unittest.TestCase):
    def test_exact_recovered_style_constants(self):
        self.assertEqual(SOURCE_CMAP_COLORS, ("#D3D3D3", "#4D4D4D"))
        self.assertEqual(SOURCE_UPPER_PERCENTILE, 99.0)
        self.assertEqual(SOURCE_POINT_SIZE, 0.30)

    def test_programs_are_read_by_zero_based_position_not_column_label(self):
        usage = pd.DataFrame(
            {
                "P1": [10.0, 11.0],
                "P2": [20.0, 21.0],
                "P3": [30.0, 31.0],
            }
        )
        vector, source_column = source_program_vector(usage, 1)
        np.testing.assert_array_equal(vector, [20.0, 21.0])
        self.assertEqual(source_column, "P2")

    def test_upper_tail_only_clipping_preserves_the_lower_end(self):
        values = np.asarray([-2.0, 0.0, 1.0, 2.0, 100.0])
        plotted, lower, upper = upper_clip_only(values, 80.0)
        self.assertEqual(lower, -2.0)
        self.assertAlmostEqual(upper, float(np.percentile(values, 80.0)))
        self.assertEqual(plotted[0], -2.0)
        self.assertEqual(plotted[-1], upper)

    def test_panels_use_the_exact_source_program_positions_and_markers(self):
        self.assertEqual(
            [(kind, feature) for kind, feature, _ in FIGURE2E_SPECS],
            [("program", 17), ("gene", "SLC8A1"), ("program", 27), ("gene", "GATA4")],
        )
        self.assertEqual(
            [(kind, feature) for kind, feature, _ in FIGURE3B_SPECS],
            [("program", 2), ("gene", "CDH5"), ("program", 34), ("gene", "ACTC1")],
        )
        self.assertEqual(
            [(kind, feature) for kind, feature, _ in FIGURE4D_SPECS],
            [("program", 15), ("gene", "FLI1"), ("program", 17), ("gene", "SLC8A1")],
        )

    def test_executed_source_p99_values_are_recorded_for_sherlock_comparison(self):
        expected = {
            ("gene", "ACTC1"): 3.287982,
            ("gene", "GATA4"): 2.463099,
            ("gene", "SLC8A1"): 2.587773,
            ("gene", "FLI1"): 3.283608,
            ("gene", "CDH5"): 3.359469,
            ("program", 2): 0.402983,
            ("program", 34): 0.184445,
            ("program", 15): 0.347766,
            ("program", 17): 0.331715,
            ("program", 27): 0.359415,
        }
        self.assertEqual(SOURCE_EXECUTED_P99, expected)


if __name__ == "__main__":
    unittest.main()
