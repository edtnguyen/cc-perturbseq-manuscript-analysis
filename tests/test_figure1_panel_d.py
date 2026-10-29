from __future__ import annotations

import unittest

import numpy as np

from cc_perturbseq.figure1_panel_d import (
    FIGURE1D_MANUSCRIPT_COLORS,
    FIGURE1D_MANUSCRIPT_PROGRAMS,
    FIGURE1D_THRESHOLD,
    assign_figure1d_programs,
)


class Figure1PanelDPaletteTests(unittest.TestCase):
    def test_manuscript_program_set_and_exact_stage_palette(self):
        self.assertEqual(len(FIGURE1D_MANUSCRIPT_PROGRAMS), 21)
        self.assertEqual(
            set(FIGURE1D_MANUSCRIPT_PROGRAMS),
            set(FIGURE1D_MANUSCRIPT_COLORS),
        )
        expected = {
            22: "#C8B4BC",
            29: "#4A4A4A",
            7: "#8A8F96",
            28: "#3E9A5F",
            30: "#1F6F4A",
            13: "#86C86A",
            48: "#4D9C72",
            15: "#F39C12",
            17: "#8E44AD",
            14: "#B794C8",
            19: "#6F78D7",
            23: "#F1C40F",
            2: "#D73027",
            3: "#FB6A4A",
            8: "#C51B8A",
            11: "#F04E37",
            42: "#8C2D04",
            24: "#6BAED6",
            33: "#9ECAE1",
            34: "#084594",
            27: "#2171B5",
        }
        self.assertEqual(FIGURE1D_MANUSCRIPT_COLORS, expected)

    def test_global_dominance_does_not_force_unselected_cells_into_palette(self):
        scores = np.zeros((3, 50), dtype=float)
        # P0 is globally dominant even though selected P2 is second-highest.
        scores[0, 0] = 0.8
        scores[0, 2] = 0.7
        # P2 is globally dominant and above the display threshold.
        scores[1, 2] = 0.8
        scores[1, 0] = 0.7
        # The source uses a strict > 0.02 threshold.
        scores[2, 2] = FIGURE1D_THRESHOLD

        dominant, maximum, displayed = assign_figure1d_programs(scores)
        self.assertEqual(dominant.tolist(), [0, 2, 2])
        np.testing.assert_allclose(maximum, [0.8, 0.8, FIGURE1D_THRESHOLD])
        self.assertEqual(displayed.tolist(), [False, True, False])

    def test_explicit_column_program_ids_are_respected(self):
        scores = np.asarray([[0.7, 0.8, 0.1]], dtype=float)
        dominant, _, displayed = assign_figure1d_programs(
            scores,
            program_ids=[2, 34, 0],
        )
        self.assertEqual(dominant.tolist(), [34])
        self.assertEqual(displayed.tolist(), [True])


if __name__ == "__main__":
    unittest.main()
