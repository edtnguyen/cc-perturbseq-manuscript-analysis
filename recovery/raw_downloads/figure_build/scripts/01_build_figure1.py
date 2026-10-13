from pathlib import Path
import matplotlib.pyplot as plt

from figure1_common import (
    ROOT, CANDIDATE, PANELS, load_screen, export_figure,
    draw_panel_1a, draw_candidate_a, draw_candidate_b, draw_candidate_c,
    draw_panel_1c, draw_panel_1d_placeholder, draw_panel_1e_placeholder,
)

screen = load_screen()

# Candidate 1b panels
candidates = [
    ("Fig1b_candidate_A_u_test_volcano_hits_only", draw_candidate_a, (4.4, 3.3)),
    ("Fig1b_candidate_B_mageck_placeholder", draw_candidate_b, (4.4, 3.3)),
    ("Fig1b_candidate_C_ranked_u_test_hits", lambda ax, df, compact=False: draw_candidate_c(ax, df, compact, selected=True), (4.4, 3.3)),
]
for name, fn, size in candidates:
    fig, ax = plt.subplots(figsize=size)
    if name.startswith("Fig1b_candidate_B"):
        fn(ax, compact=False)
    else:
        fn(ax, screen, compact=False)
    fig.tight_layout()
    export_figure(fig, CANDIDATE / name, transparent_png=True, dpi=600)
    plt.close(fig)

# Individual main panels
panel_specs = [
    ("Fig1a_screen_workflow", draw_panel_1a, (6.0, 2.7), None),
    ("Fig1b_genomewide_screen", lambda ax, compact=False: draw_candidate_c(ax, screen, compact, selected=True), (5.0, 3.0), None),
    ("Fig1c_known_regulator_recovery", lambda ax, compact=False: draw_panel_1c(ax, screen, compact), (3.0, 2.4), None),
    ("Fig1d_regulator_categories", draw_panel_1d_placeholder, (3.7, 2.4), None),
    ("Fig1e_validation", draw_panel_1e_placeholder, (3.7, 2.4), None),
]
for name, fn, size, _ in panel_specs:
    fig, ax = plt.subplots(figsize=size)
    fn(ax, compact=False)
    fig.tight_layout()
    export_figure(fig, PANELS / name, transparent_png=True, dpi=600)
    plt.close(fig)

print("Built candidate and individual Figure 1 panels.")
