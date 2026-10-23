from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from figure1_common import (
    ROOT, ASSEMBLED, load_screen, export_figure,
    draw_panel_1a, draw_candidate_c, draw_panel_1c,
    draw_panel_1d_placeholder, draw_panel_1e_placeholder, panel_letter,
)

screen = load_screen()

fig = plt.figure(figsize=(7.18, 6.15), facecolor="white")
fig.subplots_adjust(left=0.055, right=0.985, top=0.975, bottom=0.055)
outer = GridSpec(2, 1, figure=fig, height_ratios=[1.12, 0.88], hspace=0.34)

top = outer[0].subgridspec(1, 2, width_ratios=[0.42, 0.58], wspace=0.22)
bottom = outer[1].subgridspec(1, 3, width_ratios=[0.27, 0.35, 0.38], wspace=0.32)

ax_a = fig.add_subplot(top[0, 0])
ax_b = fig.add_subplot(top[0, 1])
ax_c = fig.add_subplot(bottom[0, 0])
ax_d = fig.add_subplot(bottom[0, 1])
ax_e = fig.add_subplot(bottom[0, 2])

for ax in [ax_a, ax_b, ax_c, ax_d, ax_e]:
    ax.set_facecolor("white")

# Draw panels directly so the assembled PDF/SVG remains vector.
draw_panel_1a(ax_a, compact=True)
draw_candidate_c(ax_b, screen, compact=True, selected=True)
draw_panel_1c(ax_c, screen, compact=True)
draw_panel_1d_placeholder(ax_d, compact=True)
draw_panel_1e_placeholder(ax_e, compact=True)

for ax, letter in zip([ax_a, ax_b, ax_c, ax_d, ax_e], list("abcde")):
    panel_letter(ax, letter, compact=True)

# Save white-background assembled figure.
svg = ASSEMBLED / "Figure1_balanced.svg"
pdf = ASSEMBLED / "Figure1_balanced.pdf"
png = ASSEMBLED / "Figure1_balanced.png"
fig.savefig(svg, bbox_inches="tight", facecolor="white")
fig.savefig(pdf, bbox_inches="tight", facecolor="white")
fig.savefig(png, bbox_inches="tight", facecolor="white", dpi=600)
plt.close(fig)

print(svg)
print(pdf)
print(png)
