from pathlib import Path
p=Path('src/cc_perturbseq/figure_s02_parse_mix_plots.py')
s=p.read_text()
old='    keep = cell_metrics.cell_type.isin(["K562","hiPSC"]) & cell_metrics["sample"].astype(str).isin(original.SAMPLE_ORDER)'
new='''    # The original theme hides minor grids. Explicitly suppress their automatic
    # calculation: mizani otherwise fails when the observed log range has only
    # one major tick. Major ticks, limits, values and visible source style stay intact.
    for plot in (a, b):
        plot.scales.get_scales("y").minor_breaks = []
''' + old
assert old in s
p.write_text(s.replace(old,new))
p=Path('tests/test_standalone_figure_ports.py')
p.write_text(p.read_text()+'''

def test_s2_top_guide_single_major_log_tick_is_renderable():
    from cc_perturbseq.figure_s02_parse_mix_plots import plot_panel_l
    import matplotlib.pyplot as plt
    samples = ["CS1", "DC1", "Dual", "no_guide"]
    rows = [{"sample": sample, "cell_type": kind, "top_guide_umis": value,
             "total_guide_umis": value}
            for kind in ["K562", "hiPSC"] for sample in samples
            for value in [30, 31, 34, 35, 36]]
    fig, selected = plot_panel_l(pd.DataFrame(rows))
    fig.canvas.draw()
    assert len(selected) == len(rows)
    plt.close(fig)
''')
p=Path('docs/STANDALONE_FIGURE_REPRODUCTION.md')
p.write_text(p.read_text()+'''\nA rendering-only compatibility guard disables invisible automatic log minor\nbreaks in S2l. The original theme hides these grids; major ticks, limits,\nvalues, layers, and source functions are unchanged. This prevents a mizani\nIndexError for a narrow range containing only one major tick.\n''')
