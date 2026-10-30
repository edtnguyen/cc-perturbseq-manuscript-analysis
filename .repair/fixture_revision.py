from pathlib import Path
p=Path('tests/test_standalone_figure_ports.py');s=p.read_text()
old='    matrix[48:,5:9]+=20'
new='''    matrix[48:,5:9]+=25
    # Match total marker burden between lineages, and provide a broad moderate-
    # expression lineage signal that passes the source HVG mean cutoffs. Nine
    # very highly expressed marker genes alone are intentionally not the HVGs.
    matrix[:48,10:210] += rng.poisson(1.0,size=(48,200))
    matrix[48:,210:410] += rng.poisson(1.0,size=(48,200))'''
assert old in s;s=s.replace(old,new);p.write_text(s)
