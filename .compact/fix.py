from pathlib import Path
p=Path('src/cc_perturbseq/compact_plot_inputs.py')
s=p.read_text()
old='''                frame = pd.read_csv(files[key], index_col=0, dtype={index_col: str}, float_precision="round_trip", keep_default_na=False, na_values=[""])'''
new='''                frame = pd.read_csv(files[key], dtype={index_col: str}, float_precision="round_trip", keep_default_na=False, na_values=[""])
                frame = frame.set_index(index_col)
                frame.index.name = None if index_col == "Unnamed: 0" else index_col'''
assert old in s
p.write_text(s.replace(old,new))
p=Path('tests/test_compact_plot_inputs.py')
s=p.read_text()
s+='''\n\ndef test_numeric_looking_cell_ids_remain_strings(tmp_path,inputs):
    _,_,cfg,tables=inputs
    directory=Path(cfg['compact_plot_sources']['mixed_figure_dir'])
    cells=tables['cell_metrics'].copy()
    cells.index=pd.Index([f'{i:06d}' for i in range(len(cells))])
    cells.to_csv(directory/MIX_FILES['cell_metrics'])
    spec=resolve_sources(cfg,datasets=('ipsc_k562_mix',))['ipsc_k562_mix']
    export_one(spec,tmp_path/'result','ipsc_k562_mix')
    result,_=load_plot_bundle(tmp_path/'result','ipsc_k562_mix')
    assert result['cell_metrics'].index.tolist()==cells.index.tolist()
'''
p.write_text(s)
