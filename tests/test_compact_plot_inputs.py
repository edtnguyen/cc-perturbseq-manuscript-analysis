"""Compact exports on numerical fixtures; no study-data sizes are asserted."""
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
import pytest
from cc_perturbseq.compact_plot_inputs import (resolve_sources, export_all, export_one,
    _saved_h5mu_tables, MIX_FILES, main)
from cc_perturbseq.geo_plot_data import load_plot_bundle, build_plot_bundle, prepare_hippo_dotplots, load_mixed_plot_data
from test_public_adapters import make_mixed_geo
from geo_plot_fixtures import make_hippo_geo


@pytest.fixture
def inputs(tmp_path):
    hippo = tmp_path/'submitted/hippo.h5mu';hippo.parent.mkdir()
    q,p,_=make_hippo_geo(hippo)
    mixed=tmp_path/'submitted/mixed.h5mu';bundle=make_mixed_geo(mixed)
    figures=tmp_path/'accepted/Figure_S02';figures.mkdir(parents=True)
    tables={'cell_metrics':bundle.mod['rna'].obs.copy(), 'saturation':bundle.uns['saturation'].copy(),
            'cluster_annotation':bundle.uns['cluster_annotation'].copy(),
            'median_diagnostics':bundle.uns['manuscript_median_diagnostics'].copy()}
    for key,rel in MIX_FILES.items():
        dest=figures/rel;dest.parent.mkdir(exist_ok=True,parents=True)
        tables[key].to_csv(dest,index=key=='cell_metrics')
    (figures/'Figure_S02_parse_mix_provenance.json').write_text(json.dumps({'n_cells_after_qc':len(tables['cell_metrics'])}))
    config={'compact_plot_sources':{'hippo_query':str(q),'hippo_projection':str(p),'mixed_figure_dir':str(figures)}}
    return hippo,mixed,config,tables


def assert_tables_same(left,right):
    assert set(left)==set(right)
    for k in left:
        # Stored diagnostic/annotation row IDs are metadata, not cell IDs.
        a,b=left[k].copy(),right[k].copy()
        if k not in {'cell_metrics'} and '_dot_' not in k:
            a.index=a.index.astype(str);b.index=b.index.astype(str)
        pd.testing.assert_frame_equal(a,b,check_dtype=False,check_categorical=False,check_index_type=False,rtol=1e-6,atol=1e-7)


def test_direct_sources_work_with_no_h5mu_and_match_full_route(tmp_path,inputs,monkeypatch):
    import mudata as md
    hp,mixed,config,_=inputs
    expected={}
    for key,p in [('hippo_timecourse',hp),('ipsc_k562_mix',mixed)]:
        out=tmp_path/('reference_'+key);build_plot_bundle(p,out,key)
        expected[key]=load_plot_bundle(out,key)[0]
    hp.unlink();mixed.unlink()
    monkeypatch.setattr(md,'read_h5mu',lambda *a,**k: (_ for _ in ()).throw(AssertionError('H5MU read')))
    sources=resolve_sources(config);root=tmp_path/'repo'
    dry=export_all(sources,root);assert dry['status']=='dry_run' and not root.exists()
    report=export_all(sources,root,apply=True)
    assert report['single_cell_objects_uploaded'] is False
    assert not list(root.rglob('*.h5mu')) and not list(root.rglob('*.h5ad'))
    for key in sources:
        path=root/'data/reproduction/geo_plots'/key
        data,meta=load_plot_bundle(path,key)
        assert_tables_same(data,expected[key])
        assert meta['source']['sha256_scope']=='selected uncompressed plotting tables'
        assert report['datasets'][key]['size_bytes']==sum(f.stat().st_size for f in path.iterdir())
    cells=load_mixed_plot_data(root/'data/reproduction/geo_plots/ipsc_k562_mix').cell_metrics
    assert cells.loc['cell0','total_guide_umis']==0 and np.isnan(cells.loc['cell0','guide_pool_specificity'])


def test_saved_tables_route_never_opens_count_matrices(tmp_path,inputs,monkeypatch):
    import anndata as ad
    import anndata.io
    import mudata as md
    hp,mix,_,_=inputs
    expected={}
    for key,p in [('hippo_timecourse',hp),('ipsc_k562_mix',mix)]:
        out=tmp_path/('expected_'+key);build_plot_bundle(p,out,key);expected[key]=load_plot_bundle(out,key)[0]
    def forbidden(*a,**kw):raise AssertionError('Attempted full single-cell read')
    monkeypatch.setattr(md,'read_h5mu',forbidden);monkeypatch.setattr(ad,'read_h5ad',forbidden)
    real=anndata.io.read_elem;seen=[]
    def guarded(elem):
        seen.append(elem.name)
        assert elem.name.startswith('/uns/') or elem.name=='/mod/rna/obs',elem.name
        return real(elem)
    monkeypatch.setattr(anndata.io,'read_elem',guarded)
    sources=resolve_sources({},hippo_h5mu=hp,mixed_h5mu=mix)
    export_all(sources,tmp_path/'repo',apply=True)
    assert seen
    for dataset in sources:
        tables,record=load_plot_bundle(tmp_path/'repo/data/reproduction/geo_plots'/dataset,dataset)
        assert_tables_same(tables,expected[dataset])
        assert not record['source']['recomputed_from_modalities']
        assert not record['source']['raw_matrix_integrity_verified']


def test_saved_and_direct_hippo_dotplots_match(tmp_path,inputs):
    import matplotlib.pyplot as plt
    hp,_,config,_=inputs
    from cc_perturbseq.geo_submission import _apply_source_style
    direct=tmp_path/'direct';saved=tmp_path/'saved'
    export_one(resolve_sources(config,datasets=('hippo_timecourse',))['hippo_timecourse'],direct,'hippo_timecourse')
    export_one(resolve_sources({},datasets=('hippo_timecourse',),hippo_h5mu=hp)['hippo_timecourse'],saved,'hippo_timecourse')
    _apply_source_style()
    for panel in ('Figure_3f','Figure_4c'):
        a,b=prepare_hippo_dotplots(direct,panel),prepare_hippo_dotplots(saved,panel)
        for x,y in zip(a[:2],b[:2]):
            x.make_figure();y.make_figure();x.fig.canvas.draw();y.fig.canvas.draw()
            np.testing.assert_array_equal(np.asarray(x.fig.canvas.buffer_rgba()),np.asarray(y.fig.canvas.buffer_rgba()))
        plt.close('all')


def test_required_inputs_and_mode_choice_are_explicit(tmp_path,inputs):
    hp,mix,cfg,_=inputs
    with pytest.raises(ValueError,match='not both'):resolve_sources(cfg,hippo_h5mu=hp,hippo_query='x')
    with pytest.raises(ValueError,match='not both'):resolve_sources(cfg,mixed_h5mu=mix,mixed_figure_dir='x')
    with pytest.raises(FileNotFoundError,match='No integrated H5MU'):resolve_sources({})
    mixed_dir=Path(cfg['compact_plot_sources']['mixed_figure_dir'])
    (mixed_dir/MIX_FILES['cell_metrics']).unlink()
    assert list(resolve_sources(cfg,datasets=('hippo_timecourse',)))==['hippo_timecourse']
    with pytest.raises(FileNotFoundError,match='cell_metrics'):resolve_sources(cfg)


def test_partial_s2_tables_do_not_get_published(tmp_path,inputs):
    _,_,cfg,_=inputs;spec=resolve_sources(cfg,datasets=('ipsc_k562_mix',))['ipsc_k562_mix']
    path=spec['files']['cell_metrics'];frame=pd.read_csv(path,index_col=0).iloc[1:];frame.to_csv(path)
    with pytest.raises(ValueError,match='QC cell count'):export_one(spec,tmp_path/'out','ipsc_k562_mix')
    assert not (tmp_path/'out').exists()


def test_collector_supports_compact_flag_and_keeps_other_files(tmp_path,inputs):
    from cc_perturbseq.publication_bundle import main as collect_main
    _,_,cfg,_=inputs;root=tmp_path/'repo';(root/'config').mkdir(parents=True)
    (root/'config/publication_source_files.json').write_text(json.dumps(cfg))
    marker=root/'data/reproduction/untouched.txt';marker.parent.mkdir(parents=True);marker.write_text('keep')
    result=collect_main(['--repo',str(root),'--geo-plots-only','--compact-geo-inputs','--apply'])
    assert result['status']=='collected' and marker.read_text()=='keep'
    assert len(load_plot_bundle(root/'data/reproduction/geo_plots/hippo_timecourse','hippo_timecourse')[0])==12


def test_atomic_export_does_not_touch_existing_bundle_on_failure(tmp_path,inputs,monkeypatch):
    import cc_perturbseq.compact_plot_inputs as compact
    _,_,cfg,_=inputs;root=tmp_path/'repo';sources=resolve_sources(cfg)
    export_all(sources,root,apply=True)
    old={p:p.read_bytes() for p in (root/'data').rglob('*') if p.is_file()}
    real=compact.export_one
    def broken(spec,where,ds):
        if ds=='ipsc_k562_mix':raise RuntimeError('failed second dataset')
        return real(spec,where,ds)
    monkeypatch.setattr(compact,'export_one',broken)
    with pytest.raises(RuntimeError,match='failed second'):export_all(sources,root,apply=True,overwrite=True)
    assert {p:p.read_bytes() for p in (root/'data').rglob('*') if p.is_file()}==old


def test_script_exports_only_requested_dataset_and_preserves_backup(tmp_path,inputs):
    hp,_,cfg,_=inputs;root=tmp_path/'repo';root.mkdir();cp=root/'cfg.json';cp.write_text(json.dumps(cfg))
    args=['--repo',str(root),'--source-config',str(cp),'--dataset','hippo_timecourse','--hippo-h5mu',str(hp)]
    assert main(args)['status']=='dry_run'
    assert main([*args,'--apply'])['status']=='exported'
    with pytest.raises(FileExistsError):main([*args,'--apply'])
    assert main([*args,'--apply','--overwrite'])['status']=='exported'
    assert list((root/'data').glob('hippo_timecourse.backup-*'))
    assert not (root/'data/reproduction/geo_plots/ipsc_k562_mix').exists()


def test_saved_hippo_missing_reference_fails_without_matrix_fallback(tmp_path,inputs):
    import h5py
    hp,_,_,_=inputs
    with h5py.File(hp,'a') as f:del f['uns/reference_tables/Figure_3f_genes_dot_size']
    with pytest.raises(ValueError,match='Saved plotting information absent'):
        _saved_h5mu_tables(hp,'hippo_timecourse')
