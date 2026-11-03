"""Precomputed and regenerated routes on actual HDF5 fixtures (not study data)."""
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
import pytest

from cc_perturbseq.geo_plot_data import (build_plot_bundle, load_plot_bundle,
    load_mixed_plot_data, prepare_hippo_dotplots, prepare_geo_plot_inputs, write_plot_bundle)
from cc_perturbseq.public_data import load_public_paths
from test_public_adapters import make_mixed_geo
from geo_plot_fixtures import make_hippo_geo

ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture(scope='session')
def submitted(tmp_path_factory):
    root=tmp_path_factory.mktemp('submitted')
    hippo=root/'Hippo_timecourse_integrated.h5mu';make_hippo_geo(hippo)
    mixed=root/'iPSC_K562_mixed_integrated.h5mu';make_mixed_geo(mixed)
    return hippo,mixed


def paths_at(root):
    (root/'config').mkdir(parents=True,exist_ok=True)
    shutil.copy2(ROOT/'config/paths.example.yaml',root/'config/paths.example.yaml')
    return load_public_paths(root)[1]


def test_hippo_recomputed_tables_and_rendering_match_original_producers(tmp_path,submitted):
    from cc_perturbseq import geo_submission as geo
    import mudata as md
    import matplotlib.pyplot as plt
    hp,_=submitted;out=tmp_path/'hippo'
    report=build_plot_bundle(hp,out,'hippo_timecourse')
    assert len(report['tables'])==12
    assert len(report['source']['matched_embedded_reference_tables'])==12
    data,_=load_plot_bundle(out,'hippo_timecourse')
    bundle=md.read_h5mu(hp)
    geo._apply_source_style()
    with geo._local_hippo_inputs(bundle) as (q,p):
        for panel,fn in geo._hippo_producers().items():
            direct=fn(q,p);reloaded=prepare_hippo_dotplots(out,panel)
            for i,kind in enumerate(['genes','programs']):
                np.testing.assert_allclose(direct[i].dot_color_df,reloaded[i].dot_color_df,rtol=1e-12,atol=1e-12)
                np.testing.assert_allclose(direct[i].dot_size_df,reloaded[i].dot_size_df,rtol=1e-12,atol=1e-12)
                # Real renderer parity, not just a test of serialized numbers.
                direct[i].make_figure();reloaded[i].make_figure()
                direct[i].fig.canvas.draw();reloaded[i].fig.canvas.draw()
                image1=np.asarray(direct[i].fig.canvas.buffer_rgba())
                image2=np.asarray(reloaded[i].fig.canvas.buffer_rgba())
                np.testing.assert_array_equal(image1,image2,err_msg=panel+' '+kind)
            pd.testing.assert_frame_equal(direct[2],reloaded[2],check_dtype=False,rtol=1e-12,atol=1e-12)
            pd.testing.assert_frame_equal(direct[3],reloaded[3],check_dtype=False,rtol=1e-12,atol=1e-12)
            plt.close('all')
    assert data['Figure_3f_genes']['n_cells'].min()>0


def test_mixed_labels_zero_guides_threshold_and_saturation_preserved(tmp_path,submitted):
    from cc_perturbseq.public_data import load_mixed_geo
    _,mixed=submitted
    direct=load_mixed_geo(mixed)
    out=tmp_path/'mixed';build_plot_bundle(mixed,out,'ipsc_k562_mix')
    result=load_mixed_plot_data(out)
    pd.testing.assert_frame_equal(direct.cell_metrics,result.cell_metrics,check_dtype=False,check_categorical=False)
    assert len(result.cell_metrics)==96
    assert result.cell_metrics.loc['cell0','total_guide_umis']==0
    assert np.isnan(result.cell_metrics.loc['cell0','guide_pool_specificity'])
    assert result.saturation['primer'].nunique()==4
    assert result.cell_metrics['K562_guide_umis_ge3'].le(result.cell_metrics['K562_guide_umis']).all()


def test_default_never_opens_geo_object_or_network_and_no_implicit_fallback(tmp_path,submitted,monkeypatch):
    import mudata as md
    from cc_perturbseq.geo_download import GeoClient
    paths=paths_at(tmp_path)
    build_plot_bundle(submitted[1],Path(paths['precomputed']['mixed']),'ipsc_k562_mix')
    def forbidden(*a,**k):raise AssertionError('Default tried to read a raw file or download')
    monkeypatch.setattr(md,'read_h5mu',forbidden)
    monkeypatch.setattr(GeoClient,'download',forbidden)
    info=prepare_geo_plot_inputs('ipsc_k562_mix',tmp_path,paths)
    assert info['mode']=='precomputed' and not info['requires_geo_download']
    with pytest.raises(FileNotFoundError,match='Missing published'):
        prepare_geo_plot_inputs('hippo_timecourse',tmp_path,paths)
    with pytest.raises(ValueError,match='Select mode'):
        prepare_geo_plot_inputs('ipsc_k562_mix',tmp_path,paths,accession='GSE12345')


def test_regeneration_goes_to_separate_cache_and_matches_published(tmp_path,submitted):
    from cc_perturbseq.igvf_download import digest_file
    paths=paths_at(tmp_path);pub=Path(paths['precomputed']['mixed'])
    build_plot_bundle(submitted[1],pub,'ipsc_k562_mix')
    original_hash=digest_file(pub/'manifest.json','sha256')
    report=prepare_geo_plot_inputs('ipsc_k562_mix',tmp_path,paths,mode='local_file',local_file=submitted[1])
    assert report['published_comparison']['status']=='MATCH'
    assert Path(report['directory']).is_relative_to(tmp_path/'output/cache/derived')
    assert digest_file(pub/'manifest.json','sha256')==original_hash
    mtime=(Path(report['directory'])/'manifest.json').stat().st_mtime_ns
    again=prepare_geo_plot_inputs('ipsc_k562_mix',tmp_path,paths,mode='local_file',local_file=submitted[1])
    assert Path(again['directory'])==Path(report['directory'])
    assert (Path(again['directory'])/'manifest.json').stat().st_mtime_ns==mtime


def test_geo_accession_route_downloads_then_recomputes_and_preserves_reference(tmp_path,submitted,monkeypatch):
    from cc_perturbseq.geo_download import GeoClient
    calls=[]
    def download(self,accession,cache,**kw):
        calls.append((accession,kw));return submitted[1]
    monkeypatch.setattr(GeoClient,'download',download)
    paths=paths_at(tmp_path)
    info=prepare_geo_plot_inputs('ipsc_k562_mix',tmp_path,paths,mode='geo',accession='GSE12345',filename='mix.h5mu')
    assert calls[0][0]=='GSE12345' and calls[0][1]['filename']=='mix.h5mu'
    assert info['published_comparison']['status']=='NO_PUBLISHED_REFERENCE_BUNDLE'
    data,record=load_plot_bundle(info['directory'],'ipsc_k562_mix')
    assert record['source']['recomputed_from_modalities'] and len(data['cell_metrics'])==96


def test_tampered_guide_counts_rejected_before_bundle_written(tmp_path,submitted):
    import mudata as md
    m=md.read_h5mu(submitted[1]);m.mod['crispr'].X.data[0]+=100
    p=tmp_path/'tampered.h5mu';m.write_h5mu(p)
    with pytest.raises(AssertionError):build_plot_bundle(p,tmp_path/'out','ipsc_k562_mix')
    assert not (tmp_path/'out').exists()


def test_bad_schema_or_wrong_experiment_not_interchangeable(tmp_path,submitted):
    import anndata as ad
    p=tmp_path/'expression.h5ad';ad.AnnData(np.ones((2,3))).write_h5ad(p)
    with pytest.raises(ValueError,match='expression-only'):build_plot_bundle(p,tmp_path/'out','ipsc_k562_mix')
    with pytest.raises(ValueError,match='Not the submitted'):build_plot_bundle(submitted[1],tmp_path/'out','hippo_timecourse')


def test_plot_table_integrity_and_safe_paths(tmp_path,submitted):
    p=tmp_path/'bundle';build_plot_bundle(submitted[1],p,'ipsc_k562_mix')
    manifest=json.loads((p/'manifest.json').read_text())
    f=p/manifest['tables']['cell_metrics']['file'];raw=f.read_bytes();f.write_bytes(raw+b'changed')
    with pytest.raises(ValueError,match='checksum'):load_plot_bundle(p,'ipsc_k562_mix')
    f.write_bytes(raw);manifest['tables']['cell_metrics']['file']='../outside.tsv.gz'
    (p/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='escapes'):load_plot_bundle(p,'ipsc_k562_mix')


def test_geo_only_collection_dry_run_backups_and_no_raw_upload(tmp_path,submitted):
    from cc_perturbseq.publication_bundle import collect_geo_plots_only
    root=tmp_path/'repo';root.mkdir();config=root/'sources.json'
    config.write_text(json.dumps({'geo_plot_sources':dict(zip(['hippo_timecourse','ipsc_k562_mix'],map(str,submitted)))}))
    dry=collect_geo_plots_only(root,config)
    assert dry['status']=='dry_run' and not (root/'data').exists()
    report=collect_geo_plots_only(root,config,apply=True)
    assert report['status']=='collected' and not report['raw_objects_uploaded']
    assert not list((root/'data').rglob('*.h5mu')) and not list((root/'data').rglob('*.h5ad'))
    assert len(load_plot_bundle(root/'data/reproduction/geo_plots/hippo_timecourse','hippo_timecourse')[0])==12
    with pytest.raises(FileExistsError):collect_geo_plots_only(root,config,apply=True)
    marker=root/'data/reproduction/keep_me.txt';marker.write_text('unrelated companion')
    collect_geo_plots_only(root,config,apply=True,overwrite=True)
    assert marker.read_text()=='unrelated companion' and list((root/'data').glob('geo_plots.backup-*'))
