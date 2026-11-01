from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import pytest
from cc_perturbseq.public_data import assemble_day_objects, adapt_validation_counts, extract_hippo_geo, load_mixed_geo
from cc_perturbseq.publication_bundle import validation_bridge, _read_obs_index


def validation_tables(tmp_path):
    portal=pd.DataFrame({'guide_id':['NEW_A','NEW_B','NTC'], 'gene':['A','B','NON-TARGETING'],
        'D01':[11,12,13], 'D02':[21,22,23], 'D3neg1':[31,32,33], 'D3neg2':[41,42,43],
        'D3pos1':[51,52,53], 'D3pos2':[61,62,63]})
    orig=portal.rename(columns={'guide_id':'sgRNA','gene':'Gene','D01':'D0R1','D02':'D0R2','D3neg1':'D3NR1','D3neg2':'D3NR2','D3pos1':'D3PR1','D3pos2':'D3PR2'}).iloc[[2,0,1]].copy()
    orig['sgRNA']=['NTC_old','A_old','B_old'];orig['D3TOTAL1']=[90,190,290];orig['D3TOTAL2']=[80,180,280]
    op=tmp_path/'original.tsv';pp=tmp_path/'portal.tsv.gz'
    orig.to_csv(op,sep='\t',index=False);portal.to_csv(pp,sep='\t',index=False)
    return orig,op,pp


def test_validation_bridge_preserves_labels_order_totals_and_all_six_counts(tmp_path):
    original,op,pp=validation_tables(tmp_path)
    bridge=validation_bridge(op,pp);bp=tmp_path/'bridge.tsv';bridge.to_csv(bp,sep='\t',index=False)
    output=adapt_validation_counts(pp,bp,tmp_path/'adapted.tsv')
    result=pd.read_csv(output,sep='\t')
    pd.testing.assert_frame_equal(result[original.columns].reset_index(drop=True),original.reset_index(drop=True),check_dtype=False)


def test_validation_bridge_rejects_changed_count_vectors(tmp_path):
    original,op,pp=validation_tables(tmp_path)
    original.iloc[0,original.columns.get_loc('D3PR1')]+=1;original.to_csv(op,sep='\t',index=False)
    with pytest.raises(ValueError,match='uniquely reconcile'):
        validation_bridge(op,pp)


def write_days(tmp_path, *, alter_guides=False):
    ad=pytest.importorskip('anndata');md=pytest.importorskip('mudata')
    from scipy import sparse
    paths,cnmfs=[],[]
    for i in range(4):
        ids=[f'd{i}_c0',f'd{i}_c1']
        obs=pd.DataFrame({'day':[f'd{i}']*2,'biological_sample':['r1','r2'],'leiden_res_0.5':['0','1']},index=ids)
        a=ad.AnnData(sparse.csr_matrix(np.array([[i+1,2],[3,i+4]],dtype=float)),obs=obs,var=pd.DataFrame({'gene_symbol':['G1','G2']},index=['E1','E2']))
        a.layers['raw_counts']=a.X.copy();a.obsm['X_umap']=np.array([[i,0],[i,1]],dtype=float)
        a.obsm['guide_assignment']=sparse.csr_matrix([[1,0],[0,1]])
        a.uns['guide_names']=np.array(['B','A'] if alter_guides and i==2 else ['A','B'])
        path=tmp_path/f'd{i}.h5ad';a.write_h5ad(path);paths.append(path)
        c=ad.AnnData(np.arange(100,dtype=float).reshape(2,50)+i,obs=obs,var=pd.DataFrame(index=[str(j) for j in range(50)]))
        c.obsm['guide_assignment']=a.obsm['guide_assignment'].copy();c.uns['guide_names']=a.uns['guide_names'].copy()
        c.varm['loadings']=np.arange(100).reshape(50,2);c.uns['var_names']=np.array(['E1','E2'])
        cp=tmp_path/f'd{i}.h5mu';md.MuData({'cNMF':c,'rna':a}).write_h5mu(cp);cnmfs.append(cp)
    order=[f'd{i}_c{j}' for j in (1,0) for i in (3,0,2,1)]
    op=tmp_path/'order.tsv.gz';pd.DataFrame({'cell_id':order}).to_csv(op,sep='\t',index=False)
    return paths,cnmfs,op,order


def test_day_assembly_preserves_values_identifiers_embedding_and_global_order(tmp_path):
    ad=pytest.importorskip('anndata');md=pytest.importorskip('mudata')
    paths,cnmfs,op,order=write_days(tmp_path)
    output=assemble_day_objects(paths,op,tmp_path/'combined.h5ad')
    result=ad.read_h5ad(output)
    assert list(result.obs_names)==order and list(result.uns['guide_names'])==['A','B']
    original={cell:(a.X.toarray()[j],a.obsm['X_umap'][j]) for p in paths for a in [ad.read_h5ad(p)] for j,cell in enumerate(a.obs_names)}
    for i,cell in enumerate(order):
        np.testing.assert_array_equal(result.X.toarray()[i],original[cell][0]);np.testing.assert_array_equal(result.obsm['X_umap'][i],original[cell][1])
    cp=assemble_day_objects(cnmfs,op,tmp_path/'combined.h5mu',kind='cnmf')
    c=md.read_h5mu(cp).mod['cNMF'];assert list(c.obs_names)==order
    assert c.X[0].sum()>1  # No row normalization was introduced.
    np.testing.assert_array_equal(c.varm['loadings'],np.arange(100).reshape(50,2))
    mtime=cp.stat().st_mtime_ns;assemble_day_objects(cnmfs,op,cp,kind='cnmf');assert cp.stat().st_mtime_ns==mtime


def test_changed_guide_column_order_is_not_silently_concatenated(tmp_path):
    pytest.importorskip('anndata');paths,_,op,_=write_days(tmp_path,alter_guides=True)
    with pytest.raises(ValueError,match='guide-column order'):
        assemble_day_objects(paths,op,tmp_path/'output.h5ad')


def test_missing_or_duplicate_cell_ids_are_not_invented(tmp_path):
    pytest.importorskip('anndata');paths,_,op,order=write_days(tmp_path)
    pd.DataFrame({'cell_id':order[:-1]}).to_csv(op,sep='\t',index=False)
    with pytest.raises(ValueError,match='exactly match'):
        assemble_day_objects(paths,op,tmp_path/'output.h5ad')


def test_hdf5_cell_order_reader_uses_declared_index(tmp_path):
    import h5py
    p=tmp_path/'index.h5ad'
    with h5py.File(p,'w') as f:
        g=f.create_group('obs');g.attrs['_index']='barcode';g.create_dataset('barcode',data=np.asarray(['cellB','cellA'],dtype=h5py.string_dtype()))
    assert _read_obs_index(p).tolist()==['cellB','cellA']


def make_mixed_geo(path):
    ad=pytest.importorskip('anndata');md=pytest.importorskip('mudata')
    from scipy import sparse
    from cc_perturbseq.geo_submission import SCHEMA,_metrics_from_sparse_guides
    from cc_perturbseq.figure_s02_parse_mix_data import manuscript_median_diagnostics
    n=96;rng=np.random.default_rng(12)
    obs=pd.DataFrame({'cell_type':np.repeat(['K562','hiPSC'],48),'sample':np.tile(np.repeat(['CS1','DC1','Dual','no_guide'],12),2),
        'total_counts':rng.integers(500,1500,n),'n_genes_by_counts':rng.integers(100,500,n),'pct_counts_mt':np.ones(n)},index=[f'cell{i}' for i in range(n)])
    a=ad.AnnData(sparse.csr_matrix(np.ones((n,2))),obs=obs,var=pd.DataFrame(index=['G1','G2']))
    x=rng.integers(0,40,(n,4));x[0]=0
    guides=ad.AnnData(sparse.csr_matrix(x),obs=obs.copy(),var=pd.DataFrame({'guide_pool':['K562','K562','hiPSC','hiPSC']},index=['K1','K2','I1','I2']))
    bundle=md.MuData({'rna':a,'crispr':guides});a.obs=_metrics_from_sparse_guides(bundle)
    bundle.uns.update(geo_schema=SCHEMA,dataset='ipsc_k562_mix',cluster_annotation=pd.DataFrame({'cluster':['a'],'lineage':['K562']}),manuscript_median_diagnostics=manuscript_median_diagnostics(a.obs))
    bundle.uns['saturation']=pd.DataFrame([{'primer':primer,'mean_assigned_reads_per_cell':reads,'median_transcript_umis_per_cell':reads/3,'source_file':'report.csv'} for primer in ['CS1','DC1','Dual','no_guide'] for reads in [10,20,40,80]])
    path.parent.mkdir(parents=True,exist_ok=True);bundle.write_h5mu(path)
    return bundle


def test_geo_mix_uses_frozen_labels_and_preserves_zero_guide_cells(tmp_path):
    path=tmp_path/'mix.h5mu';bundle=make_mixed_geo(path);result=load_mixed_geo(path)
    assert len(result.cell_metrics)==96 and result.cell_metrics.loc['cell0','total_guide_umis']==0
    assert result.cell_metrics.cell_type.tolist()==bundle.mod['rna'].obs.cell_type.tolist()
    assert result.provenance['original_source_files_required'] is False


def test_geo_hippo_extraction_uses_only_submitted_modalities(tmp_path):
    ad=pytest.importorskip('anndata');md=pytest.importorskip('mudata')
    from cc_perturbseq.geo_submission import SCHEMA
    rna=ad.AnnData(np.array([[1,2],[3,4]],dtype=float));cnmf=ad.AnnData(np.array([[.1,.9],[.8,.2]]))
    m=md.MuData({'rna':rna,'cnmf':cnmf});m.uns.update(geo_schema=SCHEMA,dataset='hippo_timecourse')
    p=tmp_path/'hippo.h5mu';m.write_h5mu(p)
    q,c=tmp_path/'derived/rna.h5ad',tmp_path/'derived/cnmf.h5ad';extract_hippo_geo(p,q,c)
    np.testing.assert_array_equal(ad.read_h5ad(q).X,rna.X);np.testing.assert_array_equal(ad.read_h5ad(c).X,cnmf.X)


def test_collection_is_lossless_relative_and_does_not_run_crt(tmp_path,monkeypatch):
    import h5py
    from cc_perturbseq import publication_bundle as collect_module
    from test_public_notebooks import s4_files
    root=tmp_path/'repo';root.mkdir();(root/'config').mkdir();(root/'data').mkdir()
    source=tmp_path/'source';source.mkdir()
    original,original_path,portal_path=validation_tables(source)
    expression=source/'expression.h5ad';cnmf=source/'cnmf.h5mu'
    for path,prefix in ((expression,'obs'),(cnmf,'mod/cNMF/obs')):
        with h5py.File(path,'w') as f:
            g=f.create_group(prefix);g.attrs['_index']='barcode';g.create_dataset('barcode',data=np.array(['cellB','cellA'],dtype=h5py.string_dtype()))
    numeric=source/'exact.tsv';numeric.write_text('index\tP0\ncellB\t0.123456789\ncellA\t9.8\n')
    motif=source/'modisco_compiled.h5'
    with h5py.File(motif,'w') as f:f.create_dataset('pos_patterns/0/contrib_scores',data=np.eye(4))
    _,model=s4_files(source)
    config={'files':[{'key':'test.numeric','source':str(numeric),'destination':'data/reproduction/programs/exact.tsv.gz','gzip':True},
                     {'key':'test.motif','source':str(motif),'destination':'data/reproduction/motifs/modisco_compiled.h5','gzip':False}],
            'integrated_cell_order_source':str(expression),'cnmf_cell_order_source':str(cnmf),
            'validation_guide_counts_source':str(original_path),'s4_results_source':str(model)}
    from geo_plot_fixtures import make_hippo_geo
    hippo=source/'hippo.h5mu';make_hippo_geo(hippo)
    mixed=source/'mixed.h5mu';make_mixed_geo(mixed)
    config['geo_plot_sources']={'hippo_timecourse':str(hippo),'ipsc_k562_mix':str(mixed)}
    cp=root/'config/sources.json';cp.write_text(json.dumps(config))
    (root/'config/public_data_manifest.json').write_text(json.dumps({'files':{'IGVFFI5157STFU':{}}}))
    class Client:
        def __init__(self,**kw):pass
        def download(self,*a,**kw):return portal_path
    monkeypatch.setattr(collect_module,'IgvfClient',Client)
    dry=collect_module.collect(root,cp,apply=False);assert dry['status']=='dry_run' and not (root/'data/reproduction').exists()
    result=collect_module.collect(root,cp,apply=True)
    import gzip
    with gzip.open(root/'data/reproduction/programs/exact.tsv.gz','rb') as f:assert f.read()==numeric.read_bytes()
    for record in result['files']:
        p=root/record['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==record['sha256']
    manifest=pd.read_csv(root/'data/reproduction/s4/k_50.dt_0_2/biological_replicate_effects/replicate_beta_manifest.csv')
    assert all(not Path(p).is_absolute() for p in manifest.beta_path)
    assert all((root/'data/reproduction/s4/k_50.dt_0_2/biological_replicate_effects'/p).is_file() for p in manifest.beta_path)
    assert not list((root/'data/reproduction').rglob('*run_config.json'))
    with pytest.raises(FileExistsError):collect_module.collect(root,cp,apply=True)
