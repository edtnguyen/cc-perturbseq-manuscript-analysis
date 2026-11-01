"""Portable execution tests. Fixtures are synthetic, not study results."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import numpy as np
import pandas as pd
import pytest
import yaml

ROOT=Path(__file__).resolve().parents[1]


def isolated(tmp_path):
    root=tmp_path/'checkout';root.mkdir()
    for name in ('src','config','notebooks','docs','scripts'):
        shutil.copytree(ROOT/name,root/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    (root/'data/reproduction').mkdir(parents=True)
    return root


def bundle_manifest(root, records):
    values=[]
    for key,path in records.items():
        values.append({'key':key,'path':path.relative_to(root).as_posix(),'size_bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    (root/'data/reproduction/bundle_manifest.json').write_text(json.dumps({'schema_version':1,'files':values},indent=2))


def run_notebook(root, relative):
    env=os.environ.copy();env.update(PYTHONPATH=str(root/'src'),MPLBACKEND='Agg',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
    for key in ('CCPERTURBSEQ_PUBLIC_PATHS','CCPERTURBSEQ_PATHS','IGVF_ACCESS_KEY','IGVF_SECRET_ACCESS_KEY'):
        env.pop(key,None)
    script='''
from pathlib import Path
import json, sys, socket
import IPython.display
IPython.display.display=lambda *args, **kwargs:None
# Public runtime must not attempt original filesystem paths or network calls
# for these self-contained supplementary fixtures.
def audit(event,args):
    if event=='open' and args and isinstance(args[0],(str,bytes)):
        value=str(args[0])
        if value.startswith(('/oak/','/scratch/users/','/home/users/')) or '/recovery/' in value:
            raise AssertionError('Institutional/archive dependency: '+value)
sys.addaudithook(audit)
socket.socket.connect=lambda *a,**k: (_ for _ in ()).throw(AssertionError('Unexpected network'))
nb=json.loads(Path(sys.argv[1]).read_text())
ns={'__name__':'__main__'}
for i,c in enumerate(nb['cells']):
    if c['cell_type']=='code':
        print('RUN_CELL',i,flush=True)
        exec(compile(''.join(c['source']),f'cell-{i}','exec'),ns)
print('ALL_CELLS_COMPLETED_PORTABLY',flush=True)
'''
    run=subprocess.run([sys.executable,'-c',script,relative],cwd=root,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    (root/'test-execution.log').write_text(run.stdout)
    assert run.returncode==0,run.stdout[-16000:]
    assert 'ALL_CELLS_COMPLETED_PORTABLY' in run.stdout


def s4_files(root):
    model=root/'data/reproduction/s4/k_50.dt_0_2';model.mkdir(parents=True,exist_ok=True)
    rows=[];manifest=[];counts=[];records={};rng=np.random.default_rng(5)
    for day in ['d0','d1','d2','d3']:
        for series in ['observed_target','ntc_pseudotarget']:
            for p in [.001,.02,.1,.4,.8]:
                rows.append({'day':day,'pvalue_type':'empirical','series_type':series,'pvalue':p})
        for rep in ['r1','r2']:
            path=model/'biological_replicate_effects/betas_by_replicate'/day/f'{rep}.tsv';path.parent.mkdir(parents=True,exist_ok=True)
            pd.DataFrame(rng.normal(size=(4,3)),index=['A','B','C','NTC'],columns=['P0','P1','P2']).to_csv(path,sep='\t')
            manifest.append({'day':day,'biological_sample':rep,'status':'written','beta_path':str(path.relative_to(model/'biological_replicate_effects'))})
            records[f's4.beta.{day}.{rep}']=path
            for gene in ['A','B','C']:
                counts.append({'day':day,'biological_sample':rep,'target_name':gene,'n_target_cells_with_cnmf_usage':20})
        p=model/day/f'{day}_target_pvalues_skew_wide.tsv';p.parent.mkdir(parents=True,exist_ok=True)
        pd.DataFrame(np.full((4,3),.001),index=['A','B','C','NTC'],columns=['P0','P1','P2']).to_csv(p,sep='\t');records[f's4.pvalues.{day}']=p
    p=model/'all_days_QQ_points.csv';pd.DataFrame(rows).to_csv(p,index=False);records['s4.qq']=p
    p=model/'biological_replicate_effects/replicate_beta_manifest.csv';pd.DataFrame(manifest).to_csv(p,index=False);records['s4.manifest']=p
    p=model/'biological_replicate_effects/plots/DIAGNOSTIC_target_cell_counts_by_replicate.csv';p.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(counts).to_csv(p,index=False);records['s4.counts']=p
    return records,model


def test_s4_full_notebook_no_original_paths_no_crt(tmp_path):
    root=isolated(tmp_path);records,_=s4_files(root);bundle_manifest(root,records)
    run_notebook(root,'notebooks/supplementary/Figure_S04.ipynb')
    assert len(list((root/'results/Figure_S04/panel_b').glob('*.svg')))==4
    assert len(list((root/'results/Figure_S04/panel_c').glob('*.svg')))==5
    settings=json.loads((root/'results/Figure_S04/panel_c/S4c_selection_settings.json').read_text())
    assert settings['min_target_cells_per_replicate']==20


def s6_files(root):
    import h5py
    paths=yaml.safe_load((root/'config/paths.example.yaml').read_text())
    def get(k):
        a,b=k.split('.');p=root/paths[a][b];p.parent.mkdir(parents=True,exist_ok=True);return p
    genes=[f'G{i}' for i in range(16)];records={}
    p=get('ccperturbseq.program_annotations')
    pd.DataFrame({'Program number':range(8),'Name':[f'Program {i}' for i in range(8)],'Representative genes':[', '.join([genes[(2*i+j)%16] for j in range(4)]) for i in range(8)]}).to_excel(p,index=False);records['ccperturbseq.program_annotations']=p
    summary=pd.DataFrame({'program_name':range(8),'Automatic Timepoint':np.repeat(['D0','D1','D2','D3'],2)})
    for day in ['D0','D1','D2','D3']:summary[f'Top 5 specific regulators (FDR<0.1) {day}']='R1, R2, R3'
    associations=[];motifs=[]
    for program in range(8):
        day=f'D{program//2}'
        for i,reg in enumerate(['R1','R2','R3']):associations.append({'target_name':reg,'program_name':program,'Sample':day,'approx_log2FC':(-1)**i*(program+1)/10,'adj_pval':.001*(i+1)})
        for i,motif in enumerate(['ETS','TEAD','GATA','HAND']):
            q=10.**(-i-1);motifs.append({'program':program,'day':day,'tf':motif,'qval':q,'neg_log10_qval':-np.log10(q),'AUC':.6+i/20})
    wb=root/'data/figure_02/cNMF_50_0_2_032526_6co.xlsx';wb.parent.mkdir(parents=True,exist_ok=True)
    with pd.ExcelWriter(wb) as writer:
        summary.to_excel(writer,sheet_name='Summary',index=False)
        pd.DataFrame(associations).to_excel(writer,sheet_name='Perturbation Association 1',index=False)
    weights=np.zeros((8,16))
    for i in range(8):
        weights[i,2*i:2*i+2]=10;weights[i,(2*i+2)%16]=1;weights[i,(2*i+3)%16]=.5
    p=get('ccperturbseq.cnmf_spectra');pd.DataFrame(weights,index=range(8),columns=genes).to_csv(p,sep='\t');records['ccperturbseq.cnmf_spectra']=p
    p=get('multiome.enhancer_enrichment_results');pd.DataFrame(motifs).to_csv(p,sep='\t',index=False);records['multiome.enhancer_enrichment_results']=p
    p=get('multiome.cluster_motif_map');pd.DataFrame({'cluster_id':range(4),'tf_motif_id':['ETS','TEAD','GATA','HAND']}).to_csv(p,sep='\t',index=False);records['multiome.cluster_motif_map']=p
    p=get('multiome.modisco_compiled_h5')
    with h5py.File(p,'w') as f:
        for i in range(4):f.create_dataset(f'pos_patterns/{i}/contrib_scores',data=np.eye(4)[np.arange(7)%4]*.3)
    records['multiome.modisco_compiled_h5']=p
    return records


def test_s6_full_notebook_with_bundled_motifs(tmp_path):
    root=isolated(tmp_path);records=s6_files(root);bundle_manifest(root,records)
    run_notebook(root,'notebooks/supplementary/Figure_S06.ipynb')
    assert len(list((root/'results/Figure_S06/panel_all').glob('*.svg')))==5
    tables=root/'results/Figure_S06/panel_all/_program_annotation_source'
    assert len(list(tables.glob('*_motif_enrichment_matrix.csv')))==4


def test_s2_full_notebook_reads_only_geo_integrated_file(tmp_path):
    pytest.importorskip('anndata');pytest.importorskip('mudata')
    from test_public_adapters import make_mixed_geo
    root=isolated(tmp_path);make_mixed_geo(root/'data/geo/iPSC_K562_mixed_integrated.h5mu')
    run_notebook(root,'notebooks/supplementary/Figure_S02.ipynb')
    for panel in ('i','j','l','m','n'):
        assert list((root/f'results/Figure_S02/panel_{panel}').glob('*.svg'))


def test_all_seven_sources_pass_static_audit():
    result=subprocess.run([sys.executable,str(ROOT/'scripts/audit_public_inputs.py')],capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    report=json.loads(result.stdout);assert len(report['notebooks'])==7 and report['status']=='PASS'


def test_missing_data_audit_is_not_false_success():
    result=subprocess.run([sys.executable,str(ROOT/'scripts/audit_public_inputs.py'),'--check-data'],capture_output=True,text=True)
    assert result.returncode!=0
    assert json.loads(result.stdout)['status']=='INCOMPLETE_DATA'


def test_explicit_s4_counts_do_not_probe_historical_configs(tmp_path,monkeypatch):
    from cc_perturbseq.figure_s04_postcrt import load_or_rebuild_target_counts
    records,model=s4_files(tmp_path)
    for day in ['d0','d1','d2','d3']:
        (model/day/f'{day}_run_config.json').write_text(json.dumps({'h5mu_path':'/oak/do-not-open.h5mu','usage_path':'/oak/do-not-open.txt','guide_long_csv':'/oak/do-not-open.csv'}))
    original=Path.is_file
    def guarded(p):
        assert not str(p).startswith('/oak/'),str(p)
        return original(p)
    monkeypatch.setattr(Path,'is_file',guarded)
    result=load_or_rebuild_target_counts(model,['d0','d1','d2','d3'],tmp_path/'out',counts_csv=records['s4.counts'])
    assert len(result)==24


def test_all_seven_input_preparations_with_no_institutional_inputs(tmp_path,monkeypatch):
    ad=pytest.importorskip('anndata');md=pytest.importorskip('mudata')
    from cc_perturbseq import public_data as public
    from test_public_adapters import write_days,validation_tables,make_mixed_geo
    from cc_perturbseq.publication_bundle import validation_bridge
    from cc_perturbseq.geo_submission import SCHEMA
    root=isolated(tmp_path);fixtures=root/'fixture_downloads';fixtures.mkdir()
    day_paths,cnmf_paths,order_path,order=write_days(fixtures)
    records=s6_files(root);s4,_=s4_files(root);records.update(s4)
    _,paths=public.load_public_paths(root)
    for kind in ('expression','cnmf'):
        p=root/f'data/reproduction/programs/{kind}_cell_order.tsv.gz';shutil.copy2(order_path,p);records['cell_order.'+kind]=p
    for key in ['screen.genomewide_utest_pos_vs_neg','ccperturbseq.cnmf_usage','ccperturbseq.target_metadata','ccperturbseq.figure1d_program_metadata','ccperturbseq.figure1d_program_order']:
        p=public.value(paths,key);p.parent.mkdir(parents=True,exist_ok=True);p.write_text('fixture input\n');records[key]=p
    _,original,portal=validation_tables(fixtures);bridge=validation_bridge(original,portal)
    p=root/'data/reproduction/screen/validation_guide_compatibility.tsv';bridge.to_csv(p,sep='\t',index=False);records['validation.compatibility']=p
    bundle_manifest(root,records)
    hippo=md.MuData({'rna':ad.AnnData(np.ones((2,3))),'cnmf':ad.AnnData(np.ones((2,50)))})
    hippo.uns.update(geo_schema=SCHEMA,dataset='hippo_timecourse');geo=root/'data/geo';geo.mkdir(exist_ok=True)
    hippo.write_h5mu(geo/'Hippo_timecourse_integrated.h5mu');make_mixed_geo(geo/'iPSC_K562_mixed_integrated.h5mu')
    manifest=json.loads((root/'config/public_data_manifest.json').read_text())
    downloaded={**dict(zip(manifest['groups']['expression'],day_paths)),**dict(zip(manifest['groups']['cnmf'],cnmf_paths)),
                'IGVFFI5157STFU':portal,'IGVFFI9001HIYI':portal}
    class Client:
        def download(self,record,cache):return downloaded[record['accession']]
    monkeypatch.setattr(public,'_CLIENT',Client())
    def check_path(p):
        if isinstance(p,(str,Path)):
            assert not str(p).startswith(('/oak/','/scratch/users/','/home/users/')),p
    real_open=Path.open
    def guarded(p,*a,**kw):check_path(p);return real_open(p,*a,**kw)
    monkeypatch.setattr(Path,'open',guarded)
    for figure in ['Figure_01','Figure_02','Figure_03','Figure_04','Figure_S02','Figure_S04','Figure_S06']:
        result=public.prepare_figure_inputs(figure,root,paths)
        assert result['institutional_source_files_required'] is False
    assembled=ad.read_h5ad(public.value(paths,'ccperturbseq.integrated_anndata'))
    assert list(assembled.obs_names)==order
