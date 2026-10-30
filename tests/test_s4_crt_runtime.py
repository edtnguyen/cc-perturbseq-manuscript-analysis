"""S4 RECOMPUTE_CRT=True runtime tests; the numerical programDE engine is stubbed."""
from __future__ import annotations
import ast
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
from cc_perturbseq import crt_qq_backend as qqb
from cc_perturbseq import figure_s04_runtime as runtime
ROOT = Path(__file__).resolve().parents[1]

def options(tmp_path, **overrides):
    values = dict(programde_root=tmp_path/'programde', data_dir=tmp_path/'input',
                  cnmf_dir=tmp_path/'usage', guide_long_template=str(tmp_path/'{day}_guides.csv'),
                  guide_metadata_tsv=tmp_path/'guides.tsv', output_dir=tmp_path/'out')
    values.update(overrides)
    return runtime.build_crt_argv(**values)

@pytest.mark.parametrize('flag', ['--guide-long-csv-template', '--guide-long-template'])
def test_original_flag_and_compatibility_alias(flag):
    args=qqb.parse_args([flag, '/example/{day}_guides.csv'])
    assert args.guide_long_csv_template == '/example/{day}_guides.csv'

def test_builder_validates_complete_notebook_command_not_jupyter_argv(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['ipykernel_launcher.py', '-f', 'kernel.json'])
    argv=options(tmp_path)
    args=qqb.parse_args(argv)
    assert '--guide-long-csv-template' in argv
    assert '--guide-long-template' not in argv
    assert args.B == 10000 and args.k == 50 and args.require_covariates
    assert args.ntc_label == ['SAFE', 'non-targeting', 'NTC']
    assert args.overwrite and args.no_plot
    with pytest.raises(SystemExit):
        qqb.parse_args(argv+['--not-a-real-flag'])

def test_notebook_recompute_cell_uses_actual_builder(tmp_path):
    nb=json.loads((ROOT/'notebooks/supplementary/Figure_S04.ipynb').read_text())
    cell=next(c for c in nb['cells'] if 'runtime_report = launch_crt_recomputation' in ''.join(c.get('source', [])))
    branch=ast.parse(''.join(cell['source'])).body[0]
    assert isinstance(branch,ast.If)
    captured=[]
    ns=dict(RECOMPUTE_CRT=True, CRT_STAGE="all", build_crt_argv=runtime.build_crt_argv,
            launch_crt_recomputation=lambda argv, **kw: captured.append((qqb.parse_args(argv),kw)),
            display=lambda *a:None, qqb=qqb, CRT_PYTHON=sys.executable,
            PROGRAMDE_ROOT=tmp_path/'programde', DATA_DIR=tmp_path/'input', CNMF_DIR=tmp_path/'usage',
            GUIDE_LONG_TEMPLATE=str(tmp_path/'{day}_guides.csv'),GUIDE_METADATA_TSV=tmp_path/'meta.tsv',
            RECOMPUTED_ROOT=tmp_path/'out', DAYS=['d0','d1','d2','d3'],K=50,DENSITY_THRESHOLD='0_2',
            B=10000,N_JOBS=4,NTC_LABELS=['SAFE','non-targeting','NTC'])
    exec(compile(ast.Module(body=[branch],type_ignores=[]),'S4 raw cell','exec'),ns)
    assert captured[0][0].guide_long_csv_template == ns['GUIDE_LONG_TEMPLATE']
    assert captured[0][1]['python_executable'] == sys.executable
    assert ns['INPUT_MODEL_ROOT'] == tmp_path/'out/k_50.dt_0_2'

def test_missing_python_fails_before_launch(tmp_path):
    with pytest.raises(FileNotFoundError,match='CRT Python'):
        runtime.launch_crt_recomputation(options(tmp_path),python_executable=tmp_path/'absent-python')

def test_worker_uses_one_configuration_and_propagates_overwrite(tmp_path, monkeypatch):
    from cc_perturbseq import figure_s04_compute as compute
    calls=[]
    monkeypatch.setattr(runtime,'check_crt_runtime',lambda args:{'python_executable':sys.executable})
    monkeypatch.setattr(qqb,'main',lambda argv: calls.append(('qq',qqb.parse_args(argv))) or 0)
    monkeypatch.setattr(compute,'compute_replicate_betas',lambda **kw:calls.append(('rep',kw)))
    argv=options(tmp_path,k=17,density_threshold='0.4',ntc_labels=['NTC'],overwrite=True)
    report=runtime.run_worker(argv)
    assert report['status']=='completed'
    assert [x[0] for x in calls]==['qq','rep']
    assert calls[1][1]['k']==17 and calls[1][1]['density_threshold']=='0_4'
    assert calls[1][1]['ntc_labels']==['NTC'] and calls[1][1]['overwrite'] is True

def test_preflight_failure_never_starts_crt(tmp_path,monkeypatch):
    calls=[]
    def fail(args): raise ImportError('wrong selected interpreter')
    monkeypatch.setattr(runtime,'check_crt_runtime',fail)
    monkeypatch.setattr(qqb,'main',lambda argv:calls.append(argv))
    with pytest.raises(ImportError,match='wrong selected interpreter'):
        runtime.run_worker(options(tmp_path))
    assert calls==[]

def test_failed_full_day_does_not_run_replicates(tmp_path,monkeypatch):
    from cc_perturbseq import figure_s04_compute as compute
    monkeypatch.setattr(runtime,'check_crt_runtime',lambda args:{})
    def fail(argv): raise RuntimeError('full-day failure')
    monkeypatch.setattr(qqb,'main',fail)
    calls=[]
    monkeypatch.setattr(compute,'compute_replicate_betas',lambda **kw:calls.append(kw))
    with pytest.raises(RuntimeError,match='full-day failure'):
        runtime.run_worker(options(tmp_path))
    assert calls==[]
    assert json.loads((tmp_path/'out/s4_crt_runtime.json').read_text())['status']=='failed'

STUB_ENGINE = '''
import os, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

def prepare_crt_inputs(*, adata, **kwargs):
    return adata

def run_all_genes_union_crt(*,inputs,B,calibrate_skew_normal=False,
                            return_raw_pvals=False,return_skew_normal=False,**kwargs):
    if (return_raw_pvals or return_skew_normal) and not calibrate_skew_normal:
        raise ValueError("raw/skew supplemental outputs require calibration")
    k=len(inputs.uns['program_names'])
    with open(os.environ['S4_TEST_CALLS'],'a') as f:
        f.write(json.dumps({'B':B,'n':inputs.n_obs,'k':k,'python':sys.executable,
                            'covariates':list(inputs.obsm['covar'].columns)})+'\\n')
    values=np.linspace(-1,1,3*k).reshape(3,k)
    factor=1+float(inputs.X[0,0])/100
    beta=pd.DataFrame(values*factor,index=['A','B','NTC'],columns=inputs.uns['program_names'])
    p=pd.DataFrame(np.geomspace(1e-4,0.9,3*k).reshape(3,k),index=beta.index,columns=beta.columns)
    result={'betas_df':beta,'pvals_df':p}
    if return_raw_pvals: result['pvals_raw_df']=p
    if return_skew_normal: result['pvals_skew_df']=p
    return result

def build_ntc_group_inputs(**kw): return [],{}, {}, {}
def make_ntc_groups_ensemble(**kw): return {0:[]}
def crt_pvals_for_ntc_groups_ensemble(*,inputs,**kw):
    return {0:pd.DataFrame(np.linspace(0.01,0.99,2*len(inputs.uns['program_names'])).reshape(2,-1),
                            index=['ntc1','ntc2'],columns=inputs.uns['program_names'])}
def crt_pvals_for_ntc_groups_ensemble_skew(**kw): return crt_pvals_for_ntc_groups_ensemble(**kw)
'''

def test_entire_s4_notebook_recompute_branch_fresh_python(tmp_path):
    ad=pytest.importorskip('anndata')
    md=pytest.importorskip('mudata')
    pytest.importorskip('muon')
    from test_standalone_figure_ports import _isolated_checkout, _run_notebook
    root=_isolated_checkout(tmp_path)
    inputs=root/'input';inputs.mkdir()
    usage=root/'usage';usage.mkdir()
    engine=root/'test_programde';(engine/'src/sceptre').mkdir(parents=True)
    (engine/'src/__init__.py').write_text('')
    (engine/'src/sceptre/__init__.py').write_text(STUB_ENGINE)
    metadata=root/'guide_metadata.tsv'
    pd.DataFrame({'guide_name':['A_g1','B_g1','NTC_g1'],'target_name':['A','B','NTC']}).to_csv(metadata,sep='\t',index=False)
    for day in ['d0','d1','d2','d3']:
        n=400; index=pd.Index([f'{day}_cell{i}' for i in range(n)])
        obs=pd.DataFrame({'biological_sample':['r1']*200+['r2']*200,
                          'pct_counts_mt':np.arange(n)%5+1,'total_counts':np.arange(n)+20,
                          'guides_per_cell':np.arange(n)%2+1,'guide_umi_counts':np.arange(n)%7+5,
                          'n_genes_by_counts':np.arange(n)%3+2},index=index)
        rna=ad.AnnData(np.tile(np.arange(1,4),(n,1)).astype(float),obs=obs,
                       var=pd.DataFrame(index=['G1','G2','G3']))
        guides=ad.AnnData(np.ones((n,3)),obs=pd.DataFrame(index=index),
                         var=pd.DataFrame(index=['A_g1','B_g1','NTC_g1']))
        md.MuData({'GEX':rna,'guide':guides}).write_h5mu(inputs/f'{day}_normalized_gex_and_guide.h5mu')
        pd.DataFrame(np.arange(n*50).reshape(n,50)%17+1,index=index,
                     columns=[f'P{i}' for i in range(50)]).to_csv(usage/f'{day}.usages.k_50.dt_0_2.consensus.txt',sep='\t')
        pd.DataFrame({'cell_id':index,'grna_target':['A_g1','B_g1']*(n//2),'value':1}).to_csv(inputs/f'{day}_guides.csv',index=False)
    path=root/'notebooks/supplementary/Figure_S04.ipynb';nb=json.loads(path.read_text())
    for cell in nb['cells']:
        if cell['cell_type']=='code' and 'RECOMPUTE_CRT = False' in ''.join(cell['source']):
            s=''.join(cell['source']).replace('RECOMPUTE_CRT = False','RECOMPUTE_CRT = True')
            s+=f'\nPROGRAMDE_ROOT=Path({str(engine)!r})\nDATA_DIR=Path({str(inputs)!r})\nCNMF_DIR=Path({str(usage)!r})\n'
            s+=f'GUIDE_LONG_TEMPLATE={str(inputs/"{day}_guides.csv")!r}\nN_JOBS=1\n'
            cell['source']=s.splitlines(keepends=True)
    path.write_text(json.dumps(nb))
    output=root/'outputs';call_log=root/'calls.jsonl'
    _run_notebook(root,'Figure_S04.ipynb',{'outputs':{'root':str(output)},
                  'ccperturbseq':{'target_metadata':str(metadata)}},
                  {'CCPERTURBSEQ_S4_PYTHON':sys.executable,'S4_TEST_CALLS':str(call_log)})
    calls=[json.loads(s) for s in call_log.read_text().splitlines()]
    assert len(calls)==12
    assert sorted(c['B'] for c in calls)==[127]*8+[10000]*4
    assert all(c['k']==50 for c in calls)
    assert all(Path(c['python']).resolve()==Path(sys.executable).resolve() for c in calls)
    assert all(any('biological_sample' in s for s in c['covariates']) for c in calls if c['n']==400)
    assert all(not any('biological_sample' in s for s in c['covariates']) for c in calls if c['n']==200)
    report=json.loads((output/'Figure_S04/recomputed_crt/s4_crt_runtime.json').read_text())
    assert report['status']=='completed' and report['sceptre_R_required'] is False
    assert len(list((output/'Figure_S04/panel_b').glob('*.svg')))==4
    assert len(list((output/'Figure_S04/panel_c').glob('*.svg')))==5
