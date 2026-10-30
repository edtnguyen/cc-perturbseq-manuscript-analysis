"""Exercise real programDE output semantics, metadata schema, and resume safety."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from cc_perturbseq import crt_qq_backend as qqb
from cc_perturbseq import figure_s04_contract as contract
from cc_perturbseq import figure_s04_runtime as runtime

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('calibrated', [False, True])
def test_output_flags_follow_real_api_constraint(calibrated):
    seen=[]
    def strict_engine(**kwargs):
        seen.append(kwargs)
        if (kwargs['return_raw_pvals'] or kwargs['return_skew_normal']) and not kwargs['calibrate_skew_normal']:
            raise ValueError('incompatible flags')
        return {'betas_df': 'beta'}
    result=contract.run_program_crt(SimpleNamespace(run_all_genes_union_crt=strict_engine),
                                   inputs='inputs',B=127,n_jobs=1,calibrate_skew_normal=calibrated)
    assert result['betas_df']=='beta'
    for key in ['calibrate_skew_normal','return_raw_pvals','return_skew_normal']:
        assert seen[0][key] is calibrated


def test_no_skew_outputs_are_empirical_not_mislabeled_skew():
    p=pd.DataFrame([[0.5]],index=['A'],columns=['P0'])
    empirical,skew=contract.extract_pvalue_outputs({'pvals_df':p},calibrate_skew_normal=False)
    assert empirical is p and skew is None
    raw=p/2
    empirical,skew=contract.extract_pvalue_outputs({'pvals_df':p,'pvals_raw_df':raw},calibrate_skew_normal=True)
    assert empirical is raw and skew is p
    with pytest.raises(RuntimeError,match='required empirical'):
        contract.extract_pvalue_outputs({'pvals_df':p},calibrate_skew_normal=True)


def test_v43_schema_and_missing_control_targets(tmp_path):
    p=tmp_path/'v43.tsv'
    pd.DataFrame({'guide_id':['ACVR1B_+_123','NTC_1','SAFE_1','non-targeting_1'],
                  'intended_target_name':['ACVR1B',None,np.nan,'']}).to_csv(p,sep='\t',index=False)
    names=['ACVR1B_+_123','NTC_1','SAFE_1','non-targeting_1']
    mapping=qqb.build_guide2gene(names,qqb.DEFAULT_NTC_LABELS,str(p),None,None)
    assert mapping==dict(zip(names,['ACVR1B','NTC','SAFE','non-targeting']))
    report=contract.audit_guide_metadata(qqb.parse_args(['--guide-metadata-tsv',str(p)]))
    assert report['legacy_name_difference_count']==0
    assert report['n_metadata_mapped']==1
    assert report['n_name_inferred_including_controls']==3
    with pytest.raises(KeyError,match='Explicit target column'):
        qqb.build_guide2gene(names,qqb.DEFAULT_NTC_LABELS,str(p),'guide_id','not_the_column')


def test_changed_and_conflicting_metadata_are_reported(tmp_path):
    p=tmp_path/'v43.tsv'
    pd.DataFrame({'guide_id':['OLD_guide'],'intended_target_name':['NEW']}).to_csv(p,sep='\t',index=False)
    report=contract.audit_guide_metadata(qqb.parse_args(['--guide-metadata-tsv',str(p)]))
    assert report['legacy_name_difference_count']==1
    assert report['legacy_name_differences'][0]['metadata_target']=='NEW'
    pd.DataFrame({'guide_id':['A_1','A_1'],'intended_target_name':['A','B']}).to_csv(p,sep='\t',index=False)
    with pytest.raises(ValueError,match='Conflicting'):
        contract.audit_guide_metadata(qqb.parse_args(['--guide-metadata-tsv',str(p)]))


def completed_run(tmp_path, *, legacy=True):
    root=tmp_path/'run';root.mkdir()
    meta=tmp_path/'v43.tsv'
    pd.DataFrame({'guide_id':['A_1','NTC_1'],'intended_target_name':['A',None]}).to_csv(meta,sep='\t',index=False)
    argv=runtime.build_crt_argv(programde_root=tmp_path/'programde',data_dir=tmp_path,
         cnmf_dir=tmp_path,guide_long_template=str(tmp_path/'{day}.csv'),
         guide_metadata_tsv=meta,output_dir=root,days=['d0'])
    args=qqb.parse_args(argv)
    audit=contract.audit_guide_metadata(args)
    receipt={'status':'failed','argv':argv,'error':'original incompatible flags'}
    if not legacy:receipt['guide_metadata_audit']=audit
    (root/'s4_crt_runtime.json').write_text(json.dumps(receipt))
    model=root/'k_50.dt_0_2';day=model/'d0';day.mkdir(parents=True)
    for name in ['all_days_QQ_points.csv']:(model/name).write_text('pvalue\n0.5\n')
    pd.DataFrame({'day':['d0'],'status':['completed']}).to_csv(model/'QQ_export_manifest.csv',index=False)
    outputs={'combined_QQ_points':'d0_QQ_points.csv','target_empirical_wide':'empirical.tsv','ntc_empirical_wide':'ntc.tsv'}
    for name in outputs.values():(day/name).write_text('pvalue\n0.5\n')
    cfg={'day':'d0','k':50,'B':10000,'density_threshold_token':'0_2',
         'ntc_labels':args.ntc_label,'pvalue_outputs':outputs}
    for key in ['h5mu_path','usage_path','guide_long_csv']:
        f=tmp_path/(key+'.input');f.write_text('input')
        os.utime(f,ns=(1_000_000_000,1_000_000_000));cfg[key]=str(f)
    os.utime(meta,ns=(1_000_000_000,1_000_000_000))
    (day/'d0_run_config.json').write_text(json.dumps(cfg))
    return args,audit,model


@pytest.mark.parametrize('legacy',[False,True])
def test_resume_keeps_full_day_files_untouched(tmp_path,legacy):
    args,audit,model=completed_run(tmp_path,legacy=legacy)
    before={p:p.read_bytes() for p in model.rglob('*') if p.is_file()}
    report=contract.validate_completed_full_day(args,audit)
    assert report['status']=='verified_completed_full_day'
    assert before=={p:p.read_bytes() for p in before}


@pytest.mark.parametrize('failure',['mapping','parameters','input','missing_output','incomplete_day'])
def test_resume_rejects_changed_or_incomplete_inputs(tmp_path,failure):
    args,audit,model=completed_run(tmp_path)
    if failure=='mapping':audit['legacy_name_difference_count']=1
    elif failure=='parameters':args.B=9999
    elif failure=='input':
        cfg=json.loads((model/'d0/d0_run_config.json').read_text())
        os.utime(cfg['usage_path'],ns=(9_000_000_000_000_000_000,)*2)
    elif failure=='missing_output':(model/'d0/empirical.tsv').unlink()
    else:pd.DataFrame({'day':['d0'],'status':['failed']}).to_csv(model/'QQ_export_manifest.csv',index=False)
    with pytest.raises(ValueError):contract.validate_completed_full_day(args,audit)


def test_replicate_only_worker_does_not_repeat_full_crt(tmp_path,monkeypatch):
    from cc_perturbseq import figure_s04_compute as compute
    args,audit,model=completed_run(tmp_path)
    argv=json.loads((Path(args.out_dir)/'s4_crt_runtime.json').read_text())['argv']
    monkeypatch.setattr(runtime,'check_crt_runtime',lambda args:{'guide_metadata_audit':audit})
    def forbidden(argv):raise AssertionError('full-day CRT rerun unexpectedly')
    monkeypatch.setattr(qqb,'main',forbidden)
    calls=[]
    monkeypatch.setattr(compute,'compute_replicate_betas',lambda **kwargs:calls.append(kwargs))
    report=runtime.run_worker(argv,stage='replicates')
    assert report['status']=='completed' and len(calls)==1
    assert report['full_day_resume']['status']=='verified_completed_full_day'


def test_real_programde_rejects_old_call_and_accepts_exact_replacement(tmp_path):
    root=os.environ.get('PROGRAMDE_REAL_ROOT')
    if not root:
        pytest.skip('Set PROGRAMDE_REAL_ROOT to a real programDE checkout; mandatory in repair CI.')
    script=tmp_path/'real_contract.py'
    script.write_text(r'''
from pathlib import Path
import importlib.util, json, sys
import numpy as np
import pandas as pd
from cc_perturbseq import crt_qq_backend as qqb
from cc_perturbseq.figure_s04_contract import run_program_crt, extract_pvalue_outputs, audit_guide_metadata
root=Path(sys.argv[1]).resolve()
programde=qqb.load_programde(str(root))
spec=importlib.util.spec_from_file_location('programde_test_fixtures',root/'tests/synthetic_data.py')
module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
adata,_=module.make_synthetic_adata(np.random.default_rng(91),n_cells=160,n_programs=3,
    n_genes=3,guides_per_gene=2,n_covariates=2,effect_gene='gene_0',effect_size=1.0)
inputs=programde.prepare_crt_inputs(adata)
try:
    programde.run_all_genes_union_crt(inputs=inputs,B=127,n_jobs=1,
        calibrate_skew_normal=False,return_raw_pvals=True,return_skew_normal=False)
except ValueError as exc:
    assert 'requires calibrate_skew_normal=True' in str(exc)
else:
    raise AssertionError('real API did not reject the reported invalid combination')
corrected=run_program_crt(programde,inputs=inputs,B=127,n_jobs=1,calibrate_skew_normal=False)
reference=programde.run_all_genes_union_crt(inputs=inputs,B=127,n_jobs=1,calibrate_skew_normal=False)
pd.testing.assert_frame_equal(corrected['betas_df'],reference['betas_df'])
pd.testing.assert_frame_equal(corrected['pvals_df'],reference['pvals_df'])
assert np.isfinite(corrected['betas_df'].to_numpy()).all()
empirical,skew=extract_pvalue_outputs(corrected,calibrate_skew_normal=False)
assert empirical is corrected['pvals_df'] and skew is None
calibrated=run_program_crt(programde,inputs=inputs,B=127,n_jobs=1,calibrate_skew_normal=True)
pd.testing.assert_frame_equal(calibrated['betas_df'],corrected['betas_df'])
empirical,skew=extract_pvalue_outputs(calibrated,calibrate_skew_normal=True)
assert empirical is calibrated['pvals_raw_df'] and skew is calibrated['pvals_df']
assert 'pvals_raw_df' not in corrected and 'pvals_skew_df' not in corrected
args=qqb.parse_args(['--guide-metadata-tsv',str(root/'notebooks/300genes_guide_metadata_v43.tsv')])
audit=audit_guide_metadata(args)
Path(sys.argv[2]).write_text(json.dumps({'real_programde_api':'passed',
    'old_invalid_call':'rejected','corrected_beta_matches_direct_API':True,
    'beta_unchanged_by_skew_output_flags':True,'guide_metadata_audit':audit},indent=2))
print('REAL PROGRAMDE CONTRACT PASSED')
print(json.dumps(audit,indent=2))
''')
    report=tmp_path/'real_programde_report.json'
    env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'src')
    p=subprocess.run([sys.executable,str(script),root,str(report)],cwd=tmp_path,env=env,
                     text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)
    (tmp_path/'real_programde.log').write_text(p.stdout)
    assert p.returncode==0,p.stdout
    assert json.loads(report.read_text())['real_programde_api']=='passed'
