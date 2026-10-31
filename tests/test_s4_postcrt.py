from pathlib import Path
import ast,hashlib,importlib.util,json
import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'src/cc_perturbseq/figure_s04_postcrt.py'
spec=importlib.util.spec_from_file_location('postcrt',p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def fixture_tables():
    effects=pd.DataFrame({
        'day':['d0']*6,'replicate_x':['r1']*6,'replicate_y':['r2']*6,
        'replicate_pair':['r1 vs r2']*6,
        'target_name':['A','A','B','C','D','NTC'],
        'program_name':['0','1','0','0','0','0'],
        'beta_x':[1,-2,3,0,4,999], 'beta_y':[1.5,-2.5,2,0.1,-3,-999]})
    pooled=pd.DataFrame({'day':['d0']*7,'target_name':['A','A','B','C','D','NOT_IN_BETAS','NTC'],
        'program_name':['0','1','0','0','0','0','0'],
        'pooled_pvalue':[0.001,0.6,0.001,0.04,0.9,0.8,1e-99]})
    pooled=pooled.loc[pooled.target_name!='NTC'].copy()
    pooled['pooled_qvalue']=m.bh_qvalue(pooled.pooled_pvalue)
    pooled['pooled_pvalue_source']='synthetic.tsv'
    rows=[]
    for sample in ['r1','r2']:
        for target in ['A','B','C','D']:
            n=19 if (sample,target)==('r2','B') else 20
            rows.append({'day':'d0','biological_sample':sample,'target_name':target,
                         'n_target_cells_all':1000,'n_target_cells_with_cnmf_usage':n})
    return effects,pooled,pd.DataFrame(rows)


def test_selection_count_boundary_and_daywide_bh():
    effects,pooled,counts=fixture_tables()
    selected,excluded,report=m.select_effect_pairs(effects,pooled,counts)
    # q for C = 0.04 * 6 / 3 = 0.08, despite NOT_IN_BETAS not being compared.
    assert set(zip(selected.target_name,selected.program_name))=={('A','0'),('C','0')}
    assert selected.loc[selected.target_name=='C','pooled_qvalue'].iloc[0]==pytest.approx(0.08)
    assert report.matched_finite_real_pairs.tolist()==[5]
    assert report.after_target_cell_filter.tolist()==[4]
    assert report.selected_pairs.tolist()==[2]
    assert excluded.loc[excluded.target_name=='B','passes_target_cell_filter'].tolist()==[False]
    # Zero beta is NOT removed from correlation inputs, only sign-agreement denominator.
    assert selected.loc[selected.target_name=='C','beta_x'].iloc[0]==0
    expected=effects.loc[(effects.target_name=='C') | ((effects.target_name=='A')&(effects.program_name=='0'))]
    np.testing.assert_array_equal(selected.beta_x,expected.beta_x)


def test_counts_are_unique_target_carriers_with_actual_usage():
    obs=pd.DataFrame({'biological_sample':['r1','r1','r1','r2']},index=['a','b','c','d'])
    guides=pd.DataFrame({'cell_id':['a','a','a','b','c','d','outside','d'],
        'grna_target':['A_g1','A_g2','A_g1','A_g1','B_g1','A_g1','A_g1','NTC_g1'],
        'value':[1,1,1,1,0,1,1,1]})
    counts,overlap=m.count_target_cells(obs,guides,['a','c','d'],
                                      {'A_g1':'A','A_g2':'A','B_g1':'B','NTC_g1':'NTC'},'d0')
    row=counts.set_index(['biological_sample','target_name']).loc[('r1','A')]
    assert row.n_target_cells_all==2 and row.n_target_cells_with_cnmf_usage==1
    assert set(counts.target_name)=={'A'}
    assert overlap.set_index('biological_sample').loc['r1','n_cells_with_cnmf_usage']==2


def test_q_is_computed_before_count_and_pair_filtering(tmp_path):
    folder=tmp_path/'d0';folder.mkdir()
    # Only A appears in compared pairs, but all 10 real tests enter BH.
    pd.DataFrame({'target_name':['A']+[f'X{i}' for i in range(9)]+['NTC'],
        'program_name':['0']*11,'pvalue':[0.02]+[0.9]*9+[1e-100]}).to_csv(folder/'d0_target_pvalues_skew_long.tsv',sep='\t',index=False)
    pooled,report=m.read_pooled_results(tmp_path,['d0'])
    assert len(pooled)==10 and report.n_tests_for_bh.iloc[0]==10
    assert pooled.loc[pooled.target_name=='A','pooled_qvalue'].iloc[0]==pytest.approx(0.2)


@pytest.mark.parametrize('preferred,expected',[('skew','skew'),('empirical','empirical')])
def test_pvalue_priority_and_explicit_type(tmp_path,preferred,expected):
    folder=tmp_path/'d0';folder.mkdir()
    for kind,pval in [('skew',0.001),('empirical',0.2)]:
        pd.DataFrame({'0':[pval]},index=['A']).to_csv(folder/f'd0_target_pvalues_{kind}_wide.tsv',sep='\t')
    _,report=m.read_pooled_results(tmp_path,['d0'],preferred)
    assert report.actual_pvalue_type.iloc[0]==expected


def test_pvalue_fallback_is_audited(tmp_path):
    folder=tmp_path/'d0';folder.mkdir()
    pd.DataFrame({'0':[0.01]},index=['A']).to_csv(folder/'d0_target_pvalues_empirical_wide.tsv',sep='\t')
    with pytest.warns(UserWarning,match='using empirical'):
        _,report=m.read_pooled_results(tmp_path,['d0'],'skew')
    assert report.actual_pvalue_type.iloc[0]=='empirical'


def test_missing_counts_never_silently_disable_filter(tmp_path):
    with pytest.raises(FileNotFoundError,match='NOT been silently skipped'):
        m.load_or_rebuild_target_counts(tmp_path,['d0'],tmp_path/'out')


def test_existing_source_counts_are_validated_and_recorded(tmp_path):
    _,_,counts=fixture_tables();path=tmp_path/'counts.csv';counts.to_csv(path,index=False)
    out=tmp_path/'output'
    read=m.load_or_rebuild_target_counts(tmp_path,['d0'],out,counts_csv=path)
    pd.testing.assert_frame_equal(read,counts)
    assert json.loads((out/'target_cell_count_provenance.json').read_text())['mode']=='existing_numerical_counts'


@pytest.mark.parametrize('bad',['negative','fractional','missing_day','duplicate'])
def test_invalid_count_tables_fail(tmp_path,bad):
    _,_,counts=fixture_tables()
    if bad=='negative':counts.loc[0,'n_target_cells_with_cnmf_usage']=-1
    if bad=='fractional':counts=counts.astype({'n_target_cells_with_cnmf_usage':float});counts.loc[0,'n_target_cells_with_cnmf_usage']=1.5
    if bad=='missing_day':counts['day']='d1'
    if bad=='duplicate':counts=pd.concat([counts,counts.iloc[[0]]],ignore_index=True)
    path=tmp_path/'counts.csv';counts.to_csv(path,index=False)
    with pytest.raises(ValueError):m.load_or_rebuild_target_counts(tmp_path,['d0'],tmp_path/'out',counts_csv=path)


def test_original_functions_unchanged():
    expected=json.loads((ROOT/'docs/S4_POSTCRT_SOURCE_FUNCTION_HASHES.json').read_text())
    nodes={n.name:n for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef)}
    for name,digest in expected.items():
        assert hashlib.sha256(ast.dump(nodes[name],include_attributes=False).encode()).hexdigest()==digest


def test_render_exports_actual_selected_source_values(tmp_path):
    effects,pooled,counts=fixture_tables()
    selected,_,_=m.select_effect_pairs(effects,pooled,counts)
    table,summary,plot=m.render_selected_effects(selected,tmp_path)
    assert summary.n_target_program_pairs.tolist()==[2]
    assert len(list(tmp_path.glob('*.svg')))==2
    assert len(list(tmp_path.glob('*.pdf')))==2
    assert len(list(tmp_path.glob('*.png')))==2
    saved=pd.read_csv(next(tmp_path.glob('pooled_selected_replicate_effect_pairs_*.csv')))
    assert len(saved)==2 and 'B' not in set(saved.target_name)


def test_empty_selection_is_not_top_n_fallback(tmp_path):
    effects,pooled,counts=fixture_tables()
    selected,_,report=m.select_effect_pairs(effects,pooled,counts,q_cutoff=0)
    assert selected.empty and report.selected_pairs.tolist()==[0]
    with pytest.raises(ValueError,match='No effects pass'):
        m.render_selected_effects(selected,tmp_path,q_cutoff=0)
