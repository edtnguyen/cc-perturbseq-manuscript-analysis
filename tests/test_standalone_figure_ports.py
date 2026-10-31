"""Execute re-scripted figure notebooks with only synthetic numerical inputs.

Each execution uses a fresh process and an isolated checkout with no recovery/
archive, no old notebooks, no old figures, and no accessible Oak input files.
These tests do not certify agreement with the actual study matrices.
"""
from __future__ import annotations
import ast
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
FIGURES = [*(ROOT / 'notebooks/main').glob('*.ipynb'), *(ROOT / 'notebooks/supplementary').glob('*.ipynb')]


def test_all_seven_figure_notebooks_are_not_old_notebook_launchers():
    assert len(FIGURES) == 7
    forbidden = {'execute_notebook','patch_source_notebook','resolve_s4_source_notebook',
                 'collect_s4_source_outputs','NotebookClient','ExecutePreprocessor',
                 'run_path','run_notebook','exec','eval'}
    for path in FIGURES:
        nb = json.loads(path.read_text())
        for index,cell in enumerate(nb['cells']):
            if cell['cell_type'] != 'code':
                continue
            raw = cell['source']; source = ''.join(raw) if isinstance(raw,list) else raw
            tree = ast.parse(source)
            compile(source,f'{path}:cell-{index}','exec')
            for node in ast.walk(tree):
                if isinstance(node,ast.Call):
                    name = node.func.id if isinstance(node.func,ast.Name) else getattr(node.func,'attr','')
                    assert name not in forbidden, (path,index,name)
            assert '%run ' not in source
            assert 'sha256_file(SOURCE_NOTEBOOK' not in source
    assert not (ROOT/'src/cc_perturbseq/s4_source_runner.py').exists()
    assert not (ROOT/'src/cc_perturbseq/source_notebook_runner.py').exists()


def test_ported_function_asts_match_recorded_original_functions():
    manifest = json.loads((ROOT/'docs/STANDALONE_FIGURE_SOURCE_AUDIT.json').read_text())
    for spec in manifest.values():
        if 'ported_file' not in spec:
            continue
        nodes = {n.name:n for n in ast.parse((ROOT/spec['ported_file']).read_text()).body
                 if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
        assert spec['original_function_ast_sha256']
        for name,digest in spec['original_function_ast_sha256'].items():
            assert name in nodes
            observed = hashlib.sha256(ast.dump(nodes[name],include_attributes=False).encode()).hexdigest()
            assert observed == digest, (spec['ported_file'],name)


def test_s4_numeric_pairing_and_correlations(tmp_path):
    from cc_perturbseq.figure_s04_panels import build_effect_pair_table, summarize_effect_pairs, _qq_columns_for_one_series
    rows=[]
    for sample,multiplier in [('r1',1),('r2',2)]:
        p=tmp_path/f'{sample}.tsv'
        pd.DataFrame({'P0':[1,2,900],'P1':[-1,3,-900]},index=['A','B','NTC']).mul(multiplier).to_csv(p,sep='\t')
        rows.append({'day':'d0','biological_sample':sample,'beta_path':str(p)})
    pairs=build_effect_pair_table(pd.DataFrame(rows))
    assert len(pairs)==4
    assert set(pairs.target_name)=={'A','B'}
    np.testing.assert_allclose(pairs.beta_y,2*pairs.beta_x)
    summary=summarize_effect_pairs(pairs)
    assert np.isclose(summary.pearson_r.iloc[0],1)
    assert np.isclose(summary.spearman_r.iloc[0],1)
    qq=_qq_columns_for_one_series(pd.DataFrame({'pvalue':[0.5,0.01,0.1]}))
    np.testing.assert_allclose(qq.expected_neglog10_p,-np.log10((np.arange(1,4)-0.5)/3))
    np.testing.assert_allclose(qq.observed_neglog10_p,-np.log10([0.01,0.1,0.5]))


def _isolated_checkout(tmp_path):
    root=tmp_path/'checkout'
    root.mkdir()
    for name in ['src','config','notebooks']:
        shutil.copytree(ROOT/name,root/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc','paths.local.yaml'))
    (root/'docs').mkdir()
    (root/'docs/REPRODUCTION_PLAN.md').write_text('Synthetic isolated test checkout.\n')
    assert not (root/'recovery').exists()
    return root


def _run_notebook(root,name,config,extra_env=None,motif_files=None):
    config_path=root/'config/paths.local.yaml'
    config_path.write_text(yaml.safe_dump(config))
    env=os.environ.copy()
    env.update({'PYTHONPATH':str(root/'src'),'MPLBACKEND':'Agg','OMP_NUM_THREADS':'2',
                'OPENBLAS_NUM_THREADS':'2','NUMBA_NUM_THREADS':'2','CCPERTURBSEQ_PATHS':str(config_path)})
    env.pop('CCPERTURBSEQ_S4_RESULTS',None)
    env.pop('CCPERTURBSEQ_FIGURE_S02_ROOT',None)
    env.update(extra_env or {})
    script='''
from pathlib import Path
import sys, json
import IPython.display
IPython.display.display = lambda *a, **kw: None
'''
    if motif_files:
        script+='from cc_perturbseq import motif_sources as m\n'
        for key,path in motif_files.items():
            script+=f'm.{key}=Path({str(path)!r})\n'
    script+='''
def prevent_private_sources(event,args):
    if event == 'open' and args and isinstance(args[0],(str,bytes)):
        value=str(args[0])
        if value.startswith(('/oak/','/scratch/')) or '/recovery/' in value:
            raise AssertionError('Unexpected private/archive input: '+value)
sys.addaudithook(prevent_private_sources)
'''
    script+=f'nb=json.loads(Path("notebooks/supplementary/{name}").read_text())\n'
    script+='''
ns={'__name__':'__main__'}
for index, cell in enumerate(nb['cells']):
    if cell['cell_type']=='code':
        code=cell['source']; code=''.join(code) if isinstance(code,list) else code
        print('Executing curated cell',index,flush=True)
        exec(compile(code,f'curated-cell-{index}','exec'),ns)
print('ALL_CURATED_CELLS_EXECUTED_WITHOUT_ARCHIVES',flush=True)
'''
    run=subprocess.run([sys.executable,'-c',script],cwd=root,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=300)
    (root/'execution.log').write_text(run.stdout)
    assert run.returncode==0,run.stdout[-18000:]
    assert 'ALL_CURATED_CELLS_EXECUTED_WITHOUT_ARCHIVES' in run.stdout
    assert not (root/'recovery').exists()


def test_s4_entire_notebook_from_numeric_files_without_old_notebook(tmp_path):
    root=_isolated_checkout(tmp_path)
    model=root/'numeric_inputs';model.mkdir()
    points=[];manifest=[]
    rng=np.random.default_rng(17)
    for day in ['d0','d1','d2','d3']:
        for kind in ['observed_target','ntc_pseudotarget']:
            for p in [0.001,0.02,0.1,0.4,0.8]:
                points.append({'day':day,'pvalue_type':'empirical','series_type':kind,'pvalue':p})
        for rep in ['r1','r2']:
            path=model/'biological_replicate_effects/betas_by_replicate'/day/f'{rep}.tsv'
            path.parent.mkdir(parents=True,exist_ok=True)
            pd.DataFrame(rng.normal(size=(4,3)),index=['A','B','C','NTC'],columns=['P0','P1','P2']).to_csv(path,sep='\t')
            manifest.append({'day':day,'biological_sample':rep,'status':'written',
                             'beta_path':str(path.relative_to(model/'biological_replicate_effects'))})
    pd.DataFrame(points).to_csv(model/'all_days_QQ_points.csv',index=False)
    pd.DataFrame(manifest).to_csv(model/'biological_replicate_effects/replicate_beta_manifest.csv',index=False)
    output=root/'outputs'
    _run_notebook(root,'Figure_S04.ipynb',{'outputs':{'root':str(output)}},
                  {'CCPERTURBSEQ_S4_RESULTS':str(model)})
    b=output/'Figure_S04/panel_b';c=output/'Figure_S04/panel_c'
    assert len(list(b.glob('*.svg')))==4
    assert len(list(c.glob('*.svg')))==5
    pairs=pd.read_csv(c/'biological_replicate_DE_program_effect_pairs.csv')
    assert len(pairs)==36
    assert 'NTC' not in set(pairs.target_name)


def _s6_inputs(root):
    import h5py
    input_dir=root/'synthetic_data';input_dir.mkdir()
    out={}
    genes=[f'G{i}' for i in range(16)]
    annot=pd.DataFrame({'Program number':range(8),'Name':[f'Program {i}' for i in range(8)],
                        'Representative genes':[', '.join([genes[(2*i+j)%16] for j in range(4)]) for i in range(8)]})
    out['annotations']=input_dir/'Programs_annotations.xlsx';annot.to_excel(out['annotations'],index=False)
    summary=pd.DataFrame({'program_name':range(8),'Automatic Timepoint':np.repeat(['D0','D1','D2','D3'],2)})
    for day in ['D0','D1','D2','D3']:
        summary[f'Top 5 specific regulators (FDR<0.1) {day}']='R1, R2, R3'
    associations=[];motifs=[]
    for p in range(8):
        day=f'D{p//2}'
        for i,reg in enumerate(['R1','R2','R3']):
            associations.append({'target_name':reg,'program_name':p,'Sample':day,'approx_log2FC':(-1)**i*(p+1)/10,'adj_pval':0.001*(i+1)})
        for i,motif in enumerate(['ETS','TEAD','GATA','HAND']):
            q=10.0**(-i-1)
            motifs.append({'program':p,'day':day,'tf':motif,'qval':q,'neg_log10_qval':-np.log10(q),'AUC':0.6+i/20})
    out['workbook']=root/'data/figure_02/cNMF_50_0_2_032526_6co.xlsx';out['workbook'].parent.mkdir(parents=True)
    with pd.ExcelWriter(out['workbook']) as writer:
        summary.to_excel(writer,sheet_name='Summary',index=False)
        pd.DataFrame(associations).to_excel(writer,sheet_name='Perturbation Association 1',index=False)
    weights=np.zeros((8,16),dtype=float)
    for p in range(8):
        weights[p,2*p:(2*p+2)]=10
        weights[p,(2*p+2)%16]=1
        weights[p,(2*p+3)%16]=0.5
    out['consensus']=input_dir/'consensus.tsv'
    pd.DataFrame(weights,index=range(8),columns=genes).to_csv(out['consensus'],sep='\t')
    out['enrichment']=input_dir/'enhancer_enrichment_results_min10.tsv';pd.DataFrame(motifs).to_csv(out['enrichment'],sep='\t',index=False)
    out['cluster_map']=input_dir/'cluster_to_tf_motif_id_mapping.tsv'
    pd.DataFrame({'cluster_id':range(4),'tf_motif_id':['ETS','TEAD','GATA','HAND']}).to_csv(out['cluster_map'],sep='\t',index=False)
    out['modisco']=input_dir/'modisco_compiled.h5'
    with h5py.File(out['modisco'],'w') as handle:
        for i in range(4):
            handle.create_dataset(f'pos_patterns/{i}/contrib_scores',data=np.eye(4)[np.arange(7)%4]*0.3)
    return out


def test_s6_entire_notebook_renders_from_six_data_files_without_archive(tmp_path):
    root=_isolated_checkout(tmp_path)
    paths=_s6_inputs(root);output=root/'outputs'
    _run_notebook(root,'Figure_S06.ipynb',
        {'outputs':{'root':str(output)},'ccperturbseq':{'cnmf_spectra':str(paths['consensus']),'program_annotations':str(paths['annotations'])}},
        motif_files={'CANONICAL_ENHANCER_ENRICHMENT_PATH':paths['enrichment'],
                     'CANONICAL_CLUSTER_MOTIF_MAP_PATH':paths['cluster_map'],
                     'CANONICAL_MODISCO_COMPILED_H5_PATH':paths['modisco']})
    out=output/'Figure_S06/panel_all'
    assert len(list(out.glob('*.svg'))) == 5
    tables=out/'_program_annotation_source'
    assert len(list(tables.glob('*_motif_enrichment_matrix.csv'))) == 4
    for p in tables.glob('*_motif_logo_map.csv'):
        frame=pd.read_csv(p)
        assert set(frame.status)=={'vector_cwm'}
    day0=pd.read_csv(tables/'D0_regulator_effect_matrix.csv',index_col=0)
    np.testing.assert_allclose(day0['R1'].to_numpy(),[0.1,0.2])
    with open(tables/'shared_color_scales_D0_D3.json') as handle:
        scales=json.load(handle)
    assert scales['regulator_vmin']==-scales['regulator_vmax']


def _s2_inputs(root):
    from scipy import sparse
    from scipy.io import mmwrite
    rng=np.random.default_rng(42)
    source=root/'mixed_input';dge=source/'splitpipe_out/WT/sub1/all-sample/DGE_filtered';dge.mkdir(parents=True)
    genes=['GATA1','KLF1','HBB','HBA1','HBA2','SOX2','NANOG','POU5F1','LIN28A','MT-CO1']+[f'GENE{i}' for i in range(990)]
    n=96;matrix=rng.poisson(0.3,size=(n,len(genes))).astype(np.float32)
    matrix[:48,:5]+=20
    matrix[48:,5:9]+=25
    # Match total marker burden between lineages, and provide a broad moderate-
    # expression lineage signal that passes the source HVG mean cutoffs. Nine
    # very highly expressed marker genes alone are intentionally not the HVGs.
    matrix[:48,10:210] += rng.poisson(1.0,size=(48,200))
    matrix[48:,210:410] += rng.poisson(1.0,size=(48,200))
    names=[f'cell{i}' for i in range(n)];samples=np.tile(np.repeat(['CS1','DC1','Dual','no_guide'],12),2)
    mmwrite(dge/'count_matrix.mtx',sparse.csr_matrix(matrix))
    pd.DataFrame({'gene_name':genes}).to_csv(dge/'all_genes.csv',index=False)
    pd.DataFrame({'bc_wells':names,'sample':samples}).to_csv(dge/'cell_metadata.csv',index=False)
    guide_path=source/'splitpipe_out/CRISPR_combined/grnaxcell_count.tsv';guide_path.parent.mkdir(parents=True)
    rows=[]
    for i,name in enumerate(names):
        for gid in ['K562_A','K562_B','iPSC_A','iPSC_B']:
            strong=(i<48)==gid.startswith('K562')
            rows.append({'gene':'gRNA_'+gid,'cell':name+'__s1','count':int(rng.integers(5,50) if strong else rng.integers(1,5))})
    pd.DataFrame(rows).to_csv(guide_path,sep='\t',index=False)
    reference=source/'splitpipe_out/CRISPR/sub1/all-sample/guide_RNAs_unfiltered/all_guides.csv';reference.parent.mkdir(parents=True)
    pd.DataFrame({'gene_name':['K562_A','K562_B','iPSC_A','iPSC_B']}).to_csv(reference,index=False)
    for primer in ['CS1','DC1','Dual','no_guide']:
        path=source/f'splitpipe_out/WT/sub1/{primer}/report/tscp_counts_subsampled.csv';path.parent.mkdir(parents=True)
        path.write_text(',hg38\n10,5\n20,9\n40,13\n80,17\n')
    return source


def test_s2_entire_notebook_uses_counts_and_no_original_notebook(tmp_path):
    root=_isolated_checkout(tmp_path);source=_s2_inputs(root);output=root/'outputs'
    _run_notebook(root,'Figure_S02.ipynb',{'outputs':{'root':str(output)}},
                  {'CCPERTURBSEQ_FIGURE_S02_ROOT':str(source)})
    for panel in ['i','j','l','m','n']:
        folder=output/f'Figure_S02/panel_{panel}'
        assert list(folder.glob('*.svg')),panel
        assert list(folder.glob('*.png')),panel
    assert len(list((output/'Figure_S02/panel_n').glob('*.svg')))==3
    metrics=pd.read_csv(output/'Figure_S02/Figure_S02_all_cell_metrics.csv.gz',index_col=0)
    assert len(metrics)==96
    assert set(metrics.cell_type)=={'K562','hiPSC'}


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


def test_s6_logo_rotation_stays_at_same_physical_position_across_dpi():
    import matplotlib.pyplot as plt
    from cc_perturbseq import figure_s06_annotations as s6
    fig=plt.figure(figsize=(8,6))
    ax=fig.add_axes([0.5,0.5,0.2,0.04])
    matrix=pd.DataFrame(np.eye(4)[[0,1,2,3,0,1]],columns=list('ACGT'))
    s6._draw_rotated_vector_logo(fig,ax,'TEST',{'TEST':{'type':'cwm','cwm_df':matrix}})
    assert ax.patches
    bounds=[]
    for dpi in (72,100,600):
        fig.set_dpi(dpi);fig.canvas.draw()
        bounds.append(np.asarray([p.get_window_extent().extents for p in ax.patches])/dpi)
    np.testing.assert_allclose(bounds[0],bounds[1],atol=1e-9)
    np.testing.assert_allclose(bounds[0],bounds[2],atol=1e-9)
    plt.close(fig)


def test_documented_export_adaptations_match_the_reviewed_port():
    audit=json.loads((ROOT/'docs/STANDALONE_FIGURE_SOURCE_AUDIT.json').read_text())
    spec=audit['Figure_S06']
    nodes={n.name:n for n in ast.parse((ROOT/spec['ported_file']).read_text()).body if isinstance(n,ast.FunctionDef)}
    for name,change in spec['documented_function_adaptations'].items():
        actual=hashlib.sha256(ast.dump(nodes[name],include_attributes=False).encode()).hexdigest()
        assert actual==change['ported_ast_sha256']
        assert change['source_ast_sha256'] != actual
