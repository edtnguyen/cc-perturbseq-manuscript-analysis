"""Run notebook plotting routes with private paths and external network blocked."""
from pathlib import Path
import ast
import json
import os
import shutil
import subprocess
import sys
import pytest
from test_public_notebooks import isolated
from test_public_adapters import make_mixed_geo
from geo_plot_fixtures import make_hippo_geo
from cc_perturbseq.geo_plot_data import build_plot_bundle

ROOT=Path(__file__).resolve().parents[1]


def execute_script(root,script,args=()):
    env=os.environ.copy();env.update(PYTHONPATH=str(root/'src'),MPLBACKEND='Agg',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
    for key in ('CCPERTURBSEQ_PUBLIC_PATHS','CCPERTURBSEQ_PATHS','IGVF_ACCESS_KEY','IGVF_SECRET_ACCESS_KEY'):
        env.pop(key,None)
    run=subprocess.run([sys.executable,'-c',script,*map(str,args)],cwd=root,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    (root/'test-execution.log').write_text(run.stdout)
    assert run.returncode==0,run.stdout[-16000:]

GUARD='''
from pathlib import Path
import sys, json, socket
import IPython.display
IPython.display.display=lambda *a,**k:None
def block_private(event,args):
    if event=='open' and args and isinstance(args[0],(str,bytes)):
        value=str(args[0])
        if value.startswith(('/oak/','/scratch/users/','/home/users/')) or '/recovery/' in value:
            raise AssertionError('Private source dependency: '+value)
sys.addaudithook(block_private)
socket.socket.connect=lambda *a,**k: (_ for _ in ()).throw(AssertionError('Unexpected external network'))
'''


def test_entire_s2_geo_mode_with_mock_ncbi_and_real_h5mu(tmp_path):
    root=isolated(tmp_path)
    raw=root/'served/iPSC_K562_mixed_integrated.h5mu';make_mixed_geo(raw)
    build_plot_bundle(raw,root/'data/reproduction/geo_plots/ipsc_k562_mix','ipsc_k562_mix')
    path=root/'notebooks/supplementary/Figure_S02.ipynb';nb=json.loads(path.read_text())
    for c in nb['cells']:
        if c['cell_type']=='code':
            text=''.join(c['source']).replace('DATA_SOURCE = "precomputed"','DATA_SOURCE = "geo"').replace('GEO_ACCESSION = ""','GEO_ACCESSION = "GSE12345"')
            c['source']=text.splitlines(keepends=True)
    path.write_text(json.dumps(nb))
    script=GUARD+'''
import io
from src.cc_perturbseq.geo_download import GeoClient
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12345/suppl/GSE12345_iPSC_K562_mixed_integrated.h5mu'
class Response(io.BytesIO):
    def __init__(self,raw):
        super().__init__(raw);self.headers={'Content-Length':str(len(raw))}
def serve(self,requested):
    if 'acc.cgi?' in requested:
        return Response(('^SERIES = GSE12345\\n!Series_supplementary_file = '+url+'\\n').encode())
    assert requested==url
    return Response(Path('served/iPSC_K562_mixed_integrated.h5mu').read_bytes())
GeoClient._open=serve
nb=json.loads(Path('notebooks/supplementary/Figure_S02.ipynb').read_text())
ns={'__name__':'__main__'}
for i,c in enumerate(nb['cells']):
    if c['cell_type']=='code':
        print('EXECUTING',i,flush=True)
        exec(compile(''.join(c['source']),f'S2:{i}','exec'),ns)
assert ns['ctx']['public_inputs']['geo_plotting']['published_comparison']['status']=='MATCH'
assert len(list(Path('data/geo/downloads').rglob('*.h5mu')))==1
print('ALL_S2_CELLS_GEO_ROUTE_PASSED')
'''
    execute_script(root,script)
    for panel in ('i','j','l','m','n'):
        assert list((root/f'results/Figure_S02/panel_{panel}').glob('*.svg'))


@pytest.mark.parametrize('figure,cell_index,panel',[('Figure_03',10,'f'),('Figure_04',4,'c')])
def test_actual_hippo_notebook_panel_cells_from_tables_without_geo(tmp_path,figure,cell_index,panel):
    root=isolated(tmp_path);raw=root/'private/hippo.h5mu';make_hippo_geo(raw)
    build_plot_bundle(raw,root/'data/reproduction/geo_plots/hippo_timecourse','hippo_timecourse')
    shutil.rmtree(root/'private')
    script=GUARD+'''
import ast
import pandas as pd
from src.cc_perturbseq.figure_notebook import panel_directory,write_receipt
from src.cc_perturbseq.geo_plot_data import prepare_hippo_dotplots
from src.cc_perturbseq.geo_submission import _apply_source_style
from src.cc_perturbseq.figure4_panels import FIGURE4C_CONDITION_ORDER,FIGURE4C_PROGRAMS
from IPython.display import display
_apply_source_style()
name,index,panel=sys.argv[1],int(sys.argv[2]),sys.argv[3]
nb=json.loads((Path('notebooks/main')/(name+'.ipynb')).read_text())
table=pd.DataFrame([{'panel':panel}]);table.attrs['figure_id']=name
ctx={'root':Path.cwd(),'output_dir':Path.cwd()/'results'/name,'panels':table,'git_commit':'fixture', 'config':Path('config/paths.example.yaml')}
ctx['output_dir'].mkdir(parents=True)
PLOT_DATA_DIR=Path.cwd()/'data/reproduction/geo_plots/hippo_timecourse'
DATA_SOURCE='precomputed'
# Execute the notebook's actual export function and entire Hippo panel cell.
setup=ast.parse(''.join(nb['cells'][2]['source']))
func=next(n for n in setup.body if isinstance(n,ast.FunctionDef) and n.name=='export_scanpy_dotplot')
exec(compile(ast.Module(body=[func],type_ignores=[]),'export-cell','exec'))
exec(compile(''.join(nb['cells'][index]['source']),'actual-Hippo-panel-cell','exec'))
assert len(list((ctx['output_dir']/('panel_'+panel)).glob('*.svg')))==2
print('HIPPO_NOTEBOOK_PANEL_COMPLETED_WITHOUT_RAW_INPUTS')
'''
    execute_script(root,script,(figure,cell_index,panel))
    assert not (root/'data/geo').exists()


def test_notebook_defaults_and_runtime_options_are_explicit():
    for path in (ROOT/'notebooks/main/Figure_03.ipynb',ROOT/'notebooks/main/Figure_04.ipynb',ROOT/'notebooks/supplementary/Figure_S02.ipynb'):
        n=json.loads(path.read_text());text='\n'.join(''.join(c['source']) for c in n['cells'] if c['cell_type']=='code')
        assert 'DATA_SOURCE = "precomputed"' in text
        assert 'GEO_ACCESSION = ""' in text
        assert 'initialize_figure(FIGURE, geo_options=GEO_OPTIONS)' in text
        assert 'load_mixed_geo(' not in text and 'extract_hippo_geo(' not in text
        assert 'prepare_figure3f_dotplots(' not in text and 'prepare_figure4c_dotplots(' not in text
