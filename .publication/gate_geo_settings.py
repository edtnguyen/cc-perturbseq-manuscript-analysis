"""Pass the configured accession only in GEO mode; do not loosen input guards."""
import ast
import copy
import hashlib
import json
from pathlib import Path

paths = ['notebooks/main/Figure_03.ipynb','notebooks/main/Figure_04.ipynb','notebooks/supplementary/Figure_S02.ipynb']
for name in paths:
    path=Path(name)
    before=json.loads(path.read_text())
    nb=copy.deepcopy(before)
    changed=0
    for old, cell in zip(before['cells'], nb['cells']):
        if cell['cell_type']!='code':continue
        text=''.join(cell['source'])
        pattern='GEO_OPTIONS = dict(mode=DATA_SOURCE, accession=GEO_ACCESSION, filename=GEO_FILENAME,'
        if pattern not in text:continue
        assert text.count(pattern)==1
        cell['source']=text.replace(pattern,'GEO_OPTIONS = dict(mode=DATA_SOURCE, accession=GEO_ACCESSION if DATA_SOURCE == "geo" else "", filename=GEO_FILENAME,').splitlines(keepends=True)
        left,right=ast.parse(text),ast.parse(''.join(cell['source']))
        node=next(n for n in right.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='GEO_OPTIONS' for t in n.targets))
        keyword=next(k for k in node.value.keywords if k.arg=='accession')
        expected=ast.parse('GEO_ACCESSION if DATA_SOURCE == "geo" else ""',mode='eval').body
        assert ast.dump(keyword.value,include_attributes=False)==ast.dump(expected,include_attributes=False)
        keyword.value=ast.Name(id='GEO_ACCESSION',ctx=ast.Load())
        assert ast.dump(left,include_attributes=False)==ast.dump(right,include_attributes=False)
        changed+=1
    assert changed==1,name
    for original,updated in zip(before['cells'],nb['cells']):
        original.pop('source',None);clone=copy.deepcopy(updated);clone.pop('source',None)
        assert original==clone,name
    import nbformat
    nbformat.validate(nbformat.from_dict(nb))
    for cell in nb['cells']:
        if cell['cell_type']=='code':compile(''.join(cell['source']),name,'exec')
    path.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')

path=Path('tests/test_study_geo_accessions.py')
s=path.read_text()
old="        'mode': 'precomputed', 'accession': accession, 'filename': None,"
assert s.count(old)==1
s=s.replace(old,"        'mode': 'precomputed', 'accession': '', 'filename': None,")
old="    assert namespace['DATA_SOURCE'] == 'precomputed'"
assert s.count(old)==1
s=s.replace(old,old+"\n    assert namespace['GEO_ACCESSION'] == accession")
path.write_text(s)
print('Prefilled accessions are forwarded only for DATA_SOURCE=geo; input guards and all saved notebook outputs retained.')
