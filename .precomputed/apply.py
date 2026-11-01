from pathlib import Path
import hashlib, json, subprocess, tempfile

BASE='55fe778b1fa4921e68baf2297b1d3f12f0555fd7'
def sha(data): return hashlib.sha256(data).hexdigest()
def patch_bytes(original, patch):
    with tempfile.TemporaryDirectory() as folder:
        p=Path(folder); (p/'file').write_bytes(original);(p/'change.patch').write_text(patch)
        subprocess.run(['patch','--batch','--fuzz=0','-p1','-i','change.patch'],cwd=p,check=True)
        return (p/'file').read_bytes()
def preserved(nb):
    obj=json.loads(json.dumps(nb))
    for c in obj['cells']:c['source']=None
    return sha(json.dumps(obj,sort_keys=True,ensure_ascii=False).encode())
def apply(payload):
    assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==BASE
    assert payload['baseline']==BASE
    original_notebooks={}
    for folder in ('notebooks/main','notebooks/supplementary'):
        for p in Path(folder).glob('*.ipynb'):
            original_notebooks[str(p)]=preserved(json.loads(p.read_text()))
    for spec in payload['files']:
        p=Path(spec['path']);assert not p.is_absolute() and '..' not in p.parts
        if 'patch' in spec:
            old=p.read_bytes();assert sha(old)==spec['old_sha256'],str(p)
            new=patch_bytes(old,spec['patch']);assert sha(new)==spec['new_sha256'],str(p)
        else:
            assert not p.exists(),str(p)
            new=spec['new'].encode()
        p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(new)
    for spec in payload['notebooks']:
        p=Path(spec['path']);nb=json.loads(p.read_text())
        assert len(nb['cells'])==spec['cell_count']
        for change in spec['cells']:
            cell=nb['cells'][change['index']];old=''.join(cell['source']).encode()
            assert sha(old)==change['old_sha256'],(str(p),change['index'])
            new=patch_bytes(old,change['patch']);assert sha(new)==change['new_sha256']
            cell['source']=new.decode().splitlines(keepends=True)
        assert preserved(nb)==original_notebooks[str(p)],str(p)
        p.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')
    import nbformat
    report=[]
    for path,digest in original_notebooks.items():
        nb=json.loads(Path(path).read_text());assert preserved(nb)==digest,path
        nbformat.validate(nbformat.from_dict(nb))
        for i,c in enumerate(nb['cells']):
            if c['cell_type']=='code':compile(''.join(c['source']),f'{path}:{i}','exec')
        report.append({'path':path,'unchanged_non_source_sha256':digest,
            'outputs_preserved':sum(len(c.get('outputs',[])) for c in nb['cells'])})
    for folder in ('src','scripts','tests'):
        for p in Path(folder).rglob('*.py'):compile(p.read_text(),str(p),'exec')
    Path('/tmp/geo-precomputed-output-preservation.json').write_text(json.dumps(report,indent=2)+'\n')
    Path('/tmp/geo-precomputed-changed-paths.json').write_text(json.dumps([s['path'] for s in payload['files']]+[s['path'] for s in payload['notebooks']]))
    print(json.dumps(report,indent=2))
if __name__=='__main__':
    import sys
    apply(json.loads(Path(sys.argv[1]).read_text()))
