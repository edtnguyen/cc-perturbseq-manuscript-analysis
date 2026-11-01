"""Apply the reviewed, checksum-guarded public input patch to the main baseline."""
from pathlib import Path
import base64
import hashlib
import json
import lzma
import subprocess
import tempfile

BASE = 'ebdc6a74fc7ed5296c671e249aa78eda5772f4d4'

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def patched(original, patch):
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d/'cell').write_bytes(original)
        # Rebase the unified diff paths to an isolated file.
        lines = patch.splitlines(keepends=True)
        lines[0] = '--- a/cell\n'
        lines[1] = '+++ b/cell\n'
        (d/'change.patch').write_text(''.join(lines))
        subprocess.run(['patch','--batch','--fuzz=0','-p1','-i','change.patch'],cwd=d,check=True)
        return (d/'cell').read_bytes()

def preserved(nb):
    copy = json.loads(json.dumps(nb))
    for cell in copy['cells']:
        cell['source'] = None
    return sha(json.dumps(copy,sort_keys=True,ensure_ascii=False).encode())

def apply(payload):
    observed = subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    if observed != BASE or payload['baseline'] != BASE:
        raise RuntimeError('Baseline changed; refuse to overwrite concurrent edits.')
    for spec in payload['files']:
        path=Path(spec['path'])
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('Invalid payload path')
        if 'patch' in spec:
            before=path.read_bytes()
            assert sha(before)==spec['old_sha256'],spec['path']
            after=patched(before,spec['patch'])
            assert sha(after)==spec['new_sha256'],spec['path']
        else:
            assert not path.exists(),spec['path']
            after=spec['new'].encode()
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(after)
    reports=[]
    import nbformat
    for spec in payload['notebooks']:
        path=Path(spec['path']);nb=json.loads(path.read_text())
        assert len(nb['cells'])==spec['cell_count'],str(path)
        before=preserved(nb)
        for change in spec['cells']:
            cell=nb['cells'][change['index']]
            assert cell['cell_type']==change['type']
            if change['id'] is not None:assert cell.get('id')==change['id']
            text=''.join(cell['source'])
            assert sha(text.encode())==change['old_sha256'],(str(path),change['index'])
            new=patched(text.encode(),change['patch'])
            assert sha(new)==change['new_sha256']
            cell['source']=new.decode().splitlines(keepends=True)
        assert preserved(nb)==before,str(path)
        nbformat.validate(nbformat.from_dict(nb))
        for i,cell in enumerate(nb['cells']):
            if cell['cell_type']=='code':compile(''.join(cell['source']),str(path)+':'+str(i),'exec')
        path.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')
        reread=json.loads(path.read_text())
        assert preserved(reread)==before,str(path)
        reports.append({'path':str(path),'unchanged_non_source_sha256':before,
                        'outputs_preserved':sum(len(c.get('outputs',[])) for c in nb['cells']),
                        'changed_source_cells':len(spec['cells'])})
    Path('/tmp/public-output-preservation.json').write_text(json.dumps(reports,indent=2)+'\n')
    Path('/tmp/public-changed-paths.json').write_text(json.dumps([s['path'] for s in payload['files']]+[s['path'] for s in payload['notebooks']]))
    print(json.dumps(reports,indent=2))

if __name__=='__main__':
    import sys
    apply(json.loads(Path(sys.argv[1]).read_text()))
