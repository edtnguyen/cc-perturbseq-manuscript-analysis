"""Compare against the pinned executed source; never run the archived notebook."""
from pathlib import Path
import ast
import hashlib
import json

root=Path.cwd()
source_path=root/'recovery/sherlock_notebook_text/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/iPSC_EC_QQ_full_covariates_with_replicate_correlation_bundle/iPSC_EC_QQ_full_covariates_with_replicate_correlation.ipynb.md'
raw=source_path.read_bytes()
blob=hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()
assert blob=='333850b851ddd13586f5635a3adefb33ad7e151e',blob
code=raw.decode().split('## Code cell 31\n',1)[1].split('```python\n',1)[1].split('\n```',1)[0]
original={n.name:n for n in ast.parse(code).body if isinstance(n,ast.FunctionDef)}
ported={n.name:n for n in ast.parse((root/'src/cc_perturbseq/figure_s04_postcrt.py').read_text()).body if isinstance(n,ast.FunctionDef)}

def _canonical_ast(node):
    if isinstance(node,ast.AST):
        return {'node':type(node).__name__, 'fields':{
            key:_canonical_ast(value) for key,value in ast.iter_fields(node)
            if not (key=='type_params' and not value)}}
    if isinstance(node,list):
        return [_canonical_ast(value) for value in node]
    return node

def digest(node):
    return hashlib.sha256(json.dumps(_canonical_ast(node),sort_keys=True,separators=(',',':')).encode()).hexdigest()

path=root/'docs/S4_POSTCRT_SOURCE_FUNCTION_HASHES.json'
names=json.loads(path.read_text())
for name in names:
    assert _canonical_ast(original[name])==_canonical_ast(ported[name]),name
path.write_text(json.dumps({name:digest(original[name]) for name in names},indent=2)+'\n')
p=root/'tests/test_s4_postcrt.py'
s=p.read_text()
old='assert hashlib.sha256(ast.dump(nodes[name],include_attributes=False).encode()).hexdigest()==digest'
new='assert hashlib.sha256(json.dumps(_canonical_ast(nodes[name]),sort_keys=True,separators=(",",":")).encode()).hexdigest()==digest'
assert old in s
s=s.replace(old,new)
s+='''

def _canonical_ast(node):
    # Python 3.13 changed ast.dump's treatment of empty fields. Use a stable
    # field encoding; empty generic type parameters were introduced in 3.12.
    if isinstance(node,ast.AST):
        return {'node':type(node).__name__, 'fields':{
            key:_canonical_ast(value) for key,value in ast.iter_fields(node)
            if not (key=='type_params' and not value)}}
    if isinstance(node,list):
        return [_canonical_ast(value) for value in node]
    return node
'''
p.write_text(s)
p=root/'docs/S4_POSTCRT_FILTERS.md'
p.write_text(p.read_text()+'''\nSource-helper fingerprints use sorted compact JSON of AST node names and fields,
excluding empty type_params (added in Python 3.12). They were independently
computed from code cell 31 in the executed source text, Git blob
333850b851ddd13586f5635a3adefb33ad7e151e, and compared to the port before publication.
They do not depend on Python-version-specific ast.dump formatting.\n''')
print('All seven source helpers match the pinned executed source under stable AST comparison.')
