#!/usr/bin/env python3
"""Download verified GENCODE M25 gene annotation, preserving Ensembl versions."""
import gzip,json,re,subprocess,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
REF=ROOT/'data/references/ozone';META=ROOT/'metadata/ozone_analysis'
NAME='gencode.vM25.primary_assembly.annotation.gtf.gz'
URL='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/'+NAME
MD5='c5125258a0a2c5250ddb4c192abbf4e8'

def main():
    path=REF/NAME
    if not path.exists():
        subprocess.run(['aria2c','-c','-x','4','-s','4','--file-allocation=none',
          '--auto-file-renaming=false','--checksum=md5='+MD5,'-d',str(REF),'-o',NAME,URL],check=True)
    md5,sha=hashlib.md5(),hashlib.sha256()
    with path.open('rb')as f:
        while chunk:=f.read(8*1024*1024):md5.update(chunk);sha.update(chunk)
    assert md5.hexdigest()==MD5,'Annotation checksum mismatch'
    rows=[]
    with gzip.open(path,'rt')as f:
        for line in f:
            if line.startswith('#'):continue
            fields=line.rstrip().split('\t')
            if fields[2]!='gene':continue
            a=dict(re.findall(r'(\w+) "([^"]*)"',fields[8]))
            rows.append((a['gene_id'],a.get('gene_name',''),a.get('gene_type',''),fields[0]))
    with (REF/'genes.tsv').open('w')as f:
        f.write('gene_id\tgene_symbol\tgene_type\tchromosome\n')
        for r in rows:f.write('\t'.join(r)+'\n')
    (META/'gene_annotation_source.json').write_text(json.dumps({'url':URL,'file':str(path.relative_to(ROOT)),
      'provider_md5':MD5,'md5':md5.hexdigest(),'sha256':sha.hexdigest(),'genes':len(rows),
      'note':'Gene symbols and true gene_type from gene features, not from representative transcript biotypes.'},indent=2))
    print('Verified annotation genes:',len(rows))

if __name__=='__main__':main()
