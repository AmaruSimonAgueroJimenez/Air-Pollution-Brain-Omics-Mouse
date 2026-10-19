#!/usr/bin/env python3
"""Reproducible E-MTAB-14917 FASTQ -> Salmon -> tximport/DESeq2 workflow.

Original reads are never modified. Each successful stage is resumable. Run with
the repository's .venv-analysis/bin/python; raw QC samples 200,000 pairs/library.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, os, shutil, subprocess, sys, time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / 'data/references/ozone'
META = ROOT / 'metadata/ozone_analysis'
OUT = ROOT / 'results/ozone'
BASE = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/'
FILES = {'gencode.vM25.transcripts.fa.gz': 'a821c0dde39c48b9d2c4b48d36b0180c',
         'GRCm38.primary_assembly.genome.fa.gz': '3bc591be24b77f710b6ba5d41022fc5a'}

def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False))
    temp.replace(path)

def hashes(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open('rb') as f:
        while chunk := f.read(8*1024*1024): md5.update(chunk); sha.update(chunk)
    return {'md5': md5.hexdigest(), 'sha256': sha.hexdigest(), 'bytes': path.stat().st_size}

def run(cmd, name):
    (OUT/'logs').mkdir(parents=True, exist_ok=True)
    print(time.strftime('%Y-%m-%d %H:%M:%S'), name, flush=True)
    with (OUT/'logs'/f'{name}.log').open('a') as log:
        log.write('\nCOMMAND ' + json.dumps(list(map(str, cmd)))+'\n'); log.flush()
        subprocess.run(list(map(str, cmd)), cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)

def metadata():
    groups = defaultdict(list)
    with (ROOT/'metadata/E-MTAB-14917/E-MTAB-14917.sdrf.txt').open() as f:
        for row in csv.DictReader(f, delimiter='\t'): groups[row['Source Name']].append(row)
    samples=[]
    for name, rows in sorted(groups.items(), key=lambda x:int(x[0].split()[-1])):
        assert len(rows)==2
        r=rows[0]; run_id=r['Comment[ENA_RUN]']
        assert rows[1]['Comment[ENA_RUN]']==run_id
        paths=[ROOT/f'data/E-MTAB-14917/fastq/{run_id}_{i}.fastq.gz' for i in (1,2)]
        assert all(x.exists() and x.stat().st_size>0 for x in paths)
        samples.append({'sample_id':'S'+name.split()[-1], 'source_name':name, 'run':run_id,
          'age':'M'+r['Characteristics[age]'], 'genotype':r['Characteristics[genotype]'].replace('Apo',''),
          'exposure':'FA' if r['Characteristics[environmental stress]']=='Filtered Air' else 'O3',
          'sex':'male', 'tissue':r['Characteristics[organism part]'], 'strain':r['Characteristics[strain]'],
          'fastq1':str(paths[0].relative_to(ROOT)), 'fastq2':str(paths[1].relative_to(ROOT))})
    assert len(samples)==24
    strata=defaultdict(int)
    for s in samples: strata[(s['age'],s['genotype'],s['exposure'])]+=1
    assert len(strata)==8 and set(strata.values())=={3}, strata
    META.mkdir(parents=True, exist_ok=True)
    with (META/'samples.csv').open('w') as f:
        w=csv.DictWriter(f, fieldnames=list(samples[0])); w.writeheader();w.writerows(samples)
    return samples

def prepare():
    samples=metadata(); REF.mkdir(parents=True, exist_ok=True)
    def get(item):
        name, md5=item; path=REF/name
        if not path.exists():
            part=path.with_suffix(path.suffix+'.part')
            if shutil.which('aria2c'):
                run(['aria2c','--continue=true','--max-connection-per-server=8','--split=8',
                     '--min-split-size=8M','--file-allocation=none','--auto-file-renaming=false',
                     '--check-integrity=true','--checksum=md5='+md5,'--dir',part.parent,
                     '--out',part.name,BASE+name], 'reference_'+name)
            else:
                run(['curl','--fail','--location','--retry','5','--continue-at','-','--output',part,BASE+name], 'reference_'+name)
            part.replace(path)
        rec=hashes(path)
        if rec['md5']!=md5: raise ValueError('Reference MD5 mismatch: '+name)
        return {'file':str(path.relative_to(ROOT)), 'url':BASE+name, 'provider_md5':md5, **rec}
    with ThreadPoolExecutor(max_workers=2) as pool: sources=list(pool.map(get, FILES.items()))
    dump(META/'reference_sources.json', {'assembly':'GRCm38', 'annotation':'GENCODE M25',
       'official_release':'https://www.gencodegenes.org/mouse/release_M25.html', 'sources':sources})
    txmap=REF/'tx2gene.tsv'; genes=REF/'genes.tsv'; gentrome=REF/'gentrome.fa'; decoys=REF/'decoys.txt'
    if not (gentrome.exists() and txmap.exists() and genes.exists() and decoys.exists()):
        gene_map={}; n=0
        with gentrome.with_suffix('.fa.tmp').open('w') as out, txmap.open('w') as tx:
            tx.write('transcript_id\tgene_id\n')
            with gzip.open(REF/'gencode.vM25.transcripts.fa.gz','rt') as inp:
                for line in inp:
                    if line.startswith('>'):
                        parts=line[1:].strip().split('|'); tid,gid=parts[:2]
                        symbol=parts[5]; biotype=parts[7] if len(parts)>7 else ''
                        tx.write(f'{tid}\t{gid}\n'); gene_map[gid]=(symbol,biotype); n+=1
                        out.write('>'+tid+'\n')
                    else: out.write(line)
            with gzip.open(REF/'GRCm38.primary_assembly.genome.fa.gz','rt') as inp, decoys.open('w') as d:
                for line in inp:
                    if line.startswith('>'):
                        chrom=line[1:].split()[0]; d.write(chrom+'\n'); out.write('>'+chrom+'\n')
                    else:out.write(line)
        gentrome.with_suffix('.fa.tmp').replace(gentrome)
        with genes.open('w') as f:
            f.write('gene_id\tgene_symbol\tgene_type\n')
            for gid,(sym,typ) in gene_map.items():f.write(f'{gid}\t{sym}\t{typ}\n')
        dump(META/'reference_build.json', {'transcripts':n,'genes':len(gene_map),'decoy':'entire GRCm38 primary assembly',
             'gentrome':hashes(gentrome), 'tx2gene':hashes(txmap), 'note':'Mouse reference; human APOE alleles not quantified specifically.'})
    run([sys.executable,ROOT/'scripts/ozone_reference_annotation.py'],'reference_gene_annotation')
    return samples

def index(cpus):
    target=REF/'salmon_index'
    if (target/'ozone_index_complete.json').exists():return
    run(['salmon','index','-t',REF/'gentrome.fa','-d',REF/'decoys.txt','-i',target,
         '-k','31','-p',str(cpus),'--ramLimit','48'], 'salmon_index')
    dump(target/'ozone_index_complete.json',{'version':subprocess.check_output(['salmon','--version'],text=True).strip(),'cpus':cpus})

def qc(samples):
    (OUT/'raw_qc').mkdir(parents=True, exist_ok=True); rows=[]
    for s in samples:
        path=OUT/'raw_qc'/f"{s['sample_id']}.json"
        if not path.exists():
            run(['fastp','-i',s['fastq1'],'-I',s['fastq2'],'--reads_to_process','200000','--thread','2',
                 '--detect_adapter_for_pe','--disable_quality_filtering','--disable_length_filtering',
                 '--disable_trim_poly_g','--dont_eval_duplication','--json',path,'--html',path.with_suffix('.html')], 'raw_qc_'+s['sample_id'])
        d=json.loads(path.read_text()); b=d['summary']['before_filtering']; a=d.get('adapter_cutting',{})
        frac=a.get('adapter_trimmed_reads',0)/max(1,b['total_reads'])
        rows.append({**s,'qc_sampled_reads':b['total_reads'], 'qc_sampled_pairs':b['total_reads']//2,
           'q20_rate':b['q20_rate'],'q30_rate':b['q30_rate'],'gc_content':b['gc_content'],
           'read1_mean_length':b['read1_mean_length'],'read2_mean_length':b.get('read2_mean_length'),
           'adapter_trimmed_read_fraction':frac,'trim_for_quant':frac>=0.005})
    with (OUT/'raw_qc_summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    dump(META/'qc_decisions.json',{'sampling':'first 200000 read pairs per library; not all-read FastQC',
      'trim_criterion':'adapter bases detected by paired read overlap in >=0.5% sampled reads',
      'filtering':'no quality/length filtering; no poly-G trimming; no duplicate removal',
      'decisions':[{k:s[k] for k in ('sample_id','adapter_trimmed_read_fraction','trim_for_quant')} for s in rows]})
    return rows

def quant_one(s, cpus):
    sid=s['sample_id']; target=OUT/'quant'/sid
    if (target/'ozone_quant_complete.json').exists():return
    r1,r2=ROOT/s['fastq1'],ROOT/s['fastq2']; clean=OUT/'trimmed'/sid
    if s['trim_for_quant']:
        clean.mkdir(parents=True, exist_ok=True); r1,r2=clean/'R1.fastq.gz',clean/'R2.fastq.gz'
        if not ((clean/'complete.json').exists() and r1.exists() and r2.exists()):
            run(['fastp','-i',s['fastq1'],'-I',s['fastq2'],'-o',r1,'-O',r2,'--thread',str(min(cpus,4)),
              '--compression','1','--detect_adapter_for_pe','--disable_quality_filtering','--disable_length_filtering',
              '--disable_trim_poly_g','--dont_eval_duplication','--json',clean/'fastp.json','--html',clean/'fastp.html'], 'trim_'+sid)
            dump(clean/'complete.json',{'sample_id':sid,'original_reads_untouched':True})
    run(['salmon','quant','-i',REF/'salmon_index','-l','A','-1',r1,'-2',r2,'-p',str(cpus),
        '--seqBias','--gcBias','--no-progress','-o',target], 'salmon_'+sid)
    info=json.loads((target/'aux_info/meta_info.json').read_text())
    if not (target/'quant.sf').exists():raise ValueError('No quant.sf for '+sid)
    dump(target/'ozone_quant_complete.json',{'sample_id':sid,'completed_at_unix':time.time(),
         'salmon_reported_elapsed_seconds':info.get('total_time_seconds'),
         'raw_reads':not s['trim_for_quant'],'meta_info':info,'quant_sf':hashes(target/'quant.sf')})
    # The reports remain; only derived temporary reads are removed after success.
    if clean.exists():
        for x in [clean/'R1.fastq.gz',clean/'R2.fastq.gz']:x.unlink(missing_ok=True)

def quant(rows, cpus):
    workers=2 if cpus>=8 else 1
    per_worker=max(1,cpus//workers)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(quant_one,s,per_worker) for s in rows]
        for future in futures:future.result()

def enrichment():
    import pandas as pd
    folder=OUT/'enrichment_inputs';folder.mkdir(parents=True,exist_ok=True)
    cmd=[sys.executable,ROOT/'scripts/analyze_enrichment.py','--threads','4']
    names=[]
    for age in ['M3','M17']:
        for genotype in ['E3','E4']:
            key=age+'_'+genotype;name='ozone_'+key
            d=pd.read_csv(OUT/f'de_{key}_O3_vs_FA.csv').rename(columns={'gene_symbol':'gene'})
            path=folder/f'{key}.tsv'
            d[['gene','baseMean','log2FoldChange','stat','pvalue','padj','gene_id']].to_csv(path,sep='\t',index=False)
            cmd += ['--dataset',name+'='+str(path.relative_to(ROOT))];names.append(name)
    for name in names:cmd += ['--only',name]
    run(cmd,'enrichment')

def main():
    a=argparse.ArgumentParser();a.add_argument('--stage',choices=['prepare','index','qc','quant','de','figures','enrichment','all'],default='all')
    a.add_argument('--cpus',type=int,default=4);args=a.parse_args()
    for p in [META,OUT,OUT/'logs']:p.mkdir(parents=True,exist_ok=True)
    samples=metadata()
    dump(META/'pipeline_status.json',{'stage':args.stage,'started_at':time.strftime('%Y-%m-%dT%H:%M:%S'),'pid':os.getpid(),'complete':False})
    if args.stage in ('prepare','all'):prepare()
    if args.stage in ('qc','all'):qc(samples)
    if args.stage in ('index','all'):index(args.cpus)
    if args.stage in ('quant','all'):quant(qc(samples),args.cpus)
    if args.stage in ('de','all'):run(['Rscript',ROOT/'scripts/analyze_ozone.R'], 'deseq2')
    if args.stage in ('figures','all'):run([sys.executable,ROOT/'scripts/ozone_figures.py'], 'figures')
    if args.stage in ('enrichment','all'):enrichment()
    dump(META/'pipeline_status.json',{'stage':args.stage,'finished_at':time.strftime('%Y-%m-%dT%H:%M:%S'),'pid':os.getpid(),'complete':True})

if __name__=='__main__':main()
