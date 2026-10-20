#!/usr/bin/env python3
"""Rank unchanged regional DE results and annotate gene products from frozen official records.

Annotations describe gene products/functions; no protein abundance, cell-specific exposure
response, isoform activity, causal mechanism or biological replication is inferred.
"""
from __future__ import annotations
import argparse,hashlib,json,urllib.request,urllib.parse,gzip,io
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'svg.hashsalt':'regional-priority-genes-v1','svg.fonttype':'none','pdf.fonttype':42,'font.family':'DejaVu Sans','savefig.facecolor':'white'})
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/gene_priorities'; META=ROOT/'metadata/gene_priorities'
GROUPS={
'protein_homeostasis':('Protein folding / ER homeostasis','Plegamiento / homeostasis del RE'),
'signaling':('Receptors / intracellular signaling','Receptores / señalización intracelular'),
'transcription':('Transcription / circadian regulation','Transcripción / regulación circadiana'),
'extracellular':('Extracellular matrix / adhesion','Matriz extracelular / adhesión'),
'metabolism':('Metabolism / biochemical processing','Metabolismo / procesamiento bioquímico'),
'cytoskeleton':('Cytoskeleton / cell organization','Citoesqueleto / organización celular'),
'unknown':('Incompletely characterized function','Función insuficientemente caracterizada'),
'noncoding':('Noncoding RNA / pseudogene','ARN no codificante / pseudogén')}
# Author-curated short paraphrases of the cited NCBI/UniProt mouse records.
# None describes a function experimentally measured in the pollution datasets.
CURATED={
'Hipk3':('Cinasa HIPK3, que interactúa con homeodominios','Serine/threonine kinase regulating phosphorylation and transcriptional signaling.','Cinasa de serina/treonina que regula fosforilación y señalización transcripcional.','signaling'),
'Ksr2':('Supresor de cinasa Ras 2, proteína andamio KSR2','MAPK scaffold coordinating RAF–MEK signaling; reviewed annotation reports very low in-vitro kinase activity, not an assumed active kinase in these samples.','Andamio de MAPK que coordina RAF–MEK; la anotación revisada describe actividad cinasa in vitro muy baja, sin demostrar una cinasa activa en estas muestras.','signaling'),
'Tll1':('Proteína tolloid-like 1, metaloproteasa TLL1','Extracellular protease that processes protein precursors, including matrix-associated substrates.','Proteasa extracelular que procesa precursores proteicos, incluidos sustratos asociados con la matriz.','extracellular'),
'Pcnx':('Proteína pecanex-like 1 (PCNX1)','Membrane-associated pecanex protein; the retrieved mouse records do not establish a detailed molecular function.','Proteína pecanex asociada a membrana; los registros de ratón consultados no establecen una función molecular detallada.','unknown'),
'Fzd3':('Receptor Frizzled-3','Wnt receptor with annotated roles in axon growth, guidance and planar cell polarity.','Receptor Wnt con funciones anotadas en crecimiento y guía axonal y polaridad celular planar.','signaling'),
'Gm13423':('Sin producto proteico funcional confirmado: pseudogén de Npm1','Annotated Npm1 pseudogene; its transcript signal does not establish production of functional NPM1 protein.','Pseudogén anotado de Npm1; su señal de RNA no demuestra producción de proteína NPM1 funcional.','noncoding'),
'B230312C02Rik':('ARN no codificante; sin producto proteico anotado','Noncoding RNA locus; no established molecular function is provided by the retrieved gene record.','Locus de ARN no codificante; el registro consultado no aporta una función molecular establecida.','noncoding'),
'Fstl4':('Proteína relacionada con folistatina 4','BDNF-binding regulator annotated in restraint of BDNF signaling and dendritic-spine development.','Regulador que une BDNF, anotado en la limitación de su señalización y del desarrollo de espinas dendríticas.','extracellular'),
'Ptar1':('Proteína con repeticiones de subunidad alfa de preniltransferasa 1','Predicted prenyltransferase-related component, with a predicted contribution to Rab geranylgeranylation.','Componente relacionado con preniltransferasas, con contribución predicha a la geranilgeranilación de Rab.','metabolism'),
'Mmp16':('Metaloproteinasa de matriz 16 (MMP16)','Zinc-dependent protease involved in extracellular-matrix substrate cleavage and remodeling.','Proteasa dependiente de zinc implicada en corte de sustratos y remodelación de matriz extracelular.','extracellular'),
'Lmbrd2':('Proteína LMBRD2 asociada con receptores acoplados a proteína G','Predicted regulator of adrenergic and other G-protein-coupled receptor signaling, largely supported by similarity.','Regulador predicho de señalización adrenérgica y de otros receptores acoplados a proteína G, principalmente por similitud.','signaling'),
'Hspa5':('Chaperona del retículo endoplásmico BiP / GRP78','ER chaperone involved in protein folding and unfolded-protein-response regulation.','Chaperona del RE implicada en plegamiento proteico y regulación de la respuesta a proteínas mal plegadas.','protein_homeostasis'),
'Pdia4':('Proteína disulfuro-isomerasa A4 (ERp72)','ER-associated disulfide isomerase participating in protein folding and chaperone complexes.','Isomerasa de disulfuros asociada al RE que participa en plegamiento y complejos de chaperonas.','protein_homeostasis'),
'Xbp1':('Proteína de unión a caja X 1 (XBP1)','ER-stress transcription factor; gene-level RNA does not establish the abundance of the spliced XBP1s isoform.','Factor transcripcional de estrés del RE; el RNA génico no establece la abundancia de la isoforma empalmada XBP1s.','protein_homeostasis'),
'Fosb':('Proteína FosB, subunidad del factor AP-1','AP-1 transcription-factor subunit; gene-level RNA does not distinguish FosB from ΔFosB isoform effects.','Subunidad del factor transcripcional AP-1; el RNA génico no distingue efectos de las isoformas FosB y ΔFosB.','transcription'),
'Zbtb40':('Proteína ZBTB40 con dedos de zinc y dominio BTB','Predicted DNA-binding transcription regulator with an annotated DNA-damage-response association.','Regulador transcripcional de unión a DNA predicho, con asociación anotada a respuesta al daño del DNA.','transcription'),
'Nr1d1':('Receptor nuclear NR1D1 / REV-ERBα','Nuclear transcriptional repressor coordinating circadian-clock and metabolic gene programs.','Represor transcripcional nuclear que coordina programas del reloj circadiano y genes metabólicos.','transcription'),
'Homer1':('Proteína andamio postsináptica Homer-1','Postsynaptic scaffold linking receptor and calcium signaling; gene-level counts do not isolate Homer1a.','Andamio postsináptico que conecta receptores y señalización de calcio; los conteos génicos no aíslan Homer1a.','signaling'),
'Pdia6':('Proteína disulfuro-isomerasa A6','Disulfide-processing protein with chaperone-related roles and negative regulation of ER-stress signaling.','Proteína de procesamiento de disulfuros con funciones de chaperona y regulación negativa de señalización de estrés del RE.','protein_homeostasis'),
'Egr3':('Proteína de respuesta temprana al crecimiento 3','Sequence-specific transcription regulator annotated in neural development and neuromuscular signaling.','Regulador transcripcional específico de secuencia, anotado en desarrollo neural y señalización neuromuscular.','transcription'),
'Ciart':('Represor transcripcional asociado al ritmo circadiano, CIART / CHRONO','Transcriptional repressor that contributes to the negative regulatory arm of the circadian clock.','Represor transcripcional que contribuye al componente regulador negativo del reloj circadiano.','transcription'),
'Ak4':('Adenilato-cinasa 4 mitocondrial','Mitochondrial nucleotide/energy-homeostasis protein; catalytic activity has assay-dependent reports and is not measured here.','Proteína mitocondrial asociada con nucleótidos y homeostasis energética; la actividad catalítica presenta resultados dependientes del ensayo y no se mide aquí.','metabolism'),
'Sgsm1':('Modulador de señalización de proteínas G pequeñas 1','Rab-associated regulator that modulates small-GTPase signaling and vesicle-associated processes.','Regulador asociado a Rab que modula señalización de GTPasas pequeñas y procesos vinculados a vesículas.','signaling'),
'Creld2':('Proteína disulfuro-isomerasa CRELD2','Probable disulfide isomerase associated with ER protein homeostasis; several functional assignments remain qualified.','Probable isomerasa de disulfuros asociada a homeostasis proteica del RE; varias asignaciones funcionales siguen siendo provisionales.','protein_homeostasis'),
'Fmnl1':('Proteína formin-like 1','Actin-cytoskeleton regulator involved in cell shape and migration-related organization.','Regulador del citoesqueleto de actina implicado en forma celular y organización relacionada con migración.','cytoskeleton'),
'Mlec':('Malectina','Carbohydrate-binding ER protein with a proposed role in early N-glycosylation, supported partly by similarity.','Proteína del RE que une carbohidratos, con una función propuesta en N-glicosilación temprana, parcialmente por similitud.','protein_homeostasis'),
'P4ha1':('Subunidad alfa 1 de prolil-4-hidroxilasa','Enzyme subunit that hydroxylates proline residues during collagen protein maturation.','Subunidad enzimática que hidroxila residuos de prolina durante la maduración de proteínas de colágeno.','metabolism'),
'Lyrm9':('Proteína 9 con motivo LYR','LYR-motif protein; a detailed molecular role is not established by the retrieved mouse records.','Proteína con motivo LYR; los registros de ratón consultados no establecen una función molecular detallada.','unknown'),
'Manf':('Factor neurotrófico derivado de astrocitos mesencefálicos, MANF','ER-associated factor involved in proteostasis and stress-associated survival; it also has secreted neurotrophic roles.','Factor asociado al RE implicado en proteostasis y supervivencia bajo estrés, con funciones neurotróficas extracelulares.','protein_homeostasis'),
'Cbfa2t3':('Correpresor transcripcional CBFA2T3 / ETO2','Transcriptional corepressor that recruits regulatory complexes and chromatin-modifying machinery.','Correpresor transcripcional que recluta complejos reguladores y maquinaria modificadora de cromatina.','transcription'),
'Zkscan2':('Proteína ZKSCAN2 con dedos de zinc y dominios KRAB/SCAN','Predicted sequence-specific DNA-binding regulator; detailed mouse mechanism is not established in this dataset.','Regulador predicho de unión específica a DNA; estos datos no establecen un mecanismo detallado en ratón.','transcription'),
'S1pr1':('Receptor 1 de esfingosina-1-fosfato','Sphingosine-1-phosphate GPCR involved in vascular and immune-cell trafficking biology; RNA does not identify responding cell types.','Receptor acoplado a proteína G de esfingosina-1-fosfato con funciones vasculares e inmunitarias; el RNA no identifica células respondedoras.','signaling'),
'C730002L08Rik':('ARN largo no codificante; sin producto proteico anotado','Noncoding RNA locus; no established molecular function is provided by the retrieved gene record.','Locus de ARN no codificante; el registro consultado no aporta una función molecular establecida.','noncoding'),
'Fosl2':('Antígeno Fos-like 2 / subunidad AP-1 FOSL2','AP-1-family transcription regulator involved in stimulus-responsive gene regulation.','Regulador transcripcional de la familia AP-1 implicado en regulación génica sensible a estímulos.','transcription'),
'Col19a1':('Cadena alfa 1 de colágeno XIX','Extracellular-matrix collagen chain with structural and matrix-organization roles.','Cadena de colágeno de matriz extracelular con funciones estructurales y de organización de matriz.','extracellular'),
'Pcdh8':('Protocadherina 8','Transmembrane protocadherin involved in cell-contact signaling and regulation of cadherin-mediated adhesion.','Protocadherina transmembrana implicada en señalización de contacto celular y regulación de adhesión por cadherinas.','extracellular'),
'Iqgap2':('Proteína IQGAP2 de interacción con GTPasas y calmodulina','Calmodulin- and small-GTPase-binding protein associated with cytoskeletal organization; it is not assumed to activate GTPase catalysis.','Proteína de unión a calmodulina y GTPasas pequeñas asociada al citoesqueleto; no se presupone activación catalítica de GTPasas.','cytoskeleton'),
'Usp43':('Hidrolasa 43 del extremo carboxilo de ubiquitina','Predicted deubiquitinating enzyme involved in ubiquitin processing, with mouse function partly inferred by similarity.','Enzima desubiquitinante predicha implicada en procesamiento de ubiquitina, con función en ratón parcialmente inferida por similitud.','metabolism'),
'Rorb':('Receptor nuclear huérfano relacionado con RAR, beta','DNA-binding nuclear receptor that regulates transcription, including circadian gene programs.','Receptor nuclear de unión a DNA que regula transcripción, incluidos programas circadianos.','transcription'),
'Nr4a3':('Receptor nuclear NR4A3 / NOR1','Nuclear transcription regulator annotated in neuronal survival, axonal development and stimulus-responsive transcription.','Regulador transcripcional nuclear anotado en supervivencia neuronal, desarrollo axonal y transcripción sensible a estímulos.','transcription'),
'Egr2':('Proteína EGR2 / Krox20 de respuesta temprana al crecimiento','Sequence-specific transcription factor with neural patterning and Schwann-cell myelination roles; some SUMO-ligase functions are inferred.','Factor transcripcional específico de secuencia con funciones de patrón neural y mielinización de Schwann; algunas funciones SUMO-ligasa son inferidas.','transcription'),
'Cpne9':('Copina-9','Probable calcium-dependent phospholipid-binding protein; its signaling role remains annotation-based here.','Probable proteína de unión a fosfolípidos dependiente de calcio; aquí su función de señalización procede de anotaciones.','signaling'),
'Hgf':('Factor de crecimiento de hepatocitos (HGF)','Secreted ligand of MET that regulates growth, motility and morphogenesis across several tissues.','Ligando secretado de MET que regula crecimiento, motilidad y morfogénesis en diversos tejidos.','signaling'),
'Bcl11b':('Factor transcripcional BCL11B / CTIP2','Chromatin-associated transcription regulator involved in cell differentiation, including neuronal differentiation.','Regulador transcripcional asociado a cromatina e implicado en diferenciación celular, incluida la neuronal.','transcription'),
'Car12':('Anhidrasa carbónica 12','Membrane-associated enzyme catalyzing reversible carbon-dioxide hydration and bicarbonate production.','Enzima asociada a membrana que cataliza hidratación reversible de dióxido de carbono y producción de bicarbonato.','metabolism'),
'Crym':('Ketimina-reductasa mu-cristalina','Ketimine reductase with thyroid-hormone binding and annotated roles in hormone metabolism.','Ketimina-reductasa con unión a hormona tiroidea y funciones anotadas en metabolismo hormonal.','metabolism')}

def dump(path,data):
 path.parent.mkdir(parents=True,exist_ok=True)
 path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def raw_sources(refresh=False):
 """Frozen selected records are default; refresh explicitly retrieves official public sources."""
 paths=['ncbi_gene_info_selected.json','ncbi_gene_summaries.json','uniprot_mouse_candidates.json']
 if refresh or not all((META/p).exists() for p in paths):
  groups=selections()[0];want=set(groups.gene)
  url='https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz'
  cache=META/'private_inputs/Mus_musculus.gene_info.gz';cache.parent.mkdir(parents=True,exist_ok=True)
  raw=urllib.request.urlopen(url,timeout=120).read();cache.write_bytes(raw)
  t=pd.read_csv(io.BytesIO(gzip.decompress(raw)),sep='\t',dtype=str,keep_default_na=False);rows=[]
  for gene in sorted(want):
   exact=t[t.Symbol.eq(gene)];m=exact if len(exact)==1 else t[t.Synonyms.map(lambda s:gene in s.split('|'))]
   if len(m)!=1:raise ValueError(f'Nonunique official identity: {gene}')
   rows.append({**m.iloc[0].to_dict(),'source_symbol':gene,'matching':'exact' if len(exact)==1 else 'unique_official_alias'})
  dump(META/paths[0],rows);dump(META/'ncbi_gene_info_source.json',{'url':url,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'selected_records':len(rows)})
  url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'+urllib.parse.urlencode({'db':'gene','id':','.join(r['GeneID']for r in rows),'retmode':'json','retmax':100})
  raw=urllib.request.urlopen(url,timeout=120).read();dump(META/paths[1],json.loads(raw));dump(META/'ncbi_gene_summaries_source.json',{'url':url,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)})
  query='(organism_id:10090) AND ('+' OR '.join('gene_exact:'+r['Symbol']for r in rows if r['type_of_gene']=='protein-coding')+')'
  url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode({'query':query,'format':'json','size':500,'fields':'accession,id,gene_names,protein_name,cc_function,organism_name'})
  raw=urllib.request.urlopen(url,timeout=120).read();dump(META/paths[2],json.loads(raw));dump(META/'uniprot_mouse_source.json',{'url':url,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)})
 for name,url in [('ak4_human_reviewed','https://rest.uniprot.org/uniprotkb/P27144.json'),('ak4_primary_papers','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=19130895,19766732&retmode=json')]:
  record=META/(name+('_uniprot'if name=='ak4_human_reviewed'else'')+'.json')
  if refresh or not record.exists():
   raw=urllib.request.urlopen(url,timeout=120).read();dump(record,json.loads(raw));dump(META/(name+'_source.json'),{'url':url,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'retrieved_on':datetime.now(timezone.utc).date().isoformat()})
 return [json.loads((META/p).read_text()) for p in paths]

def selections():
 h=pd.read_csv(ROOT/'results/hypothalamus/differential_expression.tsv',sep='\t');c=pd.read_csv(ROOT/'results/cortex/DE_PM_vs_PBS.tsv',sep='\t')
 panel=pd.read_csv(ROOT/'results/spatial_viewer/gene_panel_coverage.tsv',sep='\t')
 h['atlas_measured']=h.gene.isin(panel.symbol);c['atlas_measured']=c.gene_id.isin(panel.id)
 h['gene_id_ensembl']='';c['gene_id_ensembl']=c.gene_id
 h['source_gene_biotype']='not_supplied_in_original_symbol_count_matrix';c['source_gene_biotype']=c.gene_biotype
 valid=lambda t:np.isfinite(t[['padj','log2FoldChange','lfcSE','log2FoldChange_shrunk']]).all(axis=1)
 rules={
 'hypothalamus_primary_FDR':('hypothalamus',h[h.padj.lt(.05)],'All primary hypothalamus BH q<0.05 genes, sorted by q then source symbol; no atlas restriction.'),
 'hypothalamus_measured_top10':('hypothalamus',h[h.atlas_measured&valid(h)].sort_values(['padj','gene']).head(10),'Ten lowest finite primary BH q among exact-symbol measured genes with finite MLE/SE/shrunken effects. Only primary q<0.05 is significant; remaining rows are exploratory.'),
 'cortex_primary_top20':('cortex',c[c.padj.lt(.05)&c.gene_biotype.eq('protein_coding')].sort_values(['padj','gene']).head(20),'Twenty lowest primary BH q among protein-coding cortex genes with q<0.05, tie source symbol; no atlas restriction.'),
 'cortex_measured_top10':('cortex',c[c.padj.lt(.05)&c.atlas_measured&valid(c)].sort_values(['padj','gene']).head(10),'Ten lowest primary BH q among q<0.05 cortex genes with exact Ensembl measured-panel matches and finite effects; includes noncoding genes.')}
 out=[]
 for key,(region,t,rule) in rules.items():
  x=t.sort_values(['padj','gene']).copy();x['region']=region;x['selection_group']=key;x['selection_rank']=range(1,len(x)+1);x['ranking_rule']=rule;out.append(x)
 return pd.concat(out,ignore_index=True),rules

def annotations(rows,summaries,uniprot):
 result=[]
 for record in rows:
  gene=record['source_symbol'];symbol=record['Symbol'];matches=[]
  for entry in uniprot['results']:
   if not entry['entryType'].startswith('UniProtKB reviewed'):continue
   if entry['organism']['taxonId']!=10090:raise ValueError('Unexpected non-mouse UniProt record')
   primary={x.get('geneName',{}).get('value')for x in entry.get('genes',[])}
   alias={x['value']for g in entry.get('genes',[])for x in g.get('synonyms',[])}
   if symbol in primary:matches.append((0,entry))
   elif symbol in alias:matches.append((1,entry))
  best=min((a for a,_ in matches),default=None);chosen=[e for a,e in matches if a==best]
  if len(chosen)>1:raise ValueError(f'Ambiguous reviewed UniProt mapping: {gene}')
  entry=chosen[0]if chosen else None
  protein_en=entry['proteinDescription']['recommendedName']['fullName']['value']if entry else record['description']
  bio={'protein-coding':'protein_coding','ncRNA':'ncRNA','pseudo':'pseudogene'}.get(record['type_of_gene'],record['type_of_gene'])
  if bio=='ncRNA':protein_en='Noncoding RNA; no annotated protein product'
  if bio=='pseudogene':protein_en='No confirmed functional protein: Npm1 pseudogene'
  es,funen,funes,group=CURATED[gene]
  result.append({'gene':gene,'official_symbol':symbol,'gene_id_ncbi':record['GeneID'],'gene_type_ncbi':record['type_of_gene'],'biotype':bio,'chromosome':record['chromosome'],
   'protein_name_en':protein_en,'protein_name_es':es,'function_en':funen,'function_es':funes,'functional_group_en':GROUPS[group][0],'functional_group_es':GROUPS[group][1],
   'annotation_evidence':('NCBI current gene type; no established protein product'if bio!='protein_coding'else 'NCBI mouse annotation; '+('reviewed mouse UniProt record'if entry else 'no unique reviewed mouse UniProt record')+'; prediction/similarity qualifiers retained in function text where applicable')+('; mouse function inferred by similarity to human P27144, whose reviewed caution records both inactive and active in-vitro reports; no mouse catalytic assay performed here'if gene=='Ak4'else''),
   'source_url':'https://www.ncbi.nlm.nih.gov/gene/'+record['GeneID'],'uniprot_url':'https://www.uniprot.org/uniprotkb/'+entry['primaryAccession']+'/entry'if entry else '',
   'symbol_match':record['matching'],'functional_group_key':group,'protein_coding':bio=='protein_coding',
   'uniprot_match':'exact_primary_gene_name'if entry and best==0 else('unique_synonym'if entry else'not_assigned'),
   'ncbi_summary_available':bool(summaries['result'][record['GeneID']].get('summary')),
   'additional_source_urls':('https://www.uniprot.org/uniprotkb/P27144/entry;https://pubmed.ncbi.nlm.nih.gov/19130895/;https://pubmed.ncbi.nlm.nih.gov/19766732/'if gene=='Ak4'else'')})
 return pd.DataFrame(result)

def figures(t):
 palette={'protein_homeostasis':'#0072B2','signaling':'#009E73','transcription':'#CC79A7','extracellular':'#D55E00','metabolism':'#8A6B00','cytoskeleton':'#6B5B95','unknown':'#666666','noncoding':'#333333'}
 captions={};manifest={}
 for es in [False,True]:
  suffix='_es'if es else''
  for spatial,stem in [(False,'gene_priority_evidence'),(True,'spatial_candidate_evidence')]:
   groups=['hypothalamus_measured_top10','cortex_measured_top10']if spatial else['hypothalamus_primary_FDR','cortex_primary_top20']
   fig=plt.figure(figsize=(15.5,11.5 if not spatial else 8.8),layout='constrained')
   grid=fig.add_gridspec(2,2,height_ratios=[1,.13],width_ratios=[1.1,1.25])
   axes=[fig.add_subplot(grid[0,i])for i in range(2)]
   footer=fig.add_subplot(grid[1,:]);footer.axis('off')
   for i,(ax,group)in enumerate(zip(axes,groups)):
    d=t[t.selection_group.eq(group)].sort_values('selection_rank');ys=np.arange(len(d));lo=d.log2FoldChange-1.959963984540054*d.lfcSE;hi=d.log2FoldChange+1.959963984540054*d.lfcSE
    for j,(_,r)in enumerate(d.iterrows()):
     color='#0B7475' if i==0 else '#B85121'
     ax.plot([lo.iloc[j],hi.iloc[j]],[j,j],color=color,lw=1.5)
     ax.scatter(r.log2FoldChange,j,s=48,facecolors=color if r.primary_significant else'white',edgecolors=color,linewidths=1.4,zorder=4)
     ax.scatter(r.log2FoldChange_shrunk,j,s=23,marker='D',color='#111111',zorder=5)
    labels=[f"{r.gene}{'†'if r.atlas_measured else''}"+(' [ncRNA]'if r.biotype=='ncRNA'else' [pseudo]'if r.biotype=='pseudogene'else'')for r in d.itertuples()]
    ax.set_yticks(ys,labels,fontsize=10);ax.invert_yaxis();ax.axvline(0,c='#8A999C',ls='--',lw=.8);ax.grid(axis='x',c='#E7EEF0',lw=.8);ax.set_axisbelow(True)
    ax.set_xlabel('log₂FC de PM frente a control; MLE e IC puntual del 95%'if es else'PM versus control log₂FC; MLE and pointwise 95% CI',fontsize=10)
    ax.set_title(('A  Hipotálamo · '+str(len(d))+' genes'if i==0 else'B  Corteza cingulada · '+str(len(d))+' genes')if es else('A  Hypothalamus · '+str(len(d))+' genes'if i==0 else'B  Cingulate cortex · '+str(len(d))+' genes'),loc='left',fontsize=13,weight='bold')
    ax.spines[['top','right']].set_visible(False)
   fig.suptitle(('Genes medidos candidatos: evidencia regional y límites de significación'if spatial else'Genes prioritarios regionales: todos los hallazgos de hipotálamo y los 20 primeros de corteza')if es else('Measured-panel candidates: regional evidence and significance limits'if spatial else'Regional priority genes: all hypothalamic hits and the 20 leading cortex genes'),fontsize=15,weight='bold')
   legend=[Line2D([0],[0],marker='o',mfc='#4B6870',mec='#4B6870',ls='none',label='BH q<0.05'),Line2D([0],[0],marker='o',mfc='white',mec='#4B6870',ls='none',label='Exploratorio: q≥0.05'if es else'Exploratory: q≥0.05'),Line2D([0],[0],marker='D',c='#111111',ls='none',label='Efecto regularizado (sin intervalo MLE)'if es else'Shrunken effect (not the MLE interval)')]
   footer.legend(handles=legend,loc='lower center',ncol=3,fontsize=9,frameon=False)
   # Categories describe annotations, never enrichment or pathway activity.
   text=('† Gen medido en MERFISH. Colores: región analizada; las funciones se anotan en las tablas fuente.\nHipotálamo:4 FA / 3 PM; corteza:4 PBS / 4 PM. Selección por q, sin elegir animales; los intervalos no ajustan por selección.'if es else'† Gene measured by MERFISH. Colors identify regional panels; functions are annotated in the source tables.\nHypothalamus:4 FA / 3 PM; cortex:4 PBS / 4 PM. q-ranked selection without animal deletion; intervals do not adjust for selection.')
   footer.text(.5,.7,text,ha='center',va='center',fontsize=8.5,transform=footer.transAxes)
   key=stem+suffix
   for ext in ['png','svg','pdf']:
    metadata={'CreationDate':None,'ModDate':None}if ext=='pdf'else({'Date':None}if ext=='svg'else{})
    fig.savefig(OUT/'figures'/f'{key}.{ext}',dpi=300,metadata=metadata)
   plt.close(fig)
   cap=('Se muestran los diez genes medidos con menor q por región; solo Tll1 supera q<0.05 en hipotálamo y los otros nueve son exploratorios. Corteza conserva sus diez hallazgos medidos, incluido un RNA no codificante.'if spatial else'Se muestran todos los 11 hallazgos hipotalámicos y los 20 genes codificantes corticales de menor q. Los dos aumentos hipotalámicos corresponden a un pseudogén y un RNA no codificante según anotación actual; no se les asignan proteínas funcionales.')if es else('The ten lowest-q measured genes per region are displayed; only Tll1 passes q<0.05 in hypothalamus and the other nine are exploratory. Cortex retains its ten measured hits, including one noncoding RNA.'if spatial else'All 11 hypothalamic hits and the 20 lowest-q cortex protein-coding genes are displayed. The two higher hypothalamic signals correspond to a pseudogene and a noncoding RNA in current annotation; no functional proteins are assigned to them.')
   cap+=(' Círculos: MLE de RNA regional e IC puntual±1.96EE; diamantes: efectos regularizados separados. Las escalas horizontales varían entre paneles; los intervalos no ajustan por selección o multiplicidad. †:medido en MERFISH. Las categorías funcionales son anotaciones descriptivas, no mecanismos demostrados. Las células del atlas no son réplicas experimentales.'if es else' Circles: regional RNA MLE and pointwise±1.96SE interval; diamonds: separate shrunken effects. Horizontal scales vary between panels; intervals do not adjust for selection or multiplicity. †:measured by MERFISH. Functional categories are descriptive annotations, not demonstrated mechanisms. Atlas cells are not experimental exposure replicates.')
   captions[key]=cap
   manifest[key]={'caption':cap,'language':'es'if es else'en','source_table':'results/gene_priorities/prioritized_genes.tsv','groups':groups,'size_inches':[15.5,11.5 if not spatial else 8.8],'dpi':300,'files':{ext:{'path':str((OUT/'figures'/f'{key}.{ext}').relative_to(ROOT)),'sha256':sha(OUT/'figures'/f'{key}.{ext}')}for ext in['png','svg','pdf']}}
 dump(OUT/'figure_captions.json',captions);dump(OUT/'figure_manifest.json',manifest)

def artifact_context(t,a):
 """Persist interpretation, estimands and schema alongside the source-backed tables."""
 import platform
 method={'purpose':'Transparent ranked candidates and gene-product annotation; descriptive context for regional RNA and spatial reference visualization.',
 'audience':'Bilingual scientific website and downloadable research figures; no journal submission specification asserted.',
 'sample_units':{'hypothalamus':{'control_FA':4,'exposed_PM':3,'replicate':'animal'},'cortex':{'control_PBS':4,'exposed_PM':4,'replicate':'animal'}},
 'source_DE':'Unchanged primary regional results; cortex author20 reanalysis is not the ranking input.',
 'ranking':'BH q then source gene symbol; finite MLE, SE and shrunken coefficient required for measured-panel lists.',
 'interval':'Pointwise normal approximation: MLE log2FC ± 1.959963984540054 * original lfcSE. Not simultaneous or selection-adjusted; diamonds show separate shrunken coefficients.',
 'atlas_mapping':{'hypothalamus':'Exact original source symbol against MERFISH panel; source matrix does not provide original Ensembl IDs. Official current aliases are displayed separately.','cortex':'Exact source Ensembl gene ID against MERFISH panel; source IDs retained.'},
 'annotation_identity':'Current exact NCBI mouse symbol preferred; unique official alias explicitly disclosed. UniProt exact primary mouse gene name preferred to historical synonyms; only unique reviewed records assigned.',
 'functional_group_semantics':'Descriptive manual summaries from linked official records; not enrichment, pathway activation, cell identity, protein quantification or exposure-specific mechanism.',
 'cautions':['All primary hypothalamic coding hits decrease; both increases are noncoding or pseudogene signals under current annotation.','Nine of ten measured hypothalamic candidates are exploratory (q≥0.05); six show strong coefficient shrinkage.','All11 hypothalamic primary hits retain direction in7of7 leave-one-animal-out runs, but none retains q<0.05 in all7.','Cortex sex/identity metadata require verification; this study-level flag does not mean that every selected gene is sex-linked.','Regional RNA does not resolve cell composition, protein abundance, enzyme activity or transcript isoforms.','Ak4 mouse functional annotation is by similarity to human AK4; reviewed human annotation records conflicting in-vitro catalytic findings.','Ksr2 is interpreted as a scaffold; very low in-vitro kinase activity in reviewed annotation does not establish active catalysis in these data.'],
 'software':{'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__,'matplotlib':matplotlib.__version__},
 'rebuild':'./.venv-analysis/bin/python scripts/analyze_gene_priorities.py; use --refresh-annotations only to change the frozen annotation snapshot.'}
 dump(META/'methods.json',method)
 dump(OUT/'table_schema.json',{'annotations.tsv':{'rows':len(a),'key':['gene'],'columns':{k:str(v)for k,v in a.dtypes.items()}},'prioritized_genes.tsv':{'rows':len(t),'key':['region','selection_group','selection_rank'],'columns':{k:str(v)for k,v in t.dtypes.items()}},'missing_values':'TSV empty cells and JSON null denote unavailable values, never zero stability or an absent effect.','gene':'Original assay symbol retained; official_symbol is a separate annotation.','protein_names':'Official reviewed English protein name where uniquely mapped; Spanish name is a curated translation. Noncoding and pseudogene rows explicitly carry no protein assignment.','additional_source_urls':'Semicolon-delimited official records/primary papers supporting a specific qualification; presently Ak4.'})
 en='''# Regional gene priorities and encoded products

## Introduction

Regional transcriptomics can prioritize reproducible molecular candidates while a spatial reference atlas shows where their baseline RNA is detected. These two measurements answer different questions: differential expression estimates a regional exposure association, and an unexposed atlas supplies spatial context. Functional annotation identifies an encoded product or noncoding locus; it does not measure that product's abundance or activity.

## Methods

The existing primary hypothalamus and cingulate-cortex differential-expression tables were ranked without refitting models or removing animals. Four selection groups retain all 11 hypothalamic BH q<0.05 hits, the 20 lowest-q cortical protein-coding hits, the ten measured hypothalamic genes with lowest finite q, and the ten lowest-q measured cortical hits. Finite effect estimates and standard errors are required for measured-panel lists. Cortex panel availability uses exact original Ensembl IDs; the symbol-only hypothalamic source uses exact original symbols. The original Pcnx symbol is preserved and its unique current official alias Pcnx1 is disclosed.

Official NCBI mouse gene types and summaries and uniquely mapped reviewed mouse UniProt records were retrieved and frozen. Current primary names take precedence over historical aliases; this prevents the historical Fzd3 alias of Fzd9 from creating a false mapping. English protein names follow the reviewed record where available. Spanish names and both functional summaries are curated paraphrases with prediction/similarity qualifications preserved. Every gene links to its official record. Source URLs, hashes, selected records, software versions and ranking rules accompany the analysis. A default rerun uses this frozen annotation snapshot.

The forest figures show the unchanged MLE log2 fold change and pointwise normal 95% interval (MLE±1.96SE), with the shrunken estimate as a separate diamond. These are distinct estimates; the intervals are not confidence intervals for the shrunken coefficient and do not adjust for multiplicity or candidate selection. Each panel labels its own horizontal scale.

## Results and interpretation

There are 46 unique annotated loci across 51 selection rows: 43 protein-coding loci, two noncoding RNAs and one pseudogene. All nine coding hypothalamic primary hits decrease, including Hipk3, Ksr2, Tll1, Fzd3 and Fstl4. The two increases correspond to Gm13423 (an Npm1 pseudogene) and B230312C02Rik (noncoding RNA); functional NPM1 protein cannot be inferred from the former signal. All eleven directions persist in the seven leave-one-animal-out fits, but no gene passes BH q<0.05 in all seven fits.

The leading cortical genes include Hspa5/BiP, Pdia4, Xbp1, Pdia6, Creld2, Mlec and Manf, whose official annotations converge on protein folding and endoplasmic-reticulum homeostasis. This is descriptive convergence, not a new enrichment test or proof of an activated pathway. Fosb, Egr3 and Fosl2 encode transcription regulators; Homer1 encodes a postsynaptic scaffold; Nr1d1 and Ciart are circadian transcriptional repressors. Gene-level counts do not isolate ΔFosB, Homer1a or spliced XBP1s. Ak4 is described conservatively as a mitochondrial nucleotide/energy-homeostasis protein: its mouse annotation is by similarity to human AK4, whose reviewed record records both inactive and active catalytic studies. Ksr2 is treated as a MAPK scaffold, without assuming active kinase catalysis in these samples. Detailed molecular roles remain unresolved in the retrieved mouse records for Pcnx1 and Lyrm9.

Only Tll1 among the ten measured hypothalamic candidates passes primary BH q<0.05. The other nine are explicitly exploratory; six large MLE effects shrink close to zero. In cortex, all ten measured candidates pass primary BH q<0.05, including noncoding C730002L08Rik. Atlas measurement alone does not validate an exposure effect. Primary hits absent from the panel remain in the tables with atlas availability marked false.

## Limits

Small animal numbers, hypothalamic sensitivity to individual animals, and unresolved cortical sex/identity metadata limit interpretation. Missing cortical leave-one-out values are unavailable, not zero stability. The reference atlas contains unexposed cells; a regional RNA effect projected onto it remains a simulation and cannot identify exposed cell types. No protein abundance, enzyme activity, isoform-specific response or causal mechanism is established here. Statistical rank is an explicit prioritization rule rather than proof of biological importance.
'''
 es='''# Prioridades génicas regionales y productos codificados

## Introducción

La transcriptómica regional permite priorizar candidatos moleculares, mientras que un atlas espacial de referencia muestra dónde se detecta su RNA basal. La expresión diferencial estima una asociación regional con la exposición y el atlas sin exposición aporta contexto espacial. La anotación funcional identifica el producto codificado o el locus no codificante; no mide la abundancia ni la actividad de ese producto.

## Metodología

Se ordenaron los resultados principales existentes del hipotálamo y la corteza cingulada sin reajustar modelos ni eliminar animales. Cuatro grupos conservan los 11 hallazgos hipotalámicos con BH q<0.05, los 20 genes corticales codificantes con menor q, los diez genes medidos de hipotálamo con menor q finito y los diez hallazgos corticales medidos de menor q. Las listas medidas exigen efectos y errores estándar finitos. La disponibilidad en el panel se establece mediante los identificadores Ensembl originales exactos en corteza y mediante símbolos originales exactos en el hipotálamo, cuya matriz no aporta identificadores Ensembl. Se conserva Pcnx y se declara por separado su alias oficial actual único Pcnx1.

Se consultaron y congelaron los tipos génicos y resúmenes de ratón de NCBI y los registros revisados de UniProt asignables de forma única. El nombre principal actual tiene prioridad sobre los alias históricos; esto evita confundir Fzd3 con el alias histórico de Fzd9. Los nombres proteicos en inglés proceden del registro revisado cuando existe. Los nombres en español y los resúmenes funcionales bilingües son paráfrasis curadas que conservan las calificaciones de predicción o similitud. Cada gen enlaza a su registro oficial. Se entregan URLs, hashes, registros seleccionados, versiones y reglas de ordenación. La ejecución habitual reutiliza esta instantánea congelada.

Los gráficos muestran el log2 del cambio de expresión MLE original y un intervalo normal puntual del 95% (MLE±1.96EE), junto al efecto regularizado como diamante separado. Son estimaciones distintas: el intervalo no corresponde al coeficiente regularizado ni ajusta por multiplicidad o selección de candidatos. Cada panel presenta su propia escala horizontal.

## Resultados e interpretación

Se anotaron 46 loci únicos distribuidos en 51 filas de selección: 43 codificantes, dos RNA no codificantes y un pseudogén. Los nueve hallazgos hipotalámicos codificantes disminuyen, incluidos Hipk3, Ksr2, Tll1, Fzd3 y Fstl4. Los dos aumentos corresponden a Gm13423 (pseudogén de Npm1) y B230312C02Rik (RNA no codificante); la primera señal no demuestra proteína NPM1 funcional. Los once genes conservan su dirección en siete análisis que omiten un animal, pero ninguno mantiene BH q<0.05 en los siete.

Los genes corticales principales incluyen Hspa5/BiP, Pdia4, Xbp1, Pdia6, Creld2, Mlec y Manf, con anotaciones convergentes en plegamiento proteico y homeostasis del retículo endoplásmico. Esta convergencia es descriptiva: no constituye una nueva prueba de enriquecimiento ni demuestra activación de una vía. Fosb, Egr3 y Fosl2 codifican reguladores transcripcionales; Homer1, un andamio postsináptico; Nr1d1 y Ciart, represores transcripcionales circadianos. Los conteos génicos no aíslan ΔFosB, Homer1a ni XBP1s empalmado. Ak4 se describe con cautela como proteína mitocondrial asociada a nucleótidos y homeostasis energética: su anotación de ratón se apoya en similitud con AK4 humana, cuyo registro revisado recoge estudios con y sin actividad catalítica. Ksr2 se interpreta como andamio MAPK, sin presuponer actividad cinasa en estas muestras. Los registros de ratón consultados no establecen funciones moleculares detalladas para Pcnx1 y Lyrm9.

Solo Tll1 supera BH q<0.05 entre los diez candidatos hipotalámicos medidos. Los otros nueve son explícitamente exploratorios; seis efectos MLE grandes se regularizan casi hasta cero. Los diez candidatos corticales medidos superan BH q<0.05, incluido el RNA no codificante C730002L08Rik. Estar medido en el atlas no valida un efecto de exposición. Los hallazgos principales ausentes del panel se mantienen en las tablas con disponibilidad espacial falsa.

## Límites

El pequeño número de animales, la sensibilidad hipotalámica a animales individuales y las dudas pendientes sobre sexo e identidad de las muestras corticales limitan la interpretación. La estabilidad cortical faltante indica ausencia de evaluación, no estabilidad cero. El atlas contiene células sin exposición; proyectar sobre él un efecto regional sigue siendo una simulación y no identifica tipos celulares expuestos. Estos datos no establecen abundancia proteica, actividad enzimática, respuesta específica de isoformas ni mecanismos causales. El orden estadístico es una regla explícita de priorización, no una demostración de importancia biológica.
'''
 for lang,text in [('en',en),('es',es)]:
  sources='\n\n## '+('Official annotation sources'if lang=='en'else'Fuentes oficiales de anotación')+'\n\n'
  sources+='\n'.join(f"- {r.gene} ({r.official_symbol}): [NCBI Gene]({r.source_url})"+(f"; [UniProt]({r.uniprot_url})"if r.uniprot_url else'')for r in a.itertuples())
  sources+='\n\n'+('Ak4 catalytic qualification: 'if lang=='en'else'Calificación catalítica de Ak4: ')+ '[human reviewed AK4](https://www.uniprot.org/uniprotkb/P27144/entry), [Liu et al.](https://pubmed.ncbi.nlm.nih.gov/19130895/), [Panayiotou et al.](https://pubmed.ncbi.nlm.nih.gov/19766732/).\n'
  (OUT/f'narrative_{lang}.md').write_text(text+sources)

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--refresh-annotations',action='store_true');args=parser.parse_args()
 OUT.mkdir(parents=True,exist_ok=True);META.mkdir(parents=True,exist_ok=True);(OUT/'figures').mkdir(exist_ok=True)
 selected,rules=selections();rows,summaries,u=raw_sources(args.refresh_annotations);a=annotations(rows,summaries,u)
 if set(selected.gene)!=set(a.gene):raise ValueError('Frozen annotation candidates differ; refresh annotations explicitly')
 a.to_csv(OUT/'annotations.tsv',sep='\t',index=False)
 t=selected.merge(a,on='gene',how='left',validate='many_to_one')
 loo=pd.read_csv(ROOT/'results/hypothalamus/leave_one_out_primary_hits.tsv',sep='\t').set_index('gene')
 t['runs_same_lfc_direction']=[loo.at[r.gene,'runs_same_lfc_direction']if r.region=='hypothalamus'and r.gene in loo.index else np.nan for r in t.itertuples()]
 t['runs_fdr_below_0_05']=[loo.at[r.gene,'runs_fdr_below_0_05']if r.region=='hypothalamus'and r.gene in loo.index else np.nan for r in t.itertuples()]
 t['primary_significant']=t.padj.lt(.05);t['CI95_lower']=t.log2FoldChange-1.959963984540054*t.lfcSE;t['CI95_upper']=t.log2FoldChange+1.959963984540054*t.lfcSE
 t['evidence_label_en']=np.where(t.primary_significant,'Primary regional BH q<0.05','Exploratory regional RNA estimate; BH q≥0.05')
 t['evidence_label_es']=np.where(t.primary_significant,'BH regional principal q<0.05','Estimación regional exploratoria de RNA; BH q≥0.05')
 def flag(r):
  out=[]
  if r.region=='hypothalamus':out.append('individual_animal_sensitivity')
  else:out.append('study_sex_identity_metadata_requires_verification')
  if r.biotype=='pseudogene':out.append('pseudogene_read_mapping_specificity_not_evaluated')
  if r.chromosome in ['X','Y']:out.append('sex_chromosome_gene')
  if abs(r.log2FoldChange)>=1 and abs(r.log2FoldChange_shrunk)<.1*abs(r.log2FoldChange):out.append('strong_shrinkage_keep_estimands_separate')
  return';'.join(out)
 t['qc_flag']=[flag(r)for r in t.itertuples()]
 core=['region','selection_group','selection_rank']+list(a.columns)+['gene_id_ensembl','source_gene_biotype','atlas_measured','baseMean','log2FoldChange','lfcSE','CI95_lower','CI95_upper','log2FoldChange_shrunk','padj','primary_significant','evidence_label_en','evidence_label_es','runs_same_lfc_direction','runs_fdr_below_0_05','qc_flag','ranking_rule']
 t=t[core];t.to_csv(OUT/'prioritized_genes.tsv',sep='\t',index=False)
 artifact_context(t,a)
 # JSON is safe for web clients; missing cross-study stability is null, never zero.
 dump(OUT/'annotations.json',json.loads(a.to_json(orient='records',double_precision=15)))
 dump(OUT/'prioritized_genes.json',json.loads(t.to_json(orient='records',double_precision=15)))
 dump(OUT/'selection_rules.json',{k:{'region':r,'n':int(t.selection_group.eq(k).sum()),'genes':t[t.selection_group.eq(k)].sort_values('selection_rank').gene.tolist(),'rule':rule}for k,(r,_,rule)in rules.items()})
 figures(t)
 source_values=[]
 for region,path in [('hypothalamus','results/hypothalamus/differential_expression.tsv'),('cortex','results/cortex/DE_PM_vs_PBS.tsv')]:
  original=pd.read_csv(ROOT/path,sep='\t').set_index('gene')
  observed=t[t.region.eq(region)].set_index('gene')
  columns=['padj','log2FoldChange','lfcSE','log2FoldChange_shrunk','baseMean']
  source_values.append(np.array_equal(observed[columns].to_numpy(),original.loc[observed.index,columns].to_numpy(),equal_nan=True))
 checks={'all11_hypothalamus_hits_retained':len(t[t.selection_group.eq('hypothalamus_primary_FDR')])==11,'cortex_top20_rank_recomputed':len(t[t.selection_group.eq('cortex_primary_top20')])==20,'both_measured_groups_have10':all(len(t[t.selection_group.eq(k)])==10 for k in['hypothalamus_measured_top10','cortex_measured_top10']),'hypothalamus_measured_only_Tll1_primary':t[t.selection_group.eq('hypothalamus_measured_top10')&t.primary_significant].gene.tolist()==['Tll1'],'no_protein_assigned_to_noncoding_or_pseudogene':all('No 'in r.protein_name_en or'Noncoding'in r.protein_name_en for r in a[~a.protein_coding].itertuples()),'all_gene_identities_unique':not a.gene.duplicated().any(),'all_NCBI_source_links':a.source_url.str.startswith('https://www.ncbi.nlm.nih.gov/gene/').all(),'Pcnx_alias_disclosed':a.set_index('gene').at['Pcnx','official_symbol']=='Pcnx1','Fzd3_not_merged_with_Fzd9':a.set_index('gene').at['Fzd3','uniprot_url'].endswith('/Q61086/entry'),'source_statistics_unchanged':all(source_values),'selected_effects_are_finite':bool(np.isfinite(t[['padj','log2FoldChange','lfcSE','log2FoldChange_shrunk']]).all().all())}
 checks={k:bool(v)for k,v in checks.items()};dump(OUT/'validation.json',{'passed':all(checks.values()),'checks':checks,'scope':'Selection, identity, annotation-type and numeric invariants; manual source/figure review separately recorded.'})
 summary={'complete':True,'annotation_date_UTC':datetime.now(timezone.utc).isoformat(),'unique_genes':len(a),'selected_rows':len(t),'protein_coding_genes':int(a.protein_coding.sum()),'groups':{k:{'rows':len(d),'primary_FDR':int(d.primary_significant.sum()),'atlas_measured':int(d.atlas_measured.sum()),'genes':d.gene.tolist()}for k,d in t.groupby('selection_group',sort=False)},'not_protein_coding':a[~a.protein_coding][['gene','biotype','protein_name_en']].to_dict('records'),'ranking_is_not_causal_relevance':True,'no_primary_analysis_refitted':True,'no_exposed_single_cells_inferred':True}
 dump(OUT/'summary.json',summary)
 sources=[]
 for p in ['ncbi_gene_info_source.json','ncbi_gene_summaries_source.json','uniprot_mouse_source.json','ak4_human_reviewed_source.json','ak4_primary_papers_source.json']:
  sources.append(json.loads((META/p).read_text()))
 dump(META/'source_manifest.json',{'retrieved_on':'2026-10-19','sources':sources,'frozen_local_record_files':{p.name:sha(p)for p in META.glob('*.json')if p.name!='source_manifest.json'},'primary_input_hashes':{p:sha(ROOT/p)for p in ['results/hypothalamus/differential_expression.tsv','results/cortex/DE_PM_vs_PBS.tsv','results/hypothalamus/leave_one_out_primary_hits.tsv','results/spatial_viewer/gene_panel_coverage.tsv']},'annotation_identity_policy':'NCBI current exact symbol wins; a unique official alias is disclosed. UniProt exact primary gene name wins over historical synonyms; only unique reviewed mouse records are assigned.'})
 print(json.dumps(summary,indent=2))
 if not all(checks.values()):raise SystemExit('Validation check failed')
if __name__=='__main__':main()
