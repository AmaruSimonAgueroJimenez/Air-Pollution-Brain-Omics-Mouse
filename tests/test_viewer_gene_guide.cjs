// Exercise production UI functions using real field shapes and a minimal DOM.
const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');

class Element{
  constructor(tag='div'){this.tagName=tag.toUpperCase();this.children=[];this.dataset={};this.attributes={};this.listeners={};this.hidden=false;this.disabled=false;this.checked=false;this.value='all';this.ownText='';this.classList={add(){},remove(){}};}
  set textContent(value){this.ownText=String(value);this.children=[];}
  get textContent(){return this.ownText+this.children.map(c=>c.textContent).join('');}
  set innerHTML(value){this.html=value;this.children=[];this.ownText='';}
  appendChild(child){this.children.push(child);return child;}
  append(...children){this.children.push(...children);}
  replaceChildren(...children){this.children=[...children];this.ownText='';}
  setAttribute(name,value){this.attributes[name]=String(value);}
  addEventListener(name,handler){this.listeners[name]=handler;}
  querySelectorAll(selector){return this.children.flatMap(c=>[...(c.tagName===selector.toUpperCase()?[c]:[]),...c.querySelectorAll(selector)]);}
}
const clone=x=>JSON.parse(JSON.stringify(x));
function fixture({region='ACA',lang='en'}={}){
  const n=region==='ACA'?8:7,control=region==='ACA'?'PBS':'FA';
  const bulkValues=Array.from({length:n},(_,i)=>({sample:i<4?'Control'+(i+5):'PMA'+(i-3),label:i<4?'C'+(i+1):'PM'+(i-3),condition:i<4?control:'PM',count:10+i}));
  const gene=(symbol,id,significant,padj,biotype='protein_coding')=>({symbol,id,name:'Original description '+symbol,significant,estimable:true,effect:{padj,log2FoldChange_shrunk:.75,log2FoldChange:.8,lfcSE:.2},ci_low:.4,ci_high:1.2,bulk_values:clone(bulkValues),leave_one_out:null,gene_biotype:biotype});
  const genes=[gene('Fosb','ENSMUSG_Fosb',true,2.2e-19),gene('Tll1','ENSMUSG_Tll1',false,.2),gene('C730002L08Rik','ENSMUSG_Rik',true,.003,'lncRNA')];
  const nodes=new Map();
  const get=id=>{if(!nodes.has(id))nodes.set(id,new Element());return nodes.get(id);};
  get('scenario').value='central';get('focus').value='hypothalamus';
  const window={ATLAS_DATA:{region_key:region,control_label:control,n_control:4,n_PM:n-4,bulk_sample_n:n,cell_count:3,genes,region_label_en:region==='ACA'?'Anterior cingulate cortex':'Hypothalamus',region_label_es:region==='ACA'?'Corteza cingulada anterior':'Hipotálamo'},
    GENE_ANNOTATIONS:{Fosb:{biotype:'protein_coding',protein_name_en:'Protein FOSB',protein_name_es:'Proteína FOSB',function_en:'Curated function.',function_es:'Función curada.',source_url:'https://www.ncbi.nlm.nih.gov/gene/14282'},C730002L08Rik:{biotype:'lncRNA',protein_name_en:'A misleading product must not override the noncoding annotation',function_en:'A noncoding transcript.',source_url:'https://www.ncbi.nlm.nih.gov/gene/000'}},
    GENE_PRIORITIES:{[region]:[{symbol:'Fosb',padj:2.2e-19,measured:true},{symbol:'Hspa5',padj:3.1e-33,measured:false},{symbol:'C730002L08Rik',padj:.003,measured:true}]}};
  window.parent=window;
  const calls=[];
  const context=vm.createContext({window,document:{documentElement:{lang},getElementById:get,createElement:tag=>new Element(tag),querySelectorAll:()=>[],querySelector:()=>null},console,
    Plotly:{register(){},react:(id,traces,layout)=>{calls.push({id,traces,layout});return Promise.resolve();}},Map,Set,Number,Promise,JSON});
  const source=fs.readFileSync(path.join(__dirname,'../docs/assets/viewer/viewer.js'),'utf8');
  const init=source.lastIndexOf('\ninit().catch(');assert.ok(init>0);
  vm.runInContext(source.slice(0,init)+`\nglobalThis.api={bulkPanel,initializeStudyUI,proteinDescription,buildGeneGuide,updateMetrics,summarizeReference,renderComparison,measuredPriorities,exponent,notifyParentHeight,
    setCurrent:(gene,counts)=>{currentGene=gene;currentCounts=counts;},seed:(id,counts)=>cache.set(id,counts),compare:symbols=>{comparisonGenes=symbols;}};`,context);
  return {api:context.api,nodes,get,calls,genes,window};
}

test('cortex bulk plot retains all four PBS and four PM animals with source labels',()=>{
  const f=fixture();f.api.bulkPanel(f.genes[0]);
  const plot=f.calls.find(c=>c.id==='bulk');
  assert.equal(plot.traces[0].name,'PBS (n=4)');assert.equal(plot.traces[1].name,'PM2.5 (n=4)');
  assert.deepEqual(Array.from(plot.traces[0].y),[10,11,12,13]);assert.deepEqual(Array.from(plot.traces[1].y),[14,15,16,17]);
  assert.deepEqual(Array.from(plot.layout.xaxis.categoryarray),['C1','C2','C3','C4','PM1','PM2','PM3','PM4']);
});

test('hypothalamus bulk plot preserves its distinct 4FA/3PM design',()=>{
  const f=fixture({region:'HY'});f.api.bulkPanel(f.genes[0]);
  const plot=f.calls.find(c=>c.id==='bulk');
  assert.equal(plot.traces[0].name,'Filtered air (n=4)');assert.equal(plot.traces[1].name,'PM2.5 (n=3)');
  assert.equal(plot.layout.xaxis.categoryarray.length,7);
});

test('significant cortical gene remains significant and uses eight-animal context',()=>{
  const f=fixture();f.api.setCurrent(f.genes[0],[0,3,6]);f.api.updateMetrics(f.genes[0],[0,1,2],2**.75);
  assert.match(f.get('evidence').textContent,/Primary cortical bulk FDR < 0.05/);
  assert.match(f.get('stability').textContent,/8 independent animals/);
  assert.doesNotMatch(f.get('stability').textContent,/seven|Exploratory gene/);
  assert.match(f.get('coverage').textContent,/ACA cells/);
});

test('noncoding annotation does not acquire an invented protein product',()=>{
  const f=fixture();
  assert.match(f.api.proteinDescription(f.genes[2]),/Non-coding annotation; no protein product assigned/);
  assert.doesNotMatch(f.api.proteinDescription(f.genes[2]),/misleading/);
});

test('unmeasured priority stays visible with a gene-profile link and no 3D action',()=>{
  const f=fixture();f.api.buildGeneGuide();
  const table=f.get('priority-table').children[0],body=table.children[0];
  const absent=body.children.find(r=>r.children[0].textContent==='Hspa5');
  assert.ok(absent,'Retain unmeasured gene in the priorities table');
  assert.equal(absent.children[3].textContent,'Not measured');
  assert.equal(absent.children[4].querySelectorAll('button').length,0);
  assert.equal(absent.children[4].children[0].href,'../../gene_profiles.html');
  assert.equal(f.api.measuredPriorities().some(g=>g.symbol==='Hspa5'),false);
});

test('reference summaries include measured zeros in detection and mean',()=>{
  const f=fixture();const summary=f.api.summarizeReference([0,3,6]);
  assert.equal(summary.cells,3);assert.equal(summary.positive,2);assert.equal(summary.mean,3);
  assert.ok(Math.abs(summary.detection-200/3)<1e-10);
});

test('three-gene comparison retains gene-specific evidence and whole-reference summaries',async()=>{
  const f=fixture();f.api.compare(f.genes.map(g=>g.symbol));
  for(const g of f.genes)f.api.seed(g.id,[0,3,6]);
  await f.api.renderComparison();
  const body=f.get('comparison-table').children[0].children[0];
  assert.equal(body.children.length,3);
  assert.equal(body.children[0].children[1].textContent,'Primary FDR hit');
  assert.equal(body.children[1].children[1].textContent,'Exploratory · FDR ≥ 0.05');
  for(const row of body.children){assert.equal(row.children[4].textContent,'66.67%');assert.equal(row.children[5].textContent,'3');}
  assert.match(f.get('comparison-status').textContent,/every registered regional reference cell, including zeros/);
});

test('Spanish cortex labels keep PBS, regional function and sample population accurate',()=>{
  const f=fixture({lang:'es'});f.api.initializeStudyUI();
  assert.match(f.get('bulk-effect-label').textContent,/PM \/ PBS/);
  assert.match(f.get('bulk-description').textContent,/8 bibliotecas/);
  assert.match(f.get('bulk-description').textContent,/4 controles PBS y 4 animales PM2.5/);
  assert.match(f.get('regional-model').textContent,/Solo dentro de ACA/);
  assert.match(f.get('study-limit').textContent,/ocho animales corticales/);
});

test('iframe size uses intrinsic main content and does not accumulate a viewport offset',()=>{
  const f=fixture(),messages=[];
  f.window.parent={postMessage:payload=>messages.push(payload)};
  // The real document lookup is isolated, so create a dedicated minimal context.
  const source=fs.readFileSync(path.join(__dirname,'../docs/assets/viewer/viewer.js'),'utf8');
  const start=source.indexOf('function notifyParentHeight(){'),end=source.indexOf('\nfunction annotationFor',start);
  const context=vm.createContext({window:{parent:{postMessage:p=>messages.push(p)}},document:{querySelector:()=>({getBoundingClientRect:()=>({height:2375.2})})},lastReportedHeight:null,Math});
  vm.runInContext(source.slice(start,end)+'\nnotifyParentHeight();notifyParentHeight();',context);
  assert.equal(messages.length,1);assert.equal(messages[0].height,2376);
});
