/* Offline scientific viewer. Reference counts are loaded one gene at a time. */
"use strict";
const $ = id => document.getElementById(id);
const D = window.ATLAS_DATA;
const geneBySymbol = new Map(D.genes.map(g => [g.symbol.toLowerCase(), g]));
const geneById = new Map(D.genes.map(g => [g.id, g]));
const pending = new Map(), cache = new Map(), inFlight = new Map();
let currentGene = null, currentCounts = null, renderVersion = 0, sceneVersion = 0, cameraLock = false;
let lastCamera = {eye:{x:1.2,y:1.2,z:.9},up:{x:0,y:0,z:1}};
const config = {responsive:true,displaylogo:false,toImageButtonOptions:{format:"png",filename:"hypothalamus_reference_scenario",scale:2},modeBarButtonsToRemove:["lasso2d","select2d"]};
const colors = {FA:"#0d7475",PM:"#c05b23"};
const fnum = (v,d=3) => Number.isFinite(v) ? Number(v).toLocaleString("en",{maximumFractionDigits:d}) : "not estimable";
const fp = v => !Number.isFinite(v) ? "unavailable" : v < .01 ? v.toExponential(2) : v.toFixed(3);

window.receiveAtlasGene = payload => {
  const receiver = pending.get(payload.id);
  if(receiver){pending.delete(payload.id);receiver.resolve(payload);}
};

async function fetchCounts(gene){
  if(cache.has(gene.id)) return cache.get(gene.id);
  if(!("DecompressionStream" in window)) throw new Error("This offline viewer requires a recent browser with gzip decompression support. The static figure is available in the report.");
  const payload = await new Promise((resolve,reject)=>{
    pending.set(gene.id,{resolve,reject});
    const script=document.createElement("script");
    script.src="genes/"+encodeURIComponent(gene.id)+".js";
    script.onerror=()=>{pending.delete(gene.id);reject(new Error("Gene data could not be loaded: "+gene.symbol));};
    script.onload=()=>script.remove();
    document.head.appendChild(script);
  });
  const packed=Uint8Array.from(atob(payload.gzip_f32le),c=>c.charCodeAt(0));
  const stream=new Blob([packed]).stream().pipeThrough(new DecompressionStream("gzip"));
  const buffer=await new Response(stream).arrayBuffer();
  const counts=new Float32Array(buffer);
  if(counts.length!==D.cell_count||payload.n!==D.cell_count) throw new Error("Gene data and coordinate lengths differ");
  cache.set(gene.id,counts);
  if(cache.size>10) cache.delete(cache.keys().next().value);
  return counts;
}

function loadCounts(gene){
  if(cache.has(gene.id))return Promise.resolve(cache.get(gene.id));
  if(inFlight.has(gene.id))return inFlight.get(gene.id);
  const operation=fetchCounts(gene).finally(()=>inFlight.delete(gene.id));
  inFlight.set(gene.id,operation);
  return operation;
}

function meshes(){
  if(!$('anatomy').checked)return [];
  return [
    {...D.meshes.brain,type:"mesh3d",name:"Whole-brain anatomy",color:"#b6c5ca",opacity:.075,flatshading:false,hoverinfo:"skip",showscale:false,lighting:{ambient:.8,diffuse:.6}},
    {...D.meshes.hypothalamus,type:"mesh3d",name:"Hypothalamus anatomy",color:"#6a969a",opacity:.095,flatshading:false,hoverinfo:"skip",showscale:false,lighting:{ambient:.8,diffuse:.6}}
  ];
}

function ranges(){
  if($('focus').value==="brain") return {x:[0,11.4],y:[0,13.2],z:[-8.2,.1]};
  const bounds=a=>[Math.min(...a)-.3,Math.max(...a)+.3];
  return {x:bounds(D.x),y:bounds(D.y),z:bounds(D.z)};
}

function layout(){
  const r=ranges();
  const axis=(title,range)=>({title:{text:title,font:{size:10}},range,backgroundcolor:"#fff",gridcolor:"#edf2f3",zeroline:false,tickfont:{size:9},showspikes:false});
  const sizes=[r.x[1]-r.x[0],r.y[1]-r.y[0],r.z[1]-r.z[0]],largest=Math.max(...sizes);
  return {margin:{l:0,r:0,t:0,b:0},paper_bgcolor:"#fff",font:{family:"Arial,sans-serif",color:"#173239"},showlegend:false,uirevision:$('focus').value,
    scene:{camera:lastCamera,aspectmode:"manual",aspectratio:{x:sizes[0]/largest,y:sizes[1]/largest,z:sizes[2]/largest},xaxis:axis("ML (mm)",r.x),yaxis:axis("AP (mm)",r.y),zaxis:axis("−DV (mm)",r.z)}};
}

function selectedIndices(){
  const ct=$('cell-type').value,region=$('region').value,hide=$('hide-zero').checked;
  const selected=[];
  for(let i=0;i<D.cell_count;i++) if((ct==="all"||D.type_codes[i]===Number(ct))&&(region==="all"||D.region_codes[i]===Number(region))&&(!hide||currentCounts[i]>0)) selected.push(i);
  return selected;
}

function exponent(gene){
  if(!gene.estimable) return null;
  if($('scenario').value==="lower") return gene.ci_low;
  if($('scenario').value==="upper") return gene.ci_high;
  return gene.effect.log2FoldChange_shrunk;
}

function pointTrace(indices,values,scaleMax,simulated,available){
  const choose=a=>indices.map(i=>a[i]);
  const custom=indices.map(i=>[D.cell_ids[i],D.type_labels[D.type_codes[i]],D.region_labels[D.region_codes[i]],currentCounts[i],values[i]]);
  return {type:"scatter3d",mode:"markers",name:simulated?"Simulated HY cells":"Measured HY cells",x:choose(D.x),y:choose(D.y),z:choose(D.z),customdata:custom,
    marker:{size:2.3,opacity:.86,color:simulated&&!available?"#b8c3c8":choose(values).map(v=>Math.log2(1+v)),colorscale:"Viridis",cmin:0,cmax:scaleMax,showscale:!simulated||available,colorbar:{title:{text:"log₂(1 + count)",font:{size:10}},thickness:10,len:.6,tickfont:{size:9},x:.96}},
    hovertemplate:"%{customdata[1]}<br>Region: %{customdata[2]}<br>"+(simulated?(available?"Simulated value: %{customdata[4]:.3f}":"Exposure fit unavailable"):"Measured molecules: %{customdata[3]:.0f}")+"<br>Cell %{customdata[0]}<extra></extra>"};
}

function updateMetrics(gene,indices,ratio){
  $('gene-name').textContent=gene.symbol;
  $('evidence').textContent=gene.significant?"Primary bulk FDR < 0.05 · small-n sensitive":gene.estimable?"Exploratory · not FDR-significant":"Reference only · no eligible matched effect";
  $('effect').textContent=gene.estimable?fnum(gene.effect.log2FoldChange_shrunk)+" log₂FC":"Not estimable";
  $('fdr').textContent="Bulk Wald FDR: "+fp(gene.effect.padj);
  $('scenario-label').textContent='SCENARIO MULTIPLIER · '+(gene.estimable?{central:'CENTRAL',lower:'LOWER BOUND',upper:'UPPER BOUND'}[$('scenario').value]:'UNAVAILABLE');
  $('multiplier').textContent=ratio===null?"Unavailable":fnum(ratio)+"×";
  $('interval').textContent=gene.estimable?"Group MLE 95% CI: ["+fnum(gene.ci_low)+", "+fnum(gene.ci_high)+"] log₂FC":"No exposure effect is assigned";
  $('cells').textContent=indices.length.toLocaleString("en");
  $('selection-note').hidden=indices.length>0;
  const detected=indices.reduce((n,i)=>n+(currentCounts[i]>0?1:0),0);
  $('coverage').textContent="of "+D.cell_count.toLocaleString("en")+" HY cells · "+detected.toLocaleString("en")+" shown with measured count > 0";
  $('scenario').disabled=!gene.estimable;
  $('simulation-title').innerHTML=gene.estimable?'PM2.5 scenario <span class="sim-label">SIMULATED</span>':'Exposure fit unavailable <span>Reference cells in gray</span>';
  const s=gene.leave_one_out;
  $('stability').textContent=s?"Sensitivity analysis: "+gene.symbol+" retains its effect direction in "+s.runs_same_lfc_direction+"/7 leave-one-out fits and FDR < 0.05 in "+s.runs_fdr_below_0_05+"/7. No animal was selectively removed from the primary analysis.":gene.estimable?"Exploratory genes are shown without evidence of a significant bulk exposure effect. The primary hypothesis test uses seven animals, not the atlas cells.":"No reliable bulk exposure test is available for this gene. The measured reference is retained; no exposure effect is assigned.";
}

function bulkPanel(gene){
  const traces=["FA","PM"].map(condition=>{
    const rows=gene.bulk_values.map((v,i)=>({...v,index:i})).filter(v=>v.condition===condition);
    return {type:"scatter",mode:"markers",name:condition==="FA"?"Filtered air (n=4)":"PM2.5 (n=3)",x:rows.map(v=>v.label),y:rows.map(v=>v.count),marker:{color:colors[condition],size:10,symbol:condition==="FA"?"circle":"diamond"},text:rows.map(v=>v.sample),hovertemplate:"%{x}<br>%{text}<br>Normalized count: %{y:.2f}<extra></extra>"};
  });
  Plotly.react('bulk',traces,{margin:{t:15,l:65,b:45,r:15},height:280,font:{family:"Arial,sans-serif",size:11,color:"#173239"},paper_bgcolor:"white",plot_bgcolor:"white",xaxis:{title:"Independent bulk animal",categoryorder:"array",categoryarray:["FA1","FA2","FA3","FA4","PM1","PM2","PM3"]},yaxis:{title:"DESeq2-normalized counts",rangemode:"tozero",gridcolor:"#e8eff0"},legend:{orientation:"h",y:1.18}},config);
  const table=document.createElement('table');
  table.innerHTML='<thead><tr><th>Animal</th><th>Group</th><th>Normalized count</th></tr></thead>';
  const body=document.createElement('tbody');
  for(const v of gene.bulk_values){const tr=document.createElement('tr');for(const text of [v.label,v.condition,fnum(v.count,2)]){const td=document.createElement('td');td.textContent=text;tr.appendChild(td);}body.appendChild(tr);}
  table.appendChild(body);$('bulk-table').replaceChildren(table);
}

async function render(){
  if(!currentGene||!currentCounts)return;
  const revision=++sceneVersion;
  const gene=currentGene,indices=selectedIndices();
  const beta=exponent(gene),ratio=beta===null?null:2**beta;
  if(ratio!==null&&!Number.isFinite(ratio))throw new Error("The selected coefficient interval is numerically unstable. Use the central scenario.");
  const predicted=Array.from(currentCounts,v=>v*(ratio===null?1:ratio));
  const all=[];for(let i=0;i<D.cell_count;i++){all.push(Math.log2(1+currentCounts[i]));all.push(Math.log2(1+predicted[i]));}
  all.sort((a,b)=>a-b);
  const positive=all.filter(v=>v>0);
  const scaleMax=Math.max(.5,$('clip').checked?(positive.length?positive[Math.floor((positive.length-1)*.99)]:0):all[all.length-1]);
  const clipped=all.filter(v=>v>scaleMax).length;
  updateMetrics(gene,indices,ratio);
  const common=meshes();
  await Promise.all([Plotly.react('reference',[...common,pointTrace(indices,currentCounts,scaleMax,false,true)],layout(),config),Plotly.react('simulation',[...common,pointTrace(indices,predicted,scaleMax,true,gene.estimable)],layout(),config)]);
  if(gene!==currentGene||revision!==sceneVersion)return;
  $('scale-note').textContent="Common color scale: log₂(1 + imaging count), 0–"+fnum(scaleMax)+". "+($('clip').checked?clipped.toLocaleString('en')+" values across the full HY reference/scenario exceed the common positive-value 99th-percentile cap; values are retained in hover and exports.":"Full reference/scenario range; no color capping.")+" Geometry is unchanged.";
  if(!gene.estimable)$('scale-note').textContent="Reference color scale: log₂(1 + imaging count), 0–"+fnum(scaleMax)+". The right panel shows reference positions in gray; no exposure simulation is available.";
  if(window.parent!==window)window.parent.postMessage({type:"omics-atlas-height",height:document.documentElement.scrollHeight},"*");
}

async function applyGene(){
  const gene=geneBySymbol.get($('gene').value.trim().toLowerCase())||geneById.get($('gene').value.trim());
  if(!gene){$('loading').textContent="Choose a symbol or identifier from the measured panel.";$('loading').classList.add('error');return;}
  const version=++renderVersion;
  $('loading').textContent="Loading "+gene.symbol+"…";$('loading').classList.remove('error');$('apply').disabled=true;
  try{
    const counts=await loadCounts(gene);
    if(version!==renderVersion)return;
    currentGene=gene;currentCounts=counts;$('gene').value=gene.symbol;
    await render();
    if(version!==renderVersion)return;
    bulkPanel(gene);$('loading').textContent=gene.estimable?"Ready · measured panel / explicit scenario":"Ready · measured reference only";
  }catch(error){$('loading').textContent=error.message;$('loading').classList.add('error');console.error(error);}
  finally{if(version===renderVersion)$('apply').disabled=false;}
}

function linkCameras(source,target){
  $(source).on('plotly_relayout',event=>{
    if(cameraLock||!event['scene.camera'])return;
    lastCamera=event['scene.camera'];cameraLock=true;
    Plotly.relayout(target,{'scene.camera':lastCamera}).finally(()=>{cameraLock=false;});
  });
}

async function init(){
  for(const g of [...D.genes].sort((a,b)=>a.symbol.localeCompare(b.symbol))){const option=document.createElement('option');option.value=g.symbol;option.textContent=g.id+(g.significant?" · primary FDR hit":"");$('gene-list').appendChild(option);}
  for(const [id,labels] of [['cell-type',D.type_labels],['region',D.region_labels]])labels.forEach((label,i)=>{const o=document.createElement('option');o.value=i;o.textContent=label;$(id).appendChild(o);});
  $('gene').value=D.default_gene;
  $('apply').addEventListener('click',applyGene);
  $('gene').addEventListener('keydown',event=>{if(event.key==='Enter')applyGene();});
  for(const id of ['scenario','cell-type','region','clip','hide-zero','anatomy'])$(id).addEventListener('change',()=>render().catch(error=>{$('loading').textContent=error.message;$('loading').classList.add('error');}));
  $('focus').addEventListener('change',()=>{lastCamera={eye:{x:1.2,y:1.2,z:.9},up:{x:0,y:0,z:1}};render();});
  $('reset').addEventListener('click',()=>{lastCamera={eye:{x:1.2,y:1.2,z:.9},up:{x:0,y:0,z:1}};render();});
  await applyGene();
  linkCameras('reference','simulation');linkCameras('simulation','reference');
}
init().catch(error=>{$('loading').textContent=error.message;$('loading').classList.add('error');console.error(error);});
