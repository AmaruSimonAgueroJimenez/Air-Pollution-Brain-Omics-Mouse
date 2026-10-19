/* Offline scientific viewer. Reference counts are loaded one gene at a time. */
"use strict";
const $ = id => document.getElementById(id);
const spanish = /^es\b/i.test(document.documentElement.lang);
const tr = (english, spanishText) => spanish ? spanishText : english;
// Keep the source report's decimal/thousands convention in both languages.
const numberLocale = "en";
if(spanish) Plotly.register({moduleType:"locale",name:"es",dictionary:{
  "Download plot as a PNG":"Descargar gráfico en PNG","Zoom":"Acercar","Pan":"Desplazar",
  "Orbital rotation":"Rotación orbital","Turntable rotation":"Rotación sobre el eje",
  "Reset camera to default":"Restablecer cámara predeterminada","Reset camera to last save":"Restablecer cámara guardada",
  "Zoom in":"Acercar","Zoom out":"Alejar","Autoscale":"Ajustar escala","Reset axes":"Restablecer ejes",
  "Toggle spike lines":"Mostrar u ocultar líneas guía","Show closest data on hover":"Mostrar el dato más cercano al pasar el cursor",
  "Compare data on hover":"Comparar datos al pasar el cursor","Produced with Plotly.js":"Creado con Plotly.js"
},format:{decimal:".",thousands:",",grouping:[3]}});
const D = window.ATLAS_DATA;
const geneBySymbol = new Map(D.genes.map(g => [g.symbol.toLowerCase(), g]));
const geneById = new Map(D.genes.map(g => [g.id, g]));
const pending = new Map(), cache = new Map(), inFlight = new Map();
let currentGene = null, currentCounts = null, renderVersion = 0, sceneVersion = 0;
let cameraSyncRunning = false, queuedCamera = null, activeCameraTarget = null;
const cameraEchoes = new Map();
let lastCamera = {eye:{x:1.2,y:1.2,z:.9},up:{x:0,y:0,z:1}};
const config = {responsive:true,displaylogo:false,toImageButtonOptions:{format:"png",filename:"hypothalamus_reference_scenario",scale:2},modeBarButtonsToRemove:["lasso2d","select2d"]};
if(spanish) config.locale="es";
const colors = {FA:"#0d7475",PM:"#c05b23"};
const fnum = (v,d=3) => Number.isFinite(v) ? Number(v).toLocaleString(numberLocale,{maximumFractionDigits:d}) : tr("not estimable","no estimable");
const fp = v => !Number.isFinite(v) ? tr("unavailable","no disponible") : v < .01 ? v.toExponential(2) : v.toFixed(3);

window.receiveAtlasGene = payload => {
  const receiver = pending.get(payload.id);
  if(receiver){pending.delete(payload.id);receiver.resolve(payload);}
};

async function fetchCounts(gene){
  if(cache.has(gene.id)) return cache.get(gene.id);
  if(!("DecompressionStream" in window)) throw new Error(tr("This offline viewer requires a recent browser with gzip decompression support. The static figure is available in the report.","Este visor local requiere un navegador reciente con descompresión gzip. La figura estática está disponible en el informe."));
  const payload = await new Promise((resolve,reject)=>{
    pending.set(gene.id,{resolve,reject});
    const script=document.createElement("script");
    script.src="genes/"+encodeURIComponent(gene.id)+".js";
    script.onerror=()=>{pending.delete(gene.id);reject(new Error(tr("Gene data could not be loaded: ","No se pudieron cargar los datos del gen: ")+gene.symbol));};
    script.onload=()=>script.remove();
    document.head.appendChild(script);
  });
  const packed=Uint8Array.from(atob(payload.gzip_f32le),c=>c.charCodeAt(0));
  const stream=new Blob([packed]).stream().pipeThrough(new DecompressionStream("gzip"));
  const buffer=await new Response(stream).arrayBuffer();
  const counts=new Float32Array(buffer);
  if(counts.length!==D.cell_count||payload.n!==D.cell_count) throw new Error(tr("Gene data and coordinate lengths differ","Las longitudes de los datos del gen y de las coordenadas no coinciden"));
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
    {...D.meshes.brain,type:"mesh3d",name:tr("Whole-brain anatomy","Anatomía del cerebro completo"),color:"#b6c5ca",opacity:.075,flatshading:false,hoverinfo:"skip",showscale:false,lighting:{ambient:.8,diffuse:.6}},
    {...D.meshes.hypothalamus,type:"mesh3d",name:tr("Hypothalamus anatomy","Anatomía del hipotálamo"),color:"#6a969a",opacity:.095,flatshading:false,hoverinfo:"skip",showscale:false,lighting:{ambient:.8,diffuse:.6}}
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
  return {type:"scatter3d",mode:"markers",name:simulated?tr("Simulated HY cells","Células HY simuladas"):tr("Measured HY cells","Células HY medidas"),x:choose(D.x),y:choose(D.y),z:choose(D.z),customdata:custom,
    marker:{size:2.3,opacity:.86,color:simulated&&!available?"#b8c3c8":choose(values).map(v=>Math.log2(1+v)),colorscale:"Viridis",cmin:0,cmax:scaleMax,showscale:!simulated||available,colorbar:{title:{text:tr("log₂(1 + count)","log₂(1 + conteo)"),font:{size:10}},thickness:10,len:.6,tickfont:{size:9},x:.96}},
    hovertemplate:"%{customdata[1]}<br>"+tr("Region: ","Región: ")+"%{customdata[2]}<br>"+(simulated?(available?tr("Simulated value: %{customdata[4]:.3f}","Valor simulado: %{customdata[4]:.3f}"):tr("Exposure fit unavailable","Estimación de exposición no disponible")):tr("Measured molecules: %{customdata[3]:.0f}","Moléculas medidas: %{customdata[3]:.0f}"))+"<br>"+tr("Cell ","Célula ")+"%{customdata[0]}<extra></extra>"};
}

function updateMetrics(gene,indices,ratio){
  $('gene-name').textContent=gene.symbol;
  $('evidence').textContent=gene.significant?tr("Primary bulk FDR < 0.05 · small-n sensitive","FDR primaria regional < 0.05 · sensible a animales individuales"):gene.estimable?tr("Exploratory · not FDR-significant","Exploratorio · no significativo por FDR"):tr("Reference only · no eligible matched effect","Solo referencia · sin efecto coincidente válido");
  $('effect').textContent=gene.estimable?fnum(gene.effect.log2FoldChange_shrunk)+" log₂FC":tr("Not estimable","No estimable");
  $('fdr').textContent=tr("Bulk Wald FDR: ","FDR de Wald regional: ")+fp(gene.effect.padj);
  $('scenario-label').textContent=tr('SCENARIO MULTIPLIER · ','MULTIPLICADOR DEL ESCENARIO · ')+(gene.estimable?{central:tr('CENTRAL','CENTRAL'),lower:tr('LOWER BOUND','LÍMITE INFERIOR'),upper:tr('UPPER BOUND','LÍMITE SUPERIOR')}[$('scenario').value]:tr('UNAVAILABLE','NO DISPONIBLE'));
  $('multiplier').textContent=ratio===null?tr("Unavailable","No disponible"):fnum(ratio)+"×";
  $('interval').textContent=gene.estimable?tr("Group MLE 95% CI: [","IC del 95% de la MLE grupal: [")+fnum(gene.ci_low)+", "+fnum(gene.ci_high)+"] log₂FC":tr("No exposure effect is assigned","No se asigna un efecto de exposición");
  $('cells').textContent=indices.length.toLocaleString(numberLocale);
  $('selection-note').hidden=indices.length>0;
  const detected=indices.reduce((n,i)=>n+(currentCounts[i]>0?1:0),0);
  $('coverage').textContent=tr("of ","de ")+D.cell_count.toLocaleString(numberLocale)+tr(" HY cells · "," células HY · ")+detected.toLocaleString(numberLocale)+tr(" shown with measured count > 0"," mostradas con conteo medido > 0");
  $('scenario').disabled=!gene.estimable;
  $('simulation-title').innerHTML=gene.estimable?tr('PM2.5 scenario <span class="sim-label">SIMULATED</span>','Escenario PM2.5 <span class="sim-label">SIMULADO</span>'):tr('Exposure fit unavailable <span>Reference cells in gray</span>','Estimación de exposición no disponible <span>Células de referencia en gris</span>');
  const s=gene.leave_one_out;
  $('stability').textContent=s?tr("Sensitivity analysis: ","Análisis de sensibilidad: ")+gene.symbol+tr(" retains its effect direction in "," conserva la dirección del efecto en ")+s.runs_same_lfc_direction+tr("/7 leave-one-out fits and FDR < 0.05 in ","/7 ajustes omitiendo un animal cada vez y FDR < 0.05 en ")+s.runs_fdr_below_0_05+tr("/7. No animal was selectively removed from the primary analysis.","/7. No se retiró selectivamente ningún animal del análisis primario."):gene.estimable?tr("Exploratory genes are shown without evidence of a significant bulk exposure effect. The primary hypothesis test uses seven animals, not the atlas cells.","Los genes exploratorios se muestran sin evidencia de un efecto regional de exposición significativo. La prueba primaria utiliza siete animales, no las células del atlas."):tr("No reliable bulk exposure test is available for this gene. The measured reference is retained; no exposure effect is assigned.","No se dispone de una prueba regional de exposición fiable para este gen. Se conserva la referencia medida y no se asigna un efecto de exposición.");
}

function bulkPanel(gene){
  const traces=["FA","PM"].map(condition=>{
    const rows=gene.bulk_values.map((v,i)=>({...v,index:i})).filter(v=>v.condition===condition);
    return {type:"scatter",mode:"markers",name:condition==="FA"?tr("Filtered air (n=4)","Aire filtrado (n=4)"):"PM2.5 (n=3)",x:rows.map(v=>v.label),y:rows.map(v=>v.count),marker:{color:colors[condition],size:10,symbol:condition==="FA"?"circle":"diamond"},text:rows.map(v=>v.sample),hovertemplate:tr("%{x}<br>%{text}<br>Normalized count: %{y:.2f}<extra></extra>","%{x}<br>%{text}<br>Conteo normalizado: %{y:.2f}<extra></extra>")};
  });
  Plotly.react('bulk',traces,{margin:{t:15,l:65,b:45,r:15},height:280,font:{family:"Arial,sans-serif",size:11,color:"#173239"},paper_bgcolor:"white",plot_bgcolor:"white",xaxis:{title:tr("Independent bulk animal","Animal independiente"),categoryorder:"array",categoryarray:["FA1","FA2","FA3","FA4","PM1","PM2","PM3"]},yaxis:{title:tr("DESeq2-normalized counts","Conteos normalizados por DESeq2"),rangemode:"tozero",gridcolor:"#e8eff0"},legend:{orientation:"h",y:1.18}},config);
  const table=document.createElement('table');
  table.innerHTML=tr('<thead><tr><th>Animal</th><th>Group</th><th>Normalized count</th></tr></thead>','<thead><tr><th>Animal</th><th>Grupo</th><th>Conteo normalizado</th></tr></thead>');
  const body=document.createElement('tbody');
  for(const v of gene.bulk_values){const tr=document.createElement('tr');for(const text of [v.label,v.condition,fnum(v.count,2)]){const td=document.createElement('td');td.textContent=text;tr.appendChild(td);}body.appendChild(tr);}
  table.appendChild(body);$('bulk-table').replaceChildren(table);
}

async function render(){
  if(!currentGene||!currentCounts)return;
  const revision=++sceneVersion;
  const gene=currentGene,indices=selectedIndices();
  const beta=exponent(gene),ratio=beta===null?null:2**beta;
  if(ratio!==null&&!Number.isFinite(ratio))throw new Error(tr("The selected coefficient interval is numerically unstable. Use the central scenario.","El intervalo seleccionado del coeficiente es numéricamente inestable. Use el escenario central."));
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
  $('scale-note').textContent=tr("Common color scale: log₂(1 + imaging count), 0–","Escala de color común: log₂(1 + conteo de imagen), 0–")+fnum(scaleMax)+". "+($('clip').checked?clipped.toLocaleString(numberLocale)+tr(" values across the full HY reference/scenario exceed the common positive-value 99th-percentile cap; values are retained in hover and exports."," valores en la referencia y el escenario HY completos superan el límite común del percentil 99 de valores positivos; los valores se conservan al pasar el cursor y en las exportaciones."):tr("Full reference/scenario range; no color capping.","Rango completo de la referencia y el escenario; sin límite de color."))+tr(" Geometry is unchanged."," La geometría no cambia.");
  if(!gene.estimable)$('scale-note').textContent=tr("Reference color scale: log₂(1 + imaging count), 0–","Escala de color de referencia: log₂(1 + conteo de imagen), 0–")+fnum(scaleMax)+tr(". The right panel shows reference positions in gray; no exposure simulation is available.",". El panel derecho muestra las posiciones de referencia en gris; no se dispone de una simulación de exposición.");
  if(window.parent!==window)window.parent.postMessage({type:"omics-atlas-height",height:document.documentElement.scrollHeight},"*");
}

async function applyGene(){
  const gene=geneBySymbol.get($('gene').value.trim().toLowerCase())||geneById.get($('gene').value.trim());
  if(!gene){$('loading').textContent=tr("Choose a symbol or identifier from the measured panel.","Elija un símbolo o identificador del panel medido.");$('loading').classList.add('error');return;}
  const version=++renderVersion;
  $('loading').textContent=tr("Loading ","Cargando ")+gene.symbol+"…";$('loading').classList.remove('error');$('apply').disabled=true;
  try{
    const counts=await loadCounts(gene);
    if(version!==renderVersion)return;
    currentGene=gene;currentCounts=counts;$('gene').value=gene.symbol;
    await render();
    if(version!==renderVersion)return;
    bulkPanel(gene);$('loading').textContent=gene.estimable?tr("Ready · measured panel / explicit scenario","Listo · panel medido / escenario explícito"):tr("Ready · measured reference only","Listo · solo referencia medida");
  }catch(error){$('loading').textContent=error.message;$('loading').classList.add('error');console.error(error);}
  finally{if(version===renderVersion)$('apply').disabled=false;}
}

function sameCamera(a,b){
  if(!a||!b)return false;
  const defaults={eye:{x:1.25,y:1.25,z:1.25},center:{x:0,y:0,z:0},up:{x:0,y:0,z:1}};
  return Object.entries(defaults).every(([part,values])=>Object.keys(values).every(axis=>
    Math.abs((a[part]?.[axis]??values[axis])-(b[part]?.[axis]??values[axis]))<1e-12
  ))&&(a.projection?.type??'perspective')===(b.projection?.type??'perspective');
}

async function flushCameraQueue(){
  if(cameraSyncRunning)return;
  cameraSyncRunning=true;
  try{
    while(queuedCamera){
      const update=queuedCamera;queuedCamera=null;
      for(const target of update.targets){
        // A newer gesture supersedes any remaining writes from this update.
        if(queuedCamera)break;
        activeCameraTarget=target;cameraEchoes.set(target,update.camera);
        try{await Plotly.relayout(target,{'scene.camera':update.camera});}
        catch(error){console.error(error);}
        finally{cameraEchoes.delete(target);activeCameraTarget=null;}
      }
    }
  }finally{cameraSyncRunning=false;}
}

function linkCameras(source,target){
  $(source).on('plotly_relayout',event=>{
    const camera=event['scene.camera'];
    if(!camera||sameCamera(camera,cameraEchoes.get(source))||sameCamera(camera,lastCamera))return;
    lastCamera=JSON.parse(JSON.stringify(camera));
    // Keep only the latest gesture. If an older write is still targeting its
    // source, restore that source too before synchronizing the opposite view.
    const targets=activeCameraTarget===source?[source,target]:[target];
    queuedCamera={camera:lastCamera,targets};
    void flushCameraQueue();
  });
}

async function init(){
  for(const g of [...D.genes].sort((a,b)=>a.symbol.localeCompare(b.symbol))){const option=document.createElement('option');option.value=g.symbol;option.textContent=g.id+(g.significant?tr(" · primary FDR hit"," · hallazgo primario por FDR"):"");$('gene-list').appendChild(option);}
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
