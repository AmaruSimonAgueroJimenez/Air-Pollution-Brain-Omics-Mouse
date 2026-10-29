/* Compact measured-reference viewer for the bilingual HTML defense.
 * Inputs are the existing regional atlas exports. No cells are added or moved.
 * Only this document handles map controls; presentation navigation stays outside.
 */
(function () {
  'use strict';
  const params = new URLSearchParams(window.location.search);
  const spanish = params.get('lang') === 'es';
  const REGION = params.get('region') === 'ACA' ? 'ACA' : 'HY';
  const instance = params.get('instance') || '';
  const directory = REGION === 'ACA' ? '../viewer-cortex/' : '../viewer/';
  const $ = id => document.getElementById(id);
  const tr = (en, es) => spanish ? es : en;
  const num = (value, digits=3) => Number.isFinite(value) ? value.toLocaleString('en', {maximumFractionDigits:digits}) : tr('unavailable', 'no disponible');
  const pnum = value => !Number.isFinite(value) ? tr('unavailable','no disponible') : value < .01 ? value.toExponential(2) : value.toFixed(3);
  const pending = new Map(), cache = new Map(), inFlight = new Map();
  let D, annotations, priorityGroups, geneBySymbol, geneById;
  let currentGene = null, currentCounts = null, renderVersion = 0, sceneVersion = 0;
  let renderQueue = Promise.resolve(), plotsReady = false;
  const defaultCamera = () => ({eye:{x:1.2,y:1.2,z:.9},up:{x:0,y:0,z:1},projection:{type:'orthographic'}});
  let lastCamera = defaultCamera();
  let cameraSyncRunning = false, queuedCamera = null, activeCameraTarget = null;
  const cameraEchoes = new Map();
  let state = {region:REGION, gene:null, status:'loading'};

  function parentOrigin() {
    const inherited = window.location.ancestorOrigins && window.location.ancestorOrigins[0];
    if (inherited) return inherited;
    try { if (document.referrer) return new URL(document.referrer).origin; } catch (_) {}
    return window.location.origin;
  }
  const expectedParentOrigin = parentOrigin();
  function publish(status, extra={}) {
    state = {...state, status, ...extra};
    document.documentElement.dataset.mapStatus = status;
    document.documentElement.dataset.mapRegion = REGION;
    document.documentElement.dataset.mapGene = state.gene || '';
    if (window.parent !== window) {
      window.parent.postMessage({type:'omics-slide-map-state', instance, ...state}, expectedParentOrigin === 'null' ? '*' : expectedParentOrigin);
    }
  }
  function fail(error) {
    const message = error instanceof Error ? error.message : String(error);
    $('error').textContent = tr('Interactive map unavailable. ', 'Mapa interactivo no disponible. ') + message;
    $('error').hidden = false;
    publish('error', {message});
    console.error(error);
  }
  function script(url) {
    return new Promise((resolve,reject) => {
      const tag = document.createElement('script');
      tag.src = url;
      tag.onload = () => {tag.remove();resolve();};
      tag.onerror = () => {tag.remove();reject(new Error(tr('A local viewer asset could not be loaded.', 'No se pudo cargar un recurso local del visor.')));};
      document.head.appendChild(tag);
    });
  }
  window.receiveAtlasGene = payload => {
    const receiver = pending.get(payload.id);
    if (!receiver) return;
    pending.delete(payload.id);receiver.resolve(payload);
  };
  async function fetchCounts(gene) {
    if (!('DecompressionStream' in window)) throw new Error(tr('Use a recent browser with gzip decompression support.', 'Use un navegador reciente con descompresión gzip.'));
    const payload = await new Promise((resolve,reject) => {
      const tag = document.createElement('script');
      pending.set(gene.id, {resolve,reject});
      tag.src = directory+'genes/'+encodeURIComponent(gene.id)+'.js';
      tag.onerror = () => {pending.delete(gene.id);tag.remove();reject(new Error(tr('Measured gene data could not be loaded: ', 'No se pudieron cargar los datos medidos de ')+gene.symbol));};
      tag.onload = () => {
        tag.remove();
        if (pending.has(gene.id)) {pending.delete(gene.id);reject(new Error(tr('Gene payload did not match the selected gene.', 'Los datos cargados no corresponden al gen seleccionado.')));}
      };
      document.head.appendChild(tag);
    });
    const packed = Uint8Array.from(atob(payload.gzip_f32le), c => c.charCodeAt(0));
    const stream = new Blob([packed]).stream().pipeThrough(new DecompressionStream('gzip'));
    const values = new Float32Array(await new Response(stream).arrayBuffer());
    if (values.length !== D.cell_count || payload.n !== D.cell_count) throw new Error(tr('Coordinate and count lengths differ.', 'Las longitudes de coordenadas y conteos difieren.'));
    if (values.some(value => !Number.isFinite(value) || value < 0)) throw new Error(tr('Measured imaging counts must be finite and nonnegative.', 'Los conteos de imagen medidos deben ser finitos y no negativos.'));
    cache.set(gene.id, values);
    if (cache.size > 10) cache.delete(cache.keys().next().value);
    return values;
  }
  function loadCounts(gene) {
    if (cache.has(gene.id)) return Promise.resolve(cache.get(gene.id));
    if (inFlight.has(gene.id)) return inFlight.get(gene.id);
    const operation = fetchCounts(gene).finally(() => inFlight.delete(gene.id));
    inFlight.set(gene.id,operation);return operation;
  }
  function bounds(values) {
    let min=Infinity,max=-Infinity;
    for (const v of values) {if(v<min)min=v;if(v>max)max=v;}
    return [min-.3,max+.3];
  }
  function layout() {
    const r = $('brain-context').checked ? {x:[0,11.4],y:[0,13.2],z:[-8.2,.1]} : {x:bounds(D.x),y:bounds(D.y),z:bounds(D.z)};
    const sizes = [r.x[1]-r.x[0],r.y[1]-r.y[0],r.z[1]-r.z[0]], largest = Math.max(...sizes);
    const axis = (title,range) => ({title:{text:title,font:{size:13}},range,backgroundcolor:'#fff',gridcolor:'#e3ebea',zeroline:false,tickfont:{size:12},showspikes:false});
    return {margin:{l:0,r:0,t:0,b:0},paper_bgcolor:'#fff',font:{family:'Arial,sans-serif',color:'#142c35'},showlegend:false,uirevision:$('brain-context').checked?'brain':'region',
      scene:{camera:lastCamera,aspectmode:'manual',aspectratio:{x:sizes[0]/largest,y:sizes[1]/largest,z:sizes[2]/largest},
        xaxis:axis('ML (mm)',r.x),yaxis:axis('AP (mm)',r.y),zaxis:axis('−DV (mm)',r.z)}};
  }
  function meshes() {
    if (!$('brain-context').checked) return [];
    return [
      {...D.meshes.brain,type:'mesh3d',color:'#b6c5ca',opacity:.075,flatshading:false,hoverinfo:'skip',showscale:false},
      {...(D.meshes.region||D.meshes.hypothalamus),type:'mesh3d',color:'#6a969a',opacity:.095,flatshading:false,hoverinfo:'skip',showscale:false}
    ];
  }
  function trace(indices, counts, values, scaleMax, simulated, available, zeros) {
    const choose = a => indices.map(i=>a[i]);
    const custom = indices.map(i=>[D.cell_ids[i],D.type_labels[D.type_codes[i]],D.region_labels[D.region_codes[i]],counts[i],values[i]]);
    const colored = available && !zeros;
    return {type:'scatter3d',mode:'markers',x:choose(D.x),y:choose(D.y),z:choose(D.z),customdata:custom,
      marker:{size:zeros?1.5:2.8,opacity:zeros ? .18 : .9,color:colored?choose(values).map(v=>Math.log2(1+v)):'#b4bec1',colorscale:'Viridis',cmin:0,cmax:scaleMax,showscale:colored,
        colorbar:{title:{text:tr('log₂(1 + count)','log₂(1 + conteo)'),font:{size:12}},thickness:12,len:.65,tickfont:{size:11},x:.97}},
      hovertemplate:'%{customdata[1]}<br>'+tr('Anatomy: ','Anatomía: ')+'%{customdata[2]}<br>'+(simulated?(available?tr('Illustrative value: %{customdata[4]:.3f}','Valor ilustrativo: %{customdata[4]:.3f}'):tr('No eligible exposure estimate','Sin estimación de exposición válida')):tr('Measured count: %{customdata[3]:.0f}','Conteo medido: %{customdata[3]:.0f}'))+'<br>'+tr('Cell ','Célula ')+'%{customdata[0]}<extra></extra>'};
  }
  function evidence(gene) {
    if (!gene.estimable) return tr('Reference only · no eligible regional effect','Solo referencia · sin efecto regional válido');
    return gene.significant ? tr('Primary BH discovery','Hallazgo BH primario') : tr('Exploratory · not BH-significant','Exploratorio · no significativo por BH');
  }
  function product(gene) {
    const a = annotations[gene.symbol] || {};
    return (spanish?a.protein_name_es:a.protein_name_en) || a.protein_name_en || gene.name || tr('Product annotation unavailable','Producto no anotado');
  }
  function updateSummary(gene, counts, ratio, positives, shown) {
    $('current-gene').textContent=gene.symbol;
    $('product').textContent=product(gene);$('product').title=product(gene);
    $('evidence').textContent=evidence(gene);
    $('effect').textContent=gene.estimable ? 'q='+pnum(gene.effect?.padj)+' · '+num(gene.effect.log2FoldChange_shrunk)+' log₂FC · ×'+num(ratio,5) : tr('No PM effect assigned','Sin efecto PM asignado');
    $('cell-summary').textContent=positives.toLocaleString('en')+'/'+D.cell_count.toLocaleString('en')+tr(' measured positive cells',' células medidas positivas')+($('hide-zero').checked?' · '+shown.toLocaleString('en')+tr(' plotted',' mostradas'):'');
    $('scenario-title').textContent=gene.estimable ? tr('Illustrative uniform PM₂.₅ scenario','Escenario PM₂.₅ uniforme e ilustrativo')+' ×'+num(ratio,5) : tr('No exposure estimate · reference positions in gray','Sin estimación de exposición · posiciones en gris');
    for (const b of $('gene-shortcuts').querySelectorAll('button')) b.setAttribute('aria-pressed',String(b.dataset.gene===gene.symbol));
    $('summary').title = gene.estimable ? tr('The coefficient comes from a separate regional animal experiment. It is applied uniformly, not estimated per cell.','El coeficiente procede de otro experimento regional con animales. Se aplica uniformemente y no se estima por célula.') : evidence(gene);
  }
  async function renderNow(revision) {
    if (revision !== sceneVersion || !currentGene || !currentCounts) return;
    const gene=currentGene,counts=currentCounts;
    const beta=gene.estimable?gene.effect?.log2FoldChange_shrunk:null;
    const ratio=beta===null?null:2**beta;
    if (gene.estimable && (!Number.isFinite(beta)||!Number.isFinite(ratio)||ratio<=0)) throw new Error(tr('The regional scenario is numerically unavailable for this gene.','El escenario regional no está disponible numéricamente para este gen.'));
    const values=Array.from(counts,v=>v*(ratio===null?1:ratio));
    const positive=[],zeros=[],allPositive=[];
    let scaleMax=0;
    for(let i=0;i<D.cell_count;i++) {
      if (counts[i]>0) {positive.push(i);allPositive.push(Math.log2(1+counts[i]));if(gene.estimable)allPositive.push(Math.log2(1+values[i]));}
      else zeros.push(i);
    }
    if(allPositive.length) {
      if($('clip').checked){allPositive.sort((a,b)=>a-b);scaleMax=allPositive[Math.floor((allPositive.length-1)*.99)];}
      else for(const v of allPositive)if(v>scaleMax)scaleMax=v;
    }
    const displayScaleMax=scaleMax>0?scaleMax:1;
    const showZeros=!$('hide-zero').checked;
    const shared=meshes();
    const ref=[...shared],sim=[...shared];
    if(showZeros&&zeros.length){ref.push(trace(zeros,counts,counts,displayScaleMax,false,true,true));sim.push(trace(zeros,counts,values,displayScaleMax,true,gene.estimable,true));}
    if(positive.length){ref.push(trace(positive,counts,counts,displayScaleMax,false,true,false));sim.push(trace(positive,counts,values,displayScaleMax,true,gene.estimable,false));}
    // Retain explicit empty scatter traces so an all-zero or hidden-zero selection
    // still gets a rotatable scene with the original registered coordinate limits.
    if(!ref.length)ref.push({type:'scatter3d',x:[],y:[],z:[],mode:'markers'});
    if(!sim.length)sim.push({type:'scatter3d',x:[],y:[],z:[],mode:'markers'});
    updateSummary(gene,counts,ratio,positive.length,positive.length+(showZeros?zeros.length:0));
    await Promise.all([Plotly.react('reference',ref,layout(),config),Plotly.react('simulation',sim,layout(),config)]);
    if(revision!==sceneVersion||gene!==currentGene)return;
    if(!plotsReady){plotsReady=true;linkCameras('reference','simulation');linkCameras('simulation','reference');}
    $('scale-note').textContent=tr('Shared scale: ','Escala común: ')+'log₂(1 + '+tr('imaging count','conteo de imagen')+'), 0–'+num(displayScaleMax)+'. '+(allPositive.length?($('clip').checked?tr('99th-percentile color cap; values are preserved.','Límite de color al percentil 99; se conservan los valores.'):tr('Full range, no color cap.','Rango completo, sin límite de color.')):tr('All measured counts are zero.','Todos los conteos medidos son cero.'))+' '+(gene.estimable?tr('Gray dots: measured zeros.','Puntos grises: ceros medidos.'):tr('Reference gray dots: measured zeros.','Puntos grises en la referencia: ceros medidos.'))+(!gene.estimable?' '+tr('The right scene has no exposure effect.','La escena derecha no tiene efecto de exposición.'):'');
    $('error').hidden=true;
    publish('ready',{gene:gene.symbol,cellCount:D.cell_count,shownCells:positive.length+(showZeros?zeros.length:0),positiveCells:positive.length,multiplier:ratio,significant:!!gene.significant,estimable:!!gene.estimable,q:gene.effect?.padj??null,beta,scaleMax:displayScaleMax,colorCap:$('clip').checked,hideZero:$('hide-zero').checked,coordinateFrame:'ML, AP, −DV in CCF mm'});
  }
  function render() {
    const revision=++sceneVersion;
    renderQueue=renderQueue.catch(()=>{}).then(()=>renderNow(revision));
    return renderQueue;
  }
  async function applyGene() {
    const raw=$('gene').value.trim();
    const gene=geneBySymbol.get(raw.toLowerCase())||geneById.get(raw);
    if(!gene){$('error').textContent=tr('Choose a symbol from the measured gene panel.','Elija un símbolo del panel de genes medidos.');$('error').hidden=false;return;}
    const revision=++renderVersion;
    $('apply').disabled=true;$('error').hidden=true;
    $('current-gene').textContent=tr('Loading ','Cargando ')+gene.symbol+'…';
    publish('loading',{requestedGene:gene.symbol});
    try{
      const counts=await loadCounts(gene);
      if(revision!==renderVersion)return;
      currentGene=gene;currentCounts=counts;$('gene').value=gene.symbol;
      await render();
    } catch(error){if(revision===renderVersion)fail(error);}
    finally{if(revision===renderVersion)$('apply').disabled=false;}
  }
  function sameCamera(a,b) {
    if(!a||!b)return false;
    const defaults={eye:{x:1.25,y:1.25,z:1.25},center:{x:0,y:0,z:0},up:{x:0,y:0,z:1}};
    return Object.entries(defaults).every(([part,axes])=>Object.keys(axes).every(axis=>Math.abs((a[part]?.[axis]??axes[axis])-(b[part]?.[axis]??axes[axis]))<1e-12))&&(a.projection?.type??'perspective')===(b.projection?.type??'perspective');
  }
  async function flushCameraQueue() {
    if(cameraSyncRunning)return;
    cameraSyncRunning=true;
    try{while(queuedCamera){const update=queuedCamera;queuedCamera=null;for(const target of update.targets){if(queuedCamera)break;activeCameraTarget=target;cameraEchoes.set(target,update.camera);try{await Plotly.relayout(target,{'scene.camera':update.camera});}finally{cameraEchoes.delete(target);activeCameraTarget=null;}}}}
    catch(error){console.error(error);}
    finally{cameraSyncRunning=false;}
  }
  function linkCameras(source,target) {
    $(source).on('plotly_relayout',event=>{
      const camera=event['scene.camera'];
      if(!camera||sameCamera(camera,cameraEchoes.get(source))||sameCamera(camera,lastCamera))return;
      lastCamera=JSON.parse(JSON.stringify(camera));
      queuedCamera={camera:lastCamera,targets:activeCameraTarget===source?[source,target]:[target]};
      void flushCameraQueue();
    });
  }
  function resizePlots() {
    if(!plotsReady)return;
    for(const id of ['reference','simulation'])if($(id).getBoundingClientRect().width>2&&$(id).getBoundingClientRect().height>2)Plotly.Plots.resize($(id));
  }
  window.addEventListener('resize',resizePlots);
  window.addEventListener('message',event=>{
    if(event.source!==window.parent||event.origin!==expectedParentOrigin)return;
    if(event.data?.type==='omics-slide-map-resize')resizePlots();
  });
  const config={responsive:true,displaylogo:false,displayModeBar:true,toImageButtonOptions:{format:'png',filename:'brain_reference_scenario',scale:2},modeBarButtonsToRemove:['lasso2d','select2d','resetCameraLastSave3d']};
  async function init() {
    document.documentElement.lang=spanish?'es':'en';
    document.title=tr('Interactive brain reference and PM₂.₅ scenario','Referencia cerebral interactiva y escenario PM₂.₅');
    $('map-app').setAttribute('aria-label',REGION+' · '+tr('interactive measured brain maps','mapas cerebrales medidos interactivos'));
    $('gene-label').textContent=tr('Measured gene','Gen medido');
    $('apply').textContent=tr('Apply gene','Aplicar gen');$('reset').textContent=tr('Reset cameras','Restablecer cámaras');
    $('context-label').textContent=tr('Whole-brain context','Contexto cerebral');$('zero-label').textContent=tr('Hide zeros','Ocultar ceros');$('clip-label').textContent=tr('99% color cap','Límite de color 99%');
    $('priority-label').textContent=tr('Priority genes','Genes prioritarios');$('gene-shortcuts').setAttribute('aria-label',tr('Prioritized measured genes','Genes medidos priorizados'));
    $('reference-title').textContent=tr('Measured MERFISH reference','Referencia MERFISH medida');
    $('reference').setAttribute('aria-label',tr('Measured reference in three dimensions','Referencia medida en tres dimensiones'));
    $('simulation').setAttribute('aria-label',tr('Illustrative regional PM2.5 scenario in three dimensions','Escenario regional ilustrativo de PM2.5 en tres dimensiones'));
    $('panel-help').textContent=tr('One unexposed donor. Drag either map to rotate both. The uniform scenario preserves positions and detection and does not identify cellular exposure responses.','Un donante sin exposición. Arrastre un mapa para rotar ambos. El escenario uniforme conserva posiciones y detección, sin identificar respuestas celulares a la exposición.');
    await Promise.all([script('../viewer/plotly.min.js'),script(directory+'data.js'),script('../viewer/annotations.js')]);
    D=window.ATLAS_DATA;annotations=window.GENE_ANNOTATIONS||{};priorityGroups=window.GENE_PRIORITIES||{};
    if(!D||!Array.isArray(D.genes)||D.x.length!==D.cell_count||D.y.length!==D.cell_count||D.z.length!==D.cell_count)throw new Error(tr('Invalid measured atlas data.','Los datos medidos del atlas no son válidos.'));
    geneBySymbol=new Map(D.genes.map(g=>[g.symbol.toLowerCase(),g]));geneById=new Map(D.genes.map(g=>[g.id,g]));
    if(spanish){Plotly.register({moduleType:'locale',name:'es',dictionary:{'Download plot as a PNG':'Descargar gráfico en PNG','Zoom':'Acercar','Pan':'Desplazar','Orbital rotation':'Rotación orbital','Turntable rotation':'Rotación sobre el eje','Reset camera to default':'Restablecer cámara predeterminada'},format:{decimal:'.',thousands:',',grouping:[3]}});config.locale='es';}
    for(const gene of [...D.genes].sort((a,b)=>a.symbol.localeCompare(b.symbol))){const option=document.createElement('option');option.value=gene.symbol;option.textContent=gene.id;$('gene-list').appendChild(option);}
    const priorities=[],seen=new Set();
    const entries=priorityGroups[REGION]||[];
    for(const entry of entries){const g=geneBySymbol.get(String(typeof entry==='string'?entry:entry.symbol).toLowerCase());if(g&&!seen.has(g.id)){priorities.push(g);seen.add(g.id);}}
    for(const gene of [...D.genes].filter(g=>g.estimable).sort((a,b)=>(a.effect?.padj??Infinity)-(b.effect?.padj??Infinity))){if(priorities.length>=10)break;if(!seen.has(gene.id)){priorities.push(gene);seen.add(gene.id);}}
    for(const gene of priorities.slice(0,10)){const button=document.createElement('button');button.type='button';button.dataset.gene=gene.symbol;button.textContent=gene.symbol;button.title=evidence(gene)+' · q='+pnum(gene.effect?.padj);button.setAttribute('aria-pressed','false');button.addEventListener('click',()=>{$('gene').value=gene.symbol;void applyGene();});$('gene-shortcuts').appendChild(button);}
    $('controls').addEventListener('submit',event=>{event.preventDefault();void applyGene();});
    for(const id of ['brain-context','hide-zero','clip'])$(id).addEventListener('change',()=>{if(id==='brain-context')lastCamera=defaultCamera();void render().catch(fail);});
    $('reset').addEventListener('click',()=>{
      lastCamera=defaultCamera();
      if(plotsReady){queuedCamera={camera:JSON.parse(JSON.stringify(lastCamera)),targets:['reference','simulation']};void flushCameraQueue();}
      else void render().catch(fail);
    });
    const requested=params.get('defaultgene')||D.default_gene;
    $('gene').value=geneBySymbol.get(requested.toLowerCase())?.symbol||D.default_gene;
    await applyGene();
    if(window.ResizeObserver)new ResizeObserver(()=>window.requestAnimationFrame(resizePlots)).observe($('plot-grid'));
  }
  init().catch(fail);
}());
