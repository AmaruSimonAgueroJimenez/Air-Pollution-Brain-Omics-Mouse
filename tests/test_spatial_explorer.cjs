const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');
const crypto=require('node:crypto');
const C=require('../docs/assets/spatial-explorer/explorer.js');
const root=path.join(__dirname,'..');
const fixture={cell_count:6,cell_ids:['a','b','c','d','e','f'],experimental_x:[0,.02,.13,.15,.21,.3],experimental_y:[1,1.02,1,1,1,1],section_codes:[0,0,0,0,1,1],anatomy_codes:[0,0,1,1,0,1],anatomy_labels:['L2/3','L5'],type_codes:[0,1,0,1,0,1],bin_size_mm:.1};
const counts=[0,2,0,4,0,8];

test('measured detection denominator includes zeros after section, layer and class filters',()=>{
  const all=C.selectCells(fixture,0),a=C.summarize(counts,all);
  assert.deepEqual(all,[0,1,2,3]);assert.equal(a.n,4);assert.equal(a.positive,2);assert.equal(a.fraction,.5);assert.equal(a.mean,1.5);
  const layer=C.summarize(counts,C.selectCells(fixture,0,1));assert.equal(layer.n,2);assert.equal(layer.positive,1);assert.equal(layer.mean,2);
  const classZero=C.summarize(counts,C.selectCells(fixture,0,'all',0));assert.equal(classZero.n,2);assert.equal(classZero.positive,0);assert.equal(classZero.fraction,0);
  const empty=C.summarize(counts,C.selectCells(fixture,0,1,9));assert.equal(empty.n,0);assert.equal(empty.mean,null);assert.equal(empty.fraction,null);
});

test('occupied spatial bins report actual positive/all cells rather than area density',()=>{
  const bins=C.detectionBins(fixture,counts,C.selectCells(fixture,0));
  assert.equal(bins.length,2);assert.deepEqual(bins.map(b=>[b.positive,b.n,b.fraction,b.mean]),[[1,2,.5,1],[1,2,.5,2]]);
  assert.deepEqual(bins.map(b=>b.x),[.05,.15000000000000002]);
  assert.equal(bins.reduce((n,b)=>n+b.n,0),4);
});

test('anatomy summary retains both positive and zero cells',()=>{
  const groups=C.groupSummary(fixture,counts,C.selectCells(fixture,0));
  assert.deepEqual(groups.map(g=>[g.label,g.n,g.positive,g.mean]),[['L2/3',2,1,1],['L5',2,1,2]]);
});

test('Bayesian predictive lower bounds remain negative in natural-log units',()=>{
  const payload={posterior_mean:[.1],ci_lower:[-.05],ci_upper:[.25],predictive_lower:[-.9],predictive_upper:[1.1]};
  const result=C.validateModel({cell_count:1},payload);assert.equal(result.predictive_lower[0],-.9);assert.equal(result.posterior_mean[0],.1);
  assert.throws(()=>C.validateModel({cell_count:2},payload),/alignment differs/);
  assert.throws(()=>C.validateModel({cell_count:1},{...payload,ci_lower:[.2]}),/interval ordering/);
});

test('measurement alignment rejects missing cells and negative observed counts',()=>{
  assert.throws(()=>C.validateCounts(fixture,{n:5,counts}),/alignment differs/);
  assert.throws(()=>C.validateCounts(fixture,{n:6,counts:[0,2,0,4,-1,8]}),/Invalid measured/);
});

test('download CSV preserves identifiers and literal values without implicit transforms',()=>{
  const csv=C.toCSV([['id','lower'],['cell,one',-.9],['cell"two',0]]);
  assert.equal(csv,'id,lower\r\n"cell,one",-0.9\r\n"cell""two",0\r\n');
});

test('Bayesian detection accepts bounded probability intervals only for matching measured groups',()=>{
  const data={...fixture,sections:[{label:'s0'},{label:'s1'}]},section=(sc,ac,n,positive)=>({section_code:sc,label:'s'+sc,anatomy_code:ac,n,positive,observed_fraction:positive/n,p_mean:.4,p_ci_low:.03,p_ci_high:.9});
  const gene={sections:[section(0,0,2,1),section(0,1,2,1),section(1,0,1,0),section(1,1,1,1)],anatomy:[{mu_mean:.4,mu_ci_low:.1,mu_ci_high:.8,new_section_p_ci_low:0,new_section_p_ci_high:.99}]};
  assert.equal(C.validateDetection(data,gene,counts),gene);
  assert.throws(()=>C.validateDetection(data,{...gene,sections:[{...gene.sections[0],n:3},...gene.sections.slice(1)]},counts),/fraction differs|counts differ/);
  assert.throws(()=>C.validateDetection(data,{...gene,sections:[{...gene.sections[0],p_ci_low:-.1},...gene.sections.slice(1)]},counts),/not bounded/);
  assert.throws(()=>C.validateDetection(data,{...gene,sections:[...gene.sections,gene.sections[0]]},counts),/Duplicate/);
});

for(const region of ['ACA','HY'])test(region+' real assets reconcile raw-cell, section, anatomy and occupied-bin counts',()=>{
  const context={window:{}};vm.runInNewContext(fs.readFileSync(path.join(root,'docs/assets/spatial-sections/data-'+region+'.js'),'utf8'),context);
  const data=context.window['SPATIAL_SECTIONS_'+region];assert.equal(data.cell_ids.length,data.cell_count);assert.equal(new Set(data.cell_ids).size,data.cell_count);assert.equal(data.independent_reference_animals,1);
  assert.equal(data.genes.length,10);assert.equal(data.sections.reduce((n,s)=>n+s.n,0),data.cell_count);
  assert.equal(data.bins.reduce((n,b)=>n+b.n,0),data.cell_count);
  assert.equal(data.anatomy_kind,region==='ACA'?'layer':'nucleus_subregion');
  for(const gene of data.genes){let payload;const window={registerSpatialGene:(r,id,p)=>{assert.equal(r,region);assert.equal(id,gene.id);payload=p;}};
    vm.runInNewContext(fs.readFileSync(path.join(root,'docs/assets/spatial-sections/genes-'+region,gene.id+'.js'),'utf8'),{window});
    const values=C.validateCounts(data,payload),indices=Array.from({length:data.cell_count},(_,i)=>i),summary=C.summarize(values,indices);
    assert.equal(summary.positive,gene.n_positive,gene.symbol+' regional positives');
    let sectionsPositive=0,sectionsN=0,binsPositive=0,binsN=0;
    for(let section=0;section<data.sections.length;section++){
      const selected=C.selectCells(data,section),s=C.summarize(values,selected);assert.equal(s.n,data.sections[section].n);sectionsPositive+=s.positive;sectionsN+=s.n;
      const groups=C.groupSummary(data,values,selected);assert.equal(groups.reduce((n,g)=>n+g.n,0),s.n);assert.equal(groups.reduce((n,g)=>n+g.positive,0),s.positive);
      for(const b of C.detectionBins(data,values,selected)){assert.equal(b.n,data.bins[b.key].n);assert.equal(b.n,payload.bins.n[b.key]);assert.equal(b.positive,payload.bins.positive[b.key]);binsPositive+=b.positive;binsN+=b.n;}
    }
    assert.equal(sectionsN,data.cell_count);assert.equal(binsN,data.cell_count);assert.equal(sectionsPositive,summary.positive);assert.equal(binsPositive,summary.positive);
  }
});

for(const region of ['ACA','HY'])test(region+' actual baseline arrays preserve cell identity, interval ordering and reported positive-tail diagnostics',()=>{
  const context={window:{}};vm.runInNewContext(fs.readFileSync(path.join(root,'docs/assets/spatial-sections/data-'+region+'.js'),'utf8'),context);
  const data=context.window['SPATIAL_SECTIONS_'+region],directory=path.join(root,'docs/assets/spatial-bayes',region),meta=JSON.parse(fs.readFileSync(path.join(directory,'metadata.json'),'utf8'));
  const expectedHash=crypto.createHash('sha256').update(data.cell_ids.join('\n')+'\n').digest('hex');
  assert.equal(meta.cell_ids_sha256,expectedHash);assert.equal(meta.cell_count,data.cell_count);assert.equal(meta.donor_n,1);assert.match(meta.units,/ln\(1\+raw imaging count\)/);
  assert.equal(meta.cross_validation.unit,'whole physical brain section');assert.equal(meta.cross_validation.donors,1);
  for(const gene of meta.genes){const payload=JSON.parse(fs.readFileSync(path.join(directory,gene.file),'utf8')),pred=C.validateModel(data,payload);
    assert.equal(payload.id,gene.id);assert.equal(payload.region,region);assert.equal(payload.n,data.cell_count);assert.equal(payload.cell_ids_sha256,expectedHash);
    assert.ok(pred.predictive_lower.some(v=>v<0),'Gaussian support mismatch must not be silently clipped');
    assert.ok(gene.heldout.coverage_positive>=0&&gene.heldout.coverage_positive<=1);assert.ok(gene.heldout.flags.includes('nonnegative_support_mismatch'));
    assert.equal(gene.by_section.reduce((n,s)=>n+s.n_cells,0),data.cell_count);
  }
});

for(const region of ['ACA','HY'])test(region+' actual bounded detection fits reconcile all measured groups and disclose prior-only held-out anatomy',()=>{
  const context={window:{}};vm.runInNewContext(fs.readFileSync(path.join(root,'docs/assets/spatial-sections/data-'+region+'.js'),'utf8'),context);
  vm.runInNewContext(fs.readFileSync(path.join(root,'docs/assets/spatial-detection-bayes/data-'+region+'.js'),'utf8'),context);
  const data=context.window['SPATIAL_SECTIONS_'+region],model=context.window['SPATIAL_DETECTION_BAYES_'+region];
  assert.equal(model.independent_reference_animals,1);assert.equal(model.cellclass_specific,false);assert.equal(model.baseline_only,true);assert.equal(model.no_PM_effect,true);assert.equal(model.counts_unchanged,true);assert.equal(model.genes.length,10);
  for(const gene of model.genes){let payload;vm.runInNewContext(fs.readFileSync(path.join(root,'docs/assets/spatial-sections/genes-'+region,gene.id+'.js'),'utf8'),{window:{registerSpatialGene:(r,id,p)=>{payload=p;}}});
    C.validateDetection(data,gene,payload.counts);
    assert.equal(gene.sections.reduce((n,g)=>n+g.n,0),data.cell_count);
    const cv=gene.oof_summary;assert.ok(cv.count_pi_coverage>=0&&cv.count_pi_coverage<=1);assert.ok(cv.brier>=0&&cv.baseline_brier>=0);assert.ok(cv.count_pi_mean_fraction_width>=0&&cv.count_pi_mean_fraction_width<=1);assert.equal(cv.coverage_target,.95);
    assert.equal(cv.trained_anatomy_only.n_groups+cv.prior_only_groups,cv.n_groups);assert.equal(cv.trained_anatomy_only.n_cells+cv.prior_only_cells,cv.n_cells);
    assert.equal(cv.interval_target,'New physical-section positive count at its observed n');
  }
});
