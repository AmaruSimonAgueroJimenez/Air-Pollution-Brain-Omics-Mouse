// Exercise the viewer's actual handlers with deliberately unresolved Plotly writes.
const assert = require('node:assert/strict');
const {test} = require('node:test');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const copy = value => JSON.parse(JSON.stringify(value));
const camera = zoom => ({eye:{x:zoom,y:zoom,z:zoom*.75},up:{x:0,y:0,z:1}});
const tick = () => new Promise(resolve => setImmediate(resolve));

function fixture(){
  const plots = Object.fromEntries(['reference','simulation'].map(id => [id,{
    camera:camera(1.2), handlers:{}, on(name,handler){this.handlers[name]=handler;}
  }]));
  const calls=[],pending=[],errors=[];
  const context=vm.createContext({
    document:{documentElement:{lang:'en'},getElementById:id=>plots[id]},
    window:{ATLAS_DATA:{genes:[]}},console:{error:error=>errors.push(error)},
    Plotly:{relayout(target,update){
      const entry={target,camera:copy(update['scene.camera'])};calls.push(entry);
      return new Promise((resolve,reject)=>pending.push({...entry,resolve,reject}));
    }}
  });
  const source=fs.readFileSync(path.join(__dirname,'../docs/assets/viewer/viewer.js'),'utf8');
  const init=source.lastIndexOf('\ninit().catch(');
  assert.ok(init>0,'Viewer initialization entry point must be present');
  // Load every production definition, omitting only browser/data initialization.
  vm.runInContext(source.slice(0,init)+
    '\nglobalThis.testAPI={linkCameras,getLast:()=>lastCamera,isBusy:()=>cameraSyncRunning||!!queuedCamera};',context);
  context.testAPI.linkCameras('reference','simulation');
  context.testAPI.linkCameras('simulation','reference');
  function emit(id,value){
    plots[id].camera=copy(value);
    plots[id].handlers.plotly_relayout({'scene.camera':value});
  }
  async function complete({fail=false}={}){
    assert.ok(pending.length,'An asynchronous camera write should be pending');
    const entry=pending.shift();
    if(fail) entry.reject(new Error('deliberate Plotly failure'));
    else{
      // Worst case: commit after a new user gesture, with Plotly's full-camera echo.
      plots[entry.target].camera=copy(entry.camera);
      plots[entry.target].handlers.plotly_relayout({'scene.camera':{
        center:{x:0,y:0,z:0},projection:{type:'perspective'},...copy(entry.camera)
      }});
      entry.resolve();
    }
    await tick();
  }
  async function drain(){
    for(let n=0;pending.length&&n<20;n++) await complete();
    assert.equal(pending.length,0,'Programmatic camera echoes must not loop');
    assert.equal(context.testAPI.isBusy(),false,'Camera queue must reach idle');
  }
  function bothEqual(value){
    for(const plot of Object.values(plots))assert.deepEqual(plot.camera,value);
    assert.deepEqual(copy(context.testAPI.getLast()),value);
  }
  return {plots,calls,pending,errors,emit,complete,drain,bothEqual};
}

test('rapid same-source zoom events coalesce to the latest camera instead of being dropped',async()=>{
  const f=fixture();
  f.emit('reference',camera(1));
  f.emit('reference',camera(.8));
  f.emit('reference',camera(.6));
  assert.equal(f.calls.length,1,'Keep only one Plotly camera write in flight');
  await f.complete();
  assert.equal(f.calls.length,2);
  assert.deepEqual(f.calls[1].camera,camera(.6),'Skip intermediate zoom, preserve newest zoom');
  await f.drain();
  f.bothEqual(camera(.6));
  assert.equal(f.calls.length,2);
});

test('programmatic target echoes, expanded defaults and non-camera events do not bounce back',async()=>{
  const f=fixture();
  f.emit('reference',camera(.9));
  await f.drain();
  f.plots.simulation.handlers.plotly_relayout({'scene.camera':{
    projection:{type:'perspective'},center:{z:0,x:0,y:0},...camera(.9)
  }});
  f.plots.reference.handlers.plotly_relayout({autosize:true});
  await tick();
  assert.equal(f.calls.length,1);
  assert.equal(f.pending.length,0);
  f.bothEqual(camera(.9));
});

test('camera synchronization works in reverse after a completed update',async()=>{
  const f=fixture();
  f.emit('reference',camera(.8));await f.drain();
  f.emit('simulation',camera(1.6));await f.drain();
  assert.deepEqual(f.calls.map(c=>c.target),['simulation','reference']);
  f.bothEqual(camera(1.6));
});

test('reverse-direction gesture during an older target write restores both views to latest intent',async()=>{
  const f=fixture();
  f.emit('reference',camera(.8));
  // The older reference->simulation write has not committed yet.
  f.emit('simulation',camera(1.7));
  assert.equal(f.calls.length,1);
  await f.drain();
  assert.deepEqual(f.calls.map(c=>c.target),['simulation','simulation','reference']);
  f.bothEqual(camera(1.7));
});

test('later gesture supersedes an unfinished reverse-direction repair',async()=>{
  const f=fixture();
  f.emit('reference',camera(.8));f.emit('simulation',camera(1.7));
  await f.complete(); // Starts restoring simulation to 1.7.
  f.emit('reference',camera(.4));
  await f.drain();
  f.bothEqual(camera(.4));
  assert.equal(f.calls.length,3,'Obsolete pending reverse write is skipped');
});

test('a rejected write releases the queue and preserves a newer source gesture',async()=>{
  const f=fixture();
  f.emit('reference',camera(.8));f.emit('reference',camera(.3));
  await f.complete({fail:true});await f.drain();
  assert.equal(f.errors.length,1);
  f.bothEqual(camera(.3));
});
