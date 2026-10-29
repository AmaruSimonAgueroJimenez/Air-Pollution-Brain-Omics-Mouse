#!/usr/bin/env python3
"""Slide-sized 3D maps of measured atlas cells and frozen regional scenarios.

No model fitting, imputation, interpolation, mirroring or cell subsampling.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/defense_slides/figures'
MANIFEST=ROOT/'results/defense_slides/three_dimensional_figure_manifest.json'
INK,TEAL,RUST,GRAY='#142c35','#376e75','#b65531','#b4bec1'
STYLE={'font.family':'DejaVu Sans','font.size':17,'axes.titlesize':23,
       'axes.titleweight':'bold','axes.labelsize':17,'xtick.labelsize':14,
       'ytick.labelsize':14,'text.color':INK,'axes.labelcolor':INK,
       'axes.edgecolor':'#7a8f94','svg.fonttype':'none','pdf.fonttype':42,
       'figure.facecolor':'white','savefig.facecolor':'white','savefig.bbox':None}
SPEC={
 'HY': {'directory':'spatial_viewer','viewer':'viewer','cells':'registered_hypothalamus_cells.tsv.gz',
        'gene':'Tll1','stem':'spatial_3d_hy_tll1','color':TEAL,
        'name_en':'Hypothalamus','name_es':'Hipotálamo'},
 'ACA': {'directory':'spatial_cortex','viewer':'viewer-cortex','cells':'registered_cortex_cells.tsv.gz',
         'gene':'Fosb','stem':'spatial_3d_aca_fosb','color':RUST,
         'name_en':'Anterior cingulate cortex','name_es':'Corteza cingulada anterior'},
}


def sha(path):
 return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def data_for(region):
 spec=SPEC[region]
 data_path=f'docs/assets/{spec["viewer"]}/data.js'
 text=(ROOT/data_path).read_text()
 data=json.loads(text.split('=',1)[1].rstrip(';\n'))
 cell_path=f'results/{spec["directory"]}/{spec["cells"]}'
 cells=pd.read_csv(ROOT/cell_path,sep='\t',dtype={'cell_label':str})
 cache_path=f'results/{spec["directory"]}/reference_counts.npz'
 with np.load(ROOT/cache_path) as cache:
  if not np.array_equal(cache['cell_labels'],cells.cell_label.to_numpy(str)):
   raise ValueError('Cached counts and cells have different orders')
  gene=next(g for g in data['genes'] if g['symbol']==spec['gene'])
  j=int(np.flatnonzero(cache['gene_ids']==gene['id'])[0])
  counts=cache['counts'][:,j].astype(float)
 summary_path=f'results/{spec["directory"]}/summary.json'
 summary=json.loads((ROOT/summary_path).read_text())
 beta=float(gene['effect']['log2FoldChange_shrunk'])
 multiplier=float(2**beta)
 if not np.isclose(multiplier,summary['default_multiplier']):
  raise ValueError('Frozen default multiplier mismatch')
 if len(cells)!=data['cell_count'] or int((counts>0).sum())!=gene['reference_nonzero_cells']:
  raise ValueError('Frozen cell count/detection mismatch')
 if not np.isfinite(counts).all() or np.any(counts<0):
  raise ValueError('Invalid measured imaging counts')
 return dict(region=region,spec=spec,cells=cells,counts=counts,gene=gene,
             multiplier=multiplier,beta=beta,data=data,
             sources=[cell_path,cache_path,summary_path,data_path])


def frame(ax,cells,limits=None,view=(22,-58)):
 coords=cells[['z','x','y']].to_numpy(float)
 mins,maxs=coords.min(axis=0),coords.max(axis=0)
 if limits is None:
  pad=(maxs-mins)*.04
  limits=[(mins[0]-pad[0],maxs[0]+pad[0]),(mins[1]-pad[1],maxs[1]+pad[1]),
          (-maxs[2]-pad[2],-mins[2]+pad[2])]
 ax.set(xlim=limits[0],ylim=limits[1],zlim=limits[2],
        xlabel='ML (mm)',ylabel='AP (mm)',zlabel='−DV (mm)')
 ax.view_init(*view)
 ax.set_proj_type('ortho')
 ax.set_box_aspect([high-low for low,high in limits])
 ax.tick_params(pad=1,labelsize=13)
 ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
 ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
 ax.zaxis.set_major_locator(MaxNLocator(nbins=4))
 ax.xaxis.labelpad=12
 ax.yaxis.labelpad=12
 ax.zaxis.labelpad=16
 for axis in (ax.xaxis,ax.yaxis,ax.zaxis):
  axis.pane.set_facecolor((1,1,1,0))
  axis.pane.set_edgecolor((.84,.88,.88,.6))
  axis._axinfo['grid'].update({'color':(.82,.86,.86,.45),'linewidth':.5})
 return limits


def export(fig,stem,lang,caption,sources,details,manifest):
 name=stem+('_es' if lang=='es' else '')
 fig.canvas.draw()
 for ax in fig.axes:
  axes=[(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]
  if hasattr(ax,'zaxis'):axes.append((ax.zaxis,ax.get_zlim()))
  for axis,limits in axes:
   low,high=sorted(limits)
   for tick in axis.get_major_ticks():
    if tick.get_loc()<low-1e-9 or tick.get_loc()>high+1e-9:
     tick.label1.set_visible(False);tick.label2.set_visible(False)
 fig.canvas.draw()
 renderer=fig.canvas.get_renderer()
 from matplotlib.text import Text
 overflow=[]
 for artist in fig.findobj(Text):
  if artist.get_visible() and artist.get_text():
   box=artist.get_window_extent(renderer)
   w,h=fig.canvas.get_width_height()
   if box.x0 < -2 or box.y0 < -2 or box.x1>w+2 or box.y1>h+2:
    overflow.append(artist.get_text())
 if overflow:
  raise ValueError(f'{name}: text outside canvas: {overflow}')
 paths={}
 for ext in ('png','svg','pdf'):
  path=OUT/f'{name}.{ext}'
  metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else ({'Date':None} if ext=='svg' else {})
  fig.savefig(path,dpi=160,metadata=metadata)
  paths[ext]=str(path.relative_to(ROOT))
 manifest[name]={'stem':name,'language':lang,'caption':caption,'files':paths,
                 'sha256':{ext:sha(path) for ext,path in paths.items()},
                 'source_files':[{ 'path':p,'sha256':sha(p)} for p in sources],
                 'figure_size_inches':[float(v) for v in fig.get_size_inches()],
                 'png_dimensions':[int(round(v*160)) for v in fig.get_size_inches()],
                 'text_overflow':overflow,'no_model_refit':True,**details}
 plt.close(fig)


def scenario_plot(d,lang,manifest):
 es=lang=='es';spec=d['spec'];cells=d['cells'];counts=d['counts'];m=d['multiplier']
 positive=counts>0
 scenario=counts*m
 values=[np.log2(1+counts),np.log2(1+scenario)]
 vmax=float(np.max(np.concatenate(values)))
 norm=Normalize(0,vmax)
 fig=plt.figure(figsize=(16,6.2))
 axes=[fig.add_subplot(1,2,k+1,projection='3d') for k in range(2)]
 fig.subplots_adjust(left=.012,right=.90,bottom=.11,top=.92,wspace=.02)
 for ax,value in zip(axes,values):
  # Every cell is plotted, with explicit visual separation of measured zeros.
  # A shared small glyph does not encode abundance or a modeled cell size.
  ax.scatter(cells.z[~positive],cells.x[~positive],-cells.y[~positive],
             s=1.2,c=GRAY,alpha=.055,edgecolors='none',depthshade=False,rasterized=True)
  ax.scatter(cells.z[positive],cells.x[positive],-cells.y[positive],
             s=12,c=value[positive],cmap='viridis',norm=norm,
             alpha=.97,edgecolors='none',depthshade=False,rasterized=True)
  frame(ax,cells)
 axes[0].set_title('Referencia MERFISH medida' if es else 'Measured MERFISH reference',pad=9)
 axes[1].set_title(f'Escenario PM₂.₅ uniforme × {m:.3f}' if es else f'Uniform PM₂.₅ scenario × {m:.3f}',pad=9)
 cax=fig.add_axes([.918,.24,.016,.55])
 cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),cax=cax)
 cb.ax.tick_params(labelsize=14)
 cb.set_label('log₂(1 + conteo)\nEscala común sin recorte' if es else 'log₂(1 + count)\nShared scale, no clipping',labelpad=15,fontsize=16)
 n=len(cells);pos=int(positive.sum())
 legend=[Line2D([],[],marker='o',linestyle='',markerfacecolor=GRAY,markeredgecolor='none',markersize=6,
                label='Ceros medidos en gris tenue' if es else 'Measured zeros in faint gray')]
 fig.legend(handles=legend,loc='lower left',bbox_to_anchor=(.72,.008),frameon=False,fontsize=15)
 caption=(f'{spec["name_es"]}, gen {spec["gene"]}. Izquierda: conteos de imagen MERFISH realmente medidos del macho Zhuang-ABCA-2 sin exposición. Derecha: esos mismos conteos multiplicados uniformemente por 2^({d["beta"]:.6f})={m:.8f}, usando el efecto bulk regional regularizado. Se muestran las {n:,} células en ambos paneles, sin remuestreo, desplazamiento, interpolación ni espejo. Hay {pos:,} células con conteo medido positivo; el escalamiento conserva la detección y los ceros. Las cámaras, los límites de ejes y la escala de color log₂(1+conteo) son iguales; se usa el máximo conjunto completo y no se recortan valores. Las posiciones corresponden al hemisferio registrado y se grafican ML, AP, −DV en mm CCF. Estas coordenadas derivan del registro de los cortes, no de una adquisición 3D isotrópica directa. El panel derecho es una ilustración, no una medición de células expuestas ni una predicción bayesiana.' if es else
          f'{spec["name_en"]}, {spec["gene"]}. Left: actually measured MERFISH imaging counts from the unexposed male Zhuang-ABCA-2 donor. Right: those same counts uniformly multiplied by 2^({d["beta"]:.6f})={m:.8f}, using the regional shrunken bulk effect. All {n:,} cells appear in each panel, without resampling, jitter, interpolation or mirroring. There are {pos:,} measured positive cells; scaling preserves detection and measured zeros. Cameras, axis limits and the log₂(1+count) color scale are identical, with the full joint maximum and no value clipping. Positions belong to the registered hemisphere and are plotted ML, AP, −DV in CCF mm. These coordinates derive from registered sections rather than direct isotropic 3D acquisition. The right panel is an illustration, not a measurement of exposed cells or a Bayesian prediction.')
 export(fig,spec['stem'],lang,caption,d['sources'],
        {'region':d['region'],'gene':spec['gene'],'cells_per_panel':n,'measured_positive_cells':pos,
         'shrunken_bulk_log2FC':d['beta'],'scenario_multiplier':m,
         'coordinates':'ML=cell.z, AP=cell.x, vertical=-DV=-cell.y; CCF millimetres',
         'display_transformation':'log2(1+count); same full-joint-maximum scale; no clipping',
         'shared_color_min':0,'shared_color_max':vmax,'camera_elevation':22,'camera_azimuth':-58,
         'projection':'orthographic','cell_subsampling':False,'mirroring':False,'interpolation':False,
         'zero_encoding':'all measured zeros retained, faint gray glyphs','point_sizes':{'zero':1.2,'positive':12},
         'point_opacity':{'zero':.055,'positive':.97}},manifest)


def anatomy_plot(hy,aca,lang,manifest):
 es=lang=='es';fig=plt.figure(figsize=(16,7.2))
 ax=fig.add_axes([.02,.04,.67,.89],projection='3d')
 mesh=hy['data']['meshes']['brain']
 vertices=np.c_[mesh['x'],mesh['y'],mesh['z']]
 triangles=vertices[np.c_[mesh['i'],mesh['j'],mesh['k']]]
 collection=Poly3DCollection(triangles,facecolors='#b8c7cb',edgecolors='none',
                             alpha=.055,rasterized=True,zsort='average')
 ax.add_collection3d(collection)
 for d in (aca,hy):
  c=d['cells']
  ax.scatter(c.z,c.x,-c.y,s=1.2,c=d['spec']['color'],alpha=.48,
             edgecolors='none',depthshade=False,rasterized=True)
 limits=[(vertices[:,j].min(),vertices[:,j].max()) for j in range(3)]
 allcells=pd.concat([hy['cells'],aca['cells']],ignore_index=True)
 frame(ax,allcells,limits=limits,view=(20,-58))
 ax.set_title('Contexto anatómico CCF' if es else 'CCF anatomical context',pad=9)
 x=.70
 fig.text(x,.78,'Un donante sin exposición' if es else 'One unexposed donor',fontsize=22,fontweight='bold')
 fig.text(x,.69,'Zhuang-ABCA-2\nMacho adulto / MERFISH' if es else 'Zhuang-ABCA-2\nAdult male / MERFISH',fontsize=19,linespacing=1.5)
 fig.text(x,.50,'ACA  /  14,582 células' if es else 'ACA  /  14,582 cells',fontsize=22,fontweight='bold',color=RUST)
 fig.text(x,.40,'HY  /  37,811 células' if es else 'HY  /  37,811 cells',fontsize=22,fontweight='bold',color=TEAL)
 fig.text(x,.35 if es else .26,'52,393 posiciones celulares\nregistradas\nUn hemisferio registrado\nSin crear células nuevas' if es else
          '52,393 registered cell positions\nOne registered hemisphere\nNo new cells generated',fontsize=18,linespacing=1.55,va='top' if es else 'baseline')
 fig.text(.5,.024,'Contorno gris: anatomía de referencia CCF. Los puntos representan solo las regiones medidas seleccionadas.' if es else
          'Gray outline: CCF reference anatomy. Points represent only the selected measured regions.',ha='center',fontsize=16)
 caption=('Posiciones registradas de todas las 14,582 células ACA y 37,811 células HY del mismo donante Zhuang-ABCA-2. Las coordenadas CCF proceden del registro de los cortes, no de una adquisición 3D isotrópica directa. El gris es una malla anatómica de cerebro completo previamente exportada de la anotación CCF, con paso superficial de 0.15 mm. Esa superficie no representa cobertura de expresión. Los puntos coloreados son células realmente medidas del hemisferio registrado. No se rellenan huecos, no se reflejan células al hemisferio opuesto y no se añaden valores de expresión ni donantes. Los colores solo distinguen regiones, sin codificar una respuesta a PM₂.₅.' if es else
          'Registered positions of all 14,582 ACA and 37,811 HY cells from the same Zhuang-ABCA-2 donor. CCF coordinates come from registered sections rather than direct isotropic 3D acquisition. Gray is the previously exported whole-brain anatomical mesh from the CCF annotation, with 0.15-mm surface step. That surface does not represent expression coverage. Colored points are actually measured cells from the registered hemisphere. Gaps are not filled, cells are not reflected to the opposite hemisphere, and no expression values or donors are added. Colors identify regions only and do not encode a PM₂.₅ response.')
 export(fig,'spatial_3d_anatomy',lang,caption,list(dict.fromkeys(hy['sources']+aca['sources'])),
        {'total_measured_cells':52393,'ACA_cells':14582,'HY_cells':37811,
         'one_unexposed_donor':'Zhuang-ABCA-2','surface_step_mm':.15,
         'surface':'existing whole-brain CCF annotation mesh; anatomical context only',
         'cell_subsampling':False,'mirroring':False,'interpolation':False,'coordinate_frame':'ML,AP,-DV in CCF millimetres'},manifest)


def main():
 OUT.mkdir(parents=True,exist_ok=True)
 hy,aca=data_for('HY'),data_for('ACA')
 manifest={}
 with plt.rc_context(STYLE):
  for lang in ('en','es'):
   scenario_plot(hy,lang,manifest);scenario_plot(aca,lang,manifest);anatomy_plot(hy,aca,lang,manifest)
 MANIFEST.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
 print(f'Exported {len(manifest)} bilingual 3D figures, all measured cells retained.')

if __name__=='__main__':main()
