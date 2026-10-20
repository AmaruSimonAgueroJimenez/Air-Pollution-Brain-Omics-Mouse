#!/usr/bin/env python3
"""Exact conjugate Bayesian baseline spatial regression, never a PM exposure model.

All fits use the original unexposed male ABCA-2 reference cells. Whole physical
section blocks are held out, with adjacent sections guarded. Fixed priors are
not selected using held-out outcomes. Gaussian transformed-count limitations,
calibration failures and the single-donor scope are explicitly retained.
"""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import argparse,hashlib,json,time,platform,itertools
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import scipy
from scipy.linalg import cho_factor,cho_solve
from scipy.spatial import cKDTree
from scipy.stats import t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/spatial_bayes';META=ROOT/'metadata/spatial_bayes';WEB=ROOT/'docs/assets/spatial-bayes'
LENGTHS=np.array([.60,.45,.45]) # AP,DV,ML millimetres, fixed before any fit.
CONFIGS={'primary':(1.,1.),'short_length':(.5,1.),'long_length':(2.,1.),'strong_shrinkage':(1.,.25),'weak_shrinkage':(1.,4.)}
NOMINAL=np.array([.20,.50,.80,.90,.95,.99])
SEED=20261020
plt.rcParams.update({'svg.fonttype':'none','svg.hashsalt':'spatial-bayes-v1','pdf.fonttype':42,'font.family':'DejaVu Sans','savefig.facecolor':'white'})

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clean(x):
 if isinstance(x,dict):return {str(k):clean(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v)for v in x]
 if isinstance(x,np.ndarray):return clean(x.tolist())
 if isinstance(x,np.generic):return clean(x.item())
 if isinstance(x,float)and not np.isfinite(x):return None
 return x
def dump(p,x,compact=False):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(clean(x),ensure_ascii=False,indent=None if compact else 2,separators=(',',':')if compact else None,allow_nan=False)+'\n')
def id_hash(ids):return hashlib.sha256(('\n'.join(map(str,ids))+'\n').encode()).hexdigest()

def methods():
 return {'schema_version':1,'model':'Proper conjugate Normal–Inverse-Gamma spatial regression with physical-section nuisance random intercepts',
 'response':'y=ln(1+raw MERFISH imaging count); no normalization or raw-count back-transformation',
 'likelihood':'y = X beta + section_offset + epsilon; epsilon~N(0,sigma²), section_offset~N(0,0.25 sigma²); working Gaussian model is approximate for zero-heavy nonnegative data.',
 'features':['global intercept','training-observed anatomical substructure intercepts','training-centered/scaled AP,DV,ML linear coordinates','36 Gaussian RBFs with anisotropic physical lengths','training-section nuisance intercepts'],
 'celltype_predictors':False,'why_no_celltype_predictors':'Atlas expression-derived class/subclass labels could leak the target transcript; they are not used to construct features.',
 'basis':{'shape':[4,3,3],'AP_quantiles':[.05,.35,.65,.95],'DV_ML_quantiles':[.1,.5,.9],'centers':'Cartesian product of training-coordinate quantiles only','length_mm_AP_DV_ML':LENGTHS,'training_only_centering':'Coordinate means/SD and mean RBF feature values are learned only from the training cells. Outcome is not centered or scaled.','random_seed':None,'deterministic':True},
 'prior':{'beta_mean':0.,'beta_covariance':'sigma² * diagonal(V0)','intercept_variance_multiplier':25.,'anatomy_variance_multiplier':1.,'linear_coordinate_variance_multiplier':.25,'each_RBF_variance_multiplier':1/36,'section_variance_multiplier':.25,'sigma2':'InverseGamma(shape=3,scale=1); density proportional to (sigma²)^(-4) exp(-1/sigma²)','proper':True},
 'posterior':{'Vn':'inverse(X.T X + inverse(V0))','mn':'Vn X.T y','an':'3+n/2','bn':'1+0.5*(y.T y-mn.T inverse(Vn) mn)','solver':'Cholesky; normal-equation and residual-plus-prior quadratic checks recorded','marginal':'Student t with df=2an','latent_scale2':'bn/an * x.T Vn x, plus independent prior variance for unseen anatomy','predictive_scale2':'bn/an * (1+0.25+x.T Vn x), plus independent prior variance for unseen anatomy'},
 'target':'Section-adjusted spatial conditional mean; physical-section nuisance offset averaged to zero. New-cell predictive distribution integrates a NEW section offset, at the actual measured coordinate. No new cell is instantiated.',
 'unseen_anatomy':'Training vocabulary only. Unknown test anatomy has zero prior mean and an added independent anatomy-prior variance; it is never silently treated as a known group.',
 'validation':{'scheme':'Five contiguous blocks of whole physical sections, ordered by acquisition AP/experimental_z','fold_assignment':'Deterministic array_split of ordered region sections; all cells in each test section held out together. No cell/section exclusion by expression.','guard':'Nearest represented region section immediately before and after each test block removed from training; guard may exceed one acquisition cut if region is absent from intervening cuts. Exact labels and AP distances recorded.','hyperparameters':'All fixed before CV; no held-out tuning or model selection. All four sensitivity variants are reported regardless of performance.','gene_selection':'Ten measured candidates per region frozen by prior bulk-results selection; no atlas-outcome selection.','donors':1,'no_independent_animal_validation':True},
 'baselines':{'anatomy_constant':'Proper Bayesian global/anatomical intercept model, same prior and section nuisance structure; training cells only.','knn':'Unweighted average of30 nearest training cells in same anatomy when≥30available, otherwise all training anatomies; distances scaled by fixed primary AP/DV/ML lengths. No uncertainty interval assigned to deterministic interpolation.'},
 'sensitivity':{k:{'length_factor':v[0],'spatial_prior_variance_factor':v[1]}for k,v in CONFIGS.items()},
 'conditional_tail_metrics':'coverage_positive and coverage_zero condition on realized outcomes and are descriptive tail-recovery diagnostics. They need not equal the nominal marginal 95%level, even under a calibrated marginal predictive interval. The positive-tail flag is not a formal calibration test.',
 'interval_scope':'Conditional on fixed kernel lengths, prior variances, anatomy labels and Gaussian likelihood; no residual spatial-correlation, registration-error or between-animal uncertainty model. Cell-wise calibration is descriptive; sections/donors are not independent cell replicates.',
 'interpretation_limits':['Gaussian fitted means, latent intervals and predictive distributions can be negative on the physically nonnegative ln(1 + count) scale. Predictive negative mass is quantified; fitted quantities and intervals are never clipped or presented as raw-count predictions.','Latent-mean credible intervals are not new-cell predictive intervals; empirical coverage is assessed only for the latter.','Full-reference posterior arrays are in-sample fitted spatial summaries; held-out metrics come exclusively from refits excluding test/guard sections.','No PM response, differential expression, protein abundance, additional resolution, missing-gene expression, new positions, cell generation, mirroring or relocation is inferred.'],
 'computational_limits':{'basis':'36 finite-rank RBFs, no Gaussian-process or neural-network fit','all_reference_cells_retained':True,'fold_fits_total':60,'full_reference_fits_total':2,'posterior_integration':'Exact conjugate algebra; no MCMC, approximate optimization or convergence diagnostic required.','hyperparameter_uncertainty':'Length scales, anatomy/section prior variances and basis size are fixed, not integrated out.','residual_dependence':'After section nuisance effects, remaining spatial dependence and heteroscedasticity are not modeled.','resolution':'No additional spatial resolution or unmeasured-gene/cell information is created.'},
 'sources':[{'title':'Murphy: Conjugate Bayesian analysis of the Gaussian distribution','url':'https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf'},{'title':'SciPy Student t implementation and parameterization','url':'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html'},{'title':'Original MERFISH atlas','url':'https://doi.org/10.1038/s41586-023-06808-9'}]}

class Design:
 def __init__(self,coords,anatomy,sections,length_factor=1.,shrink=1.,spatial=True):
  self.anatomy=sorted(set(anatomy));self.sections=sorted(set(sections));self.spatial=spatial;self.lengths=LENGTHS*length_factor
  self.mean=coords.mean(0);self.scale=coords.std(0);self.scale[self.scale<1e-8]=1.
  qs=[[.05,.35,.65,.95],[.1,.5,.9],[.1,.5,.9]]
  self.centers=np.array(list(itertools.product(*[np.quantile(coords[:,k],q)for k,q in enumerate(qs)])))
  self.rbf_mean=self.rbf(coords).mean(0)
  v=[25.]+[1.]*len(self.anatomy)
  if spatial:v += [.25*shrink]*3+[shrink/36]*36
  v += [.25]*len(self.sections)
  self.prior=np.array(v)
 def rbf(self,c):return np.exp(-.5*np.square((c[:,None,:]-self.centers[None,:,:])/self.lengths).sum(2))
 def transform(self,c,a,s,include_section=False):
  x=[np.ones((len(c),1)),np.array([a==k for k in self.anatomy],dtype=float).T]
  if self.spatial:x += [(c-self.mean)/self.scale,self.rbf(c)-self.rbf_mean]
  x += [np.array([s==k for k in self.sections],dtype=float).T if include_section else np.zeros((len(c),len(self.sections)))]
  unknown=~np.isin(a,self.anatomy)
  return np.column_stack(x),unknown.astype(float)
 def audit(self):return {'anatomy_levels':self.anatomy,'section_levels':self.sections,'coordinate_mean_training':self.mean,'coordinate_sd_training':self.scale,'RBF_centers_training':self.centers,'RBF_mean_training':self.rbf_mean,'lengths_mm':self.lengths,'prior_variance_multipliers':self.prior,'feature_count':len(self.prior),'uses_outcomes_for_preprocessing':False}

def fit(x,y,prior):
 precision=x.T@x+np.diag(1/prior);rhs=x.T@y;ch=cho_factor(precision,lower=True,check_finite=True)
 m=cho_solve(ch,rhs);v=cho_solve(ch,np.eye(len(prior)))
 a=3+len(x)/2;b=1+.5*(np.square(y).sum(0)-(m*rhs).sum(0))
 rss=np.square(y-x@m).sum(0);quadratic=(np.square(m)/prior[:,None]).sum(0)
 b_direct=1+.5*(rss+quadratic)
 audit={'normal_equation_relative_error':float(np.linalg.norm(precision@m-rhs)/(1+np.linalg.norm(rhs))),'scale_quadratic_relative_error':float(np.max(np.abs(b-b_direct)/(1+np.abs(b_direct)))),'covariance_inverse_relative_error':float(np.linalg.norm(precision@v-np.eye(len(v)))/len(v)),'posterior_scale_min':float(b.min()),'df':2*a,'feature_count':len(prior)}
 if min(b)<=0:raise ValueError('Nonpositive posterior inverse-gamma scale')
 return {'mean':m,'V':v,'a':a,'b':b,'audit':audit}

def predict(model,x,unknown):
 h=np.einsum('ij,jk,ik->i',x,model['V'],x,optimize=True)+unknown
 if h.min()<-1e-9:raise ValueError('Negative posterior latent variance')
 h=np.maximum(h,0);mu=x@model['mean'];scale2=model['b']/model['a'];df=2*model['a']
 latent=np.sqrt(h[:,None]*scale2);pred=np.sqrt((h[:,None]+1.25)*scale2);q=student_t.ppf(.975,df)
 return {'posterior_mean':mu,'ci_lower':mu-q*latent,'ci_upper':mu+q*latent,'predictive_lower':mu-q*pred,'predictive_upper':mu+q*pred,'predictive_scale':pred,'df':df,'negative_mass':student_t.cdf(-mu/pred,df)}

def interval_score(y,lo,hi,alpha=.05):return hi-lo+2/alpha*(lo-y)*(y<lo)+2/alpha*(y-hi)*(y>hi)
def score(y,mu,lo=None,hi=None,negative=None):
 d={'n_cells':len(y),'rmse':float(np.sqrt(np.mean(np.square(y-mu)))),'mae':float(np.mean(np.abs(y-mu)))}
 if lo is not None:d.update({'predictive_95_coverage':float(np.mean((y>=lo)&(y<=hi))),'predictive_95_mean_width':float(np.mean(hi-lo)),'interval_score_95':float(np.mean(interval_score(y,lo,hi))),'coverage_positive':float(np.mean((y[y>0]>=lo[y>0])&(y[y>0]<=hi[y>0])))if(y>0).any()else None,'coverage_zero':float(np.mean((lo[y==0]<=0)&(hi[y==0]>=0)))if(y==0).any()else None,'n_positive':int((y>0).sum())})
 if negative is not None:d['negative_predictive_probability_mean']=float(np.mean(negative))
 return d

def knn_predict(train_c,test_c,train_a,test_a,y):
 out=np.empty((len(test_c),y.shape[1]));fallback=0
 for group in np.unique(test_a):
  test=np.flatnonzero(test_a==group);training=np.flatnonzero(train_a==group)
  if len(training)<30:training=np.arange(len(train_c));fallback+=len(test)
  tree=cKDTree(train_c[training]/LENGTHS);_,near=tree.query(test_c[test]/LENGTHS,k=min(30,len(training)),workers=1)
  if near.ndim==1:near=near[:,None]
  out[test]=y[training[near]].mean(1)
 return out,fallback

def load_region(region):
 folder='spatial_viewer'if region=='HY'else'spatial_cortex';name='registered_hypothalamus_cells.tsv.gz'if region=='HY'else'registered_cortex_cells.tsv.gz'
 cells=pd.read_csv(ROOT/'results'/folder/name,sep='\t');cache=np.load(ROOT/'results'/folder/'reference_counts.npz')
 if not np.array_equal(cells.cell_label.astype(str).to_numpy(),cache['cell_labels']):raise ValueError('Cell ordering differs from measured count cache')
 section_path=ROOT/'results/spatial_sections'/region/'cells_with_sections.tsv.gz'
 if section_path.exists():
  sec=pd.read_csv(section_path,sep='\t')
  if not np.array_equal(sec.cell_label.astype(str).to_numpy(),cache['cell_labels']):raise ValueError('Section-joined cell order differs')
  for col in sec.columns:
   if col not in cells:cells[col]=sec[col]
 else:
  raw=next((ROOT/'data/reference_atlases/Zhuang_2023/metadata/Zhuang-ABCA-2').glob('*/cell_metadata.csv'))
  sec=pd.read_csv(raw,usecols=['cell_label','brain_section_label','z'],dtype={'cell_label':str}).rename(columns={'z':'experimental_z'})
  cells=cells.merge(sec,on='cell_label',how='left',validate='one_to_one',sort=False)
 if 'brain_section_label'not in cells and'section_label'in cells:cells['brain_section_label']=cells.section_label
 if 'experimental_z'not in cells:
  raw=next((ROOT/'data/reference_atlases/Zhuang_2023/metadata/Zhuang-ABCA-2').glob('*/cell_metadata.csv'))
  z=pd.read_csv(raw,usecols=['cell_label','z'],dtype={'cell_label':str}).rename(columns={'z':'experimental_z'})
  cells=cells.merge(z,on='cell_label',how='left',validate='one_to_one',sort=False)
 if cells[['brain_section_label','experimental_z']].isna().any().any():raise ValueError('Missing physical section membership/ordering')
 if not np.array_equal(cells.cell_label.astype(str).to_numpy(),cache['cell_labels']):raise ValueError('Section merge reordered cells')
 p=pd.read_csv(ROOT/'results/gene_priorities/prioritized_genes.tsv',sep='\t');group='hypothalamus_measured_top10'if region=='HY'else'cortex_measured_top10'
 genes=p[p.selection_group.eq(group)].sort_values('selection_rank').copy();panel=pd.read_csv(ROOT/'results'/folder/'gene_panel_coverage.tsv',sep='\t')
 if region=='HY':genes['gene_id_ensembl']=[panel.set_index('symbol').at[g,'id']for g in genes.gene]
 positions=pd.Index(cache['gene_ids']).get_indexer(genes.gene_id_ensembl)
 if(positions<0).any()or len(genes)!=10:raise ValueError('Frozen measured candidate mismatch')
 counts=np.array(cache['counts'][:,positions],dtype=float)
 if not np.isfinite(counts).all()or(counts<0).any():raise ValueError('Invalid measured counts')
 return cells,genes,np.log1p(counts),[ROOT/'results'/folder/name,ROOT/'results'/folder/'reference_counts.npz',ROOT/'results/gene_priorities/prioritized_genes.tsv',section_path if section_path.exists()else raw]

def run_region(region):
 begin=time.monotonic();dest=OUT/region;dest.mkdir(parents=True,exist_ok=True);web=WEB/region;web.mkdir(parents=True,exist_ok=True)
 cells,genes,y,inputs=load_region(region);coords=cells[['x','y','z']].to_numpy();anatomy=cells.substructure.fillna('UNASSIGNED').astype(str).to_numpy();sections=cells.brain_section_label.astype(str).to_numpy()
 ordered=cells.groupby('brain_section_label',sort=False).experimental_z.median().sort_values(kind='stable').index.astype(str).to_numpy()
 if len(ordered)<10:raise ValueError('Insufficient physical sections for5-fold guarded validation')
 section_fold={sec:i+1 for i,b in enumerate(np.array_split(ordered,5))for sec in b};fold_ids=np.array([section_fold[s]for s in sections])
 sectiontable=pd.DataFrame({'section':ordered,'order':np.arange(len(ordered)),'fold':[section_fold[s]for s in ordered]});sectiontable['n_cells']=[int((sections==s).sum())for s in ordered];sectiontable['experimental_z']=[float(cells.loc[sections==s,'experimental_z'].median())for s in ordered]
 sectiontable.to_csv(dest/'section_folds.tsv',sep='\t',index=False)
 metrics=[];section_metrics=[];calibration=[];audits=[];fold_records=[];cv_arrays={};sensitivity=[];knn_all=np.full_like(y,np.nan);primary_test_count=np.zeros(len(y),dtype=int)
 for conf,(lf,shrink)in CONFIGS.items():
  cv={k:np.full_like(y,np.nan)for k in ['posterior_mean','predictive_lower','predictive_upper','predictive_scale','negative_mass']};cvdf=np.full(len(y),np.nan)
  models=['bayes_spatial','anatomy_constant']if conf=='primary'else['bayes_spatial']
  for modelname in models:
   for k in cv:cv[k].fill(np.nan)
   for fold in range(1,6):
    test_sections=ordered[[section_fold[s]==fold for s in ordered]];locations=np.flatnonzero(np.isin(ordered,test_sections));guard=[]
    if locations[0]>0:guard.append(ordered[locations[0]-1])
    if locations[-1]+1<len(ordered):guard.append(ordered[locations[-1]+1])
    test=np.isin(sections,test_sections);train=~(test|np.isin(sections,guard))
    if set(sections[train])&set(sections[test]):raise ValueError('Section leakage')
    design=Design(coords[train],anatomy[train],sections[train],lf,shrink,modelname=='bayes_spatial')
    x,_=design.transform(coords[train],anatomy[train],sections[train],True);fitted=fit(x,y[train],design.prior)
    xx,unknown=design.transform(coords[test],anatomy[test],sections[test]);pred=predict(fitted,xx,unknown)
    for key in cv:cv[key][test]=pred[key]
    cvdf[test]=pred['df']
    audit={'region':region,'config':conf,'model':modelname,'fold':fold,'n_train':int(train.sum()),'n_test':int(test.sum()),'test_sections':test_sections,'guard_sections':guard,'unknown_anatomy_test_cells':int(unknown.sum()),'unknown_anatomy_labels':sorted(set(anatomy[test])-set(design.anatomy)),'posterior':fitted['audit'],'design':design.audit(),'train_cell_ids_sha256':id_hash(cells.cell_label[train]),'test_cell_ids_sha256':id_hash(cells.cell_label[test])};audits.append(audit)
    if conf=='primary'and modelname=='bayes_spatial':
     primary_test_count[test]+=1
     knn,fallback=knn_predict(coords[train],coords[test],anatomy[train],anatomy[test],y[train]);knn_all[test]=knn
     fold_records.append({**{k:audit[k]for k in ['fold','n_train','n_test','test_sections','guard_sections','unknown_anatomy_test_cells']},'knn_global_fallback_cells':fallback,'minimum_acquisition_AP_gap_mm':float(np.abs(np.unique(cells.experimental_z[train])[:,None]-np.unique(cells.experimental_z[test])[None,:]).min())})
    print(f'{region} {conf} {modelname} fold{fold}:train={train.sum()},test={test.sum()}',flush=True)
   for j,r in enumerate(genes.itertuples()):
    m=score(y[:,j],cv['posterior_mean'][:,j],cv['predictive_lower'][:,j],cv['predictive_upper'][:,j],cv['negative_mass'][:,j]);m.update({'region':region,'config':conf,'model':modelname,'gene':r.gene,'gene_id':r.gene_id_ensembl,'n_sections':len(ordered)});metrics.append(m)
    for nominal in NOMINAL:
     q=student_t.ppf((1+nominal)/2,cvdf);inside=np.abs(y[:,j]-cv['posterior_mean'][:,j])<=q*cv['predictive_scale'][:,j]
     calibration.append({'region':region,'config':conf,'model':modelname,'gene':r.gene,'nominal':nominal,'observed_coverage':float(inside.mean())})
    if conf=='primary':
     for sec in ordered:
      mask=sections==sec;sm=score(y[mask,j],cv['posterior_mean'][mask,j],cv['predictive_lower'][mask,j],cv['predictive_upper'][mask,j],cv['negative_mass'][mask,j]);section_metrics.append({'region':region,'model':modelname,'gene':r.gene,'section':sec,'fold':section_fold[sec],**sm})
   if conf=='primary':cv_arrays[modelname]={k:v.copy()for k,v in cv.items()}
  if conf!='primary':sensitivity.append({'config':conf,'posterior_means':cv['posterior_mean'].copy()})
 for j,r in enumerate(genes.itertuples()):metrics.append({'region':region,'config':'primary','model':'knn30','gene':r.gene,'gene_id':r.gene_id_ensembl,'n_sections':len(ordered),**score(y[:,j],knn_all[:,j])})
 # All-reference fit is a separate descriptive posterior, never used for CV scores.
 design=Design(coords,anatomy,sections);x,_=design.transform(coords,anatomy,sections,True);fitted=fit(x,y,design.prior);xp,unknown=design.transform(coords,anatomy,sections);full=predict(fitted,xp,unknown)
 audits.append({'region':region,'config':'primary','model':'bayes_spatial','scope':'all-reference fit','posterior':fitted['audit'],'design':design.audit()})
 np.savez_compressed(dest/'posterior_and_heldout.npz',cell_ids=cells.cell_label.to_numpy(dtype=str),gene_ids=genes.gene_id_ensembl.to_numpy(dtype=str),gene_symbols=genes.gene.to_numpy(dtype=str),observed_log1p=y,fold=fold_ids,**{'full_'+k:v for k,v in full.items()if isinstance(v,np.ndarray)},**{'heldout_'+k:v for k,v in cv_arrays['bayes_spatial'].items()},knn_heldout=knn_all)
 dump(dest/'fit_audit.json',audits);dump(dest/'fold_audit.json',fold_records)
 metricframe=pd.DataFrame(metrics);metricframe.to_csv(dest/'heldout_metrics.tsv',sep='\t',index=False);pd.DataFrame(section_metrics).to_csv(dest/'heldout_section_metrics.tsv',sep='\t',index=False);pd.DataFrame(calibration).to_csv(dest/'calibration.tsv',sep='\t',index=False)
 # Independent data-level transport contains no new coordinates or synthetic cells.
 entries=[];checks=[]
 for j,r in enumerate(genes.itertuples()):
  primary=next(m for m in metrics if m['gene']==r.gene and m['config']=='primary'and m['model']=='bayes_spatial');baseline=next(m for m in metrics if m['gene']==r.gene and m['config']=='primary'and m['model']=='anatomy_constant');kn=next(m for m in metrics if m['gene']==r.gene and m['model']=='knn30')
  flags=[]
  if primary['predictive_95_coverage']<.90:flags.append('undercoverage_overall_95_interval')
  if primary['predictive_95_coverage']>.985:flags.append('overcoverage_overall_95_interval')
  if primary['coverage_positive']is not None and primary['coverage_positive']<.90:flags.append('undercoverage_positive_cells')
  if primary['rmse']>=baseline['rmse']:flags.append('no_RMSE_improvement_over_anatomy_baseline')
  if primary['rmse']>=kn['rmse']:flags.append('no_RMSE_improvement_over_knn')
  if primary['negative_predictive_probability_mean']>.05:flags.append('nonnegative_support_mismatch')
  heldout={k:primary[k]for k in ['n_cells','n_sections','rmse','mae','predictive_95_coverage','predictive_95_mean_width','negative_predictive_probability_mean','coverage_positive','coverage_zero','n_positive','interval_score_95']};heldout.update({'anatomy_baseline_rmse':baseline['rmse'],'knn_rmse':kn['rmse'],'flags':flags,'section_macro_rmse':float(np.mean([m['rmse']for m in section_metrics if m['gene']==r.gene and m['model']=='bayes_spatial']))})
  key=r.gene_id_ensembl+'.json';payload={'schema_version':1,'region':region,'gene':r.gene,'id':r.gene_id_ensembl,'n':len(cells),'array_precision_decimal_places':5,'cell_ids_sha256':id_hash(cells.cell_label),'units':'ln(1+raw imaging count)','fit_scope':'all_reference_section_adjusted_latent_mean','predictive_target':'new cell in a new physical section at the actual measured coordinate','heldout':heldout,**{k:np.round(full[k][:,j],5)for k in ['posterior_mean','ci_lower','ci_upper','predictive_lower','predictive_upper']}}
  dump(web/key,payload,compact=True);back=json.loads((web/key).read_text());checks.append(all(np.allclose(back[k],full[k][:,j],rtol=0,atol=5.01e-6)for k in ['posterior_mean','ci_lower','ci_upper','predictive_lower','predictive_upper']))
  entries.append({'id':r.gene_id_ensembl,'symbol':r.gene,'file':key,'heldout':heldout,'by_section':[{k:m[k]for k in ['section','fold','n_cells','rmse','mae','predictive_95_coverage','predictive_95_mean_width']}for m in section_metrics if m['gene']==r.gene and m['model']=='bayes_spatial']})
 train_metrics=[{'gene':r.gene,**score(y[:,j],full['posterior_mean'][:,j])}for j,r in enumerate(genes.itertuples())]
 metadata={'schema_version':1,'region':region,'dataset':'Zhuang-ABCA-2','donor_n':1,'donor_sex':'male','cell_count':len(cells),'cell_ids_sha256':id_hash(cells.cell_label),'units':'ln(1+raw imaging count)','method':'Proper conjugate Bayesian anatomy + anisotropic spatial-basis regression; section nuisance effects','method_en':'Bayesian unexposed reference baseline','method_es':'Referencia basal bayesiana sin exposición','fit_scope':'all reference cells; section-adjusted fitted spatial summaries, not held-out predictions','credible_target':'95% posterior credible interval for section-adjusted latent mean','predictive_target':'95% new-cell interval including new-section offset and observation noise','cross_validation':{'scheme':'5 contiguous whole-section blocks with 1 neighbor section guarded per side','n_folds':5,'guard_sections_per_side':1,'unit':'whole physical brain section','donors':1,'section_count':len(ordered)},'limitations':methods()['interpretation_limits'],'genes':entries,'training_descriptive_metrics':train_metrics}
 dump(web/'metadata.json',metadata)
 metaout={'region':region,'n_cells':len(cells),'n_sections':len(ordered),'anatomical_groups':len(np.unique(anatomy)),'genes':[e['symbol']for e in entries],'heldout_results':[{**e['heldout'],'gene':e['symbol']}for e in entries],'elapsed_seconds':time.monotonic()-begin,'cell_ids_sha256':id_hash(cells.cell_label)}
 dump(dest/'summary.json',metaout)
 checks={'n10_measured_genes':len(genes)==10,'all_original_cells_retained':len(cells)==len(y),'cell_order_exact':True,'unique_real_cell_ids':not cells.cell_label.duplicated().any(),'each_cell_tested_once':bool(np.all(primary_test_count==1)),'heldout_predictions_complete':bool(all(np.isfinite(v).all()for v in cv_arrays['bayes_spatial'].values())),'all_posterior_arrays_finite':bool(all(np.isfinite(v).all()for v in full.values()if isinstance(v,np.ndarray))),'posterior_vs_predictive_width_distinct':bool(np.all((full['predictive_upper']-full['predictive_lower'])>(full['ci_upper']-full['ci_lower']))),'no_test_guard_sections_in_training':all(not(set(a['test_sections'])&set(a['design']['section_levels']))and not(set(a['guard_sections'])&set(a['design']['section_levels']))for a in audits if'fold'in a),'normal_equations_verified':all(a['posterior']['normal_equation_relative_error']<1e-9 and a['posterior']['scale_quadratic_relative_error']<1e-9 and a['posterior']['covariance_inverse_relative_error']<1e-9 for a in audits),'transport_roundtrip_5decimal_precision':all(checks),'no_expression_celltype_predictors':True,'no_PM_or_bulk_effect_fitted':True,'no_generated_or_relocated_cells':True}
 return {'cells':cells,'genes':genes,'y':y,'full':full,'heldout':cv_arrays,'knn':knn_all,'metrics':pd.DataFrame(metrics),'calibration':pd.DataFrame(calibration),'section_metrics':pd.DataFrame(section_metrics),'summary':metaout,'checks':checks,'inputs':inputs,'web_metadata':metadata}

def export(fig,stem,caption_en,caption_es,manifest,language):
 for ext in ['png','svg','pdf']:
  md={'CreationDate':None,'ModDate':None}if ext=='pdf'else({'Date':None}if ext=='svg'else{})
  fig.savefig(OUT/'figures'/f'{stem}.{ext}',dpi=300,metadata=md)
 manifest[stem]={'caption_en':caption_en,'caption_es':caption_es,'caption':caption_es if language=='es'else caption_en,'language':language,'files':{ext:f'results/spatial_bayes/figures/{stem}.{ext}'for ext in ['png','svg','pdf']},'sha256':{ext:sha(OUT/'figures'/f'{stem}.{ext}')for ext in ['png','svg','pdf']}}
 plt.close(fig)

def figures(regions):
 manifest={};colors={'bayes_spatial':'#0B7475','anatomy_constant':'#B85121','knn30':'#5C507D'}
 for es in [False,True]:
  language='es'if es else'en';suffix='_es'if es else''
  fig,axes=plt.subplots(2,3,figsize=(16,10.8),layout='constrained',gridspec_kw={'width_ratios':[1.4,1,1]})
  for i,(region,obj)in enumerate(regions.items()):
   genes=obj['genes'].gene.tolist();n=len(genes);ys=np.arange(n);m=obj['metrics'];m=m[m.config.eq('primary')]
   for model,offset,marker in [('bayes_spatial',-.16,'o'),('anatomy_constant',0,'s'),('knn30',.16,'^')]:
    mm=m[m.model.eq(model)].set_index('gene').loc[genes]
    axes[i,0].scatter(mm.rmse,ys+offset,color=colors[model],marker=marker,s=30,label=model)
    if model!='knn30':
     axes[i,1].scatter(mm.predictive_95_coverage,ys+offset,color=colors[model],marker=marker,s=28)
     axes[i,2].scatter(mm.predictive_95_mean_width,ys+offset,color=colors[model],marker=marker,s=28)
   for j,ax in enumerate(axes[i]):
    ax.set_yticks(ys,genes if j==0 else[]);ax.invert_yaxis();ax.grid(axis='x',alpha=.2);ax.spines[['top','right']].set_visible(False)
    ax.set_title(f'{region} · '+(['Error fuera de sección','Cobertura predictiva 95%','Anchura predictiva 95%'][j]if es else['Out-of-section error','95% predictive coverage','95% predictive width'][j]),loc='left',fontsize=12,weight='bold')
   axes[i,0].set_xlabel('RMSE · ln(1 + conteo)'if es else'RMSE · ln(1 + count)');axes[i,0].set_xlim(left=0)
   axes[i,1].set_xlim(0,1.02);axes[i,1].set_xlabel('Fracción observada'if es else'Observed fraction')
   axes[i,2].set_xlabel('ln(1 + conteo)'if es else'ln(1 + count)');axes[i,2].set_xlim(left=0)
  names=['Bayes espacial','Media anatómica Bayes','Interpolación30 vecinos']if es else['Spatial Bayes','Anatomy-mean Bayes','30-neighbor interpolation']
  for h,label in zip(axes[0,0].collections,names):h.set_label(label)
  axes[0,0].legend(loc='lower right',fontsize=8,frameon=True)
  fig.suptitle('Validación por bloques de cortes físicos · un donante sin exposición'if es else'Whole-section block validation · one unexposed donor',fontsize=16,weight='bold')
  capen='Five contiguous physical-section blocks are held out with one neighboring represented section guarded on each side. All ten fixed measured candidates per region are scored on every original cell exactly once; no held-out outcome tunes a model. Left: RMSE for spatial Bayes, anatomy-only Bayesian mean and 30-neighbor interpolation; center: empirical coverage of nominal 95% NEW-CELL predictive intervals; right: mean width of those intervals. Units are natural ln(1+raw imaging count). kNN has no invented intervals. Differences are predictive comparisons within one donor, not independent-animal evidence; poor calibration and wide intervals are retained.'
  capes='Se omiten cinco bloques contiguos de cortes físicos y se protege un corte representado vecino por lado. Los diez candidatos medidos fijos por región se evalúan una vez en cada célula original; no se ajusta el modelo con resultados omitidos. Izquierda: RMSE de Bayes espacial, media anatómica bayesiana e interpolación de 30 vecinos; centro: cobertura empírica de intervalos predictivos de NUEVA CÉLULA del 95%; derecha: su anchura media. Unidades: ln(1+conteo de imagen). No se inventan intervalos para kNN. La comparación corresponde a un donante, no a validación entre animales; se conservan mala calibración e intervalos amplios.'
  export(fig,'spatial_bayes_performance'+suffix,capen,capes,manifest,language)
  fig,axes=plt.subplots(2,2,figsize=(13.4,10.2),layout='constrained')
  for i,(region,obj)in enumerate(regions.items()):
   cal=obj['calibration'];cal=cal[cal.config.eq('primary')];m=obj['metrics'];m=m[m.config.eq('primary')&m.model.eq('bayes_spatial')]
   for gene,cc in cal[cal.model.eq('bayes_spatial')].groupby('gene'):axes[i,0].plot(cc.nominal,cc.observed_coverage,c='#72A6A5',alpha=.35,lw=.8)
   for model,marker in [('bayes_spatial','o'),('anatomy_constant','s')]:
    cc=cal[cal.model.eq(model)].groupby('nominal').observed_coverage.mean();axes[i,0].plot(cc.index,cc,color=colors[model],marker=marker,label='Bayes espacial'if es and model=='bayes_spatial'else'Media anatómica'if es else model)
   axes[i,0].plot([0,1],[0,1],ls='--',c='#333333',lw=1);axes[i,0].set(xlim=(0,1),ylim=(0,1.01),xlabel='Cobertura nominal'if es else'Nominal coverage',ylabel='Cobertura observada'if es else'Observed coverage',title=f'{region} · '+('Calibración predictiva'if es else'Predictive calibration'));axes[i,0].legend(fontsize=8)
   ys=np.arange(len(m));axes[i,1].scatter(m.coverage_zero,ys,c='#4B6870',marker='s',label='Ceros medidos'if es else'Measured zeros');axes[i,1].scatter(m.coverage_positive,ys,c='#B85121',marker='o',label='Conteos positivos'if es else'Positive counts');axes[i,1].set_yticks(ys,m.gene);axes[i,1].invert_yaxis();axes[i,1].set(xlim=(0,1.02),xlabel='Inclusión condicional: intervalo 95%'if es else'Conditional inclusion: 95% interval',title=f'{region} · '+('Ceros frente a positivos'if es else'Zeros versus positive counts'));axes[i,1].legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.15),ncol=2,frameon=False)
   for ax in axes[i]:ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.15)
  fig.suptitle('Cobertura dominada por ceros y recuperación de la cola positiva'if es else'Zero-dominated coverage and positive-tail recovery',fontsize=15,weight='bold')
  capen='Calibration uses only held-out whole-section cells. Thin curves show individual genes; bold curves average the ten gene-specific empirical coverages without assuming independent genes/cells. The diagonal is nominal calibration. Right: inclusion in the marginal 95%interval, separately among realized zeros/positives, is a tail-recovery diagnostic; neither conditional fraction is required to equal 95% even under a calibrated marginal model. No interval is clipped at zero or recalibrated using test outcomes. No binomial confidence interval treating cells as independent donors is asserted.'
  capes='La calibración utiliza únicamente células de cortes omitidos. Líneas finas: genes individuales; líneas gruesas: promedio de diez coberturas génicas sin suponer independencia. La diagonal representa calibración nominal. Derecha: inclusión en el intervalo marginal 95%, separada entre ceros/positivos observados, como diagnóstico de recuperación de la cola. Ninguna fracción condicional debe necesariamente igualar 95% aunque el modelo marginal esté calibrado. No se truncan los intervalos ni se recalibran con datos de prueba; no se presentan intervalos binomiales que confundan células con donantes.'
  export(fig,'spatial_bayes_calibration'+suffix,capen,capes,manifest,language)
  fig,axes=plt.subplots(2,2,figsize=(14.4,10.1),layout='constrained')
  for i,(region,obj)in enumerate(regions.items()):
   m=obj['metrics'];m=m[m.model.eq('bayes_spatial')];genes=obj['genes'].gene.tolist();configs=list(CONFIGS)
   base=m[m.config.eq('primary')].set_index('gene').loc[genes].rmse.to_numpy()
   errors=np.array([m[m.config.eq(c)].set_index('gene').loc[genes].rmse.to_numpy()/base for c in configs]).T
   cover=np.array([m[m.config.eq(c)].set_index('gene').loc[genes].predictive_95_coverage.to_numpy()for c in configs]).T
   for j,(arr,title)in enumerate([(errors,'RMSE / '+('principal'if es else'primary')),(cover,'Cobertura 95%'if es else'95% coverage')]):
    im=axes[i,j].imshow(arr,aspect='auto',cmap='viridis',vmin=0 if j else min(.9,float(errors.min())),vmax=1 if j else max(1.1,float(errors.max())));axes[i,j].set_yticks(np.arange(len(genes)),genes);axes[i,j].set_xticks(np.arange(5),configs,rotation=25,ha='right',fontsize=8);axes[i,j].set_title(f'{region} · {title}',loc='left');fig.colorbar(im,ax=axes[i,j],shrink=.8)
  fig.suptitle('Sensibilidad fijada antes del ajuste · no selección del mejor resultado'if es else'Predeclared prior/length sensitivity · no best-result selection',fontsize=15,weight='bold')
  capen='All five fixed configurations are re-fitted under the same guarded whole-section folds. Short/long lengths multiply all physical RBF lengths by 0.5/2; strong/weak shrinkage multiplies spatial coefficient prior variances by 0.25/4. Primary is retained regardless of results. Anatomy and section-prior multipliers are unchanged. Left: held-out RMSE divided by primary RMSE (lower is better); right: actual 95% predictive coverage. Color scales are explicitly different for the two metrics. Exact metrics, prior values and fold-specific basis centers are downloadable.'
  capes='Las cinco configuraciones fijas se reajustan con los mismos bloques protegidos. Longitudes cortas/largas multiplican las longitudes RBF por 0.5/2; regularización fuerte/débil multiplica la varianza previa de coeficientes espaciales por 0.25/4. Se conserva el modelo principal cualquiera sea el resultado. No cambian las varianzas anatómica ni de corte. Izquierda: RMSE omitido dividido por el principal; derecha: cobertura predictiva 95% observada. Las escalas de color difieren explícitamente entre métricas; se entregan valores y centros exactos por ajuste.'
  export(fig,'spatial_bayes_sensitivity'+suffix,capen,capes,manifest,language)
  fig,axes=plt.subplots(2,4,figsize=(17,9),layout='constrained')
  for i,(region,obj)in enumerate(regions.items()):
   c=obj['cells'];order=c.groupby('brain_section_label').experimental_z.median().sort_values().index;sec=order[len(order)//2];mask=c.brain_section_label.eq(sec).to_numpy();gene='Tll1'if region=='HY'else'Fosb';j=obj['genes'].gene.tolist().index(gene);f=obj['full'];v=[obj['y'][mask,j],f['posterior_mean'][mask,j],f['ci_upper'][mask,j]-f['ci_lower'][mask,j],f['predictive_upper'][mask,j]-f['predictive_lower'][mask,j]]
   cap=max(float(np.quantile(np.r_[v[0],v[1]],.99)),.1);low=min(float(v[1].min()),0.)
   for k,ax in enumerate(axes[i]):
    im=ax.scatter(c.z[mask],-c.y[mask],c=v[k],s=7,cmap='viridis',vmin=low if k<2 else 0,vmax=cap if k<2 else max(float(v[k].max()),.01),rasterized=True)
    ax.set_aspect('equal');ax.set(xlabel='ML (mm)',ylabel='−DV (mm)');title=(['RNA medido','Media latente posterior','Anchura intervalo creíble 95%','Anchura predictiva 95%']if es else['Measured RNA','Posterior latent mean','95% credible width','95% predictive width'])[k];ax.set_title(f'{region} · {gene}\n{title}',fontsize=10);fig.colorbar(im,ax=ax,shrink=.65,label='ln(1 + conteo)'if es else'ln(1 + count)')
   axes[i,0].text(0,-.19,f'{sec} · n={mask.sum()}',transform=axes[i,0].transAxes,fontsize=8)
  fig.suptitle('Solo coordenadas celulares reales · posterior de referencia completa'if es else'Only real cell coordinates · posterior fitted to the full reference',fontsize=15,weight='bold')
  capen='Illustrative median physical section in acquisition order for each region, selected without gene outcomes; every original cell of that section is plotted at its registered ML/−DV coordinate. Tll1 and Fosb were fixed as the leading measured relevant candidates. First two panels share a gene-specific 99th-percentile color cap (values retained); uncertainty widths each have their own explicit scale. These are all-reference fitted posterior summaries, NOT held-out predictions or measured PM responses. Credible width describes the section-adjusted latent mean; predictive width additionally includes residual/new-section variation. No expression is painted between cells, and no extra resolution is claimed.'
  capes='Corte físico mediano según orden de adquisición de cada región, seleccionado sin observar expresión; se representan todas sus células originales en coordenadas ML/−DV registradas. Se fijaron Tll1 y Fosb como primeros candidatos medidos relevantes. Los dos primeros paneles comparten un límite de color por gen del percentil 99; se conservan valores originales. Las anchuras tienen escalas propias explícitas. Son resúmenes posteriores ajustados a toda la referencia, NO predicciones omitidas ni RNA expuesto medido. La anchura creíble corresponde a la media latente; la predictiva añade variación residual y de nuevo corte. No se pinta expresión entre células ni se afirma mayor resolución.'
  export(fig,'spatial_bayes_posterior_maps'+suffix,capen,capes,manifest,language)
 dump(OUT/'figure_manifest.json',manifest)

def narratives(summary):
 linesen=['# Bayesian unexposed spatial reference','\n## Methods\n','The model fits natural-log1p imaging counts with anatomical and anisotropic spatial-basis predictors and proper physical-section nuisance effects. Exact Normal–Inverse-Gamma conjugate updates produce Student-t latent-mean and new-cell intervals. All priors, kernel lengths and candidate genes are fixed; whole-section blocks and their neighboring guards are excluded before learning feature scales, knots or coefficients. Expression-derived cell labels are excluded.','\n## Results\n']
 lineses=['# Referencia espacial bayesiana sin exposición','\n## Metodología\n','Se ajusta el ln(1+conteo de imagen) con anatomía, bases espaciales anisotrópicas y efectos propios de corte con distribuciones previas propias. La actualización conjugada Normal–Gamma inversa produce intervalos t de Student distintos para media latente y nueva célula. Se fijan genes, longitudes y distribuciones previas; se omiten bloques completos de cortes y vecinos antes de aprender escalas, centros o coeficientes. No se utilizan etiquetas celulares derivadas de expresión.','\n## Resultados\n']
 for region,s in summary['regions'].items():
  h=s['heldout_results'];rmse=np.median([v['rmse']for v in h]);cov=np.median([v['predictive_95_coverage']for v in h]);positive=np.median([v['coverage_positive']for v in h]);negative=np.median([v['negative_predictive_probability_mean']for v in h]);wins=sum(v['rmse']<v['anatomy_baseline_rmse']for v in h)
  linesen.append(f"{region}: all {s['n_cells']:,} real cells across {s['n_sections']} physical sections and 10 fixed measured genes were evaluated. Median gene RMSE = {rmse:.3f}; median empirical 95% predictive coverage = {cov:.1%}, versus {positive:.1%} among positive cells. Median Gaussian predictive mass below zero = {negative:.1%}. Spatial RMSE improves over the anatomy-only baseline for {wins}/10 genes. These numerical outcomes are retained even when calibration is poor.")
  lineses.append(f"{region}: se evaluaron las {s['n_cells']:,} células reales de {s['n_sections']} cortes físicos y 10 genes medidos fijos. Mediana del RMSE por gen = {rmse:.3f}; mediana de cobertura predictiva 95% = {cov:.1%}, frente a {positive:.1%} entre células positivas. La mediana de probabilidad predictiva gaussiana bajo cero es {negative:.1%}. El RMSE espacial mejora la referencia anatómica en {wins}/10 genes. Se conservan estos resultados aunque la calibración sea deficiente.")
 linesen+=['\n## Interpretation and limits\n','The Gaussian likelihood is an inadequate working baseline for positive-cell prediction in these data, not a validated reconstruction or count-generating model. Its fitted means and latent intervals can also be negative. Impossible predictive mass below zero, weak RMSE improvement and poor positive-tail recovery together limit its usefulness, even where overall coverage looks close to 95%. Positive-only interval inclusion is an outcome-conditioned tail diagnostic, not a requirement to achieve nominal 95% coverage within that subset. Credible intervals are conditional on the specified basis/priors and do not establish uncertainty across animals or registration accuracy. Held-out sections quantify within-donor spatial transport only; they do not validate a PM response. Full-reference fitted arrays are distinct from all held-out arrays/metrics. No cell is generated, relocated or given a measured exposure response.','\n## Sources\n','[Original atlas](https://doi.org/10.1038/s41586-023-06808-9); [conjugate Gaussian analysis](https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf); [Student-t implementation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html).']
 lineses+=['\n## Interpretación y límites\n','La verosimilitud gaussiana es una referencia insuficiente para predecir células positivas en estos datos, no una reconstrucción validada ni un modelo generador de conteos. Las medias ajustadas y sus intervalos latentes también pueden ser negativos. La masa predictiva imposible bajo cero, la escasa mejora del RMSE y la recuperación deficiente de la cola positiva limitan conjuntamente su utilidad, aunque la cobertura global parezca próxima al 95%. La inclusión del intervalo entre positivos es un diagnóstico condicionado al resultado; no exige cobertura nominal 95%en ese subconjunto. Los intervalos creíbles son condicionales a bases y distribuciones previas; no describen incertidumbre entre animales ni precisión de registro. Omitir cortes evalúa transporte espacial dentro del donante; no valida una respuesta PM. Las matrices ajustadas a toda la referencia se separan de predicciones y métricas omitidas. No se generan ni desplazan células ni se les atribuye una respuesta de exposición medida.','\n## Fuentes\n','[Atlas original](https://doi.org/10.1038/s41586-023-06808-9); [análisis gaussiano conjugado](https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf); [implementación t de Student](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html).']
 for lang,lines in [('en',linesen),('es',lineses)]:(OUT/f'narrative_{lang}.md').write_text('\n\n'.join(lines)+'\n')

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--skip-figures',action='store_true');args=parser.parse_args()
 for p in [OUT,META,WEB,OUT/'figures']:p.mkdir(parents=True,exist_ok=True)
 start=time.monotonic();dump(META/'methods.json',methods());dump(OUT/'summary.json',{'complete':False,'status':'fitting','regions':{}})
 with threadpool_limits(limits=2):regions={r:run_region(r)for r in ['HY','ACA']}
 if not args.skip_figures:figures(regions)
 summary={'complete':True,'created_UTC':datetime.now(timezone.utc).isoformat(),'model':'reference_only_conjugate_spatial_bayes','donors':1,'exposure_model':False,'regions':{r:o['summary']for r,o in regions.items()},'elapsed_seconds':time.monotonic()-start,'performance_success_required_for_integrity_pass':False}
 dump(OUT/'summary.json',summary);narratives(summary)
 checks={r:o['checks']for r,o in regions.items()};passed=all(all(c.values())for c in checks.values());dump(OUT/'validation.json',{'passed':passed,'scope':'Implementation, joins, posterior algebra, leakage boundaries and transport integrity; NOT scientific-performance certification. Calibration flags are separately retained per gene.','regions':checks})
 inputs=sorted(set(p for o in regions.values()for p in o['inputs']));dump(META/'source_manifest.json',{'sources':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p),'bytes':p.stat().st_size}for p in inputs],'parent_atlas_manifest':'metadata/spatial_atlas/download_manifest.json','parent_receipts':'metadata/spatial_atlas/download_receipts.json','method_sources':methods()['sources'],'script_sha256':sha(Path(__file__)),'software':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},'selection':'Frozen10 measured candidates per region, from prior external bulk-result ranking; no current atlas-expression gene selection.'})
 dump(WEB/'manifest.json',{'schema_version':1,'units':'ln(1+raw imaging count)','regions':{r:{'metadata':r+'/metadata.json','cell_count':o['summary']['n_cells'],'cell_ids_sha256':o['summary']['cell_ids_sha256']}for r,o in regions.items()},'exposure_model':False})
 print(json.dumps(clean(summary),indent=2),flush=True)
 if not passed:raise SystemExit('Implementation integrity checks failed')
if __name__=='__main__':main()
