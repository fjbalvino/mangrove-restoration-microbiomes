from pathlib import Path
import json,hashlib,platform,shutil
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
P=Path(__file__).resolve().parent;O=Path((P/'LATEST.txt').read_text())
R=pd.read_csv(O/'eligibility.tsv',sep='\t').set_index('KO');T=pd.read_csv(O/'tests.tsv',sep='\t');E=pd.read_csv(O/'effects.tsv',sep='\t');V=pd.read_csv(O/'profile_values.tsv',sep='\t');D=pd.read_csv(O/'paired_depth_values.tsv',sep='\t');scope=pd.read_csv(P/'scope.tsv',sep='\t')
stages=['Degraded','Early restoration','Intermediate restoration','Advanced restoration'];short=['Degraded','Early','Intermediate','Advanced']
sel=[('C','K01601','rbcL · RuBisCO'),('C','K14138','acsB · acetyl-CoA synthase'),('N','K02588','nifH · nitrogenase component'),('N','K15864','nirS · nitrite/hydroxylamine transformations'),('N','K00376','nosZ · nitrous oxide reduction'),('S','K11180','dsrA · sulfite transformations'),('S','K00394','aprA · APS transformations'),('S','K17224','soxB · sulfur oxidation component'),('CH₄/C₁','K00399','mcrA · methyl-CoM reductase'),('CH₄/C₁','K10944','pmoA/amoA · shared monooxygenase'),('P','K02040','pstS · phosphate transport component'),('P','K01077','phoA · phosphomonoesterase')]
# Ecological examples fixed by process identity, not by P or effect size.
(O/'display_selection.json').write_text(json.dumps(sel,indent=2))
M=np.full((len(sel),4),np.nan);labels=[]
for i,(cy,k,name) in enumerate(sel):
 labels.append(f'{cy} | {name}'+(' †' if k in R.index and R.loc[k,'available'] and R.loc[k,'quantification']!='single_KO' else ''))
 for j,s in enumerate(stages):
  z=E[(E.variant=='primary')&(E.KO==k)&(E.target=='restoration')&(E.stage==s)]
  if len(z):M[i,j]=z.log2_ratio.iloc[0]
cmap=LinearSegmentedColormap.from_list('div',['#2166ac','#ffffff','#b2182b']);cmap.set_bad('#dddddd')
plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(12,7.6));im=ax.imshow(M,cmap=cmap,vmin=-2,vmax=2,aspect='auto');ax.set_yticks(range(len(sel)),labels);ax.set_xticks(range(4),short);ax.set_title('a  Estimated representation relative to conserved profiles',loc='left',pad=16)
for i in range(len(sel)):
 for j in range(4):ax.text(j,i,'—' if np.isnan(M[i,j]) else f'{2**M[i,j]:.2f}×',ha='center',va='center',color='white' if abs(M[i,j])>1.5 else '#202020')
c=fig.colorbar(im,ax=ax,fraction=.04,pad=.04);c.set_ticks([-2,-1,0,1,2],labels=['¼','½','1','2','4']);c.set_label('Stage / conserved (common scale)')
fig.text(.02,.04,'51 metagenomes, 17 independent profiles; locality-adjusted estimates. Conserved = 1.\n† Provisional fractional assignment. Grey: no eligible estimate. Colours do not indicate statistical support.',fontsize=10)
fig.tight_layout(rect=[0,.1,1,1]);fig.savefig(O/'Figure6a_stage_estimates.png',dpi=220);fig.savefig(O/'Figure6a_stage_estimates.pdf');plt.close(fig)
# Whole tested universe, abundance-independent order. Each KO appears once.
keys=R.index[R.eligible].tolist();whole=E[(E.variant=='primary')&(E.target=='restoration')].pivot(index='KO',columns='stage',values='log2_ratio').reindex(index=keys,columns=stages)
fig,ax=plt.subplots(figsize=(7,9));im=ax.imshow(whole,cmap=cmap,vmin=-2,vmax=2,aspect='auto');ax.set_xticks(range(4),short);ax.set_yticks([]);ax.set_ylabel(f'{len(keys)} unique eligible KOs (KO order)');ax.set_title('b  Expanded panel: every eligible KO');fig.colorbar(im,ax=ax,label='log₂(stage / conserved)',shrink=.7);fig.tight_layout();fig.savefig(O/'Figure6b_all_KOs.png',dpi=200);plt.close(fig)
# Genuine uncertainty is visible as independent profile variability; no unsupported pairwise CI.
subset=[x for x in sel if x[1] in keys][:8];fig,axes=plt.subplots(2,4,figsize=(14,7),squeeze=False);stage5=stages+['Conserved'];palette={'Carmen':'#d95f02','Cozumel':'#1b9e77','Tuxpan':'#7570b3'}
for ax,(_,k,name) in zip(axes.flat,subset):
 z=D[(D.KO==k)&(D.depth==40)]
 for j,s in enumerate(stage5):
  a=z[z.stage==s].sort_values('profile_id');ax.scatter(j+np.linspace(-.14,.14,len(a)),a.ratio,c=[palette[v] for v in a.locality],s=29);ax.plot([j-.22,j+.22],[a.ratio.median()]*2,color='black',lw=2)
 ax.axhline(1,color='gray',lw=.7);ax.set_yscale('log',base=2);ax.set_xticks(range(5),['D','E','I','A','C']);ax.set_title(name.split(' · ')[0]);ax.set_ylabel('Paired 40 / 5 cm')
fig.suptitle('Depth within stages: each point is one independent profile');fig.tight_layout();fig.savefig(O/'Depth_within_stages.png',dpi=200);plt.close(fig)
summary=T.groupby(['variant','target']).agg(tests=('p','size'),min_p=('p','min'),min_holm=('p_holm','min'),supported=('p_holm',lambda x:int((x<.05).sum()))).reset_index();summary.to_csv(O/'test_summary.tsv',sep='\t',index=False)
coverage=[]
for cy,g in scope.groupby('cycle'):
 a=R.reindex(sorted(g.KO.unique()));coverage.append(dict(cycle=cy,catalogue=len(a),available=int(a.available.sum()),eligible=int(a.eligible.sum()),exact_eligible=int(((a.quantification=='single_KO')&a.eligible).sum())))
pd.DataFrame(coverage).to_csv(O/'cycle_coverage.tsv',sep='\t',index=False)
# A subset of rule-level evidence is retained, independent from KO abundance tests.
for file in ['rule_summary.tsv','rule_sample_evidence.tsv.gz']:
 a=pd.read_csv(P.parent/'panel/10_audit_expanded_panel_results'/file,sep='\t');a=a[a.cycle.isin(['C','CH4_C1','N','S','P'])];a.to_csv(O/file,sep='\t',index=False)
# Numerical independent verification and nitrogen continuity before expanded multiplicity.
m=pd.read_csv(O/'samples.tsv',sep='\t');p=pd.read_csv(O/'profiles.tsv',sep='\t');X=np.column_stack([np.ones(len(p)),p.locality.eq('Cozumel'),p.locality.eq('Tuxpan'),*[p.restoration4.eq(s) for s in stages]]).astype(float)
for k in keys:
 y=V[V.KO==k].set_index('profile_id').loc[p.profile_id,'mean_log2'];b=np.linalg.lstsq(X,y,rcond=None)[0];e=E[(E.variant=='primary')&(E.target=='restoration')&(E.KO==k)].set_index('stage').loc[stages,'log2_ratio'];assert np.allclose(b[3:],e)
old=pd.read_csv(P.parent/'nitrogen_analysis/13_compare_nitrogen_results/tests.tsv',sep='\t');joined=old.merge(T,on=['variant','KO','target'],suffixes=('_old','_new'));assert len(joined)==len(old) and np.allclose(joined.p_old,joined.p_new)
assert T.groupby('variant').size().eq(len(keys)*3).all();assert T.p_holm.between(0,1).all()
validation=dict(unique_scope=len(R),available=int(R.available.sum()),eligible=len(keys),exact_eligible=int(((R.quantification=='single_KO')&R.eligible).sum()),family=len(keys)*3,nitrogen_p_values_identical=True,independent_OLS_verified=True,profiles=len(p),samples=len(m),bootstrap_calibration='No new expanded-panel finite-sample calibration. Exploratory inference; only 3-5 profiles per stage.',python=platform.python_version())
(O/'validation.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation));print(summary.to_string(index=False));print(pd.DataFrame(coverage).to_string(index=False))
fig,(ax,bx)=plt.subplots(1,2,figsize=(14,8),gridspec_kw={'width_ratios':[3,1]});im=ax.imshow(M,cmap=cmap,vmin=-2,vmax=2,aspect='auto');ax.set_yticks(range(len(sel)),labels,fontsize=10);ax.set_xticks(range(4),short,rotation=25,ha='right');ax.set_title('a  Ecological examples',loc='left')
for i in range(len(sel)):
 for j in range(4):ax.text(j,i,'—' if np.isnan(M[i,j]) else f'{2**M[i,j]:.2f}×',ha='center',va='center',fontsize=10,color='white' if abs(M[i,j])>1.5 else '#202020')
bx.imshow(whole,cmap=cmap,vmin=-2,vmax=2,aspect='auto');bx.set_yticks([]);bx.set_xticks(range(4),short,rotation=65,ha='right');bx.set_title(f'b  All {len(keys)} KOs',loc='left');bx.set_ylabel('Unique KOs, ordered by identifier');fig.subplots_adjust(left=.39,right=.91,bottom=.22,top=.91,wspace=.32);cax=fig.add_axes([.93,.28,.015,.5]);cb=fig.colorbar(im,cax=cax);cb.set_ticks([-2,-1,0,1,2],labels=['¼','½','1','2','4']);cb.set_label('Stage / conserved');fig.text(.02,.07,'Locality-adjusted estimates from 17 profiles (51 metagenomes). Conserved = 1; no pairwise tests.\n† Fractional quantification remains provisional. Grey: no eligible estimate. Colours saturate at ¼ and 4.\nNo primary omnibus test retained P < 0.05 after Holm correction across the expanded family.',fontsize=10);fig.savefig(O/'Figure6.png',dpi=220);fig.savefig(O/'Figure6.pdf');plt.close(fig)
# Leave-profile/locality-out sensitivity describes influence, not confidence intervals.
le=[]
for name,keep in [('profile:'+str(pid),p.profile_id.ne(pid).to_numpy()) for pid in p.profile_id]+[('locality:'+str(loc),p.locality.ne(loc).to_numpy()) for loc in p.locality.unique()]:
 xx=X[keep];pinv=np.linalg.pinv(xx)
 for k in keys:
  y=V[V.KO==k].set_index('profile_id').loc[p.profile_id,'mean_log2'].to_numpy();coef=pinv@y[keep]
  for j,s in enumerate(stages):
   contrast=np.eye(X.shape[1])[3+j];estimable=np.allclose(contrast,contrast@pinv@xx)
   le.append(dict(KO=k,omission=name,stage=s,estimable=estimable,log2_ratio=coef[3+j] if estimable else np.nan))
pd.DataFrame(le).to_csv(O/'leave_out_stage_effects.tsv',sep='\t',index=False)
# Explicit rule provenance: retain AND/OR expressions and unresolved process catalogue.
shutil.copy2(P.parent/'panel/operational_rules.tsv',O/'operational_rules.tsv');shutil.copy2(P.parent/'panel/curated_processes.tsv',O/'curated_processes.tsv');shutil.copy2(P/'scope.tsv',O/'scope.tsv')
# Compact manuscript version: full process descriptions remain in the atlas and expanded examples.
compact=[f"{cy} | {name.split(' · ')[0]}"+(' †' if k in R.index and R.loc[k,'available'] and R.loc[k,'quantification']!='single_KO' else '') for cy,k,name in sel]
fig,(ax,bx)=plt.subplots(1,2,figsize=(7,7.6),gridspec_kw={'width_ratios':[3.2,1]});im=ax.imshow(M,cmap=cmap,vmin=-2,vmax=2,aspect='auto');ax.set_yticks(range(len(sel)),compact,fontsize=10);ax.set_xticks(range(4),short,rotation=55,ha='right',fontsize=9);ax.set_title('a  Ecological examples',loc='left',fontsize=11)
for i in range(len(sel)):
 for j in range(4):ax.text(j,i,'—' if np.isnan(M[i,j]) else f'{2**M[i,j]:.2f}',ha='center',va='center',fontsize=10,color='white' if abs(M[i,j])>1.5 else '#202020')
bx.imshow(whole,cmap=cmap,vmin=-2,vmax=2,aspect='auto',interpolation='nearest');bx.set_yticks([]);bx.set_xticks(range(4),['D','E','I','A'],fontsize=9);bx.set_title('b  437 KOs',loc='left',fontsize=11);fig.subplots_adjust(left=.25,right=.87,bottom=.23,top=.93,wspace=.28);cax=fig.add_axes([.9,.3,.02,.52]);cb=fig.colorbar(im,cax=cax);cb.set_ticks([-2,-1,0,1,2],labels=['¼','½','1','2','4']);cb.ax.tick_params(labelsize=9);fig.text(.25,.03,'Stage / conserved; locality-adjusted estimates\n† Provisional quantification; grey: no estimate\nNo supported primary stage test after correction',fontsize=9);fig.savefig(O/'Figure6.png',dpi=300);fig.savefig(O/'Figure6.pdf');plt.close(fig)
