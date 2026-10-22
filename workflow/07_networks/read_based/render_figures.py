#!/usr/bin/env python3
"""Render current network configuration, uncertainty and matched consensus diagrams."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize,TwoSlopeNorm
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from scipy.linalg import svd

SC=['unadjusted','space_depth','space_depth_environment']
LABEL=['Original CLR','Locality + depth','Locality + depth + 5 axes']
COLORS=['#7e8994','#2171b5','#c46a28']
ENDS=['mean_abs_edge_weight','natural_connectivity_abs_weighted']
BLUE=LinearSegmentedColormap.from_list('HI_blue',['#deebf7','#6baed6','#2171b5','#08306b'])
DIVERGE=LinearSegmentedColormap.from_list('slope',['#2166ac','#f7f7f7','#b35806'])

def save(fig,out,name):
    for ext in ['png','pdf','svg']:fig.savefig(out/f'{name}.{ext}',dpi=450,facecolor='white')
    fig.savefig(out/f'{name}_preview.png',dpi=150,facecolor='white');plt.close(fig)

def pstr(p):return '<0.001' if p<.001 else f'{p:.3f}'

def force_layout(adjacency):
    """Fixed-seed force layout of the union of consensus edges; identical coordinates in all panels."""
    active=np.flatnonzero(adjacency.sum(1)>0);isolated=np.flatnonzero(adjacency.sum(1)==0)
    graph=adjacency[np.ix_(active,active)];n=len(active)
    rng=np.random.default_rng(809200);pos=rng.normal(size=(n,2))*.2;k=np.sqrt(1/n)
    for iteration in range(300):
        delta=pos[:,None,:]-pos[None,:,:];dist=np.maximum(np.linalg.norm(delta,axis=2),.01)
        force=(delta*(k*k/dist**2-graph*dist/k)[:,:,None]).sum(1)-.05*pos
        length=np.maximum(np.linalg.norm(force,axis=1),1e-12)
        pos+=force/length[:,None]*np.minimum(length,.08*(1-iteration/300))[:,None]
        pos-=pos.mean(0)
    radius=np.linalg.norm(pos,axis=1);scaled=radius/max(radius.max(),1e-12)
    # A monotone radial expansion improves separation in the dense centre.
    # Coordinates are a display aid, not a measure of ecological distance.
    xy=np.zeros((len(adjacency),2));xy[active]=pos/np.maximum(radius[:,None],1e-12)*(.88*scaled**.4)[:,None]
    # Isolated nodes have no attractive force: keep them visible on a fixed outer ring.
    angle=np.linspace(0,2*np.pi,len(isolated),endpoint=False)+np.pi/2
    xy[isolated]=np.c_[np.cos(angle),np.sin(angle)]*.99
    return xy

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();r=args.run;out=args.output;out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none',
        'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,'xtick.labelsize':7.5,'ytick.labelsize':7.5})
    meta=pd.read_csv(r/'tables/809_metadata_aligned.tsv',sep='\t')
    objects={s:np.load(r/'objects'/f'809_{s}_networks.npz') for s in SC}
    table=pd.read_csv(r/'tables/809_global_configuration.tsv',sep='\t')
    sens=pd.read_csv(r/'tables/810_parameter_configuration.tsv',sep='\t')
    assoc=pd.read_csv(r/'tables/809_index_associations.tsv',sep='\t')
    D=np.c_[np.ones(51),meta.depth_cm.eq(20),meta.depth_cm.eq(40),meta.locality.eq('Cozumel'),meta.locality.eq('Tuxpan')]
    Q=np.linalg.qr(D,mode='reduced')[0];norm=Normalize(meta.HI.min(),meta.HI.max());coords=[]
    fig=plt.figure(figsize=(7.3,7.2))
    for si,sid in enumerate(SC[1:]):
        ax=fig.add_axes([.09+.50*si,.61,.36,.29]);W=np.abs(objects[sid]['weights']);Z=W-Q@(Q.T@W)
        u,s,v=svd(Z,full_matrices=False);pc=u[:,:2]*s[:2];explained=s*s/(s*s).sum()*100
        # Sign orientation is deterministic, based on the largest loading, never on HI.
        for j in range(2):pc[:,j]*=1 if v[j,np.argmax(abs(v[j]))]>=0 else -1
        for depth,marker in [(5,'o'),(20,'s'),(40,'^')]:
            ix=meta.depth_cm.eq(depth).to_numpy();ax.scatter(pc[ix,0],pc[ix,1],c=meta.HI[ix],cmap=BLUE,norm=norm,
                marker=marker,s=26,edgecolors='#344656',linewidth=.35)
        ax.set_xlabel(f'PC1 ({explained[0]:.1f}%)');ax.set_ylabel(f'PC2 ({explained[1]:.1f}%)')
        ax.set_title(LABEL[si+1],loc='left',fontweight='bold',pad=12,fontsize=9)
        ax.text(-.18,1.11,'ab'[si],transform=ax.transAxes,fontweight='bold',fontsize=12)
        for j,row in meta.iterrows():coords.append(dict(construction=sid,sample_id=row.sample_id,HI=row.HI,PC1=pc[j,0],PC2=pc[j,1],PC1_percent=explained[0],PC2_percent=explained[1]))
    cax=fig.add_axes([.10,.535,.25,.012]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=BLUE),cax=cax,orientation='horizontal')
    cb.set_label('Taxonomic condition index (HI)',fontsize=8,labelpad=3);cb.ax.tick_params(labelsize=7,length=2);cb.outline.set_visible(False)
    handles=[Line2D([],[],marker=marker,linestyle='',markerfacecolor='#9ecae1',markeredgecolor='#344656',markersize=4,label=f'{d} cm') for d,marker in [(5,'o'),(20,'s'),(40,'^')]]
    fig.legend(handles=handles,loc='center',bbox_to_anchor=(.72,.54),ncol=3,frameon=False,fontsize=8)
    ax=fig.add_axes([.19,.12,.34,.285]);responses=['absolute_weight_euclidean','signed_weight_euclidean','retained_edge_sqrt_jaccard']
    for j,response in enumerate(responses):
        for k,sid in enumerate(SC[1:]):
            row=table[(table.construction==sid)&table.response.eq(response)].iloc[0];y=2*j+k*.55
            ax.plot(row.partial_R2*100,y,'o',ms=4.5,color=COLORS[k+1])
            ax.text(3.45,y,pstr(row.p_holm_six_new_configuration_tests),va='center',ha='center',fontsize=7)
            ax.text(4.65,y,pstr(row.p_maxT_six_new_configuration_tests),va='center',ha='center',fontsize=7)
    ax.set_xlim(0,5.3);ax.set_xticks([0,1,2,3]);ax.set_xlabel('Partial R² (%)',labelpad=6)
    ax.set_ylim(5.05,-.8);ax.set_yticks([.275,2.275,4.275],['Absolute\nweights','Signed\nweights','Retained\nedges']);ax.tick_params(axis='y',length=0,pad=5)
    ax.text(3.45,-.75,'Holm P',ha='center',fontsize=7);ax.text(4.65,-.75,'maxT P',ha='center',fontsize=7)
    ax.spines['left'].set_visible(False);ax.spines['bottom'].set_bounds(0,3)
    fig.text(.02,.455,'c',fontweight='bold',fontsize=12);fig.text(.08,.455,'Association with HI',fontweight='bold',fontsize=9)
    variants=['primary','ridge_005','ridge_020','density_0025','density_010','features_100','features_300','four_environment_axes','RPK_per_million']
    names=['Primary','Ridge 0.05','Ridge 0.20','Density 2.5%','Density 10%','100 KOs','300 KOs','Four axes','RPK per million']
    grid=np.full((len(variants),2),np.nan);significant=np.zeros_like(grid,dtype=bool)
    for i,variant in enumerate(variants):
        for j,sid in enumerate(SC[1:]):
            df=table if variant=='primary' else sens[sens.variant.eq(variant)]
            df=df[df.construction.eq(sid)&df.response.eq('absolute_weight_euclidean')]
            if len(df):
                grid[i,j]=100*df.partial_R2.iloc[0]
                field='p_maxT_six_new_configuration_tests' if variant=='primary' else 'p_maxT_all_sensitivity_configuration_tests'
                significant[i,j]=df[field].iloc[0]<.05
    ax=fig.add_axes([.80,.12,.15,.285]);bound=4
    im=ax.imshow(np.ma.masked_invalid(grid),cmap=BLUE,vmin=0,vmax=bound,aspect='auto')
    ax.set_yticks(range(len(names)),names,fontsize=7);ax.set_xticks([0,1],['L+D','L+D+E'],fontsize=6.8)
    ax.xaxis.tick_top();ax.tick_params(length=0,pad=4)
    for i in range(len(names)):
        for j in range(2):ax.text(j,i,'—' if np.isnan(grid[i,j]) else f'{grid[i,j]:.2f}'+('*' if significant[i,j] else ''),ha='center',va='center',fontsize=6.8,color='white' if grid[i,j]>.55*bound else '#263744')
    for sp in ax.spines.values():sp.set_visible(False)
    fig.text(.58,.455,'d',fontweight='bold',fontsize=12);fig.text(.63,.455,'Absolute-weight sensitivity',fontweight='bold',fontsize=9)
    fig.text(.63,.071,'Partial R² (%) · * maxT P < 0.05',fontsize=7)
    fig.legend(handles=[Line2D([],[],marker='o',linestyle='',color=COLORS[1],label='Locality + depth'),Line2D([],[],marker='o',linestyle='',color=COLORS[2],label='+ 5 environmental axes')],loc='lower left',bbox_to_anchor=(.07,.015),frameon=False,fontsize=7,ncol=2)
    fig.text(.09,.965,'Functional network configuration',fontsize=11,fontweight='bold')
    save(fig,out,'Functional_network_configuration')
    pd.DataFrame(coords).to_csv(out/'network_ordination_coordinates.tsv',sep='\t',index=False)

    refit_path=r/'tables/809_reconstruction_bootstrap_intervals.tsv'
    if refit_path.exists():
        refit=pd.read_csv(refit_path,sep='\t');fig=plt.figure(figsize=(7.3,5.8))
        for k,(pred,end) in enumerate([(p,e) for p in ['HI','MHI_local'] for e in ENDS]):
            col=k%2;row=k//2;ax=fig.add_axes([.25+.47*col,.58-.40*row,.24,.25])
            vals=[]
            for j,sid in enumerate(SC):
                z=assoc[(assoc.construction==sid)&assoc.predictor.eq(pred)&assoc.endpoint.eq(end)].iloc[0]
                t=refit[(refit.construction==sid)&refit.predictor.eq(pred)&refit.endpoint.eq(end)].iloc[0]
                ax.plot([z.ci_low,z.ci_high],[2-j+.11]*2,color=COLORS[j],lw=3,alpha=.32,solid_capstyle='round')
                ax.plot([t.refit_ci_low,t.refit_ci_high],[2-j-.05]*2,color=COLORS[j],lw=1.1)
                ax.plot(z.standardized_beta,2-j,'o',ms=4.5,color=COLORS[j]);vals.extend([z.ci_low,z.ci_high,t.refit_ci_low,t.refit_ci_high])
            ax.axvline(0,color='#969da4',lw=.7,ls='--');ax.set_ylim(-.45,2.45);ax.set_yticks([2,1,0],LABEL if col==0 else ['']*3)
            ax.tick_params(axis='y',length=0,pad=6);ax.set_xlim(min(-.6,min(vals)-.1),max(.8,max(vals)+.1));ax.set_xlabel('Standardized slope')
            ax.set_title(('HI' if pred=='HI' else 'MHI')+' · '+('Mean weight' if end==ENDS[0] else 'Natural connectivity'),loc='left',fontsize=9,pad=12,fontweight='bold')
            fig.text(.03+.51*col,.89-.40*row,'abcd'[k],fontweight='bold',fontsize=12)
        fig.text(.08,.055,'Thick pale line: 4,999 fixed-network bootstrap draws\nThin line: 499 complete network reconstructions (exploratory)',fontsize=8,linespacing=1.7)
        fig.text(.08,.965,'Network summaries and taxonomic or environmental condition',fontsize=11,fontweight='bold')
        save(fig,out,'Supplementary_Figure10')

    n=len(objects[SC[0]]['nodes']);a=objects[SC[0]]['node_a_index'];b=objects[SC[0]]['node_b_index']
    union=np.zeros((n,n))
    for sid in SC:
        mask=objects[sid]['stable_mask'];union[a[mask],b[mask]]=1;union[b[mask],a[mask]]=1
    xy=force_layout(union);fig,axs=plt.subplots(1,3,figsize=(9.3,4.1))
    for j,(ax,sid) in enumerate(zip(axs,SC)):
        z=objects[sid];mask=z['stable_mask'];v=z['full_weights'][mask]
        segments=np.stack([xy[a[mask]],xy[b[mask]]],axis=1)
        ax.add_collection(LineCollection(segments,colors=np.where(v>0,'#2171b5','#c46a28'),linewidths=.35+.8*abs(v)/max(abs(v)),alpha=.48))
        ax.scatter(xy[:,0],xy[:,1],s=7,c='#344c60',edgecolors='white',linewidths=.15,zorder=3)
        ax.set_xlim(-1.08,1.08);ax.set_ylim(-1.08,1.08);ax.set_aspect('equal');ax.axis('off')
        ax.set_title(LABEL[j]+f'\n{int(mask.sum())} consensus edges',fontsize=9,pad=8)
        ax.text(-.05,1.13,'abc'[j],transform=ax.transAxes,fontweight='bold',fontsize=12)
    fig.subplots_adjust(left=.025,right=.985,top=.80,bottom=.15,wspace=.08)
    fig.text(.04,.955,'Consensus coabundance networks',fontsize=11,fontweight='bold')
    fig.text(.04,.055,'200 KOs · shared node positions · edge frequency ≥50% and sign consistency ≥80%',fontsize=8)
    fig.legend(handles=[Line2D([],[],color='#2171b5',label='Positive'),Line2D([],[],color='#c46a28',label='Negative')],loc='lower right',bbox_to_anchor=(.98,.003),ncol=2,frameon=False,fontsize=8)
    save(fig,out,'Functional_network_consensus')
    pd.DataFrame({'KO':objects[SC[0]]['nodes'],'x':xy[:,0],'y':xy[:,1]}).to_csv(out/'consensus_node_coordinates.tsv',sep='\t',index=False)
    (out/'figure_specification.json').write_text(json.dumps({'configuration':'185 mm wide, unsupervised residual-edge PCA and full-dimensional tests',
        'uncertainty':'fixed-network and reconstruction intervals shown separately','consensus':'common fixed-seed layout; all 200 nodes; aggregate partial-correlation signs',
        'HI_used_for_layout_or_PCA':False,'formats':['png','pdf','svg']},indent=2)+'\n')
    print('Rendered figures in',out)

if __name__=='__main__':main()
