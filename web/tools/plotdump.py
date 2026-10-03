import json, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=json.load(open(sys.argv[1])); fr=d['frames']; n=len(fr)
fig,axs=plt.subplots(1,n,figsize=(2.6*n,3.2))
for ax,f in zip(np.atleast_1d(axs),fr):
    P=np.array(f['P'])
    for st,c in [(2,'#1f6fb2'),(1,'#e8590c'),(0,'#2f9e44')]:
        m=P[:,2]==st; ax.scatter(P[m,0],P[m,1],s=0.2,c=c,lw=0)
    if f['B']: B=np.array(f['B']); ax.scatter(B[:,0],B[:,1],s=1,c='w',edgecolors='k',lw=0.2)
    if f['D']: D=np.array(f['D']); ax.scatter(D[:,0],D[:,1],s=0.5,c='r',lw=0)
    for s in f['segs']: ax.plot([s[0],s[2]],[s[1],s[3]],'k-',lw=2)
    ax.plot(d['xN'],d['yN'],'rv')
    st=f['stats']
    ax.set_xlim(0,d['W']); ax.set_ylim(0,d['H']); ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title('t=%.2f wall=%.2f H=%.1f\nU=%.2f Ue=%.2f air=%.3f noz=%.1f'%(f['t'],st['wallFrac'],st['H']*100,st['U'],st['Ueff'],st['air'],st['nozzleHits']),fontsize=6)
fig.tight_layout(); fig.savefig(sys.argv[2],dpi=110)
