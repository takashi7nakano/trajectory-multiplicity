import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
mpl.rcParams.update({"font.size":8,"font.family":"serif"})
fig,ax=plt.subplots(figsize=(7.0,2.6)); ax.axis("off")
T=6; n=8; bs=2  # schematic: 6 layers shown, 8 positions, blocks of 2
removed={2}      # removed block index
import numpy as np
for t in range(T):
    for x in range(n):
        b=x//bs
        col = "#D55E00" if b in removed else "#DDDDDD"
        ax.add_patch(Rectangle((t,x),0.8,0.8,facecolor=col,edgecolor="white",lw=0.5))
# a surviving trajectory (avoids removed block)
traj=[1,0,1,3,4,5]
for t in range(T-1):
    ax.plot([t+0.4,t+1.4],[traj[t]+0.4,traj[t+1]+0.4],"-",color="#0072B2",lw=1.4)
ax.plot([t+1+0.4 for t in range(T)],[y+0.4 for y in traj],"o",color="#0072B2",ms=3)
ax.text(-0.5,n/2,"positions $x$",rotation=90,va="center",fontsize=8)
ax.text(T/2,-0.9,"layers $t=0,\\dots,T$",ha="center",fontsize=8)
ax.text(2.4,n+0.4,"removed block (dead-end, all layers)",color="#D55E00",fontsize=7)
ax.add_patch(Rectangle((0,-0.2),0.8,n+0.2,fill=False,edgecolor="none"))
ax.text(-0.5,n+0.4,"A ($\\mu_A$ uniform)",fontsize=7)
ax.set_xlim(-1,T+0.5); ax.set_ylim(-1.2,n+1)
fig.tight_layout(); fig.savefig("fig1_schematic.pdf",bbox_inches="tight")
fig.savefig("fig1_schematic.png",bbox_inches="tight"); print("fig1 written")
