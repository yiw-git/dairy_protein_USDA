import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Ellipse, Rectangle, Wedge

RED="#E21833"; DARK="#282828"; GREY="#636363"; ZEIN="#E8A33D"; BLUE="#2E6DA4"; PUR="#7A2E8E"; GRN="#1B7F5C"
W,H=16.0,3.78
fig=plt.figure(figsize=(9.0,2.13), dpi=400)
bg=fig.add_axes([0,0,1,1]); bg.set_xlim(0,W); bg.set_ylim(0,H); bg.axis('off')
fx=lambda x:x/W; fy=lambda y:y/H

# ---------- (1) whole particle ----------
cx,cy,R=2.05,2.05,1.28
bg.add_patch(Circle((cx,cy), R, facecolor=ZEIN, edgecolor='none', alpha=0.85))
rng=np.random.default_rng(5)
for _ in range(90):
    a=rng.random()*2*np.pi; d=rng.random()**0.5*(R-0.06)
    x0,y0=cx+d*np.cos(a), cy+d*np.sin(a); ang=rng.random()*2*np.pi
    bg.plot([x0-0.055*np.cos(ang),x0+0.055*np.cos(ang)],[y0-0.055*np.sin(ang),y0+0.055*np.sin(ang)],
            color='#C4832A', lw=0.7, solid_capstyle='round')
# 12 nm window marked on the rim (12/70 of the diameter)
half_ang=np.degrees(np.arcsin(6.0/35.0))
bg.add_patch(Wedge((cx,cy), R+0.05, 90-half_ang, 90+half_ang, width=0.10, facecolor=RED, edgecolor='none', zorder=5))
bg.annotate("", xy=(cx-R-0.06, cy), xytext=(cx+R+0.06, cy),
            arrowprops=dict(arrowstyle='<->', color=DARK, lw=0.9))
bg.text(cx, cy-0.20, "70 nm", ha='center', fontsize=8, color=DARK, fontweight='bold')
bg.text(cx, 3.58, "One whole nanoparticle", ha='center', fontsize=9.2, fontweight='bold', color=DARK)
bg.text(cx, 0.52, "~6,000 zein chains\n~3 × 10⁷ protein atoms\n~6 × 10⁷ atoms with water",
        ha='center', va='center', fontsize=7.4, color=GREY)

# ---------- arrow ----------
bg.add_patch(FancyArrowPatch((3.55,2.55),(4.55,2.55), arrowstyle='-|>', mutation_scale=14, lw=1.8, color=RED))
bg.text(4.05, 2.78, "zoom in", ha='center', fontsize=7, color=RED, fontweight='bold')

# ---------- (2) the patch ----------
px,py=6.65,2.15
bg.add_patch(FancyBboxPatch((px-1.75, py-1.05), 3.5, 2.1, boxstyle="round,pad=0,rounding_size=0.10",
                            facecolor="#FAFAFA", edgecolor="#D8D8D8", lw=1.0))
# gently curved raft of zein helices
xs=np.linspace(px-1.45, px+1.45, 9)
for i,x in enumerate(xs):
    yb=py-0.42 + 0.055*np.cos((x-px)/1.6*np.pi/2)
    bg.add_patch(FancyBboxPatch((x-0.13, yb-0.28), 0.26, 0.56,
                 boxstyle="round,pad=0,rounding_size=0.10", facecolor=ZEIN, edgecolor='#C4832A', lw=0.6))
bg.annotate("", xy=(px-1.62, py-0.90), xytext=(px+1.62, py-0.90),
            arrowprops=dict(arrowstyle='<->', color=DARK, lw=0.9))
bg.text(px, py-0.83, "12 nm", ha='center', fontsize=7.6, color=DARK, fontweight='bold')
for dx,c in [(-0.85,PUR),(0.0,GRN),(0.85,PUR)]:
    bg.add_patch(Ellipse((px+dx, py+0.52), 0.52, 0.36, facecolor=c, edgecolor='white', lw=0.7))
bg.text(px, py+0.90, "2–4 milk proteins", ha='center', fontsize=7.2, color=GREY)
bg.text(px, 3.58, "What we actually simulate", ha='center', fontsize=9.2, fontweight='bold', color=DARK)
bg.text(px, 0.52, "9 zein chains · ~2 × 10⁴ CG beads\nsurface is flat to within 0.5 nm here\n≈ 1% of the particle surface",
        ha='center', va='center', fontsize=7.4, color=GREY)

# ---------- (3) cost comparison, log scale ----------
ax=fig.add_axes([fx(10.55), fy(1.05), fx(4.55), fy(1.85)])
labels=["all-atom\nwhole particle","coarse-grained\nwhole particle","coarse-grained\n12 nm patch"]
vals=[5.8e7, 5.8e6, 2.1e4]
colors=["#B0B0B0", "#7FA6C9", RED]
y=np.arange(3)[::-1]
ax.barh(y, vals, color=colors, height=0.62)
ax.set_xscale('log'); ax.set_xlim(3e3, 3e8)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6.6)
ax.set_xticks([1e4,1e5,1e6,1e7,1e8])
ax.set_xticklabels(["10⁴","10⁵","10⁶","10⁷","10⁸"], fontsize=6.5)
ax.tick_params(length=2, pad=1.5)
for sp in ['top','right','left']: ax.spines[sp].set_visible(False)
ax.set_xlabel("particles in the simulation box", fontsize=7.4, labelpad=1.5)
txt=["6 × 10⁷ atoms","6 × 10⁶ beads","2 × 10⁴ beads"]
for yy,v,t in zip(y, vals, txt):
    ax.text(v*1.4, yy, t, fontsize=6.4, color=DARK, va='center')
bg.text(12.85, 3.58, "The cost we avoid", ha='center', fontsize=9.2, fontweight='bold', color=DARK)
bg.text(12.85, 0.30, "~2,700× fewer particles, plus a 20 fs timestep instead of 2 fs",
        ha='center', va='center', fontsize=7.2, color=GREY)

out="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures/patch_slide9.png"
plt.savefig(out, transparent=True, bbox_inches='tight', pad_inches=0.02)
print("saved")
