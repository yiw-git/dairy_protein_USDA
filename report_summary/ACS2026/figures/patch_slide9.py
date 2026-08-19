import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Ellipse, Rectangle, Wedge

RED="#E21833"; DARK="#000000"; GREY="#000000"; ZEIN="#E8A33D"; BLUE="#2E6DA4"; PUR="#7A2E8E"; GRN="#1B7F5C"
W,H=16.0,4.05
fig=plt.figure(figsize=(9.0,2.278), dpi=400)
bg=fig.add_axes([0,0,1,1]); bg.set_xlim(0,W); bg.set_ylim(0,H); bg.axis('off')
fx=lambda x:x/W; fy=lambda y:y/H

# ---------- (1) whole particle ----------
cx,cy,R=2.05,2.32,1.06
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
bg.text(cx, cy-0.24, "70 nm", ha='center', fontsize=10.5, color=DARK, fontweight='bold')
bg.text(cx, 3.88, "One whole nanoparticle", ha='center', fontsize=11.5, fontweight='bold', color=DARK)
bg.text(cx, 0.55, "~6,000 zein chains\n~3 × 10⁷ protein atoms\n~6 × 10⁷ atoms with water",
        ha='center', va='center', fontsize=9.5, color=GREY, linespacing=1.35)

# ---------- arrow ----------
bg.add_patch(FancyArrowPatch((3.50,2.78),(4.45,2.78), arrowstyle='-|>', mutation_scale=16, lw=2.0, color=RED))
bg.text(3.97, 3.02, "zoom in", ha='center', fontsize=9.0, color=RED, fontweight='bold')

# ---------- (2) the patch ----------
px,py=6.55,2.42
bg.add_patch(FancyBboxPatch((px-1.62, py-0.90), 3.24, 1.80, boxstyle="round,pad=0,rounding_size=0.10",
                            facecolor="#FAFAFA", edgecolor="#D8D8D8", lw=1.0))
# gently curved raft of zein helices
xs=np.linspace(px-1.34, px+1.34, 9)
for i,x in enumerate(xs):
    yb=py-0.37 + 0.050*np.cos((x-px)/1.6*np.pi/2)
    bg.add_patch(FancyBboxPatch((x-0.12, yb-0.24), 0.24, 0.48,
                 boxstyle="round,pad=0,rounding_size=0.10", facecolor=ZEIN, edgecolor='#C4832A', lw=0.6))
bg.annotate("", xy=(px-1.50, py-0.76), xytext=(px+1.50, py-0.76),
            arrowprops=dict(arrowstyle='<->', color=DARK, lw=0.9))
bg.text(px, py-0.82, "12 nm", ha='center', fontsize=10.5, color=DARK, fontweight='bold')
for dx,c in [(-0.78,PUR),(0.0,GRN),(0.78,PUR)]:
    bg.add_patch(Ellipse((px+dx, py+0.44), 0.48, 0.33, facecolor=c, edgecolor='white', lw=0.7))
bg.text(px, py+0.78, "2–4 milk proteins", ha='center', fontsize=9.5, color=GREY)
bg.text(px, 3.88, "What we actually simulate", ha='center', fontsize=11.5, fontweight='bold', color=DARK)
bg.text(px, 0.55, "9 zein chains · ~2 × 10⁴ CG beads\nsurface is flat to within 0.5 nm here\n≈ 1% of the particle surface",
        ha='center', va='center', fontsize=9.5, color=GREY, linespacing=1.35)

# ---------- (3) cost comparison, log scale ----------
ax=fig.add_axes([fx(10.85), fy(1.50), fx(4.25), fy(1.70)])
labels=["all-atom\nwhole particle","coarse-grained\nwhole particle","coarse-grained\n12 nm patch"]
vals=[5.8e7, 5.8e6, 2.1e4]
colors=["#B0B0B0", "#7FA6C9", RED]
y=np.arange(3)[::-1]
ax.barh(y, vals, color=colors, height=0.62)
ax.set_xscale('log'); ax.set_xlim(3e3, 3e8)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.5)
ax.set_xticks([1e4,1e5,1e6,1e7,1e8])
ax.set_xticklabels(["10⁴","10⁵","10⁶","10⁷","10⁸"], fontsize=8.5)
ax.tick_params(length=2, pad=1.5)
for sp in ['top','right','left']: ax.spines[sp].set_visible(False)
ax.set_xlabel("particles in the simulation box", fontsize=9.5, labelpad=2.0)
txt=["6 × 10⁷ atoms","6 × 10⁶ beads","2 × 10⁴ beads"]
for yy,v,t in zip(y, vals, txt):
    ax.text(v*1.5, yy, t, fontsize=8.5, color=DARK, va='center')
bg.text(12.98, 3.88, "The cost we avoid", ha='center', fontsize=11.5, fontweight='bold', color=DARK)
bg.text(12.98, 0.55, "~2,700× fewer particles, plus a\n20 fs timestep instead of 2 fs",
        ha='center', va='center', fontsize=9.5, color=GREY, linespacing=1.35)

out="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures/patch_slide9.png"
plt.savefig(out, transparent=True, bbox_inches='tight', pad_inches=0.02)
print("saved")
