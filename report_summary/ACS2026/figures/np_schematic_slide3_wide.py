import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Ellipse, Polygon, Rectangle

RED="#E21833"; GOLD="#FFD200"; DARK="#282828"; GREY="#636363"; ZEIN="#E8A33D"; CAS="#4A7FB5"
W,H = 16.0, 4.0
fig, ax = plt.subplots(figsize=(9.0,2.25), dpi=400)
ax.set_xlim(0,W); ax.set_ylim(0,H); ax.axis('off')
CY, R = 1.85, 0.80

def shell(cx, cy, r, n=24, amp=0.08):
    for i in range(n):
        th=2*np.pi*i/n; t=np.linspace(0,1,36); rr=r+0.28*t; wob=amp*np.sin(6*np.pi*t)
        ax.plot(cx+rr*np.cos(th)-wob*np.sin(th), cy+rr*np.sin(th)+wob*np.cos(th),
                color=CAS, lw=1.3, solid_capstyle='round', zorder=2)

def core(cx, cy, r, seed=3):
    ax.add_patch(Circle((cx,cy), r, facecolor=ZEIN, edgecolor='none', zorder=3))
    rng=np.random.default_rng(seed)
    for _ in range(26):
        a=rng.random()*2*np.pi; d=rng.random()**0.5*(r-0.15)
        x0,y0=cx+d*np.cos(a), cy+d*np.sin(a); ang=rng.random()*2*np.pi
        ax.plot([x0-0.12*np.cos(ang), x0+0.12*np.cos(ang)],[y0-0.12*np.sin(ang), y0+0.12*np.sin(ang)],
                color='#C4832A', lw=1.9, solid_capstyle='round', zorder=4)

# ---- Panel 1 ----
cx1 = 2.05
shell(cx1, CY, R); core(cx1, CY, R)
ax.text(cx1, 3.70, "Zein–caseinate nanoparticle", ha='center', fontsize=11, fontweight='bold', color=DARK)
ax.annotate("packed zein chains\n(hydrophobic core)", xy=(cx1-0.34,CY-0.26), xytext=(0.02,0.34),
            fontsize=8.5, color=DARK, ha='left', arrowprops=dict(arrowstyle='-', color=GREY, lw=0.9))
ax.annotate("caseinate shell\n(charge + steric)", xy=(cx1+0.70,CY+0.72), xytext=(3.05,2.70),
            fontsize=8.5, color=DARK, ha='left', arrowprops=dict(arrowstyle='-', color=GREY, lw=0.9))
ax.text(cx1, 3.22, "d ≈ 150 nm · ζ ≈ −11 mV (pH 7)", ha='center', fontsize=8, color=GREY)

# ---- arrow 1 ----
ax.add_patch(FancyArrowPatch((4.55,CY),(5.55,CY), arrowstyle='-|>', mutation_scale=17, lw=2.0, color=RED))
ax.text(5.05, CY+0.42, "+ dairy proteins", ha='center', fontsize=8.5, color=RED, fontweight='bold')

# ---- Panel 2 ----
cx2 = 7.55
shell(cx2, CY, R); core(cx2, CY, R)
for ang,c in [(0,'#7A2E8E'),(52,'#1B7F5C'),(112,'#7A2E8E'),(168,'#C1121F'),(215,'#1B7F5C'),(268,'#C1121F'),(318,'#7A2E8E')]:
    th=np.radians(ang); rr=R+0.40
    ax.add_patch(Ellipse((cx2+rr*np.cos(th), CY+rr*np.sin(th)), 0.40, 0.28, angle=ang,
                         facecolor=c, edgecolor='white', lw=0.8, zorder=5))
ax.text(cx2, 3.70, "Protein corona", ha='center', fontsize=11, fontweight='bold', color=DARK)
ax.text(cx2, 3.22, "size ↑ 180–250 nm · |ζ| ↓ toward 0", ha='center', fontsize=8, color=GREY)
for i,(lab,c) in enumerate([("β-lactoglobulin",'#7A2E8E'),("α-lactalbumin",'#1B7F5C'),("caseins",'#C1121F')]):
    x=5.55+1.75*i
    ax.add_patch(Ellipse((x,0.22), 0.17, 0.125, facecolor=c, edgecolor='none'))
    ax.text(x+0.14, 0.22, lab, va='center', fontsize=8, color=DARK)

# ---- arrow 2 ----
ax.add_patch(FancyArrowPatch((10.05,CY),(11.05,CY), arrowstyle='-|>', mutation_scale=17, lw=2.0, color=RED))

# ---- Panel 3: food matrices as photos ----
import matplotlib.image as mpimg
PH="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures/photos"
ax.text(13.60, 3.70, "Real food matrices", ha='center', fontsize=11, fontweight='bold', color=DARK)
DIA = 1.45   # data units
for i,(name,cond,fn) in enumerate([("Milk","pH 6.6 · 80 mM","milk"),("Yogurt","pH 4.5 · 20 mM","yogurt"),("Bread","pH 5.3 · 50 mM","bread")]):
    x = 11.85+1.70*i
    img = mpimg.imread(f"{PH}/{fn}_circle.png")
    ax.imshow(img, extent=[x-DIA/2, x+DIA/2, 2.10-DIA/2, 2.10+DIA/2], zorder=5, interpolation='lanczos')
    ax.text(x, 1.15, name, ha='center', fontsize=9.5, fontweight='bold', color=DARK)
    ax.text(x, 0.85, cond, ha='center', fontsize=7.8, color=GREY)

plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
out="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures/np_schematic_slide3_wide.png"
plt.savefig(out, transparent=True, bbox_inches='tight', pad_inches=0.02)
print("saved", out)
