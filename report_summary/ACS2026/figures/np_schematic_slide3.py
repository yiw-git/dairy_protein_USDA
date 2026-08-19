import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Ellipse, Polygon, Rectangle

RED="#E21833"; GOLD="#FFD200"; DARK="#282828"; GREY="#636363"; ZEIN="#E8A33D"; CAS="#4A7FB5"
fig, ax = plt.subplots(figsize=(12,4.9), dpi=300)
ax.set_xlim(0,12); ax.set_ylim(0,4.9); ax.axis('off')

def caseinate_shell(cx, cy, r, n=26, amp=0.09, color=CAS, lw=1.4):
    for i in range(n):
        th = 2*np.pi*i/n
        t = np.linspace(0,1,40)
        rr = r + 0.30*t
        wob = amp*np.sin(6*np.pi*t)
        x = cx + rr*np.cos(th) - wob*np.sin(th)
        y = cy + rr*np.sin(th) + wob*np.cos(th)
        ax.plot(x, y, color=color, lw=lw, solid_capstyle='round', zorder=2)

# ---------------- Panel 1: the particle ----------------
cx, cy, r = 1.95, 2.55, 0.85
caseinate_shell(cx, cy, r)
ax.add_patch(Circle((cx,cy), r, facecolor=ZEIN, edgecolor='none', zorder=3))
np.random.seed(3)
for _ in range(28):  # hint of packed zein chains in the core
    a=np.random.rand()*2*np.pi; d=np.random.rand()**0.5*(r-0.16)
    x0,y0=cx+d*np.cos(a), cy+d*np.sin(a); ang=np.random.rand()*2*np.pi
    ax.plot([x0-0.13*np.cos(ang), x0+0.13*np.cos(ang)],[y0-0.13*np.sin(ang), y0+0.13*np.sin(ang)],
            color='#C4832A', lw=2.0, solid_capstyle='round', zorder=4)
ax.text(cx, 4.55, "Zein–caseinate nanoparticle", ha='center', fontsize=11, fontweight='bold', color=DARK)
ax.annotate("packed zein chains\n(hydrophobic core)", xy=(cx-0.40,cy-0.30), xytext=(0.05,1.02),
            fontsize=8.5, color=DARK, ha='left',
            arrowprops=dict(arrowstyle='-', color=GREY, lw=0.9))
ax.annotate("sodium caseinate shell\n(charge + steric layer)", xy=(cx+0.72,cy+0.80), xytext=(2.65,4.02),
            fontsize=8.5, color=DARK, ha='left',
            arrowprops=dict(arrowstyle='-', color=GREY, lw=0.9))
ax.text(cx, 0.35, "d ≈ 150 nm   ζ ≈ −11 mV (pH 7)", ha='center', fontsize=8.5, color=GREY)

# ---------------- arrow 1 ----------------
ax.add_patch(FancyArrowPatch((4.05,2.55),(5.15,2.55), arrowstyle='-|>', mutation_scale=18, lw=2.0, color=RED))
ax.text(4.6, 2.80, "+ dairy\nproteins", ha='center', fontsize=8.5, color=RED, fontweight='bold')

# ---------------- Panel 2: corona ----------------
cx2 = 6.4
caseinate_shell(cx2, cy, r)
ax.add_patch(Circle((cx2,cy), r, facecolor=ZEIN, edgecolor='none', zorder=3))
corona = [(0,'#7A2E8E'),(48,'#1B7F5C'),(105,'#7A2E8E'),(160,'#C1121F'),(210,'#1B7F5C'),(262,'#C1121F'),(310,'#7A2E8E')]
for ang,c in corona:
    th=np.radians(ang); rr=r+0.42
    ax.add_patch(Ellipse((cx2+rr*np.cos(th), cy+rr*np.sin(th)), 0.42, 0.30,
                         angle=ang, facecolor=c, edgecolor='white', lw=0.8, zorder=5))
ax.text(cx2, 4.55, "Protein corona", ha='center', fontsize=11, fontweight='bold', color=DARK)
for i,(lab,c) in enumerate([("β-lactoglobulin",'#7A2E8E'),("α-lactalbumin",'#1B7F5C'),("casein fractions",'#C1121F')]):
    y=0.95-0.26*i
    ax.add_patch(Ellipse((5.30,y), 0.22, 0.16, facecolor=c, edgecolor='none'))
    ax.text(5.47, y, lab, va='center', fontsize=8.2, color=DARK)
ax.text(7.70, 4.02, "size ↑ 180–250 nm\n|ζ| ↓ toward 0", fontsize=8.5, color=GREY, ha='left')

# ---------------- arrow 2 ----------------
ax.add_patch(FancyArrowPatch((8.55,2.55),(9.25,2.55), arrowstyle='-|>', mutation_scale=18, lw=2.0, color=RED))

# ---------------- Panel 3: food matrices ----------------
ax.text(10.6, 4.55, "Real food matrices", ha='center', fontsize=11, fontweight='bold', color=DARK)
def milk(x,y,s=1.0):
    ax.add_patch(Rectangle((x-0.26*s,y-0.42*s), 0.52*s, 0.62*s, facecolor='#EDEDED', edgecolor=GREY, lw=1.1))
    ax.add_patch(Polygon([[x-0.26*s,y+0.20*s],[x+0.26*s,y+0.20*s],[x,y+0.50*s]], facecolor='#EDEDED', edgecolor=GREY, lw=1.1))
def yogurt(x,y,s=1.0):
    ax.add_patch(Polygon([[x-0.28*s,y+0.22*s],[x+0.28*s,y+0.22*s],[x+0.20*s,y-0.34*s],[x-0.20*s,y-0.34*s]],
                         facecolor='#F7F2E8', edgecolor=GREY, lw=1.1))
    ax.add_patch(Rectangle((x-0.30*s,y+0.20*s), 0.60*s, 0.09*s, facecolor=GOLD, edgecolor=GREY, lw=0.9))
def bread(x,y,s=1.0):
    ax.add_patch(FancyBboxPatch((x-0.28*s,y-0.34*s), 0.56*s, 0.50*s, boxstyle="round,pad=0.02,rounding_size=0.10",
                                facecolor='#E2B87A', edgecolor=GREY, lw=1.1))
    ax.add_patch(Ellipse((x,y+0.18*s), 0.56*s, 0.26*s, facecolor='#E2B87A', edgecolor=GREY, lw=1.1))
items=[("Milk","pH 6.6 · 80 mM", milk),("Yogurt","pH 4.5 · 20 mM", yogurt),("Bread","pH 5.3 · 50 mM", bread)]
for i,(name,cond,fn) in enumerate(items):
    y = 3.75-1.15*i
    fn(9.90, y, 0.90)
    ax.text(10.32, y+0.11, name, fontsize=9.5, color=DARK, fontweight='bold', va='center')
    ax.text(10.32, y-0.15, cond, fontsize=8.2, color=GREY, va='center')

plt.tight_layout(pad=0.2)
out="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures/np_schematic_slide3.png"
plt.savefig(out, transparent=True, bbox_inches='tight')
print("saved", out)
