import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Circle, FancyBboxPatch

HYDC="#E8A33D"; POLC="#C9C9C9"; CHGC="#2E6DA4"; DARK="#282828"; GREY="#636363"; RED="#E21833"
F="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures"
comp={d['name']:d for d in json.load(open('/tmp/comp.json'))}

fig=plt.figure(figsize=(9.0,2.60), dpi=400)
bg=fig.add_axes([0,0,1,1]); bg.set_xlim(0,16); bg.set_ylim(0,4.62); bg.axis('off')

W=16.0; H=4.62
def fx(x): return x/W
def fy(y): return y/H

cols=[2.05, 6.05, 10.05, 14.05]
IMG_H=1.55
IMG_CY=3.52

def put_img(path, cx, max_w=3.5):
    im=mpimg.imread(path)
    ar=im.shape[1]/im.shape[0]
    h=IMG_H; w=h*ar
    if w>max_w: w=max_w; h=w/ar
    bg.imshow(im, extent=[cx-w/2, cx+w/2, IMG_CY-h/2, IMG_CY+h/2], zorder=3, interpolation='lanczos')

put_img(f"{F}/structures/beta_lactoglobulin.png", cols[0])
put_img(f"{F}/structures/alpha_lactalbumin.png", cols[1])
put_img(f"{F}/structures/alpha_zein.png", cols[3])

# --- casein micelle schematic ---
cx=cols[2]
rng=np.random.default_rng(7)
bg.add_patch(Circle((cx, IMG_CY), 0.62, facecolor="#F0F0F0", edgecolor="#CCCCCC", lw=1.0, zorder=2))
for _ in range(26):                                   # disordered casein chains inside
    a=rng.random()*2*np.pi; d=rng.random()**0.5*0.48
    x0,y0=cx+d*np.cos(a), IMG_CY+d*np.sin(a)
    t=np.linspace(0,1,20); ang=rng.random()*2*np.pi
    xs=x0+0.15*t*np.cos(ang)+0.030*np.sin(9*np.pi*t)*np.sin(ang)
    ys=y0+0.15*t*np.sin(ang)-0.030*np.sin(9*np.pi*t)*np.cos(ang)
    bg.plot(xs, ys, color=HYDC, lw=1.0, zorder=3, solid_capstyle='round')
for _ in range(11):                                   # CCP nanoclusters
    a=rng.random()*2*np.pi; d=0.15+rng.random()*0.36
    bg.add_patch(Circle((cx+d*np.cos(a), IMG_CY+d*np.sin(a)), 0.062,
                        facecolor="#8899A6", edgecolor='none', zorder=4))
for i in range(30):                                   # kappa-casein brush
    th=2*np.pi*i/30; t=np.linspace(0,1,26)
    rr=0.62+0.24*t; wob=0.045*np.sin(5*np.pi*t)
    bg.plot(cx+rr*np.cos(th)-wob*np.sin(th), IMG_CY+rr*np.sin(th)+wob*np.cos(th),
            color=CHGC, lw=1.0, zorder=2, solid_capstyle='round')
bg.annotate("κ-casein\nbrush", xy=(cx+0.80, IMG_CY+0.32), xytext=(cx+1.15, IMG_CY+0.75),
            fontsize=6.8, color=CHGC, ha='left', va='center',
            arrowprops=dict(arrowstyle='-', color=CHGC, lw=0.8))

# --- captions ---
caps=[("β-lactoglobulin","18 kDa · ~50% of whey\ncompact β-barrel, soluble"),
      ("α-lactalbumin","14 kDa · ~20% of whey\ncompact α/β fold, soluble"),
      ("Caseins (αs1, αs2, β, κ)","~80% of milk protein · no fixed fold\nassembled into ~150 nm micelles"),
      ("α-zein","~22 kDa · maize prolamin\nwater-insoluble, helical")]
for cx,(n,d) in zip(cols, caps):
    bg.text(cx, 2.52, n, fontsize=9.3, fontweight='bold', color=DARK, ha='center', va='top')
    bg.text(cx, 2.20, d, fontsize=7.4, color=GREY, ha='center', va='top')

# --- composition bars ---
ax=fig.add_axes([fx(1.45), fy(0.28), fx(13.2), fy(1.15)])
names=["β-lactoglobulin","α-lactalbumin","αs1-casein","β-casein","κ-casein","α-zein"]
short=["β-lg","α-la","αs1-casein","β-casein","κ-casein","α-zein"]
hyd=[comp[n]['hyd'] for n in names]; pol=[comp[n]['pol'] for n in names]; chg=[comp[n]['chg'] for n in names]
y=np.arange(len(names))[::-1]
ax.barh(y, hyd, color=HYDC, height=0.66, label="hydrophobic")
ax.barh(y, pol, left=hyd, color=POLC, height=0.66, label="polar")
ax.barh(y, chg, left=np.array(hyd)+np.array(pol), color=CHGC, height=0.66, label="charged")
for yy,h,p,c in zip(y,hyd,pol,chg):
    ax.text(min(h+p+c-0.6, 99.4), yy, f"{c:.0f}%", fontsize=6.6, color='white',
            va='center', ha='right', fontweight='bold')
ax.set_yticks(y); ax.set_yticklabels(short, fontsize=6.4)
ax.set_xlim(0,100); ax.set_xticks([0,25,50,75,100])
ax.set_xticklabels(["0","25","50","75","100%"], fontsize=6.5)
ax.tick_params(length=2, pad=1.5)
for sp in ['top','right','left']: ax.spines[sp].set_visible(False)
# label the three classes inside the top bar
ytop=y[0]
ax.text(hyd[0]/2, ytop, "hydrophobic", fontsize=6.4, color='white', ha='center', va='center', fontweight='bold')
ax.text(hyd[0]+pol[0]/2, ytop, "polar", fontsize=6.4, color='#555555', ha='center', va='center', fontweight='bold')
ax.text(hyd[0]+pol[0]+chg[0]/2, ytop, "charged", fontsize=6.4, color='white', ha='center', va='center', fontweight='bold')
ax.text(93.0, y[-1], "only 2% charged", fontsize=7.0, color=RED, ha='right', va='center', fontweight='bold')

out=f"{F}/proteins_slide7.png"
plt.savefig(out, transparent=True, bbox_inches='tight', pad_inches=0.02)
print("saved")
