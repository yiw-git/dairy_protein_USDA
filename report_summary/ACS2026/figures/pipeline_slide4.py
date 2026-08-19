import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

RED="#E21833"; DARK="#000000"; GREY="#000000"; LIGHT="#F2F2F2"
fig, ax = plt.subplots(figsize=(9.0,2.32), dpi=400)
ax.set_xlim(0,16); ax.set_ylim(0,4.12); ax.axis('off')

cards = [
 ("1","DLVO / extended DLVO","whole particle · ~100 nm","Does a stability\nbarrier exist?","V(D) barrier across\npH × ionic strength"),
 ("2","Molecular docking","residue level · Å","Which contacts\nform?","protein ranking +\ninterface hypothesis"),
 ("3","Coarse-grained MD","surface patch · ~10 nm, µs","Does the interface\nhold?","contact persistence\n+ steric term"),
]
CW, CH = 4.80, 2.76
CY0, CY1 = 0.98, 0.98+CH
centers = [2.62, 8.00, 13.38]

for (num,title,scale,q,out), cx in zip(cards, centers):
    ax.add_patch(FancyBboxPatch((cx-CW/2, CY0), CW, CH,
                 boxstyle="round,pad=0.0,rounding_size=0.18",
                 facecolor=LIGHT, edgecolor="#D8D8D8", lw=1.0, zorder=1))
    ax.add_patch(Circle((cx-CW/2+0.44, CY1-0.42), 0.255, facecolor=RED, edgecolor='none', zorder=3))
    ax.text(cx-CW/2+0.44, CY1-0.425, num, ha='center', va='center', fontsize=10.5,
            color='white', fontweight='bold', zorder=4)
    ax.text(cx-CW/2+0.82, CY1-0.42, title, ha='left', va='center', fontsize=12,
            fontweight='bold', color=DARK, zorder=4)
    ax.text(cx, CY1-0.92, scale, ha='center', va='center', fontsize=10.5, color=GREY, zorder=4)
    ax.text(cx, CY1-1.58, q, ha='center', va='center', fontsize=11.5, style='italic', color=DARK, zorder=4, linespacing=1.35)
    ax.text(cx, CY1-2.34, out, ha='center', va='center', fontsize=10.5, color=GREY, zorder=4, linespacing=1.35)

# forward arrows
for x0,x1 in [(5.10,5.52),(10.48,10.90)]:
    ax.add_patch(FancyArrowPatch((x0, CY0+CH/2),(x1, CY0+CH/2), arrowstyle='-|>',
                 mutation_scale=17, lw=2.4, color=RED, zorder=5))

# feedback arrow: CG-MD -> extended DLVO
ax.add_patch(FancyArrowPatch((13.38, CY0-0.04),(2.62, CY0-0.04),
             connectionstyle="arc3,rad=-0.10", arrowstyle='-|>', mutation_scale=14,
             lw=1.6, color=GREY, linestyle=(0,(5,3)), zorder=0))
ax.text(8.0, 0.04, "effective steric / adhesion term feeds back into extended DLVO",
        ha='center', va='bottom', fontsize=10.5, color=GREY, zorder=6)

plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
out="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures/pipeline_slide4.png"
plt.savefig(out, transparent=True, bbox_inches='tight', pad_inches=0.02)
print("saved")
