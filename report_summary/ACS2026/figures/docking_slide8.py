import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.image as mpimg

RED="#E21833"; DARK="#000000"; GREY="#000000"
F="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures"
D="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/scripts_and_methods/Molecular_docking/docking_trials/Megadocking/Megadock_result/combined_casein_analysis"
df=pd.read_csv(f"{D}/all5_zein_partners_summary.csv")
nice={"alpha_whey":"α-lactalbumin","beta_whey":"β-lactoglobulin","alphaS1_casein":"αs1-casein",
      "beta_casein":"β-casein","kappa_casein":"κ-casein"}
df["label"]=df.ligand.map(nice)
df=df.sort_values("ppi_e_score")

W,H=16.0,4.38
fig=plt.figure(figsize=(9.0,2.464), dpi=400)
bg=fig.add_axes([0,0,1,1]); bg.set_xlim(0,W); bg.set_ylim(0,H); bg.axis('off')
fx=lambda x: x/W; fy=lambda y: y/H

# ---- (a) scores ----
ax=fig.add_axes([fx(2.45), fy(1.52), fx(3.15), fy(1.88)])
y=np.arange(len(df))
ax.barh(y, df.ppi_e_score, color=RED, height=0.62)
for yy,v in zip(y, df.ppi_e_score):
    ax.text(v+0.14, yy, f"{v:.2f}", fontsize=8.5, color=DARK, va='center')
ax.set_yticks(y); ax.set_yticklabels(df.label, fontsize=9.0)
ax.set_xlim(0, 7.6); ax.set_xticks([0,2,4,6])
ax.tick_params(labelsize=8.5, length=2.5, pad=1.5)
ax.set_xlabel("docking score", fontsize=9.5, labelpad=2.0)
for sp in ['top','right','left']: ax.spines[sp].set_visible(False)
bg.text(3.85, 4.28, "No protein wins", fontsize=11.5, fontweight='bold', color=DARK, ha='center', va='top', linespacing=1.22)
bg.text(3.85, 0.58, "spread 0.66 across five proteins;\nrun-to-run noise ±0.004", fontsize=9.5,
        color=GREY, ha='center', va='center', linespacing=1.35)

# ---- (b) five models ----
im=mpimg.imread(f"{F}/structures/zein_5models.png")
ar=im.shape[1]/im.shape[0]; h=1.92; w=h*ar
bg.imshow(im, extent=[8.55-w/2, 8.55+w/2, 2.48-h/2, 2.48+h/2], zorder=3, interpolation='lanczos')
bg.text(8.55, 4.28, "The zein model\nis not reproducible", fontsize=11.5, fontweight='bold',
        color=DARK, ha='center', va='top', linespacing=1.22)
bg.text(8.55, 0.58, "five independent predictions\nof the same sequence\nCα RMSD 13.8–21.9 Å",
        fontsize=9.5, color=GREY, ha='center', va='center', linespacing=1.35)

# ---- (c) helix vs ideal ----
im2=mpimg.imread(f"{F}/zein_A_cartoon_plddt.png")
ar2=im2.shape[1]/im2.shape[0]; h2=1.92; w2=h2*ar2
bg.imshow(im2, extent=[13.55-w2/2, 13.55+w2/2, 2.48-h2/2, 2.48+h2/2], zorder=3, interpolation='lanczos')
bg.text(13.55, 4.28, "Only one region\nis confident", fontsize=11.5, fontweight='bold',
        color=DARK, ha='center', va='top', linespacing=1.22)
bg.text(13.55, 0.58, "coloured by prediction confidence\n(red low → blue high)\nconfident part = residues 84–115,\na plain α-helix",
        fontsize=9.5, color=GREY, ha='center', va='center', linespacing=1.35)

out=f"{F}/docking_slide8.png"
plt.savefig(out, transparent=True, bbox_inches='tight', pad_inches=0.02)
print("saved")
