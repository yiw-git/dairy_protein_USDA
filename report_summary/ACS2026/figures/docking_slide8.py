import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.image as mpimg

RED="#E21833"; DARK="#282828"; GREY="#636363"
F="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures"
D="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/scripts_and_methods/Molecular_docking/docking_trials/Megadocking/Megadock_result/combined_casein_analysis"
df=pd.read_csv(f"{D}/all5_zein_partners_summary.csv")
nice={"alpha_whey":"α-lactalbumin","beta_whey":"β-lactoglobulin","alphaS1_casein":"αs1-casein",
      "beta_casein":"β-casein","kappa_casein":"κ-casein"}
df["label"]=df.ligand.map(nice)
df=df.sort_values("ppi_e_score")

W,H=16.0,4.18
fig=plt.figure(figsize=(9.0,2.35), dpi=400)
bg=fig.add_axes([0,0,1,1]); bg.set_xlim(0,W); bg.set_ylim(0,H); bg.axis('off')
fx=lambda x: x/W; fy=lambda y: y/H

# ---- (a) scores ----
ax=fig.add_axes([fx(2.35), fy(1.30), fx(3.05), fy(2.10)])
y=np.arange(len(df))
ax.barh(y, df.ppi_e_score, color=RED, height=0.62)
for yy,v in zip(y, df.ppi_e_score):
    ax.text(v+0.12, yy, f"{v:.2f}", fontsize=6.4, color=DARK, va='center')
ax.set_yticks(y); ax.set_yticklabels(df.label, fontsize=6.8)
ax.set_xlim(0, 7.6); ax.set_xticks([0,2,4,6])
ax.tick_params(labelsize=6.5, length=2, pad=1.5)
ax.set_xlabel("docking score", fontsize=7.5, labelpad=1.5)
for sp in ['top','right','left']: ax.spines[sp].set_visible(False)
bg.text(3.85, 3.95, "No protein wins", fontsize=9.5, fontweight='bold', color=DARK, ha='center')
bg.text(3.85, 0.30, "spread 0.66 across five proteins;\nrun-to-run noise ±0.004", fontsize=7,
        color=GREY, ha='center', va='center')

# ---- (b) five models ----
im=mpimg.imread(f"{F}/structures/zein_5models.png")
ar=im.shape[1]/im.shape[0]; h=2.55; w=h*ar
bg.imshow(im, extent=[8.55-w/2, 8.55+w/2, 2.22-h/2, 2.22+h/2], zorder=3, interpolation='lanczos')
bg.text(8.55, 3.95, "The zein model is not reproducible", fontsize=9.5, fontweight='bold',
        color=DARK, ha='center')
bg.text(8.55, 0.30, "five independent predictions of the same sequence\nCα RMSD 13.8–21.9 Å",
        fontsize=7, color=GREY, ha='center', va='center')

# ---- (c) helix vs ideal ----
im2=mpimg.imread(f"{F}/zein_A_cartoon_plddt.png")
ar2=im2.shape[1]/im2.shape[0]; h2=2.55; w2=h2*ar2
bg.imshow(im2, extent=[13.55-w2/2, 13.55+w2/2, 2.22-h2/2, 2.22+h2/2], zorder=3, interpolation='lanczos')
bg.text(13.55, 3.95, "Only one region is confident", fontsize=9.5, fontweight='bold',
        color=DARK, ha='center')
bg.text(13.55, 0.30, "coloured by prediction confidence (red low → blue high)\nthe confident part, residues 84–115, is a plain α-helix",
        fontsize=7, color=GREY, ha='center', va='center')

out=f"{F}/docking_slide8.png"
plt.savefig(out, transparent=True, bbox_inches='tight', pad_inches=0.02)
print("saved")
