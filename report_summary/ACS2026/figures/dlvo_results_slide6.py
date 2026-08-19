import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RED="#E21833"; DARK="#282828"; GREY="#636363"; GOLD="#E9A200"; GREEN="#2E7D32"
ROOT="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/scripts_and_methods/DLVO"
OUT=f"{ROOT}/outputs/real_data_08102026"
raw=pd.read_csv(f"{ROOT}/data_DLVO/DLS_data_zein_caseinate_08102026.csv")
cs=pd.read_csv(f"{OUT}/condition_summary.csv")
tg=pd.read_csv(f"{OUT}/target_condition_predictions.csv")

fig, axes = plt.subplots(1, 3, figsize=(9.0, 2.45), dpi=400,
                         gridspec_kw=dict(width_ratios=[1.0,0.85,1.25], wspace=0.42))
for ax in axes:
    for sp in ['top','right']: ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=7, length=2.5, pad=1.5)

# ---- (a) zeta vs pH with SD ----
ax=axes[0]
colors={10:'#8ECAE6',50:'#219EBC',150:'#023047'}
for I,g in cs.groupby('ionic_strength_mM'):
    g=g.sort_values('pH')
    ax.errorbar(g.pH, g.zeta_mV_mean, yerr=g.zeta_mV_sd, fmt='o-', ms=4, lw=1.1,
                capsize=2.0, elinewidth=0.9, color=colors[I], label=f"{int(I)} mM")
ax.axhline(0, color='#BBBBBB', lw=0.8)
ax.axvline(4.94, color=RED, ls='--', lw=1.2)
ax.text(5.03, 6.5, "IEP ≈ 4.94", fontsize=7, color=RED, va='center', ha='left')
ax.axvline(5.8, color=GREY, ls=':', lw=1.1)
ax.text(5.90, -41.4, "bare zein ≈ 5.8 (lit.)", fontsize=7, color=GREY, va='bottom', ha='left')
ax.set_xlabel("pH", fontsize=8, labelpad=1.5); ax.set_ylabel("ζ (mV)", fontsize=8, labelpad=1.5)
ax.set_xlim(3.5, 7.5); ax.set_ylim(-42, 14)
ax.legend(fontsize=6.5, frameon=False, loc='upper right', handlelength=1.2, borderaxespad=0.2)
ax.set_title("Caseinate shifts the IEP", fontsize=9, fontweight='bold', color=DARK, pad=4)

# ---- (b) radius vs ionic strength ----
ax=axes[1]
for pH,c,mk in [(4.0,RED,'o'),(5.5,'#7B4FA8','^'),(7.0,'#023047','s')]:
    g=cs[cs.pH==pH].sort_values('ionic_strength_mM')
    ax.errorbar(g.ionic_strength_mM, g.radius_nm_mean, yerr=g.radius_nm_sd, fmt=mk+'-',
                ms=4.5, lw=1.3, capsize=2.0, elinewidth=0.9, color=c, label=f"pH {pH:g}")
ax.set_yscale('log'); ax.set_xscale('log'); ax.set_xticks([10,50,150])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel("ionic strength (mM)", fontsize=8, labelpad=1.5)
ax.set_ylabel("z-average radius (nm)", fontsize=8, labelpad=1.5)
ax.annotate("8000 nm", xy=(150, 8000), xytext=(52, 3000), fontsize=7, color=RED,
            arrowprops=dict(arrowstyle='->', color=RED, lw=0.9))
ax.axvspan(50, 150, color=RED, alpha=0.07)
ax.text(87, 200, "CCC in\n50–150 mM", fontsize=7, color=RED, ha='center')
ax.set_ylim(50, 30000); ax.set_xlim(7.5, 210)
ax.legend(fontsize=6.5, frameon=False, loc='upper left', handlelength=1.2, borderaxespad=0.2)
ax.set_title("Aggregation onset at pH 4", fontsize=9, fontweight='bold', color=DARK, pad=4)

# ---- (c) V_max bars ----
ax=axes[2]
lab=[f"pH {r.pH:g}, {int(r.ionic_strength_mM)} mM" for r in cs.itertuples()]+list(tg.condition)
val=list(cs.V_max_kT)+list(tg.V_max_kT)
kind=['grid']*len(cs)+['food']*len(tg)
# key off the model's barrier_present flag; targets have no column, so a
# genuine interior maximum is equivalent to V_max > 0 for these conditions.
has=[bool(b) for b in cs.barrier_present]+[v > 0 for v in tg.V_max_kT]

# no-barrier conditions first (bottom), then real barriers ascending
o=sorted(range(len(val)), key=lambda i: (has[i], val[i] if has[i] else 0.0))
lab=[lab[i] for i in o]; val=[val[i] for i in o]
kind=[kind[i] for i in o]; has=[has[i] for i in o]

y=np.arange(len(val))
bars=[v if h else 0.0 for v,h in zip(val,has)]
ax.barh(y, bars, color=[(RED if k=='grid' else '#1F6FB2') for k in kind], height=0.72)
ax.axvline(0, color='#888888', lw=0.8)
ax.axvline(10, color=GOLD, ls='--', lw=1.1); ax.axvline(15, color=GREEN, ls='--', lw=1.1)
ax.text(10, len(val)-0.15, "10 kT", fontsize=6.5, color=GOLD, ha='center', va='bottom')
ax.text(15.6, len(val)-0.15, "15 kT", fontsize=6.5, color=GREEN, ha='left', va='bottom')

n_no = sum(1 for h in has if not h)
if n_no:
    # shade the no-barrier block and label it once, rather than plotting
    # boundary values that are artifacts of the D_max search window
    ax.axhspan(-0.62, n_no-0.38, color='#9AA0A6', alpha=0.10, zorder=0)
    for i in range(n_no):
        ax.plot(0, i, marker='x', ms=3.4, mew=1.0, color=GREY, zorder=3)
    ax.text(21.0, (n_no-1)/2.0, "no barrier — V(D) attractive at all separations",
            fontsize=6.3, color=GREY, va='center', ha='left', style='italic')

for i,(v,h) in enumerate(zip(val,has)):
    if h and v > 3:
        ax.text(v+1.0, i, f"{v:.1f}" if v < 10 else f"{v:.0f}",
                fontsize=6.3, color=DARK, va='center', ha='left')

ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=6.2)
ax.set_xlim(-1.5, 60); ax.set_ylim(-0.7, len(val)+0.3)
ax.set_xlabel("energy barrier $V_{max}$ ($k_BT$)", fontsize=8, labelpad=1.5)
ax.set_title("Stable only at low ionic strength", fontsize=9, fontweight='bold', color=DARK, pad=4)
ax.text(59, n_no+0.2, "red = measured grid\nblue = food condition", fontsize=6.3, color=GREY,
        ha='right', va='bottom')

plt.subplots_adjust(left=0.055, right=0.995, top=0.86, bottom=0.20)
out="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures/dlvo_results_slide6.png"
plt.savefig(out, transparent=True, bbox_inches='tight', pad_inches=0.02)
print("saved")
