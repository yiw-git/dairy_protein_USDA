import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

RED="#E21833"; DARK="#282828"; GREY="#636363"; BLUE="#3B6FA8"
fig = plt.figure(figsize=(9.0,2.45), dpi=400)
bg = fig.add_axes([0,0,1,1]); bg.set_xlim(0,16); bg.set_ylim(0,4.36); bg.axis('off')

# ---------------- left: equations ----------------
bg.text(0.05, 4.00, "Two competing forces", fontsize=10.5, fontweight='bold', color=DARK, va='top')
bg.text(0.05, 3.35, r"$V(D)\;=\;V_{EDL}(D)\;+\;V_{vdW}(D)$", fontsize=11, color=DARK, va='center')
bg.text(0.05, 2.90, "electrostatic double-layer repulsion", fontsize=7.5, color=BLUE, va='center', style='italic')
bg.text(0.05, 2.50, r"$V_{EDL}=2\pi\varepsilon_0\varepsilon_r\,a\,\zeta^{2}\,\ln\!\left(1+e^{-\kappa D}\right)$",
        fontsize=9.5, color=DARK, va='center')
bg.text(0.05, 2.02, "van der Waals attraction", fontsize=7.5, color=RED, va='center', style='italic')
bg.text(0.05, 1.62, r"$V_{vdW}=-\frac{A}{6}\left[\frac{2a^{2}}{D(D+4a)}+\frac{2a^{2}}{(D+2a)^{2}}+\ln\frac{D(D+4a)}{(D+2a)^{2}}\right]$",
        fontsize=8, color=DARK, va='center')
bg.text(0.05, 0.95, "inputs we already measure", fontsize=8.5, fontweight='bold', color=DARK, va='center')
for i,t in enumerate([r"$a$, $\zeta$  —  DLS / ELS",
                      r"$\kappa^{-1}$  —  Debye length from $I$ (0.79–3.04 nm at 150–10 mM)",
                      r"$A$  —  literature, $1\times10^{-20}$ J (zein)"]):
    bg.text(0.30, 0.60-0.28*i, "· "+t, fontsize=8, color=GREY, va='center')

# ---------------- middle: schematic V(D) ----------------
kB=1.380649e-23; T=298.15; e0=8.8541878128e-12; er=78.5; NA=6.02214076e23; qe=1.602176634e-19
def kappa(I_M): return np.sqrt(2*NA*qe**2*(I_M*1000.0)/(e0*er*kB*T))  # I: mol/L -> mol/m^3
def V(D_nm, a_nm, zeta_mV, I_M, A=1e-20):
    D=D_nm*1e-9; a=a_nm*1e-9; z=zeta_mV*1e-3; k=kappa(I_M)
    edl=2*np.pi*e0*er*a*z**2*np.log(1+np.exp(-k*D))
    vdw=-(A/6)*(2*a**2/(D*(D+4*a))+2*a**2/(D+2*a)**2+np.log(D*(D+4*a)/(D+2*a)**2))
    return edl/(kB*T), vdw/(kB*T)

axc = fig.add_axes([0.415, 0.20, 0.245, 0.66])
D=np.linspace(0.15, 8, 3000)
edl, vdw = V(D, 100, 45, 0.050)
tot = edl+vdw
axc.axhline(0, color='#BBBBBB', lw=0.8)
axc.plot(D, edl, lw=1.0, ls='--', color=BLUE, label=r'$V_{EDL}$')
axc.plot(D, vdw, lw=1.0, ls='--', color=RED, label=r'$V_{vdW}$')
axc.plot(D, tot, lw=2.0, color=DARK, label=r'$V_{total}$')
i=np.argmax(tot); axc.plot(D[i], tot[i], 'o', ms=4, color=DARK)
axc.annotate(r"$V_{max}$", xy=(D[i], tot[i]), xytext=(D[i]+2.2, tot[i]+12),
             fontsize=8.5, color=DARK, arrowprops=dict(arrowstyle='->', color=DARK, lw=0.9))
axc.axhline(15, color='#2E7D32', ls=':', lw=1.0); axc.axhline(10, color='#E9A200', ls=':', lw=1.0)
axc.text(7.85, 22, "15 kT stable / 10 kT marginal", fontsize=6.3, color=GREY, ha='right', va='bottom')

axc.set_xlim(0,8); axc.set_ylim(-60, float(np.nanmax(tot))*1.45)
axc.set_xlabel("separation $D$ (nm)", fontsize=8, labelpad=1.5)
axc.set_ylabel("$V$ ($k_BT$)", fontsize=8, labelpad=1.5)
axc.tick_params(labelsize=7, length=2.5, pad=1.5)
for sp in ['top','right']: axc.spines[sp].set_visible(False)
axc.legend(fontsize=6.5, frameon=False, loc='upper right', handlelength=1.4, borderaxespad=0.2)
bg.text(9.15, 4.00, "Barrier height sets stability", fontsize=10.5, fontweight='bold',
        color=DARK, va='top', ha='right')
bg.text(8.55, 0.02, "illustrative curve: a = 100 nm, ζ = 45 mV, I = 50 mM",
        fontsize=7, color=GREY, ha='center', va='bottom')

# ---------------- right: vial photo placeholder ----------------
bg.add_patch(FancyBboxPatch((11.25, 0.62), 4.55, 2.85,
             boxstyle="round,pad=0.0,rounding_size=0.14",
             facecolor="#DCE7F2", edgecolor=BLUE, lw=1.4, ls='--'))
bg.text(13.52, 2.20, "vial photos", ha='center', va='center', fontsize=10.5,
        color=BLUE, fontweight='bold')
bg.text(13.52, 1.78, "(placeholder — to be added)", ha='center', va='center', fontsize=8, color=BLUE)
bg.text(13.52, 4.00, "Measured grid: pH 4–7 × 10–150 mM", fontsize=10.5, fontweight='bold',
        color=DARK, va='top', ha='center')
bg.text(13.52, 0.30, "9 conditions · 3 replicates", ha='center', fontsize=8, color=GREY)

out="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures/dlvo_method_slide5.png"
plt.savefig(out, transparent=True, bbox_inches='tight', pad_inches=0.02)
print("saved")
