import pymol, os
from pymol import cmd
pymol.finish_launching(['pymol','-qc'])
PS="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/scripts_and_methods/Molecular_docking/protein_structures/zein_alpha_22597"
OUT="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures/structures"
cmd.bg_color("white"); cmd.set("ray_opaque_background",0); cmd.set("ray_shadows",0); cmd.set("antialias",2)

# ---- five ColabFold models superimposed ----
ranks=["001_alphafold2_ptm_model_3","002_alphafold2_ptm_model_4","003_alphafold2_ptm_model_2",
       "004_alphafold2_ptm_model_5","005_alphafold2_ptm_model_1"]
cols=["0xE21833","0x2E6DA4","0x2E7D32","0xE8A33D","0x7A2E8E"]
cmd.delete("all")
for i,r in enumerate(ranks):
    f=f"{PS}/zein_alpha_22597_unrelaxed_rank_{r}_seed_000.pdb"
    cmd.load(f, f"m{i}")
rms=[]
for i in range(1,5):
    r=cmd.align(f"m{i}", "m0", cycles=0)
    rms.append(round(r[0],1))
print("pairwise RMSD to model 1:", rms)
cmd.hide("everything"); cmd.show("cartoon")
for i,c in enumerate(cols): cmd.color(c, f"m{i}")
cmd.set("cartoon_transparency", 0.15)
cmd.orient("m0"); cmd.zoom("all", 2.0)
cmd.png(f"{OUT}/zein_5models.png", width=1200, height=1000, dpi=300, ray=1)
print("wrote zein_5models.png")

# ---- helix 84-115 vs ideal poly-Ala helix ----
cmd.delete("all")
cmd.load(f"{PS}/zein_alpha_22597_relaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb", "zein")
cmd.create("helix", "zein and resi 84-115 and polymer")
cmd.delete("zein")
cmd.fab("A"*32, "ideal", ss=1)
r=cmd.align("ideal////CA", "helix////CA", cycles=0)
print("ideal-helix RMSD:", round(r[0],2), "over", r[1], "atoms")
cmd.hide("everything"); cmd.show("cartoon")
cmd.color("0xE21833", "helix"); cmd.color("0x9E9E9E", "ideal")
cmd.set("cartoon_transparency", 0.0)
cmd.orient("helix"); cmd.turn("z", 90); cmd.zoom("all", 2.5)
cmd.png(f"{OUT}/zein_helix_vs_ideal.png", width=1400, height=700, dpi=300, ray=1)
print("wrote zein_helix_vs_ideal.png")
