import pymol, os
from pymol import cmd
pymol.finish_launching(['pymol','-qc'])
PS="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/scripts_and_methods/Molecular_docking/protein_structures"
OUT="/sessions/wizardly-funny-sagan/mnt/dairy_protein_USDA/report_summary/ACS2026/figures/structures"
os.makedirs(OUT, exist_ok=True)
HYD="ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO+GLY"
POL="SER+THR+ASN+GLN+CYS+TYR"
CHG="ASP+GLU+LYS+ARG+HIS"
cmd.bg_color("white"); cmd.set("ray_opaque_background",0); cmd.set("ray_shadows",0)
cmd.set("antialias",2); cmd.set("surface_quality",1); cmd.set("specular",0.15)

def render(path, name, out, sele="polymer", chain=None):
    cmd.delete("all")
    cmd.load(path, name)
    cmd.remove("solvent or hetatm")
    if chain: cmd.remove(f"not chain {chain}")
    cmd.hide("everything")
    cmd.show("surface", name)
    cmd.color("0xE8A33D", f"{name} and resn {HYD}")
    cmd.color("0xC9C9C9", f"{name} and resn {POL}")
    cmd.color("0x2E6DA4", f"{name} and resn {CHG}")
    cmd.orient(name)
    cmd.zoom(name, 2.0)
    cmd.png(os.path.join(OUT,out), width=1100, height=1100, dpi=300, ray=1)
    print("wrote", out)

render(f"{PS}/whey/beta_whey_clean.pdb", "blg", "beta_lactoglobulin.png", chain="A")
render(f"{PS}/whey/alpha_whey_clean.pdb", "ala", "alpha_lactalbumin.png", chain="A")
render(f"{PS}/zein_alpha_22597/zein_alpha_22597_relaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb",
       "zein", "alpha_zein.png")
