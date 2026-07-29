#!/usr/bin/env python3
"""
prepare_structures.py -- structure prep helper for Section 3.1.3 docking (3B).

Two jobs, both boring-but-important bookkeeping around the actual docking
(which happens on the HDOCK / ClusPro web servers):

  1. DOWNLOAD + CLEAN a PDB structure so the docking server sees exactly
     one clean protein chain (no water, no salts, no duplicate copies,
     no bound ligands). Raw PDB files are messy and servers can dock to
     the wrong thing.

        python prepare_structures.py --pdb-id 1BEB --keep-chain A
        python prepare_structures.py --pdb-id 1F6S --keep-chain A

  2. SURFACE HINT for the nanoparticle-aware binding-site restriction
     (README Step 4 / Section 2c). Given your predicted zein model, it
     estimates which residues face OUTWARD (water-liking + solvent
     exposed) -- i.e. the ones a milk protein could actually reach on a
     real NP surface. Paste these residue numbers into the HDOCK
     binding-site box so docking doesn't bury the ligand in zein's core.

        python prepare_structures.py --local zein_model.pdb --surface-hint

Requires: biopython  (pip install -r requirements.txt)
Docking itself is NOT done here -- see README Section 5.
"""
import argparse
import os
import sys
import urllib.request

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")

# Kyte-Doolittle hydropathy. POSITIVE = hydrophobic (buried in NP core);
# NEGATIVE = hydrophilic (faces water / outward on the NP surface).
KD = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
    "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
    "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
    "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}


def _pdb_url(pdb_id):
    return f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"


def download_pdb(pdb_id, dest):
    url = _pdb_url(pdb_id)
    print(f"  downloading {url}")
    urllib.request.urlretrieve(url, dest)
    print(f"  saved raw file -> {dest}")


def load_structure(in_path):
    """Parse a .pdb or .cif file, picking the right biopython parser."""
    from Bio.PDB import PDBParser, MMCIFParser

    ext = os.path.splitext(in_path)[1].lower()
    if ext == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("s", in_path)


def clean_structure(in_path, out_path, keep_chain):
    """Keep one chain, drop waters/hetero-atoms/ligands, keep protein only.

    Works for both .pdb and .cif input; always writes a clean .pdb out
    (so a downloaded .cif is converted to docking-ready .pdb in one step).
    """
    from Bio.PDB import PDBIO, Select

    class CleanSelect(Select):
        def accept_chain(self, chain):
            return keep_chain is None or chain.id == keep_chain

        def accept_residue(self, residue):
            # hetflag " " => standard amino acid; "W"/"H_..." => water/ligand
            return residue.id[0] == " "

    structure = load_structure(in_path)
    # use first model only
    model = next(structure.get_models())
    chains = [c.id for c in model.get_chains()]
    if keep_chain is not None and keep_chain not in chains:
        print(f"  WARNING: chain '{keep_chain}' not found. Chains present: "
              f"{chains}. Keeping all chains instead.")
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, select=CleanSelect())

    # count kept residues for a sanity message
    n = 0
    for c in model.get_chains():
        if keep_chain is None or c.id == keep_chain:
            n += sum(1 for r in c if r.id[0] == " ")
    print(f"  cleaned -> {out_path}  ({n} amino-acid residues kept, "
          f"chain={keep_chain or 'ALL'})")
    return out_path


def surface_hint(pdb_path, sasa_cutoff=0.20, kd_cutoff=0.0, top=None):
    """Print residues likely exposed on the NP surface.

    A residue is flagged 'surface' if it is BOTH:
      * hydrophilic  (Kyte-Doolittle < kd_cutoff), and
      * solvent-exposed (relative SASA > sasa_cutoff), when SASA is
        available; if the biopython SASA module isn't present we fall
        back to hydrophilicity alone and say so.
    """
    from Bio.PDB.Polypeptide import is_aa

    structure = load_structure(pdb_path)
    model = next(structure.get_models())

    # --- try to compute per-atom SASA, aggregate to residues ---
    have_sasa = False
    try:
        from Bio.PDB.SASA import ShrakeRupley
        ShrakeRupley().compute(model, level="R")
        have_sasa = True
    except Exception as e:  # module missing or failed -> degrade gracefully
        print(f"  (note: SASA not computed [{type(e).__name__}]; "
              f"using hydrophilicity only)")

    # rough max SASA per residue (A^2) for relative accessibility
    MAXASA = {
        "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
        "GLN": 225, "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197,
        "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
        "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174,
    }

    rows = []
    for chain in model.get_chains():
        for res in chain:
            if not is_aa(res, standard=True):
                continue
            rn = res.resname
            kd = KD.get(rn, 0.0)
            rel = None
            if have_sasa and hasattr(res, "sasa"):
                rel = res.sasa / MAXASA.get(rn, 200.0)
            hydrophilic = kd < kd_cutoff
            exposed = True if rel is None else (rel > sasa_cutoff)
            is_surface = hydrophilic and exposed
            rows.append({
                "chain": chain.id,
                "resnum": res.id[1],
                "resname": rn,
                "kd": kd,
                "rel_sasa": rel,
                "surface": is_surface,
            })

    surf = [r for r in rows if r["surface"]]
    if top:
        # rank by most exposed then most hydrophilic
        surf.sort(key=lambda r: (-(r["rel_sasa"] or 0), r["kd"]))
        surf = surf[:top]

    print(f"\n  Surface-exposed (outward-facing) residues in {pdb_path}:")
    print(f"  criteria: hydrophilic (KD<{kd_cutoff}) "
          + (f"AND rel.SASA>{sasa_cutoff}" if have_sasa
             else "(SASA unavailable -> hydrophilicity only)"))
    if not surf:
        print("  (none flagged -- loosen cutoffs or check the model)")
        return surf
    ids = ",".join(f"{r['chain']}{r['resnum']}" for r in surf)
    print(f"\n  --> paste these into the HDOCK binding-site box:\n\n  {ids}\n")
    print(f"  {len(surf)} residues flagged. "
          f"Full per-residue table would list chain/resnum/resname/KD/SASA.")
    return surf


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb-id", help="RCSB PDB ID to download + clean (e.g. 1BEB)")
    ap.add_argument("--local", help="path to a local .pdb (e.g. your predicted zein_model.pdb)")
    ap.add_argument("--keep-chain", default="A",
                    help="chain to keep when cleaning (default A; use ALL to keep every chain)")
    ap.add_argument("--surface-hint", action="store_true",
                    help="print likely surface-exposed residues for the NP binding-site restriction")
    ap.add_argument("--top", type=int, default=None,
                    help="limit surface hint to the N most-exposed residues")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    keep_chain = None if str(args.keep_chain).upper() == "ALL" else args.keep_chain

    if args.pdb_id:
        raw = os.path.join(OUTDIR, f"{args.pdb_id.upper()}_raw.pdb")
        clean = os.path.join(OUTDIR, f"{args.pdb_id.upper()}_clean.pdb")
        try:
            download_pdb(args.pdb_id, raw)
        except Exception as e:
            print(f"  ERROR downloading {args.pdb_id}: {e}")
            print("  (no internet in this environment? download manually from "
                  "rcsb.org and re-run with --local)")
            sys.exit(1)
        clean_structure(raw, clean, keep_chain)
        if args.surface_hint:
            surface_hint(clean, top=args.top)

    elif args.local:
        if not os.path.exists(args.local):
            print(f"  ERROR: file not found: {args.local}")
            sys.exit(1)
        if args.surface_hint:
            surface_hint(args.local, top=args.top)
        else:
            clean = os.path.join(
                OUTDIR,
                os.path.splitext(os.path.basename(args.local))[0] + "_clean.pdb")
            clean_structure(args.local, clean, keep_chain)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
