#!/usr/bin/env python3
"""
analyze_docking_results.py -- post-docking analysis for Section 3.1.3 (3B).

After you run the docking on HDOCK / ClusPro (web) and download the
results, this script does the two things you actually want out of docking
(see README Sections 1 and 5, Step 6):

  1. RANK the milk proteins by how strongly they bind the zein surface.
     Reads a scores CSV (scores_template.csv) that you fill in with the
     numbers the servers report. Produces a ranked table + a bar plot.
     REMEMBER (README Section 2b/2c): treat this as a RANKING, not an
     absolute affinity -- the nanoparticle avidity effect is not captured,
     so only the ORDER is trustworthy.

  2. MAP THE INTERFACE for each docked complex: which residues on the
     zein (receptor) side and the milk-protein (ligand) side are actually
     touching. This contact list is the mechanism hypothesis you hand to
     the coarse-grained MD step (3C). Reads each complex .pdb and lists
     residue pairs within a distance cutoff.

Usage:
  # ranking only:
  python analyze_docking_results.py --scores scores_template.csv

  # ranking + interface maps (needs the downloaded complex .pdb files):
  python analyze_docking_results.py --scores scores_template.csv \
      --complexes complexes/ --receptor-chain A --ligand-chain B

Score-direction note: different servers use different sign/scale
conventions (HDOCK: more NEGATIVE = stronger; some report positive
"confidence"). Tell the script which with --better {low,high}. Default is
'low' (HDOCK-style, lower is better).

Requires: pandas, matplotlib, biopython (interface maps only).
"""
import argparse
import glob
import os
import sys

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")


def rank_scores(scores_csv, better="low"):
    import pandas as pd

    df = pd.read_csv(scores_csv)
    df.columns = [c.strip() for c in df.columns]
    if "score" not in df.columns:
        print("  ERROR: scores CSV needs a 'score' column. Columns found: "
              f"{list(df.columns)}")
        sys.exit(1)

    df = df.dropna(subset=["score"]).copy()
    if df.empty:
        print("  No numeric scores filled in yet -- fill scores_template.csv "
              "with the HDOCK/ClusPro numbers, then re-run.")
        return None

    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score"])

    ascending = (better == "low")  # lower score is stronger for HDOCK
    # rank within each server so different score scales aren't mixed
    server_col = "server" if "server" in df.columns else None
    if server_col:
        df["rank_in_server"] = (
            df.groupby(server_col)["score"]
            .rank(ascending=ascending, method="min").astype(int))
        df = df.sort_values([server_col, "rank_in_server"])
    else:
        df["rank"] = df["score"].rank(ascending=ascending, method="min").astype(int)
        df = df.sort_values("rank")

    os.makedirs(OUTDIR, exist_ok=True)
    out_csv = os.path.join(OUTDIR, "docking_ranking.csv")
    df.to_csv(out_csv, index=False)
    print(f"  ranking table -> {out_csv}")
    print("\n" + df.to_string(index=False))

    # --- bar plot ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        label_col = "ligand_name" if "ligand_name" in df.columns else df.columns[0]
        fig, ax = plt.subplots(figsize=(8, 4.5))
        if server_col:
            for srv, g in df.groupby(server_col):
                ax.bar([f"{l}\n({srv})" for l in g[label_col]], g["score"], label=srv)
            ax.legend(title="server")
        else:
            ax.bar(df[label_col].astype(str), df["score"])
        ax.set_ylabel("docking score")
        ax.set_title("Docking score by milk protein (RANK only, not affinity)\n"
                     f"{'lower' if ascending else 'higher'} = stronger binding")
        ax.axhline(0, color="k", lw=0.6)
        plt.xticks(rotation=30, ha="right", fontsize=8)
        plt.tight_layout()
        out_png = os.path.join(OUTDIR, "docking_ranking.png")
        plt.savefig(out_png, dpi=150)
        print(f"  ranking plot  -> {out_png}")
    except Exception as e:
        print(f"  (plot skipped: {type(e).__name__}: {e})")

    return df


def interface_residues(complex_pdb, receptor_chain, ligand_chain, cutoff=5.0):
    """Return contacting residues between two chains (any-atom within cutoff)."""
    from Bio.PDB import PDBParser, NeighborSearch
    from Bio.PDB.Polypeptide import is_aa

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("c", complex_pdb)
    model = next(structure.get_models())

    chain_ids = [c.id for c in model.get_chains()]
    if receptor_chain not in chain_ids or ligand_chain not in chain_ids:
        return None, (f"chains {receptor_chain}/{ligand_chain} not both present "
                      f"(found {chain_ids})")

    rec_atoms = [a for a in model[receptor_chain].get_atoms()]
    lig_atoms = [a for a in model[ligand_chain].get_atoms()]
    ns = NeighborSearch(rec_atoms + lig_atoms)

    rec_res, lig_res = set(), set()
    for a, b in ns.search_all(cutoff, level="A"):
        ca, cb = a.get_parent().get_parent(), b.get_parent().get_parent()
        if {ca.id, cb.id} != {receptor_chain, ligand_chain}:
            continue  # only cross-chain contacts
        for atom in (a, b):
            res = atom.get_parent()
            ch = res.get_parent().id
            if not is_aa(res, standard=True):
                continue
            tag = (res.resname, res.id[1])
            (rec_res if ch == receptor_chain else lig_res).add(tag)

    return {"receptor": sorted(rec_res, key=lambda t: t[1]),
            "ligand": sorted(lig_res, key=lambda t: t[1])}, None


def analyze_complexes(complex_dir, receptor_chain, ligand_chain, cutoff):
    import pandas as pd

    files = sorted(glob.glob(os.path.join(complex_dir, "*.pdb")))
    if not files:
        print(f"  no .pdb files found in {complex_dir}")
        return
    rows = []
    for f in files:
        result, err = interface_residues(f, receptor_chain, ligand_chain, cutoff)
        name = os.path.basename(f)
        if err:
            print(f"  {name}: SKIPPED ({err})")
            continue
        rec = "; ".join(f"{r}{n}" for r, n in result["receptor"])
        lig = "; ".join(f"{r}{n}" for r, n in result["ligand"])
        print(f"\n  {name}")
        print(f"    zein(receptor) interface  : {rec or '(none)'}")
        print(f"    milk-protein(ligand) iface: {lig or '(none)'}")
        rows.append({"complex": name,
                     "n_receptor_contacts": len(result["receptor"]),
                     "n_ligand_contacts": len(result["ligand"]),
                     "receptor_interface_residues": rec,
                     "ligand_interface_residues": lig})
    if rows:
        os.makedirs(OUTDIR, exist_ok=True)
        out = os.path.join(OUTDIR, "interface_residues.csv")
        pd.DataFrame(rows).to_csv(out, index=False)
        print(f"\n  interface table -> {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", help="CSV of docking scores (see scores_template.csv)")
    ap.add_argument("--better", choices=["low", "high"], default="low",
                    help="does a LOWER or HIGHER score mean stronger binding? "
                         "(HDOCK=low, default)")
    ap.add_argument("--complexes", help="directory of downloaded docked complex .pdb files")
    ap.add_argument("--receptor-chain", default="A", help="zein/surface chain id (default A)")
    ap.add_argument("--ligand-chain", default="B", help="milk-protein chain id (default B)")
    ap.add_argument("--cutoff", type=float, default=5.0,
                    help="contact distance cutoff in Angstrom (default 5.0)")
    args = ap.parse_args()

    if not args.scores and not args.complexes:
        ap.print_help()
        return

    if args.scores:
        print("== RANKING ==")
        rank_scores(args.scores, better=args.better)
    if args.complexes:
        print("\n== INTERFACE RESIDUES ==")
        analyze_complexes(args.complexes, args.receptor_chain,
                          args.ligand_chain, args.cutoff)


if __name__ == "__main__":
    main()
