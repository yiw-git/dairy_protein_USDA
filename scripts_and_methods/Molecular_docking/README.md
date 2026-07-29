# Molecular Docking — Section 3.1.3 (Method 3B)

A plain-language, step-by-step guide to the molecular docking leg of the
Section 3.1.3 dry-lab strategy, written for someone new to docking.

It implements **"B. Molecular docking — corrected guidance"** from
*Dry-Lab Strategy for Section 3.1.3 — v3*, for the zein-caseinate
nanoparticle (NP) interacting with whey (β-lactoglobulin,
α-lactalbumin) and casein proteins.

**Before anything else, read Section 2 below ("The nanoparticle
concern").** It is the single most important thing to understand about
what docking can and cannot tell you for *your* system, and it changes
how you should present the results.

---

## 1. What molecular docking actually is (30-second version)

Docking is a computer method that takes the 3D shapes of two molecules
and asks: *if these two molecules stuck together, how would they fit,
and how tightly?* The computer tries thousands of ways of pressing the
two shapes together, scores each one (roughly, "how happy are the two
molecules in this arrangement"), and hands you back the best few poses
plus a score.

Two things it gives you that you care about:

- **A ranking.** Of β-lactoglobulin, α-lactalbumin, and casein, which
  one binds a zein surface most strongly? This tells you which milk
  protein is most likely to dominate the coating ("corona") on your NP.
- **An interface map.** *Which specific parts* (which amino acid
  residues) of each protein touch each other. This is a testable
  hypothesis about the mechanism, and it is the input the coarse-grained
  MD step (3C) and any later mutation/wet-lab work build on.

**What it is not:** it is not a wet-lab result and not a free-energy
measurement you can quote in kcal/mol as if it were experimental. The
score is a *relative* number for ranking, not an absolute binding
affinity. Treat docking as a **hypothesis generator**, not a
measurement.

### Why protein–protein docking, not AutoDock Vina

This is the correction flagged in the memo. There are two different
docking problems:

| Problem | Tool | What binds what |
|---|---|---|
| Small molecule into a protein pocket | **AutoDock Vina** | curcumin / resveratrol → ERα (this is what Objective 2 / Sec. 3.2.1 uses) |
| One whole protein onto another whole protein | **HDOCK / ClusPro** | whey/casein protein → zein surface (this is *your* Section 3.1.3 problem) |

Your problem is a protein sticking to another protein (corona
formation), so you want **HDOCK or ClusPro**, both free web servers.
Vina is the right tool for a *different* objective in the same grant, which
is why it appears in the project narrative — don't let that confuse you
into using it here.

---

## 2. THE NANOPARTICLE CONCERN (read this carefully)

Your instinct is correct and important: **docking was designed for two
free-floating, well-folded, dissolved proteins. You are not working with
dissolved zein — you are working with a zein nanoparticle**, which is a
solid-ish clump of many zein chains packed together, with only the
outer surface exposed to the milk proteins. This mismatch is real, and
being upfront about it is what makes your writeup credible instead of
naive.

Here is exactly where the mismatch bites, and what to do about each part.

### 2a. What the mismatch actually is

Standard protein–protein docking assumes:

1. **Two single molecules**, each a well-defined folded shape.
2. **The whole surface of each is available** to bind.
3. **A 1-to-1 interaction** — one protein meets one protein.
4. **Rigid-ish, water-soluble** partners.

Your real system violates several of these:

- Zein in your NP is **not one dissolved chain** — it's a packed
  aggregate. A milk protein can only touch the **outer surface patch**,
  not a whole free zein molecule.
- The NP surface is **curved and crowded**. Many zein chains sit side by
  side; the neighbors and the curvature change what's accessible.
- The real interaction is **many-to-one and multivalent** — a milk
  protein may contact several zein chains at once, and many milk
  proteins pile onto one particle. That "grip from multiple contacts at
  once" effect (called **avidity**) makes real binding stronger and
  stickier than any single one-to-one docking pose suggests.

So a raw "one zein chain + one β-lactoglobulin" docking run is a
**simplified proxy**, not a literal picture of your NP.

### 2b. Why docking is still worth doing anyway

It is still genuinely useful, for three honest reasons:

1. **The chemistry of a single contact is local.** When a milk protein
   touches the NP surface, the thing that decides whether it sticks *at
   that spot* is short-range chemistry — which residues meet, charge,
   hydrophobic patches. Docking one zein chain against one milk protein
   captures exactly that local contact chemistry, which is the part
   that transfers from "dissolved" to "surface."
2. **Ranking is more robust than absolute numbers.** Even if docking
   over- or under-estimates the true stickiness, the *order*
   (which milk protein binds zein best) tends to survive the
   simplification, because all three candidates are simplified the same
   way. And a ranking is exactly what you need for the corona story.
3. **It produces the interface hypothesis** that the coarse-grained MD
   step (3C) then tests under more realistic, crowded, surface-like
   conditions. Docking is the cheap, fast first pass (hours); CG-MD is
   the expensive, more-realistic second pass (days–weeks). Doing docking
   first tells CG-MD where to look.

**One-sentence framing for your memo/report:**
> "Docking is used as a fast, residue-level screen of which milk
> proteins bind the zein surface and through which contacts; it
> deliberately simplifies the nanoparticle to a single representative
> surface chain, and the multivalent, crowded, curved reality of the NP
> surface is handled downstream by the coarse-grained MD step. Docking
> ranks and hypothesizes; it does not measure."

### 2c. Concrete things to do so docking respects the NP reality

These are small choices that make your docking a fair proxy for a
*surface* rather than a free protein. Each is explained in plain terms.

1. **Dock milk protein → a single zein chain, and read it as
   "milk protein → zein surface patch," never "zein → zein."** You are
   modeling a milk protein landing on the NP surface. The zein chain is
   standing in for "the surface," so treat *zein as the big fixed thing*
   (the "receptor") and the *milk protein as the thing landing on it*
   (the "ligand"). HDOCK and ClusPro both let you label which is which.

2. **Restrict binding to the zein surface that would actually be
   exposed.** In a real NP, the water-hating (hydrophobic) parts of zein
   are buried in the particle core and the water-liking (hydrophilic)
   parts face outward. So the *outward-facing* residues are the ones a
   milk protein can reach. If you have any information about which face
   is exposed (from the CG-MD surface patch in 3C, or simply from which
   residues are hydrophilic), tell the docking server to focus the milk
   protein onto that face using its "binding site restriction" option.
   This prevents the unrealistic result where a milk protein docks into
   the buried core of zein — a pose that could never happen on a real
   particle. **This one restriction is the most important NP-aware
   step.** (See Section 5, step 4.)

3. **Don't quote the score as an affinity. Rank only.** Because the NP
   avidity effect is not captured, the *absolute* stickiness is wrong
   (too weak, usually). Only compare scores *between* the milk proteins,
   under the same settings. Report "β-lg > α-la > casein-fragment," not
   "β-lg binds at −X kcal/mol."

4. **For casein, dock fragments, not a folded blob** — but this is not
   an NP workaround, it's just casein being casein (next section).

5. **Hand the interface off to CG-MD to add the missing physics.** The
   things docking can't do — many proteins crowding one particle,
   multivalent grip, surface curvature — are exactly what the
   coarse-grained MD in Section 3C is for. So the docking deliverable
   isn't "the answer"; it's "the ranked shortlist + the contact map that
   tells CG-MD which pairs and which residues to simulate."

If you do those five things, you've turned a method built for dissolved
proteins into a defensible screen for a nanoparticle surface, *and*
you've documented the assumptions honestly — which is what a reviewer or
your mentor will actually want to see.

---

## 3. The players in your system (and their structure situation)

Docking needs a 3D structure file (a `.pdb` file) for each partner.
Here's where each one comes from and why.

| Molecule | Role here | 3D structure source | Why |
|---|---|---|---|
| **β-lactoglobulin** (β-lg) | major whey protein; likely main corona former | **PDB 1BEB** — download directly | Already solved by crystallography; nothing to predict. |
| **α-lactalbumin** (α-la) | second whey protein | **PDB 1F6S** — download directly | Also solved; download and clean. |
| **Zein** | the NP material (the "surface") | **Predict** from its UniProt sequence with **AlphaFold2/ColabFold** or **ESMFold** | Zein has *no* crystal structure — it's hard to crystallize. Prediction is standard and free. |
| **Casein** (α_s1, α_s2, β, κ) | major milk protein family | **Dock short surface-active fragments**, not a whole folded model | Casein is *intrinsically disordered* — it has no single stable fold, so "predict one structure and dock it" is scientifically wrong. Docking known surface-active peptide/phosphopeptide fragments is the standard, literature-accepted approach. |

Two things worth understanding here:

- **Why zein must be predicted:** crystallography needs the protein to
  pack into an orderly crystal. Zein is greasy and floppy at the ends
  and refuses to crystallize, so no one has an experimental structure.
  AlphaFold/ESMFold predict a plausible fold from the amino-acid
  sequence alone. Use the prediction, but remember it's a model — note
  the confidence score (pLDDT); low-confidence floppy regions are
  genuinely floppy, not errors to hide.
- **Why casein is different:** most proteins fold into one shape.
  Casein doesn't — it's a "wet noodle" that only takes shape when it
  touches something. So docking a single made-up casein fold would be
  modeling a shape that doesn't exist. Instead you dock the specific
  short stretches known to be sticky/surface-active (e.g. the
  phosphorylated regions of β- and α_s1-casein). This is why the memo
  says "representative peptide/phosphopeptide fragments."

---

## 4. The tools you'll use (all free, all web-based)

You do **not** install docking software. You upload structures to a web
server, it runs on their computers, and emails/streams you results in
a few hours. You are new to this, so start with the web servers.

| Tool | What it's for | Link | Notes |
|---|---|---|---|
| **RCSB PDB** | download β-lg (1BEB), α-la (1F6S) | rcsb.org | Free, no account. |
| **ColabFold** (AlphaFold2 in Google Colab) | predict the zein structure | the ColabFold notebook on Google Colab | Free; runs in your browser; needs the zein amino-acid sequence (from UniProt). |
| **ESMFold** | faster alternative for predicting zein | esmatlas.com/resources?action=fold | One-click; good for a quick first zein model. |
| **UniProt** | get the zein amino-acid sequence | uniprot.org | Search "zein" + *Zea mays*; copy the sequence. |
| **HDOCK** | the actual protein–protein docking | hdock.phys.hust.edu.cn | Upload two structures, get ranked poses + scores. Handles the whole job in one step; good default for beginners. |
| **ClusPro** | second docking server (cross-check) | cluspro.bu.edu | Free account required. Running the same pair on both servers and seeing if they agree is a cheap, powerful sanity check. |

The Python scripts in this folder do **not** do the docking (that
happens on the servers). They do the boring-but-important bookkeeping
around it: preparing/cleaning the structure files before you upload, and
analyzing/ranking the results after you download them. That's where
scripting actually saves you time and mistakes.

---

## 5. Step-by-step protocol

Reasons are given for every step so you understand *why*, not just
*what*.

### Step 1 — Get the whey structures (10 minutes)
Download **1BEB** (β-lg) and **1F6S** (α-la) from rcsb.org as `.pdb`.
*Why:* these are your two whey docking partners, already solved.

Then clean them (remove water molecules, extra copies, stray ligands) so
the docking server sees only the protein. Use:
```bash
python prepare_structures.py --pdb-id 1BEB --keep-chain A
python prepare_structures.py --pdb-id 1F6S --keep-chain A
```
*Why clean them:* raw PDB files contain water, salts, and sometimes
several copies of the protein. Docking servers can choke on these or
dock to the wrong thing. Keeping one clean protein chain avoids that.

### Step 2 — Get the zein sequence and predict its structure (1–3 hours, mostly waiting)
1. On UniProt, find a zein entry for *Zea mays* (maize) and copy the
   amino-acid sequence. (Zein is a family — α-zein is the most abundant;
   start there. Note which specific entry you used.)
2. Paste it into **ESMFold** (fast, one click) for a first model, and/or
   **ColabFold** (slower, usually higher quality) for a better one.
3. Download the predicted `.pdb`. **Write down the confidence (pLDDT).**
   *Why:* the predicted structure is a hypothesis. High-confidence
   regions are trustworthy; low-confidence regions are genuinely
   flexible/disordered and you should say so, not pretend they're solid.

### Step 3 — Prepare the casein fragments (30 minutes)
Instead of one casein structure, pick 2–3 short surface-active fragments
(the `docking_targets_template.csv` in this folder lists suggested
starting fragments and where they come from). You can build a simple
extended-peptide `.pdb` for each with ESMFold too (short peptides fold
fast).
*Why:* explained in Section 3 — casein has no single fold, so fragments
are the scientifically correct unit to dock.

### Step 4 — Decide which zein face is "the surface" (important NP step)
Before docking, look at your predicted zein model and identify the
**hydrophilic (water-liking) residues** — these are the ones that face
*outward* on a real NP and that a milk protein can actually reach. Note
their residue numbers. `prepare_structures.py --surface-hint` prints a
simple hydrophilicity-based list to get you started.
*Why:* this is the core fix for the nanoparticle concern (Section 2c,
point 2). It stops the docking server from burying the milk protein into
zein's greasy core — a pose that's impossible on a real particle.

### Step 5 — Run the docking (a few hours per pair, on the server)
For each milk protein (β-lg, α-la, each casein fragment):
1. Go to **HDOCK**.
2. Upload **zein as "Receptor 1"** (the big fixed surface) and the
   **milk protein as "Ligand 2"** (the thing landing on it).
3. In the optional binding-site box, enter the **exposed zein residues**
   from Step 4 to focus binding onto the realistic outer face.
4. Submit. Repeat on **ClusPro** for the same pairs as a cross-check.
*Why receptor/ligand labeling matters:* it encodes "milk protein lands
on zein surface," which is the physically correct direction for corona
formation (Section 2c, point 1).

### Step 6 — Analyze and rank the results (30 minutes)
Download each server's top-pose `.pdb` and its score. Then:
```bash
python analyze_docking_results.py --scores scores.csv \
    --complexes complexes/ --receptor-chain A --ligand-chain B
```
This produces a **ranking table** (which milk protein binds zein best)
and, for each complex, the **interface residues** (which residues on
each side are actually touching). Outputs land in `outputs/`.
*Why:* the ranking is your corona-order result; the interface residues
are the mechanism hypothesis you hand to CG-MD (3C) and, later, to any
wet-lab mutant/validation work.

### Step 7 — Sanity-check and write it up honestly
- Do HDOCK and ClusPro agree on the *order*? If yes, the ranking is
  robust. If they disagree, say so and lean on CG-MD to break the tie.
- Does the ranking make chemical sense (e.g. the most surface-active,
  most abundant protein winning)?
- State the assumptions from Section 2 explicitly in your writeup.

---

## 6. Files in this folder

| File | Purpose |
|---|---|
| `README.md` | This guide. |
| `docking_targets_template.csv` | The list of docking pairs to run (whey structures, zein prediction, casein fragments), with PDB/UniProt IDs, receptor/ligand role, and a status column to track your progress. Start here. |
| `prepare_structures.py` | Downloads a PDB by ID, strips waters/extra chains/ligands to leave one clean protein, and (with `--surface-hint`) prints the likely surface-exposed residues for the NP-aware binding-site restriction (Step 4). |
| `analyze_docking_results.py` | After docking: reads your scores CSV → builds a ranked table; reads the docked complex `.pdb` files → lists the interface (contacting) residues on each side. Saves CSVs + a ranking plot to `outputs/`. |
| `scores_template.csv` | Blank template for typing in the HDOCK/ClusPro scores as you get them. |
| `requirements.txt` | Python dependencies (only needed for the two helper scripts). |
| `outputs/` | Created automatically; holds the ranking table, interface tables, and plot. |

## 7. Quick start

```bash
pip install -r requirements.txt

# Prep the two whey structures (downloads + cleans):
python prepare_structures.py --pdb-id 1BEB --keep-chain A
python prepare_structures.py --pdb-id 1F6S --keep-chain A

# See the surface-exposed residues of your predicted zein model
# (feed these into the HDOCK binding-site box in Step 5):
python prepare_structures.py --local zein_model.pdb --surface-hint

# ...do the docking on HDOCK / ClusPro (web) ...

# Then rank the results and map the interfaces:
python analyze_docking_results.py --scores scores_template.csv \
    --complexes complexes/ --receptor-chain A --ligand-chain B
```

## 8. How this connects to the rest of Section 3.1.3

Docking is one of four legs, and it's deliberately the *cheap, fast,
hypothesis-generating* one. The four triangulate on the same question
(which milk proteins coat the zein NP, and how stable is it):

- **DLVO (3A):** physics of whether two coated particles repel or
  aggregate — the colloidal-stability angle.
- **Docking (3B, this folder):** which milk protein sticks to zein and
  through which residues — the *ranking + interface hypothesis*.
- **Coarse-grained MD (3C):** takes docking's interface hypothesis and
  tests it under realistic crowded, curved, multivalent surface
  conditions — this is where the "nanoparticle, not dissolved protein"
  physics actually lives.
- **ML / curcumin manuscript (Section 4):** literature-scale predictions
  of zeta potential across conditions.

**So the honest role of docking in the story:** it's the fast screen
that narrows "which pairs and which contacts are worth the expensive
CG-MD run," and provides a residue-level mechanism picture. It is not,
by itself, the load-bearing answer for a nanoparticle system — and your
writeup should say exactly that.

## 9. Key assumptions and limitations (state these when presenting)

- **Rigid / simplified partners.** Docking treats the shapes as mostly
  rigid; real proteins flex when they bind. HDOCK/ClusPro allow some
  give, but induced-fit changes are only approximated.
- **Single-chain proxy for the NP surface.** The biggest one — see
  Section 2. One zein chain stands in for a crowded, curved,
  multi-chain surface. Avidity/multivalency is *not* captured; absolute
  scores underestimate real stickiness. Use rankings, not absolute
  numbers.
- **Predicted (not experimental) zein structure.** The zein model is an
  AlphaFold/ESMFold prediction; low-pLDDT regions are uncertain. Report
  the confidence.
- **Casein represented by fragments.** Fragments capture local
  surface-active chemistry but not any larger-scale casein behavior.
- **Score is relative, not a measured affinity.** Docking scores rank;
  they are not kcal/mol you can compare to an experiment. A binding
  free energy would require the CG-MD umbrella-sampling / PMF step (3C)
  or experimental measurement.
- **Solvent, ions, and pH are implicit at best.** Web docking servers
  don't let you dial in "pH 4.5, 20 mM NaCl" the way your DLS matrix or
  CG-MD can. So docking gives a *pH-agnostic* interface hypothesis; the
  pH/ionic-strength dependence lives in DLVO (3A), CG-MD (3C), and the
  ML model (Section 4), not here. Don't over-claim condition-specific
  docking results.
- **SDS-PAGE is a later validation, not an input.** If you get gel data
  later, you can check whether docking's affinity *ranking* matches band
  intensities. It is a nice-to-have cross-check, never a prerequisite.

## 10. Key reference papers (see also Section 5 of the memo)

- Lee, H. *Recent Advances in Simulation Studies on the Protein Corona.*
  Pharmaceutics **2024**, 16(11), 1419. — Current review of docking +
  MD for corona prediction; the best orientation piece for a newcomer.
- Yan, Y.; Zhang, D.; Zhou, P.; Li, B.; Huang, S.-Y. *HDOCK: a web
  server for protein–protein and protein–DNA/RNA docking based on a
  hybrid strategy.* Nucleic Acids Res. **2017**, 45(W1), W365–W373. —
  The HDOCK method paper (cite when you use HDOCK).
- Kozakov, D.; Hall, D. R.; Xia, B.; et al. *The ClusPro web server for
  protein–protein docking.* Nat. Protoc. **2017**, 12(2), 255–278. —
  The ClusPro method paper (cite when you use ClusPro).
- Jumper, J.; et al. *Highly accurate protein structure prediction with
  AlphaFold.* Nature **2021**, 596, 583–589. — For the predicted zein
  structure (cite AlphaFold/ColabFold).
- Lin, Z.; et al. *Evolutionary-scale prediction of atomic-level protein
  structure with a language model.* Science **2023**, 379, 1123–1130. —
  The ESMFold paper (cite if you use ESMFold for zein/peptides).
- (Casein disorder context) Holt, C.; Carver, J. A.; Ecroyd, H.; Thorn,
  D. C. *Caseins and the casein micelle: their biological functions,
  structures, and behavior in foods.* J. Dairy Sci. **2013**, 96(10),
  6127–6146. — Supports the "casein is intrinsically disordered, dock
  fragments" choice.
