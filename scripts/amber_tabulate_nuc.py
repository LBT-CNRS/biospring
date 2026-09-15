#!/usr/bin/env python3
"""The nucleic half of amber_tabulate: AMBER's torsion terms for DNA and RNA,
enumerated from the residue topology and tabulated the same way.

Same principle, same file format, different chemistry: a nucleotide chain is
linked O3'(i) - P(i+1) rather than C(i) - N(i+1), and the families split
backbone / glycosidic / sugar the way BioSpring's own enum does -- the four
furanose ring bonds get their own family because the pucker is the single lever
choosing the A or B helical form, so it has to be switchable on its own.

  python amber_tabulate_nuc.py DNA out.bi.ff
  python amber_tabulate_nuc.py RNA out.bi.ff
"""
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_nucleic_forcefield as N

KIND = (sys.argv[1] if len(sys.argv) > 1 else "DNA").upper()
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    REPO_ROOT, f"data/reducerules/{KIND}AtomTorsion.bi.ff")
BINS = int(os.environ.get("BIOSPRING_TORSION_BINS", "512"))
PHI = np.linspace(-np.pi, np.pi, BINS + 1)

XML = {"DNA": ("amber14", "DNA.OL15.xml"), "RNA": ("amber14", "RNA.OL3.xml")}[KIND]
RESIDUES = {"DNA": ["DA", "DC", "DG", "DT"], "RNA": ["A", "C", "G", "U"]}[KIND]
tables_ff = N.ForceFieldTables(os.path.join(N.DATA, *XML))

# The furanose ring, whose four bonds carry the pucker.
SUGAR_BONDS = {frozenset(b) for b in (("C1'", "C2'"), ("C2'", "C3'"), ("C3'", "C4'"),
                                      ("C4'", "O4'"), ("O4'", "C1'"))}
GLYCOSIDIC = {frozenset(("C1'", N.GLYCOSIDIC_N[r])) for r in RESIDUES if r in N.GLYCOSIDIC_N}


def family_of(n2, n3, same_residue):
    b = frozenset((n2, n3))
    if b in SUGAR_BONDS:
        return "NUCLEIC_SUGAR"
    if b in GLYCOSIDIC:
        return "NUCLEIC_CHI"
    # The 2'-OH rotor hangs off the sugar, not the chain.
    if b == frozenset(("C2'", "O2'")):
        return "NUCLEIC_CHI"
    backbone = {"P", "O5'", "C5'", "C4'", "C3'", "O3'", "OP1", "OP2"}
    if n2 in backbone and n3 in backbone:
        return "NUCLEIC_BACKBONE"
    return "NUCLEIC_CHI"


def tabulate(terms):
    e = np.zeros_like(PHI)
    t = np.zeros_like(PHI)
    for (n, k, phase) in terms:
        e += k * (1.0 + np.cos(n * PHI - phase))
        t += n * k * np.sin(n * PHI - phase)
    return e, t


def lookup(cls):
    """Most specific first, then the axis wildcard -- AMBER's own order."""
    tp = tables_ff.torsion_params
    for key in (cls, tuple(reversed(cls)), ("", cls[1], cls[2], ""), ("", cls[2], cls[1], "")):
        if key in tp:
            return tp[key]
    return None


tables, table_of, records, seen = [], {}, [], set()
per_family = defaultdict(int)

for resname in RESIDUES:
    classes = tables_ff.atoms_of(resname)
    internal, _external = tables_ff.bonds_of(resname)
    if not classes:
        print(f"  {resname}: absent du fichier, ignore")
        continue

    # Three residues of the same kind, linked O3'(i) - P(i+1): enough for every
    # torsion the central one owns, including the two that cross a link.
    atoms = {}       # (offset, name) -> class
    for off in (-1, 0, 1):
        for n, c in classes.items():
            atoms[(off, n)] = c
    adj = defaultdict(list)
    for off in (-1, 0, 1):
        for a, b in internal:
            adj[(off, a)].append((off, b))
            adj[(off, b)].append((off, a))
    for off in (-1, 0):
        adj[(off, "O3'")].append((off + 1, "P"))
        adj[(off + 1, "P")].append((off, "O3'"))

    def spell(k):
        off, n = k
        return ("+" if off > 0 else "-") * abs(off) + n

    for a2 in [(0, n) for n in classes]:
        for a3 in adj[a2]:
            if (a3, a2) < (a2, a3):
                continue
            for a1 in adj[a2]:
                if a1 == a3:
                    continue
                for a4 in adj[a3]:
                    if a4 == a2 or a4 == a1:
                        continue
                    cls = tuple(atoms[k] for k in (a1, a2, a3, a4))
                    terms = lookup(cls)
                    if not terms:
                        continue
                    names = tuple(spell(k) for k in (a1, a2, a3, a4))
                    key = (resname,) + names
                    if key in seen:
                        continue
                    seen.add(key)
                    tkey = tuple(sorted((int(n), round(float(k), 6), round(float(p), 6)) for (n, k, p) in terms))
                    if tkey not in table_of:
                        table_of[tkey] = len(tables)
                        tables.append(tabulate(terms))
                    fam = family_of(a2[1], a3[1], a2[0] == a3[0])
                    records.append(f"TORSION\t{resname}\t{fam}\t" + "\t".join(names) + f"\t{table_of[tkey]}")
                    per_family[fam] += 1

out = [
    f"# AMBER's {KIND} torsion terms, tabulated, one record per four-atom",
    "# quadruplet enumerated from the residue topology.",
    "#",
    "# TORSIONTABLE <id> <bins> <E_0> <T_0> ... over phi in [-180, 180]",
    "# TORSION <res> <FAMILY> <a1> <a2> <a3> <a4> <table>",
]
for i, (e, t) in enumerate(tables):
    out.append(f"TORSIONTABLE {i} {BINS} " + " ".join(f"{a:.5g} {b:.5g}" for a, b in zip(e, t)))
out += records
open(OUT, "w").write("\n".join(out) + "\n")

print(f"\n  {KIND} : {len(records)} torsions, {len(tables)} tables de {BINS} intervalles")
for fam, n in sorted(per_family.items()):
    print(f"    {fam:<20s} {n:5d}")
print(f"  -> {OUT}")
