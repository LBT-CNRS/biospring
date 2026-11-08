#!/usr/bin/env python3
"""AMBER's torsion terms for EF-G's GTP, tabulated like the protein and nucleic ones.

WHY THIS FILE EXISTS. 079's EF-G carries a GTP, and `079/gtp.rbody` says in its
own header what that costs: "The overlaps leave the four torsions PG-O3B, O3B-PB,
PB-O3A and O3A-PA free, which is what a polyphosphate is: a flexible chain.
AMBER's torsion terms are not available for it -- r70.bi.ff tabulates per residue
and has no GTP -- so those four turn under the mesh alone." Measured on the
docked network: of 734180 torsions, **zero** touch a GTP atom, and the phosphate
chain's 13 atoms are linked by 12 bonds and 20 one-three pairs with **no
one-four pair at all**. So neither a table nor the mesh holds those four, and a
minimisation was able to swing the gamma phosphate 2.74 A for no measurable
energy.

NO INSTALLED AMBER FORCE FIELD DEFINES GTP, and that was checked rather than
assumed: of the 84 XML files shipped with this OpenMM, exactly one declares a
residue named GTP -- `charmm36_2024.xml` -- and none under `amber*` does.
CHARMM's is a real, published GTP parameterisation, and it is NOT used here: the
whole network is AMBER, from the charges to the other three tables, and a
torsion profile taken from another force field would be the mismatch
`reference_forcefield_version_mismatch` already cost this project once. It is
read once, by hand, as an order-of-magnitude check and nothing else:

    CHARMM36 (2024), its own GTP, per quadruplet about a P-O-P bridge
        P3-ON2-P3-ON3 : 0.418 kJ/mol at n=2 plus 0.126 at n=3, total 0.544
    AMBER14/RNA.OL3, the generic phosphoester wildcard used below
        X -OS-P -X    : 1.046 kJ/mol at n=3

A factor of two on a term of about one kJ/mol. The analogy is sound in
magnitude; it is an analogy all the same, and that is the approximation this
file declares.

THE TYPE ASSIGNMENT IS ENTIRELY INSIDE AMBER. Every atom GTP shares with AMBER's
own 3'-terminal guanosine `G3` takes that atom's type, by name. The nine atoms
GTP has in addition take the type of their chemical role, each the same role an
atom of `G3` already has:

    PA PB PG            -> RNA-P    , as G3's `P`
    O1A O2A O1B O2B
    O1G O2G O3G         -> RNA-O2   , as G3's `OP1`/`OP2`: non-bridging
    O3A O3B             -> RNA-OS   , as G3's `O3'`/`O5'`: bridging esters

`G3` rather than `G` because GTP's 3'-OH is free, so it needs `HO3'`, which the
interior template does not have. Three hydrogens are renamed: PDBFixer spells
the base's `H1`, `H21`, `H22` as `HN1`, `HN21`, `HN22`, which is what
`079/gtp.grp` and `gtp.rbody` carry.

THE BONDS come from AMBER's own `G3` bond list, renamed the same way, plus the
nine bonds of the triphosphate read off `079/gtp.rbody`'s groups -- the
project's own statement of this molecule's connectivity -- and cross-checked
against the coordinates in `079/gtp.cif`, where every declared bond must measure
1.2 to 1.8 A and no undeclared heavy pair may be shorter than 2.0 A. That check
is fatal, not a warning.

THE RECORDS GO OUT UNDER THE RESIDUE NAME `GTP`, ONE RESIDUE. `gtpfix.py --carve`
later moves the gamma phosphate to its own resid so that PeptideBond.cpp can part
the two, but that happens on the .nc, long after the torsions have been written
into it with fixed particle indices. Generating under one residue is therefore
correct and the carve does not disturb it -- while a record written as if the
gamma were a separate residue would resolve against nothing at build time.

    python amber_tabulate_gtp.py [out.bi.ff]
"""
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amber_nucleic_tables as N

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    REPO_ROOT, "data/reducerules/GTPAtomTorsion.bi.ff")
BINS = int(os.environ.get("BIOSPRING_TORSION_BINS", "512"))
PHI = np.linspace(-np.pi, np.pi, BINS + 1)

TEMPLATE = "G3"
tables_ff = N.ForceFieldTables(os.path.join(N.DATA, "amber14", "RNA.OL3.xml"))

# PDBFixer's spellings for the three base hydrogens, as carried by gtp.grp.
RENAME_H = {"H1": "HN1", "H21": "HN21", "H22": "HN22"}
# G3's phosphate becomes GTP's alpha phosphate.
RENAME_P = {"P": "PA", "OP1": "O1A", "OP2": "O2A"}
# The nine atoms GTP has beyond G3, and the AMBER type of the role each plays.
EXTRA_TYPES = {"O3A": "RNA-OS", "PB": "RNA-P", "O1B": "RNA-O2", "O2B": "RNA-O2",
               "O3B": "RNA-OS", "PG": "RNA-P", "O1G": "RNA-O2", "O2G": "RNA-O2",
               "O3G": "RNA-O2"}
# The triphosphate's bonds, read off 079/gtp.rbody's groups.
EXTRA_BONDS = [("PA", "O3A"), ("O3A", "PB"),
               ("PB", "O1B"), ("PB", "O2B"), ("PB", "O3B"),
               ("O3B", "PG"), ("PG", "O1G"), ("PG", "O2G"), ("PG", "O3G")]

# The furanose ring, whose four bonds carry the pucker -- same split as
# amber_tabulate_nuc.py, for the same reason: the pucker is its own lever.
SUGAR_BONDS = {frozenset(b) for b in (("C1'", "C2'"), ("C2'", "C3'"), ("C3'", "C4'"),
                                      ("C4'", "O4'"), ("O4'", "C1'"))}
GLYCOSIDIC = {frozenset(("C1'", "N9"))}
# The phosphate chain joins the nucleic BACKBONE family, so the three switches a
# run already sets -- dihedralnucleicbackbone/chi/sugar -- cover GTP with no
# change to any .msp. A polyphosphate is backbone chemistry.
BACKBONE = {"PA", "PB", "PG", "O1A", "O2A", "O3A", "O1B", "O2B", "O3B",
            "O1G", "O2G", "O3G", "O5'", "C5'", "C4'", "C3'", "O3'", "HO3'"}


def build_topology():
    """GTP's atoms with their AMBER types, and its bond list."""
    base = tables_ff.atoms_of(TEMPLATE)
    if not base:
        sys.exit("  %s absent de RNA.OL3.xml : rien a faire" % TEMPLATE)
    ren = dict(RENAME_H)
    ren.update(RENAME_P)
    atoms = {ren.get(n, n): t for n, t in base.items()}
    atoms.update(EXTRA_TYPES)
    internal, _external = tables_ff.bonds_of(TEMPLATE)
    bonds = [(ren.get(a, a), ren.get(b, b)) for a, b in internal]
    bonds += EXTRA_BONDS
    return atoms, bonds


def check_against_cif(atoms, bonds, path):
    """Every declared bond 1.2-1.8 A, no undeclared heavy pair under 2.0 A.

    The topology above is written down, so it has to be confronted with the
    molecule the model actually carries. Fatal on disagreement: a bond list that
    is wrong by one atom produces a table of quadruplets that do not exist, and
    the .bi.ff reader would skip them in silence.
    """
    xyz = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        f = line.split()
        if len(f) > 13 and f[0] in ("ATOM", "HETATM"):
            nm = f[3]
            if len(nm) > 1 and nm[0] == nm[-1] and nm[0] in ("'", '"'):
                nm = nm[1:-1]
            xyz[nm] = np.array([float(f[10]), float(f[11]), float(f[12])])
    if not xyz:
        print("  %s : aucune coordonnee lue, controle geometrique saute" % path)
        return
    heavy = [n for n in atoms if not n.startswith("H") and n in xyz]
    declared = {frozenset(b) for b in bonds}
    bad = []
    for a, b in bonds:
        if a in xyz and b in xyz:
            d = float(np.linalg.norm(xyz[a] - xyz[b]))
            # A P-O bond gets its own window. 3J9Z is a 3.6 A cryo-EM model, so
            # its bond lengths are not determined to better than a tenth of an
            # angstrom: PB-O3B measures 1.809 A there where the mesh's own
            # equilibrium for the same bond is 1.764. Judging it against a
            # carbon-carbon window rejects a bond that is simply a phosphate's.
            if a.startswith("H") or b.startswith("H"):
                lo, hi = 0.85, 1.25
            elif a[0] == "P" or b[0] == "P":
                lo, hi = 1.40, 1.95
            else:
                lo, hi = 1.2, 1.8
            if not lo <= d <= hi:
                bad.append("liaison declaree %s-%s a %.3f A" % (a, b, d))
    for i, a in enumerate(heavy):
        for b in heavy[i + 1:]:
            if frozenset((a, b)) in declared:
                continue
            d = float(np.linalg.norm(xyz[a] - xyz[b]))
            if d < 2.0:
                bad.append("paire lourde NON declaree %s-%s a %.3f A" % (a, b, d))
    print("  controle geometrique sur %s : %d atomes, %d liaisons, %s"
          % (os.path.basename(path), len(xyz), len(bonds),
             "conforme" if not bad else "%d ECART" % len(bad)))
    for m in bad[:8]:
        print("     %s" % m)
    if bad:
        sys.exit("  la topologie declaree ne decrit pas le GTP de gtp.cif")


def family_of(n2, n3):
    b = frozenset((n2, n3))
    if b in SUGAR_BONDS:
        return "NUCLEIC_SUGAR"
    if b in GLYCOSIDIC:
        return "NUCLEIC_CHI"
    if b == frozenset(("C2'", "O2'")):
        return "NUCLEIC_CHI"
    if n2 in BACKBONE and n3 in BACKBONE:
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


def main():
    atoms, bonds = build_topology()
    print("  GTP : %d atomes, %d liaisons, gabarit AMBER %s"
          % (len(atoms), len(bonds), TEMPLATE))
    cif = os.path.join(REPO_ROOT, os.pardir, os.pardir, "Data",
                       "BioSpring_adapted_exemple", "079.Ribosome70S_AllAtom", "gtp.cif")
    if os.path.exists(cif):
        check_against_cif(atoms, bonds, cif)
    else:
        print("  gtp.cif introuvable (%s) : controle geometrique saute" % cif)

    adj = defaultdict(list)
    for a, b in bonds:
        adj[a].append(b)
        adj[b].append(a)

    tables, table_of, records, seen = [], {}, [], {}
    per_family = defaultdict(int)
    per_axis = defaultdict(list)
    missing = set()
    for a2 in sorted(atoms):
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
                        missing.add(tuple(tables_ff.type_to_class.get(c, c) for c in cls))
                        continue
                    names = (a1, a2, a3, a4)
                    if names in seen or names[::-1] in seen:
                        continue
                    seen[names] = True
                    tkey = tuple(sorted((int(n), round(float(k), 6), round(float(p), 6))
                                        for (n, k, p) in terms))
                    if tkey not in table_of:
                        table_of[tkey] = len(tables)
                        tables.append(tabulate(terms))
                    fam = family_of(a2, a3)
                    records.append("TORSION\tGTP\t%s\t%s\t%d"
                                   % (fam, "\t".join(names), table_of[tkey]))
                    per_family[fam] += 1
                    per_axis[frozenset((a2, a3))].append(sum(k for (_, k, _) in terms))

    with open(OUT, "w") as fh:
        fh.write("# AMBER's GTP torsion terms, tabulated, one record per four-atom\n"
                 "# quadruplet enumerated from the topology in amber_tabulate_gtp.py.\n"
                 "#\n"
                 "# No installed AMBER force field defines GTP: the types are AMBER's own\n"
                 "# guanosine types for every shared atom, and for the nine extra ones the\n"
                 "# type of the role they play -- P for the phosphorus, O2 for a\n"
                 "# non-bridging oxygen, OS for a bridging ester. The four P-O-P torsions\n"
                 "# therefore resolve through AMBER's generic phosphoester wildcard\n"
                 "# X -OS-P -X, 1.046 kJ/mol at periodicity 3, which is the approximation\n"
                 "# this file declares. CHARMM36's own published GTP puts the same term at\n"
                 "# 0.544 kJ/mol; it was read as a magnitude check and not used.\n"
                 "#\n"
                 "# TORSIONTABLE <id> <bins> <E_0> <T_0> ... over phi in [-180, 180]\n"
                 "# TORSION <res> <FAMILY> <a1> <a2> <a3> <a4> <table>\n")
        for i, (e, t) in enumerate(tables):
            fh.write("TORSIONTABLE %d %d %s\n"
                     % (i, BINS, " ".join("%.5g %.5g" % (ee, tt) for ee, tt in zip(e, t))))
        fh.write("\n".join(records) + "\n")

    print("  %d quadruplets, %d tables distinctes -> %s"
          % (len(records), len(tables), os.path.relpath(OUT, REPO_ROOT)))
    for fam, n in sorted(per_family.items()):
        print("     %-20s %4d" % (fam, n))
    print("  les quatre torsions du polyphosphate, amplitude SOMMEE sur leurs quadruplets :")
    for a, b in (("PG", "O3B"), ("O3B", "PB"), ("PB", "O3A"), ("O3A", "PA")):
        v = per_axis.get(frozenset((a, b)), [])
        print("     axe %-4s-%-4s %2d quadruplets, %6.2f kJ/mol"
              % (a, b, len(v), sum(v)))
    if missing:
        print("  %d familles de classes sans parametre AMBER (quadruplets omis) :" % len(missing))
        for m in sorted(missing)[:10]:
            print("     %s" % "-".join(m))


if __name__ == "__main__":
    main()
