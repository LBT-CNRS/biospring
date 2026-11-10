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
import amber_nucleic_tables as N

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
    # The two TERMINAL hydroxyl rotors are not given that exception, and the
    # difference is the axis rather than the hydrogen: HO5' turns about O5'-C5'
    # and HO3' about C3'-O3', which are the beta and epsilon axes themselves,
    # where the 2'-OH turns about C2'-O2', a bond the chain never passes
    # through. So they fall through to the test below and come out
    # NUCLEIC_BACKBONE -- which also means -dihedralbackbone alone holds them,
    # and no run can enable the backbone and leave a chain end spinning.
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


# AMBER's terminal variants of a residue -- A5/A3/AN beside A, DA5/DA3/DAN
# beside DA -- and they are not spelling variants: the 5' template drops the
# phosphate and caps O5' with HO5', the 3' template caps O3' with HO3'.
# Enumerating the interior template alone is what left both terminal hydroxyl
# rotors with no torsion at all, measured on 075's hairpin against
# amber14/RNA.OL3: 6 quadruplets short, 6.2133 kJ/mol, the 5' and 3' hydroxyl
# hydrogens free to spin while the mesh held their bond and their angle.
#
# The records go out under the INTERIOR name all the same. The .bi.ff reader
# matches a record on the PDB residue name, and a PDB calls its terminal
# guanosine "G", never "G3" -- and it skips any record whose atoms it cannot
# resolve, so the two hydroxyl rules cost an interior residue nothing.
TERMINAL_SUFFIXES = ("5", "3", "N")


def templates_of(resname):
    """The interior template, then whichever terminal variants the file has.

    Interior first, so it is what assigns the table ids and what wins a
    disagreement -- see the conflict check in the enumeration below, which is
    fatal rather than silent.
    """
    return [resname] + [resname + s for s in TERMINAL_SUFFIXES
                        if tables_ff.atoms_of(resname + s)]


def chain_of(template, interior):
    """Three residues linked O3'(i) - P(i+1), the middle one `template`.

    The two neighbours are always the INTERIOR template, whatever the middle
    one is: no quadruplet crossing a link can reach a terminal hydroxyl, since
    a residue that HAS a neighbour on that side carries a phosphate there and
    no -OH at all.

    The links themselves are read off the template's own ExternalBond list
    rather than assumed. A 5' template has no P, so nothing joins it to the
    residue before it; a 3' template's O3' bonds to no next P. Adding a link
    the template denies would invent quadruplets AMBER never parameterised.
    """
    mid, (mid_internal, external) = tables_ff.atoms_of(template), tables_ff.bonds_of(template)
    nei, (nei_internal, _) = tables_ff.atoms_of(interior), tables_ff.bonds_of(interior)
    atoms, adj = {}, defaultdict(list)
    for off, classes, internal in ((-1, nei, nei_internal), (0, mid, mid_internal),
                                   (1, nei, nei_internal)):
        for n, c in classes.items():
            atoms[(off, n)] = c
        for a, b in internal:
            adj[(off, a)].append((off, b))
            adj[(off, b)].append((off, a))
    if "P" in external:
        adj[(-1, "O3'")].append((0, "P"))
        adj[(0, "P")].append((-1, "O3'"))
    if "O3'" in external:
        adj[(0, "O3'")].append((1, "P"))
        adj[(1, "P")].append((0, "O3'"))
    return atoms, adj


def spell(k):
    """The .bi.ff spelling of an atom, relative to the central residue."""
    off, n = k
    return ("+" if off > 0 else "-") * abs(off) + n


tables, table_of, records, seen = [], {}, [], {}
per_family = defaultdict(int)


def remember(key, tkey, template, names):
    """True the first time a quadruplet is met, False when it repeats.

    Nearly every quadruplet of a terminal template is one the interior
    template already carried, so repeats are the rule here and not a mistake.
    A repeat carrying DIFFERENT parameters would be one, and an unanswerable
    one: a record names its residue the way a PDB does, so there is no way to
    say "this one, but only at a 3' end". Neither amber14 nucleic file does
    it -- the terminal O5'/O3' is retyped OS -> OH, but every quadruplet it
    sits in resolves through the same wildcard key either way -- so this fires
    on a force field that has changed under the generator, not on these two.
    """
    if key in seen:
        if seen[key] != tkey:
            sys.exit(f"  {template} : {' '.join(names)} serait emis sous {key[0]} avec des "
                     f"parametres que le gabarit interieur ne partage pas -- le format .bi.ff "
                     f"n'a pas de quoi distinguer un terminus")
        return False
    seen[key] = tkey
    return True


for resname in RESIDUES:
    if not tables_ff.atoms_of(resname):
        print(f"  {resname}: absent du fichier, ignore")
        continue

    for template in templates_of(resname):
        atoms, adj = chain_of(template, resname)
        for a2 in [(0, n) for n in tables_ff.atoms_of(template)]:
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
                        tkey = tuple(sorted((int(n), round(float(k), 6), round(float(p), 6)) for (n, k, p) in terms))
                        if not remember((resname,) + names, tkey, template, names):
                            continue
                        if tkey not in table_of:
                            table_of[tkey] = len(tables)
                            tables.append(tabulate(terms))
                        fam = family_of(a2[1], a3[1], a2[0] == a3[0])
                        records.append(f"TORSION\t{resname}\t{fam}\t" + "\t".join(names) + f"\t{table_of[tkey]}")
                        per_family[fam] += 1

# ---- impropers -------------------------------------------------------------
#
# Same convention as the protein side, and it is the one thing here that cannot
# be guessed: the FIRST field of the entry is the central atom, and the torsion
# is emitted as (2, 3, 1, 4) -- central third. Nucleic files key their
# parameters by type rather than by class, which atoms_of already returns.
improper_entries = []
for _e in tables_ff.root.find("PeriodicTorsionForce").findall("Improper"):
    key = tuple((_e.get(f"type{i}") if tables_ff.keyed_by_type else _e.get(f"class{i}")) or "" for i in (1, 2, 3, 4))
    terms, i = [], 1
    while _e.get(f"periodicity{i}") is not None:
        k = float(_e.get(f"k{i}"))
        if k != 0.0:
            terms.append((int(_e.get(f"periodicity{i}")), k, float(_e.get(f"phase{i}"))))
        i += 1
    if terms:
        improper_entries.append((key, terms, sum(1 for c in key if c == "")))

import itertools

for resname in RESIDUES:
    if not tables_ff.atoms_of(resname):
        continue
    # Over the terminal templates too, for the same reason as the propers, even
    # though neither hydroxyl adds one: a hub needs exactly three substituents
    # and a capped O5'/O3' has two, just as the interior one does. Enumerating
    # them anyway is what makes that a measurement rather than an assumption.
    for template in templates_of(resname):
        atoms, adj = chain_of(template, resname)
        for hub in [(0, n) for n in tables_ff.atoms_of(template)]:
            subs = adj[hub]
            if len(subs) != 3:
                continue
            best = None
            for (key, terms, nwild) in improper_entries:
                if key[0] not in ("", atoms[hub]):
                    continue
                for perm in itertools.permutations(subs):
                    pc = [atoms[k] for k in perm]
                    if all(key[j] in ("", pc[k]) for j, k in ((1, 0), (2, 1), (3, 2))):
                        if best is None or nwild < best[0]:
                            best = (nwild, perm, terms)
                        break
            if best is None:
                continue
            _, perm, terms = best
            names = tuple(spell(k) for k in (perm[0], perm[1], hub, perm[2]))
            tkey = tuple(sorted((int(n), round(float(k), 6), round(float(p), 6)) for (n, k, p) in terms))
            if not remember((resname, "IMP") + names, tkey, template, names):
                continue
            if tkey not in table_of:
                table_of[tkey] = len(tables)
                tables.append(tabulate(terms))
            # Impropers share the PLANARITY family with the protein side: it is
            # what the enum calls an improper, whatever the chemistry.
            records.append(f"TORSION\t{resname}\tPLANARITY\t" + "\t".join(names) + f"\t{table_of[tkey]}")
            per_family["PLANARITY"] += 1

# ---- bond stiffnesses ------------------------------------------------------
#
# NOT A NEW TERM. The --rigidbody mesh already springs every one of AMBER's
# bonds and every one of its valence angles -- measured on 075's hairpin, 419
# of 419 and 746 of 746, with no surplus pair. What it lacks is the VALUE:
# every spring carries pdb2spn's uniform --stiffness where AMBER's bonds run
# from 1925 to 4770 kJ.mol-1.A-2. These records say what each bond's constant
# is, and the reader only ever retunes a spring that already exists.
#
# The angles are left out on purpose: a 1-3 distance spring is not an angle
# spring, converting one into the other gets the soft modes wrong (see
# 7660b34), and at the uniform value it already delivers 0.60 of AMBER's
# angular stiffness.
#
# TERMINAL TEMPLATES DISAGREE HERE, unlike for the torsions. A 5' end retypes
# O5' from OS to OH and a 3' end does the same to O3', so the C5'-O5' bond
# genuinely has a different constant at a free end -- and the nucleic .bi.ff
# writes every record under the interior residue name, because that is what a
# PDB calls a terminal nucleotide, so there is nowhere to put the variant. The
# interior value is emitted and the disagreement is REPORTED rather than
# silently dropped.
STIFF_CONVERT = 0.01  # kJ.mol-1.nm-2 -> kJ.mol-1.A-2

stiffness_records, stiff_seen, stiff_missing, stiff_clash = [], {}, 0, []
for resname in RESIDUES:
    if not tables_ff.atoms_of(resname):
        continue
    for template in templates_of(resname):
        atoms, adj = chain_of(template, resname)
        for a in list(adj):
            for b in adj[a]:
                if (b, a) < (a, b):
                    continue
                if a[0] != 0 and b[0] != 0:
                    continue  # neither end is in the residue the record belongs to
                if abs(a[0]) > 1 or abs(b[0]) > 1:
                    continue
                par = tables_ff.bond_params.get(frozenset((atoms[a], atoms[b])))
                if par is None:
                    stiff_missing += 1
                    continue
                k = par[1] * STIFF_CONVERT
                names = (spell(a), spell(b))
                key = (resname,) + names
                if key in stiff_seen:
                    if abs(stiff_seen[key] - k) > 1e-3:
                        stiff_clash.append((resname, template, names, stiff_seen[key], k))
                    continue
                stiff_seen[key] = k
                stiffness_records.append(f"STIFFNESS\t{resname}\t" + "\t".join(names) + f"\t{k:.6g}")
if stiff_missing:
    print(f"  {stiff_missing} liaison(s) sans parametre AMBER, ignoree(s)")
for resname, template, names, k0, k1 in stiff_clash:
    print(f"  {resname} {' '.join(names)} : {template} donne {k1:.1f} contre {k0:.1f} kJ/mol/A2 "
          f"a l'interieur ({100.0 * abs(k1 - k0) / k0:.1f} %) -- l'interieur est emis")

out = [
    f"# AMBER's {KIND} torsion terms, tabulated, one record per four-atom",
    "# quadruplet enumerated from the residue topology -- the interior templates",
    "# and AMBER's 5'/3' terminal ones, whose records are written under the",
    "# interior name because that is what a PDB calls a terminal residue.",
    "#",
    "# TORSIONTABLE <id> <bins> <E_0> <T_0> ... over phi in [-180, 180]",
    "# TORSION <res> <FAMILY> <a1> <a2> <a3> <a4> <table>",
    "# STIFFNESS <res> <a1> <a2> <k>   k in kJ.mol-1.A-2",
    "#",
    "# A STIFFNESS record adds NOTHING: the rigid-body mesh already springs every",
    "# bond and every valence angle, and this only says what AMBER's constant is",
    "# for the bond it names, so -bondstiffness can retune the spring that is",
    "# already there. The angles are left out on purpose -- a 1-3 distance spring",
    "# is not an angle spring. The value is the INTERIOR template's: a 5'/3' end",
    "# retypes its own O5'/O3' and the format has no terminal qualifier here.",
]
for i, (e, t) in enumerate(tables):
    out.append(f"TORSIONTABLE {i} {BINS} " + " ".join(f"{a:.5g} {b:.5g}" for a, b in zip(e, t)))
out += records
out += stiffness_records
open(OUT, "w").write("\n".join(out) + "\n")

print(f"\n  {KIND} : {len(records)} torsions, {len(tables)} tables de {BINS} intervalles")
print(f"  {len(stiffness_records)} raideurs de liaison")
for fam, n in sorted(per_family.items()):
    print(f"    {fam:<20s} {n:5d}")
print(f"  -> {OUT}")
