#!/usr/bin/env python3
"""Tabulate AMBER's torsion terms, enumerated from the residue topology itself.

Every four-atom quadruplet of every bond, taken from the chain the generator
already builds, with the terms AMBER assigns to it and AMBER's own phases. The
ring generator is not consulted at all: harvesting quadruplets from its
internals gave two pairs per backbone axis instead of six, and leaked its
internal atom spelling into the file. Enumerating from atoms and bonds is
complete by construction and has no ordering to respect.

The generator does every trigonometric evaluation, once. The runtime measures a
dihedral, indexes a table and interpolates -- it knows nothing about harmonics
or phases.

  TORSIONTABLE <id> <bins> <E_0> <T_0> ...   over phi in [-180, 180]
  TORSION <res> <FAMILY> <a1> <a2> <a3> <a4> <id>

E in kJ/mol, T = -dV/dphi in kJ.mol-1.rad-1.
"""
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bonded_axes as BA
import amber_protein_tables as G

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO_ROOT,
                                                        "data/reducerules/ProteinAtomTorsion.bi.ff")
BINS = int(os.environ.get("BIOSPRING_TORSION_BINS", "512"))
PHI = np.linspace(-np.pi, np.pi, BINS + 1)

# Which residues get records. The chain brackets every residue of interest with
# ALA spacers, so an interior residue has real neighbours on both sides; the
# spacers themselves would only duplicate ALA.
WANTED = sorted({r.name for r in G.residues})

adj = defaultdict(list)
for a, b in G.bonds:
    adj[a.index].append(b)
    adj[b.index].append(a)
byIndex = {a.index: a for r in G.residues for a in r.atoms}


def name_for(atom, central):
    """The .bi.ff spelling, relative to the residue the record belongs to."""
    d = atom.residue.index - central
    return ("+" if d > 0 else "-") * abs(d) + atom.name if d else atom.name


def family_of(b, c):
    nb, nc = b.name, c.name
    same = b.residue.index == c.residue.index
    if same and {nb, nc} == {"N", "CA"}:
        return "PHI"
    if same and {nb, nc} == {"CA", "C"}:
        return "PSI"
    if not same and {nb, nc} == {"C", "N"}:
        return "OMEGA"
    return "SIDECHAIN"


def tabulate(terms):
    e = np.zeros_like(PHI)
    t = np.zeros_like(PHI)
    for (n, k, phase) in terms:
        e += k * (1.0 + np.cos(n * PHI - phase))
        t += n * k * np.sin(n * PHI - phase)
    return e, t


tables, table_of = [], {}
records, seen = [], set()
per_family = defaultdict(int)
# Torsions kept under an owner that still names an unusual neighbour atom
# (see the UNIVERSAL filter below). Expected to stay at zero; reported so it
# cannot become a silent exception.
unconstrained = 0

for res in G.residues:
    if res.name not in WANTED:
        continue
    central = res.index
    for a2 in res.atoms:
        for a3 in adj[a2.index]:
            # One record per bond, not two: the reversed quadruplet is the
            # same torsion.
            if (a3.index, a2.index) < (a2.index, a3.index):
                continue
            for a1 in adj[a2.index]:
                if a1.index == a3.index:
                    continue
                for a4 in adj[a3.index]:
                    if a4.index == a2.index or a4.index == a1.index:
                        continue
                    cls = tuple(G.resolved_class(x.residue.name, x.name, False, False) for x in (a1, a2, a3, a4))
                    # Both axis orders: the table stores X-a-b-X under one
                    # spelling only, and GLU's CT-CT-C-O2 is the case that
                    # showed it -- AMBER matches it, a single-order lookup
                    # does not.
                    terms = BA.lookup_torsion_specific(*cls)
                    if terms is None:
                        terms = BA.lookup_torsion_specific(*reversed(cls))
                    if terms is None:
                        terms = BA.lookup_torsion_wildcard(cls[1], cls[2])
                    if terms is None:
                        terms = BA.lookup_torsion_wildcard(cls[2], cls[1])
                    if not terms:
                        continue
                    # A torsion spanning two residues is written under each of
                    # them, spelled relative to that one. Attributing it to a2's
                    # residue alone makes it fire only where that residue type
                    # happens to precede -- PRO's CA-C-N-CD is the case that
                    # showed it: recorded under the chain's ALA spacer, it never
                    # fired on a proline that follows anything else. The runtime
                    # drops the duplicate once both resolve to the same atoms.
                    # sorted, not a bare set: a set of residue objects iterates
                    # in id()-based hash order, which changes between runs, and
                    # the emitted file was not byte-reproducible because of it.
                    owners = sorted({a.residue for a in (a1, a2, a3, a4)},
                                    key=lambda r: r.index)
                    spellings = []
                    for owner in owners:
                        names = tuple(name_for(x, owner.index) for x in (a1, a2, a3, a4))
                        if any(len(n) - len(n.lstrip("+-")) > 1 for n in names):
                            continue  # more than one residue away: not addressable
                        spellings.append((owner, names))

                    # A record carries the residue name of its OWNER only: the
                    # format cannot say what the neighbour must be. So a "+X"
                    # or "-X" naming an atom that only some residue types have
                    # fires on every neighbour that happens to have that name.
                    # Measured on 072: ALA's copy of proline's own
                    # C(i-1)-N-CD-CG spelled "C +N +CD +CG" and fired on every
                    # Glu, Gln and Arg that followed an alanine -- seven of
                    # them, each a torsion across atoms 3.9 to 5.0 A apart,
                    # plus the matching planarity improper and two omegas.
                    # Restricting the neighbour reference to the four atoms
                    # EVERY residue has drops exactly those copies, and never
                    # the torsion: the residue that owns the unusual atoms
                    # spells them locally and keeps its own copy (proline's
                    # "-C N CD CG" survives, and it is the correct one).
                    UNIVERSAL = {"N", "CA", "C", "O"}
                    safe = [(o, n) for (o, n) in spellings
                            if all(x.lstrip("+-") in UNIVERSAL
                                   for x in n if x[0] in "+-")]
                    # A torsion whose every spelling reaches for an unusual
                    # neighbour atom would vanish entirely; keep the one with
                    # the most local atoms rather than lose the term.
                    if not safe and spellings:
                        safe = [max(spellings,
                                    key=lambda s: sum(1 for x in s[1] if x[0] not in "+-"))]
                        unconstrained += 1
                    for owner, names in safe:
                        key = (owner.name,) + names
                        if key in seen:
                            continue
                        seen.add(key)
                        tkey = tuple(sorted((int(n), round(float(k), 6), round(float(p), 6))
                                            for (n, k, p) in terms))
                        if tkey not in table_of:
                            table_of[tkey] = len(tables)
                            tables.append(tabulate(terms))
                        fam = family_of(a2, a3)
                        records.append(f"TORSION\t{owner.name}\t{fam}\t" + "\t".join(names) +
                                       f"\t{table_of[tkey]}")
                        per_family[fam] += 1

# ---- impropers -------------------------------------------------------------
#
# AMBER writes an improper with the CENTRAL atom third, and matches it against
# an sp2 hub's three substituents in any order, preferring the entry with the
# fewest wildcards. What it holds is planarity, which the mesh also holds -- but
# only as well as its stiffness, so the term is not redundant.
improper_entries = []
for _t in G.root.find("PeriodicTorsionForce").findall("Improper"):
    _cls = tuple(_t.get(f"class{i}") or "" for i in (1, 2, 3, 4))
    _terms, _i = [], 1
    while _t.get(f"periodicity{_i}") is not None:
        _k = float(_t.get(f"k{_i}"))
        if _k != 0.0:
            _terms.append((int(_t.get(f"periodicity{_i}")), _k, float(_t.get(f"phase{_i}"))))
        _i += 1
    if _terms:
        improper_entries.append((_cls, _terms, sum(1 for c in _cls if c == "")))

import itertools

n_improper = 0
for res in G.residues:
    for hub in res.atoms:
        subs = adj[hub.index]
        if len(subs) != 3:
            continue
        hub_cls = G.resolved_class(hub.residue.name, hub.name, False, False)
        best = None
        for (cls, terms, nwild) in improper_entries:
            # class1 is the CENTRAL atom of the entry, and the torsion is
            # emitted as (class2, class3, class1, class4) -- central third.
            # Reading the entry as if the central were already third found the
            # aromatic rings, whose wildcards forgive it, and missed every
            # backbone amide: AMBER's N-C-CT-H, 4.6 kJ/mol, is what keeps the
            # peptide's N-H in its plane.
            if cls[0] not in ("", hub_cls):
                continue
            for perm in itertools.permutations(subs):
                pc = [G.resolved_class(a.residue.name, a.name, False, False) for a in perm]
                if all(cls[j] in ("", pc[k]) for j, k in ((1, 0), (2, 1), (3, 2))):
                    if best is None or nwild < best[0]:
                        best = (nwild, perm, terms)
                    break
        if best is None:
            continue
        _, perm, terms = best
        quad = (perm[0], perm[1], hub, perm[2])
        owners = sorted({a.residue for a in quad}, key=lambda r: r.index)
        spellings = []
        for owner in owners:
            names = tuple(name_for(x, owner.index) for x in quad)
            if any(len(n) - len(n.lstrip("+-")) > 1 for n in names):
                continue
            spellings.append((owner, names))
        # Same neighbour-reference rule as the proper torsions above: proline's
        # ring-nitrogen improper spelled from the previous residue reads
        # "C +CD +N +CA" and fires on any Glu, Gln or Arg that follows.
        safe = [(o, n) for (o, n) in spellings
                if all(x.lstrip("+-") in {"N", "CA", "C", "O"}
                       for x in n if x[0] in "+-")]
        if not safe and spellings:
            safe = [max(spellings,
                        key=lambda s: sum(1 for x in s[1] if x[0] not in "+-"))]
            unconstrained += 1
        for owner, names in safe:
            key = ("IMP", owner.name) + names
            if key in seen:
                continue
            seen.add(key)
            tkey = tuple(sorted((int(n), round(float(k), 6), round(float(p), 6)) for (n, k, p) in terms))
            if tkey not in table_of:
                table_of[tkey] = len(tables)
                tables.append(tabulate(terms))
            records.append(f"TORSION\t{owner.name}\tPLANARITY\t" + "\t".join(names) + f"\t{table_of[tkey]}")
            per_family["PLANARITY"] += 1
            n_improper += 1

worst = 0.0
fine = np.linspace(-np.pi, np.pi, 8 * BINS + 1)
for (e, _t), tkey in zip(tables, table_of):
    ef = np.zeros_like(fine)
    for (n, k, p) in tkey:
        ef += k * (1.0 + np.cos(n * fine - p))
    worst = max(worst, float(np.abs(np.interp(fine, PHI, e) - ef).max()))

out = [
    "# AMBER's torsion terms, tabulated, one record per four-atom quadruplet",
    "# enumerated from the residue topology. Nothing is summed over pairs, so no",
    "# phase origin is chosen and no sign convention is guessed.",
    "#",
    "# TORSIONTABLE <id> <bins> <E_0> <T_0> ... over phi in [-180, 180]",
    "# TORSION <res> <FAMILY> <a1> <a2> <a3> <a4> <table>",
    "# E in kJ/mol, T = -dV/dphi in kJ.mol-1.rad-1.",
]
for i, (e, t) in enumerate(tables):
    out.append(f"TORSIONTABLE {i} {BINS} " + " ".join(f"{a:.5g} {b:.5g}" for a, b in zip(e, t)))
out += records
open(OUT, "w").write("\n".join(out) + "\n")

print(f"\n  {len(records)} torsions, {len(tables)} tables de {BINS} intervalles")
for fam, n in sorted(per_family.items()):
    print(f"    {fam:<12s} {n:5d}")
if unconstrained:
    print(f"    {'(voisin non contraint)':<12s} {unconstrained:5d}")
print(f"\n  erreur d'interpolation la pire : {worst:.4f} kJ/mol")
print(f"  -> {OUT}")
