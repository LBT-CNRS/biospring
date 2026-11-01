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
  TORSION <res>[:NTER|:CTER] <FAMILY> <a1> <a2> <a3> <a4> <id>

A free chain end is written under its own qualifier, because AMBER types it
from its own residue template: NALA's nitrogen is an ammonium N3 and CALA's
carbonyl oxygen a carboxylate O2, and amber99sb assigns no torsion at all to
N3-CT-C-N, HP-CT-C-O or H1-CT-C-O2 while assigning one to the identical
interior quadruplets. The residue name stays AMBER's own, as in amber.grp and
the .rbody files; the qualifier carries the position.

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


# The four atoms EVERY residue has. A record carries the residue name of its
# OWNER only: the format cannot say what the neighbour must be, so a "+X" or
# "-X" naming an atom that only some residue types have fires on every
# neighbour that happens to have that name. Measured on 072: ALA's copy of
# proline's own C(i-1)-N-CD-CG spelled "C +N +CD +CG" and fired on every Glu,
# Gln and Arg that followed an alanine -- seven of them, each a torsion across
# atoms 3.9 to 5.0 A apart, plus the matching planarity improper and two
# omegas. Restricting the neighbour reference to these four drops exactly those
# copies, and never the torsion: the residue that owns the unusual atoms spells
# them locally and keeps its own copy (proline's "-C N CD CG" survives, and it
# is the correct one).
UNIVERSAL = {"N", "CA", "C", "O"}

# Of those four, the ones whose own AMBER class depends on the residue they
# belong to being a free terminus. Read off the templates rather than listed,
# so a force field that retyped another backbone atom would be followed: "-"
# looks back, and the previous residue may be an N-terminus (its N becomes N3);
# "+" looks forward, and the next one may be a C-terminus (its O becomes O2).
_HAS_TERMINAL_TEMPLATE = sorted({r.name for r in G.residues if r.is_nterm or r.is_cterm})


def _neighbour_sensitive(prefix):
    out = set()
    for name in UNIVERSAL:
        for base in _HAS_TERMINAL_TEMPLATE:
            if G.resolved_class(base, name, prefix == "-", prefix == "+") != \
                    G.resolved_class(base, name, False, False):
                out.add(name)
    return out


SENSITIVE = {"-": _neighbour_sensitive("-"), "+": _neighbour_sensitive("+")}


# PDB v3 names the three hydrogens of a free N-terminus H1/H2/H3, which is
# what AMBER's NXXX templates use and so what these records name. OpenMM
# writes the first one "H" -- both all-atom protein inputs used here, 071's
# helix_loop_sheet.pdb and 072's gkinase.pdb, carry "CREATED WITH OPENMM" and
# spell it that way -- and the reader resolves atoms by name, so with H1 alone
# three of the nine ammonium rotor paths would silently never fire (measured on
# 071: 0.0920 of 0.2305 kJ.mol-1). So the record is written under both
# spellings, the way amber.grp carries ACE's methyl hydrogens under H1/H2/H3
# and HH31/HH32/HH33 at once. The reader keeps one torsion per quadruplet, so
# a structure using either spelling gets the term exactly once.
ALIASES = {"H1": ("H",)}


def spell(names, owner):
    """`names` plus the alternative spellings of the same quadruplet."""
    out = [names]
    if owner.is_nterm:
        for i, n in enumerate(names):
            for alias in ALIASES.get(n, ()):
                out.append(names[:i] + (alias,) + names[i + 1:])
    return out


def position_of(res):
    return "NTER" if res.is_nterm else ("CTER" if res.is_cterm else "interieur")


def qualified(res):
    """The record's residue column: AMBER's residue name, and the terminal
    qualifier when the record only applies at a free chain end."""
    if res.is_nterm:
        return res.name + ":NTER"
    if res.is_cterm:
        return res.name + ":CTER"
    return res.name


def owners_for(quad):
    """Which residues a quadruplet's record is written under, and how spelled.

    Returns (spellings, unconstrained): a list of (owner, four names relative
    to that owner), and whether the UNIVERSAL rule had to be waived.

    A torsion spanning two residues can be spelled from either of them. It is
    written under every spelling that reaches no further than the atoms
    UNIVERSAL names -- but a spelling reaching for a neighbour atom that is
    itself terminal-sensitive is dropped whenever another one survives. That is
    what lets a terminal residue's own record set decide its terminal-sensitive
    torsions. Psi, N-CA-C-+N, is the case: amber99sb gives N-CT-C-N three
    harmonics and N3-CT-C-N none, and the quadruplet can be spelled "N CA C +N"
    under the residue whose psi it is, or "-N -CA -C N" under the next one.
    Keeping the second would re-apply the interior psi to a residue whose own
    :NTER set deliberately omits it -- from a record that residue does not own,
    and that the reader has no way to qualify, since the format carries one
    residue name per record. Dropping it costs nothing: the surviving spelling
    resolves on exactly the same quadruplets.
    """
    owners = sorted({a.residue for a in quad}, key=lambda r: r.index)
    spellings = []
    for owner in owners:
        names = tuple(name_for(x, owner.index) for x in quad)
        if any(len(n) - len(n.lstrip("+-")) > 1 for n in names):
            continue  # more than one residue away: not addressable
        spellings.append((owner, names))

    safe = [(o, n) for (o, n) in spellings
            if all(x.lstrip("+-") in UNIVERSAL for x in n if x[0] in "+-")]
    # A torsion whose every spelling reaches for an unusual neighbour atom
    # would vanish entirely; keep the one with the most local atoms rather
    # than lose the term.
    if not safe and spellings:
        return [max(spellings, key=lambda s: sum(1 for x in s[1] if x[0] not in "+-"))], True

    insensitive = [(o, n) for (o, n) in safe
                   if not any(x[0] in "+-" and x.lstrip("+-") in SENSITIVE[x[0]] for x in n)]
    return (insensitive or safe), False


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
per_position = defaultdict(int)
# Torsions kept under an owner that still names an unusual neighbour atom --
# owners_for had to waive the UNIVERSAL rule because no spelling satisfied it.
# Reported rather than silent. One, at the time of writing: ACE's own
# CH3-C-+N-+H, whose only other spelling reads "-CH3 -C N H" and is no better.
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
                    # Each atom from its OWN residue's template: a free
                    # terminus is typed NXXX/CXXX, not XXX, which is what
                    # retypes N to N3, HA to HP and O to O2.
                    cls = tuple(G.resolved_class(x.residue.name, x.name,
                                                 x.residue.is_nterm, x.residue.is_cterm)
                                for x in (a1, a2, a3, a4))
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
                    # A torsion spanning two residues is written under each
                    # spelling owners_for keeps, not under a2's residue alone:
                    # that made it fire only where that residue type happens to
                    # precede -- PRO's CA-C-N-CD is the case that showed it,
                    # recorded under the chain's ALA spacer and never firing on
                    # a proline that follows anything else. The runtime drops
                    # the duplicate once both resolve to the same atoms.
                    safe, waived = owners_for((a1, a2, a3, a4))
                    unconstrained += waived
                    for owner, spelling in safe:
                        for names in spell(spelling, owner):
                            key = (qualified(owner),) + names
                            if key in seen:
                                continue
                            seen.add(key)
                            tkey = tuple(sorted((int(n), round(float(k), 6), round(float(p), 6))
                                                for (n, k, p) in terms))
                            if tkey not in table_of:
                                table_of[tkey] = len(tables)
                                tables.append(tabulate(terms))
                            fam = family_of(a2, a3)
                            records.append(f"TORSION\t{qualified(owner)}\t{fam}\t" + "\t".join(names) +
                                           f"\t{table_of[tkey]}")
                            per_family[fam] += 1
                            per_position[position_of(owner)] += 1

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
        hub_cls = G.resolved_class(hub.residue.name, hub.name,
                                   hub.residue.is_nterm, hub.residue.is_cterm)
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
                pc = [G.resolved_class(a.residue.name, a.name,
                                       a.residue.is_nterm, a.residue.is_cterm) for a in perm]
                if all(cls[j] in ("", pc[k]) for j, k in ((1, 0), (2, 1), (3, 2))):
                    if best is None or nwild < best[0]:
                        best = (nwild, perm, terms)
                    break
        if best is None:
            continue
        _, perm, terms = best
        quad = (perm[0], perm[1], hub, perm[2])
        # Same owner rule as the proper torsions above: proline's ring-nitrogen
        # improper spelled from the previous residue reads "C +CD +N +CA" and
        # fires on any Glu, Gln or Arg that follows.
        safe, waived = owners_for(quad)
        unconstrained += waived
        for owner, spelling in safe:
            for names in spell(spelling, owner):
                key = ("IMP", qualified(owner)) + names
                if key in seen:
                    continue
                seen.add(key)
                tkey = tuple(sorted((int(n), round(float(k), 6), round(float(p), 6)) for (n, k, p) in terms))
                if tkey not in table_of:
                    table_of[tkey] = len(tables)
                    tables.append(tabulate(terms))
                records.append(f"TORSION\t{qualified(owner)}\tPLANARITY\t" + "\t".join(names) +
                               f"\t{table_of[tkey]}")
                per_family["PLANARITY"] += 1
                per_position[position_of(owner)] += 1
                n_improper += 1

# ---- the invariant the reader relies on ------------------------------------
#
# A terminal residue takes its :NTER/:CTER records INSTEAD of the unqualified
# ones, so nothing else may re-create a torsion the terminal set deliberately
# omits. Only a record owned by the neighbour can: it reaches into the terminus
# with a "+"/"-" reference, and no record says what its neighbour must be. So
# an unqualified record must never name a neighbour atom whose own class
# depends on that neighbour being a free terminus -- which is what owners_for()
# arranges, and what this checks. Loud rather than silent: getting it wrong
# puts back exactly the 17.8 kJ.mol-1 terminal psi well the qualifier removes.
leaks = []
for record in records:
    fields = record.split("\t")
    if ":" in fields[1]:
        continue
    for name in fields[3:7]:
        if name[0] in "+-" and name.lstrip("+-") in SENSITIVE[name[0]]:
            leaks.append(record)
            break
if leaks:
    raise SystemExit("ERREUR : %d enregistrement(s) d'interieur nomment un atome voisin "
                     "sensible a sa position terminale, donc un terminus en heriterait :\n  %s"
                     % (len(leaks), "\n  ".join(leaks[:5])))

# And a residue AMBER gives a terminal template must get both qualified sets,
# or the reader falls back to the interior records there -- correct for the
# nucleic files, which have no qualifier at all, but a silent regression here.
_qualified_names = {r.split("\t")[1] for r in records if ":" in r.split("\t")[1]}
_missing = sorted(name + ":" + q for name in _HAS_TERMINAL_TEMPLATE for q in ("NTER", "CTER")
                  if name + ":" + q not in _qualified_names)
if _missing:
    raise SystemExit("ERREUR : aucun enregistrement pour %s" % ", ".join(_missing))

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
    "# TORSION <res>[:NTER|:CTER] <FAMILY> <a1> <a2> <a3> <a4> <table>",
    "# E in kJ/mol, T = -dV/dphi in kJ.mol-1.rad-1.",
    "#",
    "# A record qualified :NTER or :CTER applies ONLY to a free chain end, and",
    "# an unqualified one only to a residue with neighbours on both sides: a",
    "# terminus is typed from AMBER's own NXXX/CXXX template (N -> N3, HA -> HP,",
    "# O -> O2), and amber99sb assigns no torsion at all to N3-CT-C-N,",
    "# HP-CT-C-O or H1-CT-C-O2 while assigning one to the identical interior",
    "# quadruplets. The ammonium rotor (X-N3-CT-X) and the carboxylate improper",
    "# exist only at a terminus for the same reason.",
]
for i, (e, t) in enumerate(tables):
    out.append(f"TORSIONTABLE {i} {BINS} " + " ".join(f"{a:.5g} {b:.5g}" for a, b in zip(e, t)))
out += records
open(OUT, "w").write("\n".join(out) + "\n")

print(f"\n  {len(records)} torsions, {len(tables)} tables de {BINS} intervalles")
for fam, n in sorted(per_family.items()):
    print(f"    {fam:<12s} {n:5d}")
print()
for pos in ("interieur", "NTER", "CTER"):
    print(f"    {pos:<12s} {per_position[pos]:5d}")
if unconstrained:
    print(f"    {'(voisin non contraint)':<12s} {unconstrained:5d}")
print(f"\n  erreur d'interpolation la pire : {worst:.4f} kJ/mol")
print(f"  -> {OUT}")
