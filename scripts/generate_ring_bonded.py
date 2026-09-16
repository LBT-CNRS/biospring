#!/usr/bin/env python3
"""AMBER's real bond and angle constants for the furanose ring.

Everywhere else the --rigidbody mesh is enough: it already carries a spring for
every one of AMBER's bonds and angles (measured: 100 % of a DNA duplex's 1368
bonds and 2484 angles, 99.8 % of ubiquitin's with no extra spring at all), and
a uniform stiffness reproduces the torsion wells to about a degree. The
furanose is the exception, for a structural reason: a five-ring's pucker is not
a torsion about a bond, it is a collective coordinate fixed by the ring's own
bonds and angles, so no torsion term can hold it and the mesh is its only
support.

A 1-3 distance spring cannot stand in for the angle term there. Both give the
same minimum -- on an isolated ring, bonds agree to 15 mA and angles to
0.59 deg -- but not the same curvature: at equilibrium each term contributes
k.(dq/dx)(dq/dx)^T to the Hessian, and dtheta/dr is perpendicular to the bond
while dd/dr lies along the 1-3 line. Two rank-one matrices pointing differently.
On one isolated angle that is invisible; summed around a closed ring the soft
collective modes come out wrong -- +73 % on the softest, which is the pucker,
and -52 % on the next.

  STRETCH <res> <a1> <a2> <r0> <k>            r0 in A, k in kJ.mol-1.A-2
  ANGLE   <res> <a1> <a2> <a3> <theta0> <k>   theta0 in deg, k in kJ.mol-1.rad-2

Both in the convention BioSpring and OpenMM share, E = 1/2.k.(q-q0)^2, so the
XML's constants are used exactly as read. STRETCH needs no new force -- a 1-2
mesh spring already is AMBER's bond term, it only lacks its constant. ANGLE is
the one genuinely new thing, and confining it to the sugar is what keeps it at
25 terms per nucleotide instead of one per valence angle in the system: on the
40-mer duplex, 1000 angle terms against the 2484 the molecule has.

  generate_ring_bonded.py DNA|RNA [out.bi.ff]
"""
import os
import sys

sys.path.insert(0, "scripts")
import generate_nucleic_forcefield as N  # noqa: E402

KIND = (sys.argv[1] if len(sys.argv) > 1 else "DNA").upper()
OUT = sys.argv[2] if len(sys.argv) > 2 else f"{KIND}AtomRing.bi.ff"
XML = {"DNA": ("amber14", "DNA.OL15.xml"), "RNA": ("amber14", "RNA.OL3.xml")}[KIND]
RESIDUES = {"DNA": ["DA", "DC", "DG", "DT"], "RNA": ["A", "C", "G", "U"]}[KIND]

# The furanose, as a cycle: consecutive pairs are its bonds.
RING = ["C1'", "C2'", "C3'", "C4'", "O4'"]
BONDS = [(RING[i], RING[(i + 1) % 5]) for i in range(5)]

# ANGLES are every angle whose VERTEX is a ring atom, exocyclic arms included,
# not only the five angles that lie entirely inside the ring. Measured, on the
# B-DNA duplex at k = 500 kJ.mol-1.A-2: the five ring angles alone take the
# median torsion error from 2.48 to 0.99 deg, and the vertex rule takes it to
# 0.48. The difference is delta, which is a readout of the pucker AND of where
# C5' and O3' sit -- 15.67 deg with the mesh alone, 6.89 with the ring angles,
# 0.49 with the vertex rule.
#
# BONDS deliberately stay inside the ring. AMBER's ring bonds are heavy-heavy
# (2594 kJ.mol-1.A-2 at a reduced mass of 6 Da, period 30 fs, not limiting),
# but giving an X-H bond its real constant -- 4628 at 0.95 Da, period 9 fs --
# would cap the timestep near 1.8 fs, which is the whole thing being avoided.

tables = N.ForceFieldTables(os.path.join(N.DATA, *XML))

records, seen, missing = [], set(), []
for resname in RESIDUES:
    ids = tables.atoms_of(resname)
    if not all(n in ids for n in RING):
        print(f"  {resname}: furanose incomplet, ignore")
        continue

    internal, _external = tables.bonds_of(resname)
    adj = {}
    for (x, y) in internal:
        adj.setdefault(x, set()).add(y)
        adj.setdefault(y, set()).add(x)
    angles = []
    for v in RING:
        ns = sorted(adj.get(v, ()))
        for i in range(len(ns)):
            for j in range(i + 1, len(ns)):
                angles.append((ns[i], v, ns[j]))

    nb = na = 0
    for (a, b) in BONDS:
        p = tables.bond_params.get(frozenset((ids[a], ids[b])))
        if p is None:
            missing.append(f"STRETCH {resname} {a}-{b}")
            continue
        length_nm, k_nm = p
        key = (resname,) + tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        # nm -> A, kJ.mol-1.nm-2 -> kJ.mol-1.A-2
        records.append(f"STRETCH\t{resname}\t{a}\t{b}\t{length_nm * 10.0:.4f}\t{k_nm / 100.0:.2f}")
        nb += 1
    for (a, v, c) in angles:
        p = tables.angle_params.get((ids[v], frozenset((ids[a], ids[c]))))
        if p is None:
            missing.append(f"ANGLE {resname} {a}-{v}-{c}")
            continue
        theta_deg, k_rad = p
        ends = tuple(sorted((a, c)))
        key = (resname, ends[0], v, ends[1])
        if key in seen:
            continue
        seen.add(key)
        records.append(f"ANGLE\t{resname}\t{a}\t{v}\t{c}\t{theta_deg:.4f}\t{k_rad:.2f}")
        na += 1
    print(f"  {resname}: {nb} liaisons, {na} angles")

header = [
    f"# AMBER's {KIND} furanose ring, as real bonded terms. Generated by",
    "# generate_ring_bonded.py -- do not edit by hand.",
    "#",
    "# STRETCH <res> <a1> <a2> <r0 A> <k kJ.mol-1.A-2>",
    "# ANGLE   <res> <a1> <a2> <a3> <theta0 deg> <k kJ.mol-1.rad-2>",
    "#",
    "# E = 1/2.k.(q-q0)^2, the convention BioSpring and OpenMM share.",
    "#",
    "# The ring's five bonds, and every angle centred on a ring atom. Every other",
    "# bond and angle in the molecule keeps its mesh spring at the uniform",
    "# stiffness -- the pucker is the one coordinate no torsion term can hold.",
]
open(OUT, "w").write("\n".join(header + records) + "\n")
if missing:
    print("  sans parametre AMBER : " + ", ".join(missing))
print(f"\n  {len(records)} enregistrements -> {OUT}")
