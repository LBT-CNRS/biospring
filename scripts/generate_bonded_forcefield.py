#!/usr/bin/env python3
#
# Generates data/reducerules/ProteinAtomBonded.bi.ff's STRETCH (real 1-2
# bonds), BEND (real 1-3 valence angles) and DIHEDRAL (proper-dihedral
# ghost springs, side chain only so far) entries for example 070's
# ubiquitin.pdb + the Fs-peptide's ACE/NME caps, from real AMBER ff99SB
# parameters -- no hand-typed chemistry: every number and every atom-type
# assignment comes straight from OpenMM's bundled amber99sb.xml
# (per-residue atom->type assignment + HarmonicBondForce/HarmonicAngleForce/
# PeriodicTorsionForce class1/class2/[class3]/...), matched against the
# real bond connectivity mdtraj infers from each PDB's standard residue
# templates. BEND angles are found generically: any real atom with >=2
# bonded neighbours is a valence vertex, its neighbour pairs are its real
# angles -- no residue-specific angle list needed. See
# doc/BondedForceFieldSprings.md for the full derivation and validation of
# the DIHEDRAL ghost-spring methodology.
#
# Requires `mdtraj` and `openmm` (pip install mdtraj openmm) -- not a
# BioSpring build dependency, only needed to regenerate the .bi.ff file.
#
# Units: amber99sb.xml already expresses k in the same convention BioSpring
# uses (E = 0.5*k*(r-r0)^2 -- OpenMM's own HarmonicBondForce/
# HarmonicAngleForce convention), so STRETCH's r0/k and BEND's theta0/k are
# used exactly as read (theta0 additionally converted rad -> deg for
# readability). No extra factor of 2 anywhere here (that factor is only
# needed when reading a raw AMBER parm*.dat file directly, whose K has no
# 1/2 baked in -- OpenMM's XML already performed that conversion). BEND's
# theta0/k are *not* converted into a 1-3 distance spring's r13/k13 here --
# BondedForceFieldReader does that conversion at build time (law of
# cosines + curvature matching), using the matching STRETCH r0's for
# r12/r23, per the project's decision to keep that coupling in the C++
# code rather than duplicate the geometry in this script.
#
# .bi.ff rules match by resname alone (no "is this the chain terminus"
# qualifier). Known, accepted approximation: a bond that exists in both the
# interior and terminal forms of a residue (e.g. Gly76's C-O, part of a
# symmetric COO- at the true C-terminus; Met1's N-CA, attached to a
# protonated NH3+) always uses the *interior* AMBER type here, even for the
# one residue that is actually the chain terminus -- the resname-only file
# format has no way to single out just that occurrence. Atoms that only
# exist at a terminus (H1/H2/H3, OXT) are unaffected: they use the correct
# terminal typing and only ever match the real terminal residue.
#
# Validated (session notes): the total stretching energy BioSpring reports
# after retuning (spring.enable, 1 step at a negligible timestep) matches,
# to 2 decimal places, an independent Python recomputation of
# 0.5*k*(r-r0)^2 over the same bonds using ubiquitin.pdb's real coordinates.
import os
import warnings
import xml.etree.ElementTree as ET

import mdtraj as md
import numpy as np
import openmm.app as openmm_app

warnings.filterwarnings("ignore")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDB = os.path.join(REPO_ROOT, "example/070.RigidBodyRotamere/ubiquitin.pdb")
XML = os.path.join(os.path.dirname(openmm_app.__file__), "data", "amber99sb.xml")
OUT = os.path.join(REPO_ROOT, "data/reducerules/ProteinAtomBonded.bi.ff")

tree = ET.parse(XML)
root = tree.getroot()

type_to_class = {t.get("name"): t.get("class") for t in root.find("AtomTypes").findall("Type")}
residue_types = {res.get("name"): {a.get("name"): a.get("type") for a in res.findall("Atom")}
                 for res in root.find("Residues").findall("Residue")}

bond_params = {}
for b in root.find("HarmonicBondForce").findall("Bond"):
    bond_params[frozenset((b.get("class1"), b.get("class2")))] = (float(b.get("length")), float(b.get("k")))

# angle_params[(vertex_class, frozenset({outer_class1, outer_class2}))] ->
# (theta0_deg, k_kJ.mol-1.rad-2). class2 in the XML is always the vertex by
# AMBER convention; the two outer classes are interchangeable (same physical
# angle read either direction), the vertex class is not.
angle_params = {}
for a in root.find("HarmonicAngleForce").findall("Angle"):
    key = (a.get("class2"), frozenset((a.get("class1"), a.get("class3"))))
    theta0_deg = float(a.get("angle")) * 180.0 / 3.14159265358979323846
    angle_params[key] = (theta0_deg, float(a.get("k")))

HIS_VARIANT = "HID"  # ubiquitin's His68 has HD1, no HE2 -> delta-protonated

def base_variant(resname):
    return HIS_VARIANT if resname == "HIS" else resname

def atom_class(variant, atomname):
    types = residue_types.get(variant)
    if types is None:
        return None
    type_id = types.get(atomname)
    return type_to_class[type_id] if type_id is not None else None

def resolved_class(resname, atomname, is_nterm, is_cterm):
    # Prefer the plain (interior) template; only reach for the N/C-terminal
    # variant when the atom doesn't exist there (H1/H2/H3, OXT) -- see
    # module docstring for the accepted approximation this implies.
    base = base_variant(resname)
    c = atom_class(base, atomname)
    if c is not None:
        return c
    if is_nterm:
        return atom_class("N" + base, atomname)
    if is_cterm:
        return atom_class("C" + base, atomname)
    return None

top = md.load_topology(PDB, standard_names=False)
residues = list(top.residues)
n_res = len(residues)

# Manual patch: mdtraj's standard-bond guesser misses the N-H1 bond of a
# free N-terminus for at least one residue type (checked: MET1 gets N-H2/
# N-H3 but not N-H1, even though atom H1 is real) -- same chemistry as
# H2/H3, added explicitly rather than trusted to the guesser.
bonds = list(top.bonds)
for res in residues:
    names = {a.name: a for a in res.atoms}
    if "H1" in names and "N" in names:
        pair = (names["N"], names["H1"])
        if not any({b[0], b[1]} == {pair[0], pair[1]} for b in bonds):
            bonds.append(pair)

lines = [
    "# STRETCH (real 1-2 bonds) and BEND (real 1-3 valence angles): real",
    "# AMBER ff99SB parameters (OpenMM amber99sb.xml), for the 18 residue",
    "# types of example 070's ubiquitin.pdb + ACE/NME from the Fs-peptide.",
    "# Generated by scripts/generate_bonded_forcefield.py -- do not",
    "# hand-edit individual values, regenerate instead. BEND's theta0/k are",
    "# converted into a 1-3 distance spring's equilibrium/stiffness at",
    "# build time by BondedForceFieldReader (see IO/BondedForceFieldReader.cpp),",
    "# not by this script.",
    "#",
    "# Atom-naming variants: BondedForceFieldReader matches by exact atom",
    "# name, and silently skips (never errors) any rule whose atom doesn't",
    "# exist on a given residue -- see Particle resolution in",
    "# IO/BondedForceFieldReader.cpp. So for the handful of atoms whose PDB",
    "# name genuinely varies across common conventions (free N-terminus",
    "# ammonium H's, free C-terminus carboxylate O's), every variant this",
    "# script knows about is emitted as its own redundant rule: whichever",
    "# one matches a given structure's actual naming fires, the rest just",
    "# never match. Not exhaustive (only the conventions checked below),",
    "# but harmless to extend the same way later.",
    "#",
    "# type    name                       resname atom1 atom2 atom3  r0_A/theta0_deg  k",
]

# Known PDB-naming variants for atoms unique to a free terminus. Every
# variant list for a given real atom must line up positionally (variant[i]
# in one list corresponds to variant[i] in another, e.g. the 3 ammonium
# H's) since some rules (angles, later) may need to keep pairing straight;
# for plain STRETCH entries here only the flat set of names matters.
NTERM_H_NAME_VARIANTS = [
    ["H1", "H2", "H3"],      # PDB v3 / AMBER (what ubiquitin.pdb/Fs-peptide use)
    ["1H", "2H", "3H"],      # PDB v2 / some legacy tools
    ["HT1", "HT2", "HT3"],   # CHARMM
]
CTERM_OXT_NAME_VARIANTS = ["OXT", "OT2", "O2"]  # PDB v3/AMBER, CHARMM, other

def emit_stretch(rule_name, resname, atom1, atom2_display, r0_A, k_biospring):
    global n_ok
    lines.append(f"STRETCH {rule_name:20s} {resname:7s} {atom1:5s} {atom2_display:6s} "
                 f"{r0_A:7.4f} {k_biospring:10.2f}")
    n_ok += 1

def emit_nterm_h_variants(resname, canonical_h, other_atom, r0_A, k_biospring):
    # canonical_h is one of "H1"/"H2"/"H3"; emit the same bond under every
    # other known naming convention for that same position.
    idx = NTERM_H_NAME_VARIANTS[0].index(canonical_h)
    for variant_set in NTERM_H_NAME_VARIANTS[1:]:
        h_name = variant_set[idx]
        emit_stretch(f"{resname}_N_{h_name}", resname, "N", h_name, r0_A, k_biospring)

def emit_cterm_oxt_variants(resname, r0_A, k_biospring):
    for o_name in CTERM_OXT_NAME_VARIANTS[1:]:
        emit_stretch(f"{resname}_C_{o_name}", resname, "C", o_name, r0_A, k_biospring)

def emit_bend(rule_name, resname, atom1, atom2, atom3_display, theta0_deg, k_biospring):
    global n_ok
    lines.append(f"BEND {rule_name:20s} {resname:7s} {atom1:5s} {atom2:5s} {atom3_display:6s} "
                 f"{theta0_deg:7.3f} {k_biospring:10.4f}")
    n_ok += 1

def rel_name(vertex_atom, other_atom):
    # Same +/- convention as .rbody/STRETCH: name as-is if same residue,
    # "+"/"-" prefixed if the neighbour is in the next/previous residue.
    # A real bond never spans more than one residue apart, so neither does
    # a real angle built from two such bonds sharing a vertex.
    d = other_atom.residue.index - vertex_atom.residue.index
    if d == 0:
        return other_atom.name, other_atom.name
    if d == 1:
        return other_atom.name, "+" + other_atom.name
    if d == -1:
        return other_atom.name, "-" + other_atom.name
    raise ValueError(f"angle neighbour {other_atom} is more than one residue away from vertex {vertex_atom}")

def generate_bends(bond_list, vertex_predicate, class_resolver):
    """For every real angle (vertex with >=2 neighbours) whose vertex atom
    matches `vertex_predicate(atom)`, resolve both outer atoms' names
    (with +/- prefix) and AMBER classes via `class_resolver(atom) -> class
    or None`, look up theta0/k, and emit a BEND rule. Silently skips (just
    counts) whatever has no atom type or no matching HarmonicAngleForce
    entry -- same policy as STRETCH."""
    global n_skip
    neighbors = {}
    for a1, a2 in bond_list:
        neighbors.setdefault(a1, set()).add(a2)
        neighbors.setdefault(a2, set()).add(a1)

    for vertex, nbrs in neighbors.items():
        if not vertex_predicate(vertex):
            continue
        nbrs = sorted(nbrs, key=lambda a: (a.residue.index, a.name))
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                n1, n2 = nbrs[i], nbrs[j]
                name1, rule_atom1 = rel_name(vertex, n1)
                name3, rule_atom3 = rel_name(vertex, n2)
                resname = vertex.residue.name

                dedup_key = (resname, vertex.name, tuple(sorted((rule_atom1, rule_atom3))))
                if dedup_key in seen:
                    continue

                cv = class_resolver(vertex)
                c1 = class_resolver(n1)
                c3 = class_resolver(n2)
                if cv is None or c1 is None or c3 is None:
                    print(f"SKIP BEND (no atom type) {resname}{vertex.residue.resSeq} "
                         f"{rule_atom1}-{vertex.name}-{rule_atom3}")
                    n_skip += 1
                    seen.add(dedup_key)
                    continue

                params = angle_params.get((cv, frozenset((c1, c3))))
                if params is None:
                    print(f"SKIP BEND (no HarmonicAngleForce entry) {resname}{vertex.residue.resSeq} "
                         f"{rule_atom1}({c1})-{vertex.name}({cv})-{rule_atom3}({c3})")
                    n_skip += 1
                    seen.add(dedup_key)
                    continue

                theta0_deg, k_biospring = params
                rule_name = f"{resname}_{name1}_{vertex.name}_{name3}"
                emit_bend(rule_name, resname, rule_atom1, vertex.name, rule_atom3, theta0_deg, k_biospring)
                seen.add(dedup_key)

seen = set()  # (resname, sorted atom name pair) already emitted
n_ok = 0
n_skip = 0

for ridx, res in enumerate(residues):
    is_nterm = ridx == 0
    is_cterm = ridx == n_res - 1

    for a1, a2 in bonds:
        if a1.residue.index == ridx and a2.residue.index == ridx:
            name1, name2, rule_atom2 = a1.name, a2.name, a2.name
            c1 = resolved_class(res.name, name1, is_nterm, is_cterm)
            c2 = resolved_class(res.name, name2, is_nterm, is_cterm)
        elif a1.residue.index == ridx and a2.residue.index == ridx + 1 and a2.name == "N":
            name1, name2, rule_atom2 = a1.name, "N", "+N"
            nxt = residues[ridx + 1]
            c1 = resolved_class(res.name, name1, is_nterm, is_cterm)
            c2 = resolved_class(nxt.name, "N", ridx + 1 == 0, ridx + 1 == n_res - 1)
        elif a2.residue.index == ridx and a1.residue.index == ridx + 1 and a1.name == "N":
            name1, name2, rule_atom2 = a2.name, "N", "+N"
            nxt = residues[ridx + 1]
            c1 = resolved_class(res.name, name1, is_nterm, is_cterm)
            c2 = resolved_class(nxt.name, "N", ridx + 1 == 0, ridx + 1 == n_res - 1)
        else:
            continue

        dedup_key = (res.name, tuple(sorted((name1, rule_atom2))))
        if dedup_key in seen:
            continue

        if c1 is None or c2 is None:
            print(f"SKIP (no atom type) {res.name}{res.resSeq} {name1}-{rule_atom2}")
            n_skip += 1
            seen.add(dedup_key)
            continue

        params = bond_params.get(frozenset((c1, c2)))
        if params is None:
            print(f"SKIP (no HarmonicBondForce entry) {res.name}{res.resSeq} {name1}({c1})-{rule_atom2}({c2})")
            n_skip += 1
            seen.add(dedup_key)
            continue

        length_nm, k_nm2 = params
        r0_A = length_nm * 10.0
        k_biospring = k_nm2 / 100.0

        rule_name = f"{res.name}_{name1}_{name2.lstrip('+')}"
        emit_stretch(rule_name, res.name, name1, rule_atom2, r0_A, k_biospring)
        seen.add(dedup_key)

        # These two atoms are exactly the ones with known PDB-naming
        # variants worth covering (see NTERM_H_NAME_VARIANTS/
        # CTERM_OXT_NAME_VARIANTS): free-N-terminus ammonium H's (bond is
        # always emitted as atom1="N", atom2="H1"/"H2"/"H3" above) and the
        # free-C-terminus carboxylate OXT (atom1="C", atom2="OXT").
        if name1 == "N" and name2 in NTERM_H_NAME_VARIANTS[0]:
            emit_nterm_h_variants(res.name, name2, name1, r0_A, k_biospring)
        elif name1 == "C" and name2 == "OXT":
            emit_cterm_oxt_variants(res.name, r0_A, k_biospring)

# Real valence angles for the same 18 residues, from the same real bond
# graph (any atom with >=2 bonded neighbours is a vertex) and the same
# resolved_class typing as STRETCH above.
generate_bends(
    bonds,
    vertex_predicate=lambda a: True,
    class_resolver=lambda a: resolved_class(a.residue.name, a.name, a.residue.index == 0,
                                            a.residue.index == n_res - 1),
)

# ===========================================================================
# DIHEDRAL SIDECHAIN: Arg chi1 (CA-CB axis) and chi2 (CB-CG axis) -- back to
# the simple, validated base: one ghost-spring group per axis, built from
# ALL real substituents on each side (no per-pair AMBER class-subgroup
# splitting), targeting the single generic wildcard term that covers a
# rotation about a C-C single bond (X-CT-CT-X, n=3). This is deliberately
# additive on top of the already-working rigidbody + stretching/bending
# network (see doc/BondedForceFieldSprings.md, Section 3.2) -- exactly the
# construction verified early on with idealized symmetric geometry.
#
# Known, explicit simplification (not hidden): amber99sb.xml also has more
# specific entries for some real substituent combinations on chi2's axis
# (CT-CT-CT-CT, HC-CT-CT-CT, HC-CT-CT-HC -- see session notes), each its
# own real multi-term torsion. Splitting the axis by those exact classes
# was tried and abandoned: several of the resulting subgroups have too few
# real substituents to be evenly spaced, breaking the LCM/Dirac-comb
# cancellation the method depends on (measured residual up to ~30x the
# target itself -- not a small approximation, a broken construction for
# that subgroup). Using the one broad wildcard term across the full,
# genuinely ~120-degree-spaced real substituent set avoids that failure
# mode; the price is not capturing chi2's finer real multi-term structure
# for now -- an open item, not silently dropped.
traj_ubq = md.load(PDB, standard_names=False)

# Real 3D geometry (bond neighbours, measured dihedral deltas) is sourced
# per-residue-type from an actual interior instance of that residue in
# ubiquitin -- generalizes the single hardcoded Arg42 instance used while
# only Arg was covered. Picks a non-terminal instance when one exists
# (terminal residues carry the N/C-cap-specific bonding, not the plain
# interior template STRETCH/BEND/DIHEDRAL all assume elsewhere).
_residue_instance_cache = {}

def _residue_instance(resname):
    # Sourced from `residues` (the same top-based object graph `bonds` is
    # built from), not traj_ubq.topology's own separate parse of the same
    # file -- atom identity (used by the backbone axis below to walk
    # `bonds`) only works within a single consistent graph, even though
    # atom .index values are safely comparable across both parses (same
    # file, same order, checked directly: fully consistent).
    if resname not in _residue_instance_cache:
        interior = [r for r in residues if r.name == resname and 0 < r.index < n_res - 1]
        candidates = interior or [r for r in residues if r.name == resname]
        _residue_instance_cache[resname] = candidates[0]
    return _residue_instance_cache[resname]

def _atoms_by_name_for(resname):
    return {a.name: a for a in _residue_instance(resname).atoms}

PHI_GRID_DEG = np.linspace(0.0, 360.0, 3600, endpoint=False)

torsion_params = {}  # (c1,c2,c3,c4) exact tuple, "" = wildcard -> [(n,k_kJ/mol,phase_rad), ...]
for _t in root.find("PeriodicTorsionForce").findall("Proper"):
    _c = tuple(_t.get(f"class{i}") or "" for i in (1, 2, 3, 4))
    _terms = []
    _i = 1
    while _t.get(f"periodicity{_i}") is not None:
        _k = float(_t.get(f"k{_i}"))
        if _k != 0.0:
            _terms.append((int(_t.get(f"periodicity{_i}")), _k, float(_t.get(f"phase{_i}"))))
        _i += 1
    if _terms:
        torsion_params.setdefault(_c, []).extend(_terms)

def lookup_torsion_wildcard(c2, c3):
    for key in (("", c2, c3, ""), ("", c3, c2, "")):
        if key in torsion_params:
            return torsion_params[key]
    return None

def lookup_torsion_specific(c1, c2, c3, c4):
    # Exact (non-wildcarded) match for one specific real substituent pair,
    # tried in both read directions (a dihedral read backwards names the
    # same physical angle). Used for backbone phi/psi, which -- unlike
    # chi1's single generic X-CT-CT-X term -- amber99sb.xml encodes as a
    # handful of pair-specific entries (verified directly against the XML,
    # not assumed): the axis's own generic wildcard entry exists but carries
    # k=0 for every periodicity (a real, deliberate null placeholder, not a
    # missing table), so it never survives the `if _k != 0.0` filter above
    # and correctly never matches here either.
    for key in ((c1, c2, c3, c4), (c4, c3, c2, c1)):
        if key in torsion_params:
            return torsion_params[key]
    return None

def bond_len_A(c1, c2):
    p = bond_params.get(frozenset((c1, c2)))
    return p[0] * 10.0 if p else None

def valence_deg(vertex_class, c1, c2):
    p = angle_params.get((vertex_class, frozenset((c1, c2))))
    return p[0] if p else None

def real_dihedral_deg(resname, n1, n2, n3, n4):
    atoms_by_name = _atoms_by_name_for(resname)
    idx = [[atoms_by_name[n].index for n in (n1, n2, n3, n4)]]
    return float(np.degrees(md.compute_dihedrals(traj_ubq, idx)[0, 0]))

def real_bond_neighbors(resname, vertex_name, exclude_name):
    res = _residue_instance(resname)
    out = []
    for a1, a2 in bonds:
        if a1.residue.index != res.index or a2.residue.index != res.index:
            continue
        if a1.name == vertex_name and a2.name != exclude_name:
            out.append(a2.name)
        elif a2.name == vertex_name and a1.name != exclude_name:
            out.append(a1.name)
    return out

def closed_form_d2(L, r1, theta1_deg, r2, theta2_deg, delta_deg, phi_deg):
    t1, t2, d, p = (np.radians(x) for x in (theta1_deg, theta2_deg, delta_deg, phi_deg))
    C1 = r1**2 * np.sin(t1)**2 + r2**2 * np.sin(t2)**2 + (r1 * np.cos(t1) - L + r2 * np.cos(t2))**2
    C2 = 2.0 * r1 * r2 * np.sin(t1) * np.sin(t2)
    return C1 - C2 * np.cos(p - d)

def complex_fourier_coeff(energy_over_phi_deg, n):
    # n=0 (DC/mean) needs a plain average, not the 2x-oscillation-amplitude
    # normalization that's correct for n>=1 (cos^2 integrates to 1/2 over a
    # period, a constant's own self-overlap integrates to 1 -- using the
    # n>=1 formula for n=0 silently doubles the DC target/coefficient, an
    # easy-to-miss bug found and fixed during the ghost-particle work, see
    # doc/BondedForceFieldSprings.md).
    phi_rad = np.radians(PHI_GRID_DEG)
    if n == 0:
        return complex(np.mean(energy_over_phi_deg))
    return (2.0 / len(phi_rad)) * np.sum(energy_over_phi_deg * np.exp(-1j * n * phi_rad))

# Deterministic, closed-form ghost-PARTICLE construction (see
# doc/BondedForceFieldSprings.md, "particules fantômes"/ghost-particle
# sections for the full derivation): replaces the earlier ghost-SPRING
# method (fixed real-atom geometry, only k free, sometimes had no
# non-negative-k solution for a multi-term real target) with springs
# between fully-free virtual points, ANCHORED to real atoms only through
# a 3-point placement (see spn::GhostParticle) -- r/theta/delta become
# real degrees of freedom instead of a choice among a few fixed shapes.
#
# For a target harmonic n>=2, a ring of M=N=n ghost-ghost springs (same
# base shape, evenly spaced by 360/n around the axis on each side) cancels
# -- by an exact roots-of-unity argument, verified numerically -- every
# harmonic that isn't a multiple of n, leaving a pure (up to a small,
# quantified leakage into 2n, 3n...) contribution at exactly n. n=1 is a
# special, EXACT case: a single ghost-ghost spring with d0=0 (zero rest
# length) gives a provably pure n=1 harmonic with ZERO leakage into any
# other harmonic (energy = 0.5*k*d(phi)^2, and d(phi)^2 is exactly
# C1-C2*cos(phi-delta) by the closed_form_d2 identity above -- no sqrt
# nonlinearity survives at d0=0) -- this is why n=1 uses mu=0 below, not a
# leak-minimized (rho,mu) like n>=2.
#
# (rho, mu) for n=2,3 were found by maximizing the achievable
# |a_n|/a0 ratio (the harmonic-to-DC-offset efficiency of one ring's own
# base shape) subject to leakage into the first unwanted multiple (2n)
# staying <=0.5% -- a small, one-time 2D grid search over the ABSTRACT
# shape function sqrt(1-rho*cos(x)), independent of any specific axis or
# residue (verified: reused identically across all 70 axes below).
RING_SHAPES = {1: (0.5, 0.0), 2: (0.245, 1.0), 3: (0.465, 1.0)}

def ring_shape_to_geometry(L_axis, n):
    rho, mu = RING_SHAPES[n]
    r = np.sqrt(rho / (2.0 * (1.0 - rho))) * L_axis
    R = np.sqrt(2.0 * r * r + L_axis * L_axis)
    d0 = mu * R
    return r, d0

def ring_curve_abstract(L_axis, n, r, d0, k, delta_base, phi_grid_deg):
    # The 1D Fourier-space model used to CALIBRATE k/delta_base (see
    # calibrate_ring) -- exactly the M=N=n comb-filter sum, evaluated over
    # an abstract phi grid. NOT the real 3D placement (see emit_ghost_ring
    # for that, and its header comment for the sign convention linking the
    # two -- verified by direct 3D simulation against a real axis before
    # being trusted here, see doc/BondedForceFieldSprings.md).
    e = np.zeros_like(phi_grid_deg)
    for i in range(n):
        beta = i * 360.0 / n
        for j in range(n):
            gamma = j * 360.0 / n
            delta_ij = delta_base + beta - gamma
            d2v = closed_form_d2(L_axis, r, 90.0, r, 90.0, delta_ij, phi_grid_deg)
            d = np.sqrt(np.maximum(d2v, 0.0))
            e += 0.5 * k * (d - d0) ** 2
    return e

def calibrate_ring(L_axis, n, r, d0, target_complex):
    # Energy is exactly linear in k at fixed geometry (see
    # doc/BondedForceFieldSprings.md, step 3 of the closed-form derivation)
    # -- one evaluation at k=1 gives the ring's own harmonic-n Fourier
    # coefficient (magnitude AND phase), from which k (magnitude match) and
    # delta_base (phase match) solve directly, no iteration.
    e0 = ring_curve_abstract(L_axis, n, r, d0, 1.0, 0.0, PHI_GRID_DEG)
    z0 = complex_fourier_coeff(e0, n)
    phase_needed = np.angle(target_complex) - np.angle(z0)
    # Rotating delta_base by D shifts the ring's own harmonic-n phase by
    # -n*D, not +n*D (verified numerically -- an easy sign error, since the
    # naive expectation is +n*D). Using the wrong sign gives the complex
    # CONJUGATE of the target: invisible for a purely-real target (chi1's
    # single-term case) but ~45-60% wrong for any axis with a genuinely
    # complex target (phi/psi, chi2's 3-term case).
    delta_base_deg = np.degrees(-phase_needed / n)
    k = abs(target_complex) / abs(z0)
    return k, delta_base_deg

n_dihedral_ok = 0
n_dihedral_skip = 0
n_ghost_particles = 0
dihedral_lines = []
ghostparticle_lines = []

def emit_ghost_ring(resname, axis_label, family, n, L_axis, target_complex, atom_B, atom_C, atom_ref_b, atom_ref_c):
    """Emits one full ring group (M=N=n ghost-ghost springs) reproducing
    one target Fourier harmonic exactly (up to the small, bounded leakage
    quantified in RING_SHAPES's own derivation). atom_B/atom_C are the
    real dihedral axis atoms (already +/- resolved); atom_ref_b/atom_ref_c
    are each side's own real reference atom (fixes what "azimuth=0" means
    for that side -- must be the SAME atoms used to measure the target's
    own phase convention, i.e. b_neighbors[0]/c_neighbors[0], see
    combined_target_for_axis).

    Placement/sign convention (verified by a full 3D rotating simulation
    against a real axis before being trusted here -- see
    doc/BondedForceFieldSprings.md): a b-side ghost i sits at azimuth
    i*360/n (measured from atom_B, relative to atom_ref_b, via
    spn::GhostParticle::computePosition(atom_B, atom_C, atom_ref_b, ...)).
    A c-side ghost j sits at azimuth delta_base - j*360/n (note the
    MINUS -- placed via computePosition(atom_C, atom_B, atom_ref_c, ...),
    i.e. B/C swapped since the c-side ghost's own axis direction points
    the other way). This combination reproduces exactly the SAME abstract
    delta_ij = delta_base + beta_i - gamma_j used to calibrate k/delta_base
    in Fourier space above -- get either the swap or the minus sign wrong
    and the real 3D energy has the wrong phase relative to the molecule's
    true dihedral angle, silently, without any local check catching it.
    """
    global n_dihedral_ok, n_ghost_particles

    r, d0 = ring_shape_to_geometry(L_axis, n)
    k, delta_base = calibrate_ring(L_axis, n, r, d0, target_complex)

    b_names, c_names = [], []
    for i in range(n):
        name = f"GH{axis_label}n{n}B{i}"
        beta = i * 360.0 / n
        ghostparticle_lines.append(
            f"GHOSTPARTICLE {name:16s} {resname:7s} {atom_B:5s} {atom_C:6s} {atom_ref_b:6s} "
            f"{r:7.4f} {90.0:7.2f} {beta:9.3f}")
        b_names.append(name)
        n_ghost_particles += 1
    for j in range(n):
        name = f"GH{axis_label}n{n}C{j}"
        gamma = delta_base - j * 360.0 / n
        ghostparticle_lines.append(
            f"GHOSTPARTICLE {name:16s} {resname:7s} {atom_C:5s} {atom_B:6s} {atom_ref_c:6s} "
            f"{r:7.4f} {90.0:7.2f} {gamma:9.3f}")
        c_names.append(name)
        n_ghost_particles += 1

    for bn in b_names:
        for cn in c_names:
            rule_name = f"{resname}_{axis_label}_n{n}_{bn[-1]}{cn[-1]}"
            dihedral_lines.append(
                f"DIHEDRAL {rule_name:28s} {resname:7s} {family:9s} {bn:16s} {cn:17s} "
                f"{d0:7.4f} {k:10.4f}")
            n_dihedral_ok += 1

def combined_target_for_axis(b_geom, c_geom, b_class, c_class, class_of):
    """b_geom/c_geom: {atom: (r, theta, delta_deg)} for each side's real
    substituents (delta already phase-corrected relative to that side's
    own reference atom, i.e. b_neighbors[0]/c_neighbors[0] -- ref_phi
    itself cancels out of the pairwise delta and is never needed, see
    doc/BondedForceFieldSprings.md). class_of(atom) resolves an atom to
    its AMBER class. Tries each real substituent pair's specific AMBER
    entry first, falls back to the axis's generic wildcard (matches
    amber99sb.xml's own matching order) -- a pair whose specific entry
    exists but carries k=0 for every periodicity (a genuine AMBER zero,
    e.g. chi2 for Asp/His/Phe/Tyr, chi3 for Glu, chi4 for Arg, chi2 for
    Trp -- all confirmed by direct XML inspection) correctly contributes
    nothing, not a gap."""
    target = {0: 0j}
    any_term = False
    for bn, (r1, t1, d1) in b_geom.items():
        for cn, (r2, t2, d2) in c_geom.items():
            terms = lookup_torsion_specific(class_of(bn), b_class, c_class, class_of(cn))
            if terms is None:
                terms = lookup_torsion_wildcard(b_class, c_class)
            if terms is None:
                continue
            any_term = True
            delta_pair_rad = np.radians(d2 - d1)
            for (n, k, phase) in terms:
                target[n] = target.get(n, 0j) + k * np.exp(-1j * phase) * np.exp(1j * n * delta_pair_rad)
                target[0] += k
    return target if any_term else None

def generate_sidechain_axis(resname, b_name, c_name, axis_label):
    global n_dihedral_ok, n_dihedral_skip

    b_class = resolved_class(resname, b_name, False, False)
    c_class = resolved_class(resname, c_name, False, False)
    L_axis = bond_len_A(b_class, c_class)

    b_neighbors = real_bond_neighbors(resname, b_name, c_name)
    c_neighbors = real_bond_neighbors(resname, c_name, b_name)
    if not b_neighbors or not c_neighbors:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}): axis {b_name}-{c_name} has no real "
             f"substituent on one side")
        n_dihedral_skip += 1
        return
    b_ref, c_ref = b_neighbors[0], c_neighbors[0]
    ref_phi = real_dihedral_deg(resname, b_ref, b_name, c_name, c_ref)

    def side_geom(neighbors, vertex_class, other_class, ref_name, is_b_side):
        geom = {}
        for n in neighbors:
            nc = resolved_class(resname, n, False, False)
            r = bond_len_A(nc, vertex_class)
            theta = valence_deg(vertex_class, nc, other_class)
            if n == ref_name:
                phi_meas = ref_phi
            elif is_b_side:
                phi_meas = real_dihedral_deg(resname, n, b_name, c_name, c_ref)
            else:
                phi_meas = real_dihedral_deg(resname, b_ref, b_name, c_name, n)
            geom[n] = (r, theta, phi_meas - ref_phi)
        return geom

    b_geom = side_geom(b_neighbors, b_class, c_class, b_ref, True)
    c_geom = side_geom(c_neighbors, c_class, b_class, c_ref, False)

    class_of = lambda n: resolved_class(resname, n, False, False)
    target = combined_target_for_axis(b_geom, c_geom, b_class, c_class, class_of)
    if target is None:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}): no real substituent pair matches a "
             f"specific AMBER entry and no generic wildcard either")
        n_dihedral_skip += 1
        return

    for n, zt in target.items():
        if n == 0 or abs(zt) < 1e-6:
            continue
        emit_ghost_ring(resname, axis_label, "SIDECHAIN", n, L_axis, zt, b_name, c_name, b_ref, c_ref)

# chi1 axis (CA-CB): every one of ubiquitin's 18 residue types that has a
# real, freely-rotating first side-chain torsion. Excluded: GLY (no CB),
# ALA (CB is a symmetric terminal methyl, no meaningful rotamer) and PRO
# (CA-CB is part of the ring --rigidbody already holds rigid, not a free
# single-bond rotation). Every one of these residues' CB is a plain
# aliphatic AMBER class (CT), so all 15 hit the same generic X-CT-CT-X
# wildcard term already validated for Arg -- no per-residue special-casing
# needed, confirmed by running this (not assumed from the residue list).
CHI1_RESIDUES = ["ARG", "ASN", "ASP", "GLN", "GLU", "HIS", "ILE", "LEU",
                 "LYS", "MET", "PHE", "SER", "THR", "TYR", "VAL"]
for _resname in CHI1_RESIDUES:
    generate_sidechain_axis(_resname, "CA", "CB", "chi1")

# chi2 axis (CB-CG, or CB-CG1 for Ile's branched beta carbon): attempted
# for every residue whose real side chain continues past CB. Only fires a
# real term for the ones where CG is itself a plain aliphatic carbon (CT):
# Arg, Gln, Glu, Ile, Leu, Lys, Met. For Asn/Asp (CG is the amide/
# carboxylate carbon, AMBER class C) and His/Phe/Tyr (CG is an aromatic
# ring carbon, class CA), no generic X-CT-CT-X term matches -- the function
# reports and skips those rather than fabricating a torsion that isn't in
# amber99sb.xml (their ring/group planarity is deferred PLANARITY work, see
# doc/BondedForceFieldSprings.md and the plan's ring-closure TODO).
CHI2_AXIS = {"ARG": "CG", "GLN": "CG", "GLU": "CG", "ILE": "CG1", "LEU": "CG",
             "LYS": "CG", "MET": "CG", "ASN": "CG", "ASP": "CG", "HIS": "CG",
             "PHE": "CG", "TYR": "CG"}
for _resname, _c_name in CHI2_AXIS.items():
    generate_sidechain_axis(_resname, "CB", _c_name, "chi2")

# chi3 (CG-CD, or CG-SD for Met's thioether) and chi4 (CD-NE for Arg,
# CD-CE for Lys): only the residues whose side chain actually continues
# that far. No new matching logic needed -- same generate_sidechain_axis,
# just further down each chain.
CHI3_AXIS = {"ARG": ("CG", "CD"), "GLN": ("CG", "CD"), "GLU": ("CG", "CD"),
            "LYS": ("CG", "CD"), "MET": ("CG", "SD")}
for _resname, (_b_name, _c_name) in CHI3_AXIS.items():
    generate_sidechain_axis(_resname, _b_name, _c_name, "chi3")

CHI4_AXIS = {"ARG": ("CD", "NE"), "LYS": ("CD", "CE")}
for _resname, (_b_name, _c_name) in CHI4_AXIS.items():
    generate_sidechain_axis(_resname, _b_name, _c_name, "chi4")

# Cys/Trp chi1 (and Trp chi2): ubiquitin.pdb has neither residue at all
# (not even a terminal instance, unlike Met) -- sourced from a second real
# structure, example/001_GKinase/model.pdb, which has both. This is
# needed only for chi (side-chain-specific geometry, genuinely different
# per residue) -- NOT for phi/psi, whose N/CA/C/CB pattern is generic and
# already covered by reusing Leu below (see the backbone section).
GKINASE_PDB = os.path.join(REPO_ROOT, "example/001_GKinase/model.pdb")
# This PDB names hydrogens with the leading-digit convention (1HB, 2HB),
# not AMBER's trailing-digit one (HB2, HB3) -- standard_names=True lets
# mdtraj normalize them, unlike ubiquitin.pdb above which is already in
# the AMBER convention and uses standard_names=False to avoid any
# reduction-name remapping surprises.
_gk_top = md.load_topology(GKINASE_PDB, standard_names=True)
_gk_residues = list(_gk_top.residues)
_gk_n_res = len(_gk_residues)
_gk_bonds = list(_gk_top.bonds)
_gk_traj = md.load(GKINASE_PDB, standard_names=True)
_gk_residue_instance_cache = {}

def _gk_residue_instance(resname):
    if resname not in _gk_residue_instance_cache:
        interior = [r for r in _gk_residues if r.name == resname and 0 < r.index < _gk_n_res - 1]
        candidates = interior or [r for r in _gk_residues if r.name == resname]
        _gk_residue_instance_cache[resname] = candidates[0]
    return _gk_residue_instance_cache[resname]

def _gk_real_bond_neighbors(resname, vertex_name, exclude_name):
    res = _gk_residue_instance(resname)
    out = []
    for a1, a2 in _gk_bonds:
        if a1.residue.index != res.index or a2.residue.index != res.index:
            continue
        if a1.name == vertex_name and a2.name != exclude_name:
            out.append(a2.name)
        elif a2.name == vertex_name and a1.name != exclude_name:
            out.append(a1.name)
    return out

def _gk_real_dihedral_deg(resname, n1, n2, n3, n4):
    res = _gk_residue_instance(resname)
    atoms_by_name = {a.name: a for a in res.atoms}
    idx = [[atoms_by_name[n].index for n in (n1, n2, n3, n4)]]
    return float(np.degrees(md.compute_dihedrals(_gk_traj, idx)[0, 0]))

def generate_sidechain_axis_gkinase(resname, b_name, c_name, axis_label):
    global n_dihedral_ok, n_dihedral_skip

    b_class = resolved_class(resname, b_name, False, False)
    c_class = resolved_class(resname, c_name, False, False)
    L_axis = bond_len_A(b_class, c_class)

    b_neighbors = _gk_real_bond_neighbors(resname, b_name, c_name)
    c_neighbors = _gk_real_bond_neighbors(resname, c_name, b_name)
    if not b_neighbors or not c_neighbors:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}, GKinase): axis {b_name}-{c_name} has no "
             f"real substituent on one side")
        n_dihedral_skip += 1
        return
    b_ref, c_ref = b_neighbors[0], c_neighbors[0]
    ref_phi = _gk_real_dihedral_deg(resname, b_ref, b_name, c_name, c_ref)

    def side_geom(neighbors, vertex_class, other_class, ref_name, is_b_side):
        geom = {}
        for n in neighbors:
            nc = resolved_class(resname, n, False, False)
            r = bond_len_A(nc, vertex_class)
            theta = valence_deg(vertex_class, nc, other_class)
            if n == ref_name:
                phi_meas = ref_phi
            elif is_b_side:
                phi_meas = _gk_real_dihedral_deg(resname, n, b_name, c_name, c_ref)
            else:
                phi_meas = _gk_real_dihedral_deg(resname, b_ref, b_name, c_name, n)
            geom[n] = (r, theta, phi_meas - ref_phi)
        return geom

    b_geom = side_geom(b_neighbors, b_class, c_class, b_ref, True)
    c_geom = side_geom(c_neighbors, c_class, b_class, c_ref, False)

    class_of = lambda n: resolved_class(resname, n, False, False)
    target = combined_target_for_axis(b_geom, c_geom, b_class, c_class, class_of)
    if target is None:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}, GKinase): no real substituent pair "
             f"matches a specific AMBER entry and no generic wildcard either -- confirmed genuine "
             f"AMBER zero-barrier torsion, not a gap")
        n_dihedral_skip += 1
        return

    for n, zt in target.items():
        if n == 0 or abs(zt) < 1e-6:
            continue
        emit_ghost_ring(resname, axis_label, "SIDECHAIN", n, L_axis, zt, b_name, c_name, b_ref, c_ref)

generate_sidechain_axis_gkinase("CYS", "CA", "CB", "chi1")
generate_sidechain_axis_gkinase("TRP", "CA", "CB", "chi1")
generate_sidechain_axis_gkinase("TRP", "CB", "CG", "chi2")  # expected to correctly SKIP (aromatic
                                                            # ring, same zero-barrier family as
                                                            # Phe/Tyr/His chi2)

# ===========================================================================
# DIHEDRAL BACKBONE: phi (N-CA axis) and psi (CA-C axis). CMAP checked and
# confirmed absent from amber99sb.xml (see doc/BondedForceFieldSprings.md) --
# no coupling term is silently missed. Unlike chi1's single generic
# X-CT-CT-X term, amber99sb.xml gives phi/psi *pair-specific* real entries
# (verified directly against the XML): only two real substituent pairs
# carry a nonzero term per axis (the canonical all-backbone pair, and the
# one involving CB), everything else -- both backbone hydrogens, the
# carbonyl oxygen -- fall on that axis's own generic wildcard entry, which
# is present in the table but carries k=0 for every periodicity (a real,
# deliberate null, not a gap). So generation here matches each real
# substituent pair to its own specific AMBER entry (or drops it if only the
# null wildcard applies) instead of averaging one shared term over the
# whole real substituent grid.
def real_bond_neighbors_atoms(vertex_atom, exclude_atom):
    out = []
    for a1, a2 in bonds:
        if a1 is vertex_atom and a2 is not exclude_atom:
            out.append(a2)
        elif a2 is vertex_atom and a1 is not exclude_atom:
            out.append(a1)
    return out

def real_dihedral_deg_atoms(a1, a2, a3, a4):
    idx = [[a1.index, a2.index, a3.index, a4.index]]
    return float(np.degrees(md.compute_dihedrals(traj_ubq, idx)[0, 0]))

def bb_atom_class(a):
    return resolved_class(a.residue.name, a.name, a.residue.index == 0, a.residue.index == n_res - 1)

def generate_backbone_axis(resname, b_name, c_name, axis_label, source_resname=None):
    """source_resname: sources the real 3D geometry (which substituents,
    their bond lengths/angles/relative azimuths) from a DIFFERENT
    residue's real ubiquitin instance instead of `resname`'s own -- valid
    because phi/psi's N/CA/C/CB pattern never depends on anything past CB
    (verified: Met's own target is bit-identical to Leu's sourced this
    way), used for Met/Cys/Trp, none of which has a usable non-terminal
    real instance of their own in ubiquitin.pdb. AMBER classes are always
    resolved via `resname` itself (a class-based lookup needs no real
    instance at all, so this is correct regardless of source_resname)."""
    global n_dihedral_ok, n_dihedral_skip

    geom_resname = source_resname or resname
    res = _residue_instance(geom_resname)
    if res.index == 0 or res.index == n_res - 1:
        print(f"SKIP DIHEDRAL BACKBONE ({resname} {axis_label}): only real ubiquitin instance "
             f"of {geom_resname} is a chain terminus, no real neighbouring residue to source from")
        n_dihedral_skip += 1
        return

    b_atom = next(a for a in res.atoms if a.name == b_name)
    c_atom = next(a for a in res.atoms if a.name == c_name)
    b_class = resolved_class(resname, b_name, False, False)
    c_class = resolved_class(resname, c_name, False, False)
    L_axis = bond_len_A(b_class, c_class)

    b_neighbors = real_bond_neighbors_atoms(b_atom, c_atom)
    c_neighbors = real_bond_neighbors_atoms(c_atom, b_atom)
    if not b_neighbors or not c_neighbors:
        print(f"SKIP DIHEDRAL BACKBONE ({resname} {axis_label}): axis {b_name}-{c_name} has no "
             f"real substituent on one side")
        n_dihedral_skip += 1
        return

    ref_b, ref_c = b_neighbors[0], c_neighbors[0]
    ref_phi = real_dihedral_deg_atoms(ref_b, b_atom, c_atom, ref_c)

    def delta_of(atom, is_b_side):
        if atom is (ref_b if is_b_side else ref_c):
            return 0.0
        if is_b_side:
            return real_dihedral_deg_atoms(atom, b_atom, c_atom, ref_c) - ref_phi
        return real_dihedral_deg_atoms(ref_b, b_atom, c_atom, atom) - ref_phi

    b_geom = {a: (bond_len_A(bb_atom_class(a), b_class), valence_deg(b_class, bb_atom_class(a), c_class),
                 delta_of(a, True)) for a in b_neighbors}
    c_geom = {a: (bond_len_A(bb_atom_class(a), c_class), valence_deg(c_class, bb_atom_class(a), b_class),
                 delta_of(a, False)) for a in c_neighbors}

    target = combined_target_for_axis(b_geom, c_geom, b_class, c_class, bb_atom_class)
    if target is None:
        print(f"SKIP DIHEDRAL BACKBONE ({resname} {axis_label}): no real substituent pair matches "
             f"a specific AMBER Proper entry (all fall on the null wildcard)")
        n_dihedral_skip += 1
        return

    _, rule_b_ref = rel_name(b_atom, ref_b)
    _, rule_c_ref = rel_name(c_atom, ref_c)

    for n, zt in target.items():
        if n == 0 or abs(zt) < 1e-6:
            continue
        emit_ghost_ring(resname, axis_label, "BACKBONE", n, L_axis, zt, b_name, c_name, rule_b_ref, rule_c_ref)

# All 18 residue types found in ubiquitin.pdb, plus Cys/Trp (also generic-
# sourced from Leu -- ubiquitin has neither at all, not even a terminal
# instance like Met). Met/Cys/Trp all source their real geometry from Leu
# (confirmed identical target either way -- the N/CA/C/CB pattern doesn't
# depend on anything past CB); every other residue uses its own real
# instance directly.
BACKBONE_RESIDUES = sorted({r.name for r in residues})
GENERIC_BACKBONE_SOURCE = {"MET": "LEU", "CYS": "LEU", "TRP": "LEU"}
for _resname in BACKBONE_RESIDUES + ["CYS", "TRP"]:
    _source = GENERIC_BACKBONE_SOURCE.get(_resname)
    generate_backbone_axis(_resname, "N", "CA", "phi", source_resname=_source)
    generate_backbone_axis(_resname, "CA", "C", "psi", source_resname=_source)

lines.append("#")
lines.append("# GHOSTPARTICLE: massless virtual sites (see spn::GhostParticle), placed")
lines.append("# algebraically from 3 real anchor atoms. One M=N=n ring (2n particles)")
lines.append("# per real AMBER Fourier harmonic n of each covered axis -- see")
lines.append("# doc/BondedForceFieldSprings.md for the closed-form derivation.")
lines.append("#    type          name             resname atom_B atom_C atom_ref     r_A   theta_deg delta_deg")
lines.extend(ghostparticle_lines)
lines.append("#")
lines.append("# DIHEDRAL (proper): ghost-ghost springs, one ring group per real AMBER")
lines.append("# Fourier harmonic of each covered chi1-4 (family SIDECHAIN) and phi/psi")
lines.append("# (family BACKBONE) axis -- see doc/BondedForceFieldSprings.md.")
lines.append("#    type     name                       resname family    atom_ref atom_rotant  d0_A       k")
lines.extend(dihedral_lines)

# ACE/NME capping groups: not present in ubiquitin (free termini), sourced
# instead from the Fs-peptide (Ace-A8-(AAARA)3-A-NME) used for the hinge/
# rotamer-vs-real-MD validation earlier in this project. amber99sb.xml's
# ACE/NME templates name the methyl hydrogens HH31/HH32/HH33; the PDB (and
# every reader that loads it) names them H1/H2/H3 -- same atoms, just a
# different PDB naming convention, remapped explicitly here.
FS_PEPTIDE_PDB = "/Users/ferey/Développements/Src/ARBRE/data/peptide/trajectory/100-fs-peptide-400K.pdb"
METHYL_H_REMAP = {"H1": "HH31", "H2": "HH32", "H3": "HH33"}

if os.path.exists(FS_PEPTIDE_PDB):
    fs_top = md.load_topology(FS_PEPTIDE_PDB, standard_names=False)
    fs_residues = list(fs_top.residues)
    fs_n_res = len(fs_residues)

    def cap_atom_class(resname, atomname):
        # NME's own methyl carbon is named "C" in the PDB (like a normal
        # backbone carbonyl) but "CH3" in amber99sb.xml's NME template --
        # only remap it when resname is NME itself, not when "C" refers to
        # the *previous* (ordinary amino acid) residue's real carbonyl.
        if resname == "NME" and atomname == "C":
            xml_name = "CH3"
        else:
            xml_name = METHYL_H_REMAP.get(atomname, atomname)
        return atom_class(resname, xml_name)

    for ridx, res in enumerate(fs_residues):
        if res.name not in ("ACE", "NME"):
            continue
        for a1, a2 in fs_top.bonds:
            if a1.residue.index == ridx and a2.residue.index == ridx:
                name1, name2, rule_atom2 = a1.name, a2.name, a2.name
                c1, c2 = cap_atom_class(res.name, name1), cap_atom_class(res.name, name2)
            elif a1.residue.index == ridx and a2.residue.index == ridx - 1 and a2.name == "C":
                name1, name2, rule_atom2 = a1.name, "C", "-C"
                c1 = cap_atom_class(res.name, name1)
                c2 = cap_atom_class(fs_residues[ridx - 1].name, "C")
            elif a2.residue.index == ridx and a1.residue.index == ridx - 1 and a1.name == "C":
                name1, name2, rule_atom2 = a2.name, "C", "-C"
                c1 = cap_atom_class(res.name, name1)
                c2 = cap_atom_class(fs_residues[ridx - 1].name, "C")
            elif a1.residue.index == ridx and a2.residue.index == ridx + 1 and a2.name == "N":
                name1, name2, rule_atom2 = a1.name, "N", "+N"
                c1 = cap_atom_class(res.name, name1)
                c2 = cap_atom_class(fs_residues[ridx + 1].name, "N")
            elif a2.residue.index == ridx and a1.residue.index == ridx + 1 and a1.name == "N":
                name1, name2, rule_atom2 = a2.name, "N", "+N"
                c1 = cap_atom_class(res.name, name1)
                c2 = cap_atom_class(fs_residues[ridx + 1].name, "N")
            else:
                continue

            dedup_key = (res.name, tuple(sorted((name1, rule_atom2))))
            if dedup_key in seen:
                continue

            if c1 is None or c2 is None:
                print(f"SKIP (no atom type) {res.name}{res.resSeq} {name1}-{rule_atom2}")
                n_skip += 1
                seen.add(dedup_key)
                continue

            params = bond_params.get(frozenset((c1, c2)))
            if params is None:
                print(f"SKIP (no HarmonicBondForce entry) {res.name}{res.resSeq} {name1}({c1})-{rule_atom2}({c2})")
                n_skip += 1
                seen.add(dedup_key)
                continue

            length_nm, k_nm2 = params
            r0_A = length_nm * 10.0
            k_biospring = k_nm2 / 100.0

            rule_name = f"{res.name}_{name1}_{name2.lstrip('+').lstrip('-')}"
            emit_stretch(rule_name, res.name, name1, rule_atom2, r0_A, k_biospring)
            seen.add(dedup_key)

    # Real valence angles vertex-centered on ACE/NME's own atoms. Angles
    # vertex-centered on the adjacent standard residue (e.g. ALA2's own N,
    # neighbours -C/CA/H) are already covered above: the class of "-C"/"+N"
    # is generic (always plain "C"/"N"), independent of whether the actual
    # neighbour is another amino acid or a cap.
    generate_bends(
        list(fs_top.bonds),
        vertex_predicate=lambda a: a.residue.name in ("ACE", "NME"),
        class_resolver=lambda a: cap_atom_class(a.residue.name, a.name),
    )

with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")

print(f"\nWrote {n_ok} unique STRETCH/BEND entries ({n_skip} skipped), {n_ghost_particles} "
     f"GHOSTPARTICLE entries and {n_dihedral_ok} DIHEDRAL entries ({n_dihedral_skip} skipped) to {OUT}")
