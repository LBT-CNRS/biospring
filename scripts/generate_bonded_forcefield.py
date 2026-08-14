#!/usr/bin/env python3
#
# Generates data/reducerules/ProteinAtomBonded.bi.ff's STRETCH (real 1-2
# bonds), BEND (real 1-3 valence angles, converted to a ghost-particle
# distance spring at build time -- see below), GHOSTPARTICLE (virtual
# sites, spn::GhostParticle) and DIHEDRAL (proper-dihedral ghost-particle
# rings: chi1-4, phi/psi, and omega -- the peptide-bond torsion, the only
# axis whose own two axis atoms span a residue boundary) entries, from real
# AMBER ff99SB parameters -- no hand-typed chemistry: every number and
# every atom-type assignment comes straight from OpenMM's bundled
# amber99sb.xml (per-residue atom->type assignment +
# HarmonicBondForce/HarmonicAngleForce/PeriodicTorsionForce
# class1/class2/[class3]/...), matched against the real bond connectivity
# mdtraj infers from each PDB's standard residue templates (ubiquitin.pdb
# for most residues, example/001_GKinase/model.pdb for Cys/Trp chi1, which
# ubiquitin lacks entirely). BEND angles are found generically: any real
# atom with >=2 bonded neighbours is a valence vertex, its neighbour pairs
# are its real angles -- no residue-specific angle list needed. See
# doc/BondedForceFieldSprings.md for the full derivation and validation of
# the ghost-particle methodology (why ghost PARTICLES rather than springs
# directly between real atoms: a real-atom spring for BEND leaks bond-
# stretch noise, and for DIHEDRAL cannot always reach a real AMBER
# multi-term target without an unphysical negative stiffness).
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
import re
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

def idealized_azimuth_deg(vertex_class, ref_class, x_class, y_class):
    """Azimuthal offset of substituent y relative to x, viewed down the
    vertex->ref axis, derived PURELY from AMBER's own HarmonicAngleForce
    equilibrium angles at `vertex_class` -- no PDB structure of any kind.
    Spherical law of cosines: for 3 unit vectors from a common vertex at
    known pairwise angles (ref-x, ref-y, x-y), the azimuth phi of y around
    the ref axis (x's own azimuth taken as zero) satisfies
    cos(a_rx)*cos(a_ry) + sin(a_rx)*sin(a_ry)*cos(phi) = cos(a_xy).
    Assumes the vertex is planar (a_rx+a_ry+a_xy ~= 360 deg, true for any
    real sp2 center to within AMBER's own small parameterization
    residual) -- verified algebraically that exact planarity forces
    cos(phi)=-1 (phi=180 deg exactly); AMBER's real angle values give
    176-177 deg for omega's C(i)/N(i+1) vertices (~3 deg short of ideal
    360, a genuine feature of the real parameterization, not an error).
    Returns a value in [0,180] (sign/handedness needs one more real-world
    convention to resolve, irrelevant here since 180 deg is its own
    mirror image -- omega's case only, not a general-purpose solver)."""
    a_rx = np.radians(valence_deg(vertex_class, ref_class, x_class))
    a_ry = np.radians(valence_deg(vertex_class, ref_class, y_class))
    a_xy = np.radians(valence_deg(vertex_class, x_class, y_class))
    cos_phi = (np.cos(a_xy) - np.cos(a_rx) * np.cos(a_ry)) / (np.sin(a_rx) * np.sin(a_ry))
    return np.degrees(np.arccos(np.clip(cos_phi, -1.0, 1.0)))

def _dihedral_from_points(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x, y = np.dot(v, w), np.dot(np.cross(b1n, v), w)
    return np.degrees(np.arctan2(y, x))

def _build_two_vectors(a_rx_deg, a_ry_deg, a_xy_deg, sign):
    """ref along +z; x at angle a_rx from ref (azimuth 0); y at angle a_ry
    from ref, azimuth +-arccos(...) from x (sign picks which of the 2
    mirror-image solutions -- chirality-blind on its own, see callers for
    how each resolves it)."""
    a_rx, a_ry, a_xy = np.radians(a_rx_deg), np.radians(a_ry_deg), np.radians(a_xy_deg)
    ref = np.array([0.0, 0.0, 1.0])
    x = np.array([np.sin(a_rx), 0.0, np.cos(a_rx)])
    cos_phi = np.clip((np.cos(a_xy) - np.cos(a_rx) * np.cos(a_ry)) / (np.sin(a_rx) * np.sin(a_ry)), -1.0, 1.0)
    phi = sign * np.arccos(cos_phi)
    y = np.array([np.sin(a_ry) * np.cos(phi), np.sin(a_ry) * np.sin(phi), np.cos(a_ry)])
    return ref, x, y

_symmetric_tetrahedral_cache = {}

def _build_symmetric_tetrahedral_frame(vertex_class, axis_class, branch_class, h_class):
    """A tetrahedral vertex with exactly 3 non-axis substituents, 2 of
    which SHARE an AMBER class (e.g. CB's CG + HB2 + HB3 -- one real
    branch, two chemically-identical H's) -- the common case for every
    aliphatic chi vertex except the backbone C-alpha and Thr's C-beta
    (see formula_derived_deltas). Unlike the backbone stereocenter, the
    mirror-image ambiguity that 3 pairwise angles alone cannot resolve is
    PROVEN irrelevant here, not fixed by a chirality rule: swapping which
    of the 2 identical substituents is placed via the free sign choice
    (h1) vs solved from the other 3 (h2) is exactly a relabelling of two
    atoms that combined_target_for_axis sums over identically (their
    contribution is 1 + 2*cos(n*X), an even function of the mirror
    parameter X -- verified algebraically and, for a concrete case,
    numerically against combined_target_for_axis itself). sign=+1 is used
    arbitrarily; a genuine 3-distinct-class stereocenter (only the
    backbone Cα and Thr's Cβ, among every axis in this file) is NOT this
    case and must not be routed through this function."""
    key = (vertex_class, axis_class, branch_class, h_class)
    if key in _symmetric_tetrahedral_cache:
        return _symmetric_tetrahedral_cache[key]
    a_axis_branch = valence_deg(vertex_class, axis_class, branch_class)
    a_axis_h = valence_deg(vertex_class, axis_class, h_class)
    a_branch_h = valence_deg(vertex_class, branch_class, h_class)
    a_hh = valence_deg(vertex_class, h_class, h_class)
    axis_v, branch_v, h1_v = _build_two_vectors(a_axis_branch, a_axis_h, a_branch_h, +1.0)
    A = np.array([axis_v, branch_v, h1_v])
    b = np.array([np.cos(np.radians(a_axis_h)), np.cos(np.radians(a_branch_h)), np.cos(np.radians(a_hh))])
    h2_v = np.linalg.solve(A, b)
    h2_v = h2_v / np.linalg.norm(h2_v)
    frame = {"axis": axis_v, "branch": branch_v, "h1": h1_v, "h2": h2_v}
    _symmetric_tetrahedral_cache[key] = frame
    return frame

def formula_derived_deltas(vertex_class, axis_class, atoms_with_classes, ref_atom):
    """Attempts to derive delta (deg) for EVERY atom in atoms_with_classes
    (a list of (atom, class) -- every real non-axis substituent at this
    vertex) purely from AMBER's own angle tables, relative to ref_atom
    (matching whatever the caller's own ref_b/ref_c convention already
    is -- no PDB structure of any kind. Returns None if this vertex isn't
    one of the solvable cases:
    - exactly 2 substituents (trigonal/planar vertex): idealized_azimuth_deg,
      forced unambiguous by planarity (see that function).
    - exactly 3 substituents, all 3 sharing one class (a methyl's 3 H's,
      no branch atom at all -- see the "all-same-class" block below):
      pure 3-fold symmetry, no chirality question, no frame to build.
    - exactly 3 substituents with >=2 sharing a class (tetrahedral, one
      real branch + 2 chemically-identical substituents):
      _build_symmetric_tetrahedral_frame, sign-independent (see its own
      docstring).
    A genuine 3-distinct-class tetrahedral stereocenter (only the
    backbone Cα -- handled separately by ca_tetrahedral_delta -- and
    Thr's C-beta, among every axis in this file) returns None here; the
    caller falls back to real-structure measurement for that one case."""
    if len(atoms_with_classes) == 2:
        (a1, c1), (a2, c2) = atoms_with_classes
        delta_a2 = idealized_azimuth_deg(vertex_class, axis_class, c1, c2)
        deltas = {a1: 0.0, a2: delta_a2}
        shift = deltas[ref_atom]
        return {a: d - shift for a, d in deltas.items()}
    if len(atoms_with_classes) == 3:
        classes = [c for _, c in atoms_with_classes]
        if len(set(classes)) == 1:
            # A methyl's 3 H's (or any vertex with 3 non-axis substituents
            # that all share one class): an even MORE symmetric case than
            # the 1-branch+2-identical one below, with no branch atom to
            # build a frame from at all. combined_target_for_axis sums
            # this vertex's contribution over all 3 real substituents, and
            # every one of them looks up the identical (n,k,phase) AMBER
            # term when paired with any given other-side atom (same
            # class) -- so the sum only depends on the SET of 3 azimuths,
            # never on which specific atom gets which value. Placing them
            # at {0, +120, -120} (any assignment) is therefore exact, not
            # an approximation: algebraically, summing
            # exp(i*n*{0,+120,-120}) gives 1+2*cos(n*120deg), which is
            # exactly 3 when n is a multiple of 3 and exactly 0 otherwise
            # -- the vertex's own 3-fold symmetry forces every non-3x
            # harmonic to cancel on its own, matching amber99sb.xml's own
            # choice to parametrize this term as pure n=3 with no n=1/n=2
            # companion (confirmed by direct XML inspection, not assumed).
            a0, a1, a2 = (a for a, _ in atoms_with_classes)
            deltas = {a0: 0.0, a1: 120.0, a2: -120.0}
            shift = deltas[ref_atom]
            return {a: d - shift for a, d in deltas.items()}
        if set(classes) == {"OH", "CT", "H1"}:
            # Threonine's C-beta -- the one genuine 3-distinct-class
            # tetrahedral stereocenter among every chi vertex in this file
            # besides the backbone Cα (see thr_cb_tetrahedral_delta's own
            # docstring for the CIP derivation). Class-matched, not
            # resname-matched, so it fires wherever this exact substituent
            # pattern occurs -- today that is only Thr's CB.
            by_class = {c: a for a, c in atoms_with_classes}
            deltas = {by_class["CT"]: thr_cb_tetrahedral_delta("OG1", "CG2"),
                     by_class["OH"]: 0.0,
                     by_class["H1"]: thr_cb_tetrahedral_delta("OG1", "HB")}
            shift = deltas[ref_atom]
            return {a: d - shift for a, d in deltas.items()}
        common_class = next((c for c in classes if classes.count(c) >= 2), None)
        h_atoms = [a for a, c in atoms_with_classes if c == common_class]
        branch_atoms = [a for a, c in atoms_with_classes if c != common_class]
        if common_class is None or len(h_atoms) != 2 or len(branch_atoms) != 1:
            return None
        branch_atom = branch_atoms[0]
        h1_atom, h2_atom = h_atoms
        branch_class = next(c for a, c in atoms_with_classes if a is branch_atom)
        frame = _build_symmetric_tetrahedral_frame(vertex_class, axis_class, branch_class, common_class)
        origin = np.zeros(3)
        pts = {"origin": origin, "axis": origin + frame["axis"], "branch": origin + frame["branch"],
              "h1": origin + frame["h1"], "h2": origin + frame["h2"]}
        deltas = {
            branch_atom: 0.0,
            h1_atom: _dihedral_from_points(pts["branch"], pts["origin"], pts["axis"], pts["h1"]),
            h2_atom: _dihedral_from_points(pts["branch"], pts["origin"], pts["axis"], pts["h2"]),
        }
        shift = deltas[ref_atom]
        return {a: d - shift for a, d in deltas.items()}
    return None

def _build_ca_stereocenter_frame():
    """The backbone C-alpha's 4 substituents (N, C, CB, HA -- generic
    classes N/C/CT/H1) as unit vectors from CA, derived ONCE from AMBER's
    own real HarmonicAngleForce angles (chirality-blind: any 3 known
    pairwise angles fix a tetrahedral arrangement up to mirror
    reflection) PLUS the universal L-amino-acid (S)-configuration
    convention (CIP/CORN rule: viewed with HA pointing away from the
    observer, N -> C(carbonyl) -> CB winds in the sense fixed below --
    true without exception for all 19 chiral proteinogenic amino acids,
    including Cys, whose CIP *label* differs (R not S) only because
    sulfur reorders substituent priority, not because its physical
    arrangement differs).

    This is NOT measured from any PDB structure -- the sign ambiguity
    that bond angles alone cannot resolve (found and abandoned as
    unsolvable 2026-08-01 when attempted via pairwise-numeric-consistency
    instead of the actual documented chirality rule) is resolved here by
    the textbook convention itself, a fixed discrete fact, not a
    continuously-variable structural measurement. Verified (2026-08-08)
    against real ubiquitin.pdb (ARG42): the 4 dihedrals this frame is
    actually used for (see ca_tetrahedral_delta) reproduce the real
    measured values to 0.15-2.3 deg -- ordinary real-bond-angle-vs-AMBER-
    equilibrium variation, not a sign or convention error."""
    a_nc = valence_deg("CT", "N", "C")
    a_ncb = valence_deg("CT", "N", "CT")
    a_ccb = valence_deg("CT", "C", "CT")
    a_nha = valence_deg("CT", "N", "H1")
    a_cha = valence_deg("CT", "C", "H1")
    a_cbha = valence_deg("CT", "CT", "H1")

    # sign=+1 here (not -1) is the verified L-amino-acid choice -- confirmed
    # by direct comparison against real ubiquitin.pdb, see docstring above.
    N, C, CB = _build_two_vectors(a_nc, a_ncb, a_ccb, +1.0)

    A = np.array([N, C, CB])
    b = np.array([np.cos(np.radians(a_nha)), np.cos(np.radians(a_cha)), np.cos(np.radians(a_cbha))])
    HA = np.linalg.solve(A, b)
    HA = HA / np.linalg.norm(HA)
    return {"N": N, "C": C, "CB": CB, "HA": HA}

_CA_STEREOCENTER = _build_ca_stereocenter_frame()

def ca_tetrahedral_delta(ref_name, b_name, c_name, other_name):
    """dihedral(ref, b_name, c_name, other), read directly off the
    resolved CA stereocenter frame (_CA_STEREOCENTER) instead of measured
    from any real PDB structure. Matches group_target's own delta_of
    formula exactly: verified numerically against real ubiquitin.pdb
    (ARG42) that BOTH branches of delta_of (is_b_side True or False)
    reduce to this SAME dihedral(ref_atom, b_atom, c_atom, other_atom)
    form (mod 360 deg, immaterial -- exp(i*n*x) is exactly 360-periodic
    for any integer n, so a +-360 deg relabeling changes nothing
    downstream in combined_target_for_axis). Names in
    {"N","C","CB","HA","CA"}."""
    CA = np.zeros(3)
    pts = {k: CA + v for k, v in _CA_STEREOCENTER.items()}
    pts["CA"] = CA
    return _dihedral_from_points(pts[ref_name], pts[b_name], pts[c_name], pts[other_name])

def _build_thr_cb_stereocenter_frame():
    """Threonine's C-beta -- the one genuine SECOND tetrahedral
    stereocenter among every chi axis in this file (besides the backbone
    Cα itself): real proteinogenic L-threonine is (2S,3R)-2-amino-3-
    hydroxybutanoic acid (Thr and Ile are the two standard amino acids
    with 2 stereocenters; the (2S,3S) diastereomer is a distinct real
    molecule, allo-threonine, never used in proteins).

    CIP priority at C-beta: OG1 (O, highest) > CA (attached to N/C/H) >
    CG2 (attached to H/H/H) > HB (lowest, H). For R configuration (the
    real one), viewed with HB pointing away from the observer, OG1 -> CA
    -> CG2 winds CLOCKWISE -- resolved here the same way as the backbone
    Cα (chirality-blind bond angles + the documented, universal
    configuration fact, never measured from a PDB structure).

    All 6 pairwise angles at this vertex are AMBER's generic tetrahedral
    109.5 deg (amber99sb.xml has no CT-OH/CT-CT/CT-H1-specific
    differentiation here) -- a perfectly regular tetrahedron.

    Verified (2026-08-08) against real ubiquitin.pdb (THR7): reproduces
    the real dihedral(CG2,CA,CB,OG1)/dihedral(CG2,CA,CB,HB) to ~1.3 deg
    (ordinary real-vs-ideal angle variation, not a sign error) -- same
    validation standard as every other formula-derived vertex in this
    file."""
    a = 109.5
    # sign=+1 verified R/real-Thr (winding formula calibrated against the
    # already-verified backbone case: negative winding = S/CCW there,
    # positive winding = R/CW here, both cross-checked against real
    # ubiquitin.pdb before trusting this).
    CA, OG1, CG2 = _build_two_vectors(a, a, a, +1.0)
    M = np.array([CA, OG1, CG2])
    b = np.array([np.cos(np.radians(a))] * 3)
    HB = np.linalg.solve(M, b)
    HB = HB / np.linalg.norm(HB)
    return {"CA": CA, "OG1": OG1, "CG2": CG2, "HB": HB}

_THR_CB_STEREOCENTER = _build_thr_cb_stereocenter_frame()

def thr_cb_tetrahedral_delta(ref_name, other_name):
    """dihedral(ref, CA, CB, other), read directly off the resolved Thr
    C-beta frame (_THR_CB_STEREOCENTER) instead of measured from any real
    PDB structure. Names in {"CA","OG1","CG2","HB"}."""
    CB = np.zeros(3)
    pts = {k: CB + v for k, v in _THR_CB_STEREOCENTER.items()}
    pts["CB"] = CB
    return _dihedral_from_points(pts[ref_name], pts["CA"], pts["CB"], pts[other_name])

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

# Overall size of an n>=2 ring, as a multiple of the axis bond length.
#
# A ghost ring's energy is NOT a pure function of the torsion angle: with
# ghost planes separated by Lam and ring radius rho_r, the ghost-ghost
# distance is d^2 = Lam^2 + 2*rho_r^2*(1-cos psi) (exactly what
# closed_form_d2 computes), so d(d^2)/dL = 2*Lam -- the energy also depends
# on the REAL distance L between the two axis atoms. Expanding
# E = sum_i 0.5k(d_i-d0)^2 over a ring of M>=2 springs, the sum_i cos(psi_i)
# term vanishes and at the torsion minimum everything collapses to
#     E_min = 0.5 * k * M * (sqrt(2*rho_r^2 + Lam^2) - d0)^2
# a positive-definite quadratic in L. Since mu=1.0 sets d0 to exactly that
# square root evaluated at AMBER's IDEAL bond length, E_min is zero only when
# the real bond happens to be ideal; any real deviation adds spurious energy
# (found 2026-08-09 on real Fs-peptide MD frames, where guanidinium C-N bonds
# deviate by up to 0.06 A: +32 kJ/mol of pure artifact, see
# doc/BondedForceFieldSprings.md).
#
# The fix, and why it is free: the harmonic content depends only on
# eps = 2*rho_r^2/(Lam^2 + 2*rho_r^2), which is exactly RING_SHAPES' own rho.
# Scaling rho_r and Lam TOGETHER leaves eps -- hence the whole torsion curve --
# untouched while multiplying the artifact by 1/scale^2 (verified numerically:
# scale 1..6 leaves the curve identical to 1.7e-13 kJ/mol while the artifact
# falls 14.87 -> 1.58 kJ/mol). Lam is decoupled from the real bond length by
# tilting the ghosts off their anchor's perpendicular plane (theta != 90 deg,
# already supported by the .bi.ff format and by GhostParticle::computePosition,
# and already used by the n=1 case) -- no format or C++ change needed.
#
# n=1 needs no scaling: its exact d0=0/theta=60 construction already puts
# Lam = L - 2*r*cos(60deg) = L - L_axis = 0 at the ideal bond length, so it is
# already at the optimum this constant exists to reach for n>=2.
#
# 4.0 buys a 16x artifact reduction while keeping ghosts within ~3 A of their
# anchors; going further has diminishing returns (1/scale^2) against steadily
# longer force-redistribution lever arms on the real anchor atoms.
RING_SCALE = 4.0

def ring_shape_to_geometry(L_axis, n):
    """n=1 is special-cased to an EXACT closed-form solution (found
    2026-08-05, see combined_target_for_axis/emit_ghost_ring's own
    comments): with theta free (not fixed to 90 deg), a single ghost-ghost
    spring (d0=0) can be made to reproduce ANY AMBER k*(1+cos(phi-phase))
    term with ZERO residual, DC included -- not just harmonic-pure (that
    part was already exact) but DC/amplitude-ratio-exact too, by choosing
    cos(theta) = L_axis/(2r) for any r >= L_axis/2. r=L_axis (a plain,
    unremarkable choice on the same scale as the axis bond itself) gives
    theta=60 deg exactly. Verified numerically: max|E_ring-E_amber|=0 to
    machine precision, for any target k/phase. n>=2 keeps the RING_SHAPES
    (rho, mu) table (theta=90 deg) -- the M=N=n comb filter's own base
    function is no longer purely first-harmonic once d0!=0 (needed to get
    ANY power at n>=2 at all, see doc/BondedForceFieldSprings.md), so the
    same one-shot exact trick does not directly apply there; that ring's
    own residual (harmonic leakage into 2n, 3n... plus whatever DC/
    amplitude deviation RING_SHAPES leaves) is smaller (already close to a
    1:1 ratio) but not yet proven exact -- open follow-up, not silently
    assumed equival to n=1's case."""
    if n == 1:
        r = L_axis
        theta = np.degrees(np.arccos(L_axis / (2.0 * r)))
        return r, theta, 0.0
    rho, mu = RING_SHAPES[n]
    # Unscaled shape: ghost planes Lam=L_axis apart (theta=90 deg), ring
    # radius rho_r chosen so that 2*rho_r^2/(Lam^2+2*rho_r^2) == rho.
    rho_r = np.sqrt(rho / (2.0 * (1.0 - rho))) * L_axis
    # Scale radius and plane separation together: eps -- hence the harmonic
    # content and the whole torsion curve -- is invariant, only the spurious
    # bond-length sensitivity shrinks (as 1/RING_SCALE^2). See RING_SCALE.
    rho_r *= RING_SCALE
    lam = RING_SCALE * L_axis
    # Realize that plane separation from anchors that are only L_axis apart by
    # tilting each ghost off the perpendicular: closed_form_d2's own axial term
    # is (r1*cos(t1) - L + r2*cos(t2)), so with both sides equal the separation
    # is Lam = L_axis - 2*r*cos(theta) -- solve for the axial offset.
    axial = (L_axis - lam) / 2.0                     # negative: tilted outward
    r = float(np.hypot(rho_r, axial))
    theta = float(np.degrees(np.arctan2(rho_r, axial)))
    d0 = mu * np.sqrt(2.0 * rho_r * rho_r + lam * lam)
    return r, theta, d0

def ring_curve_abstract(L_axis, n, r, theta_deg, d0, k, delta_base, phi_grid_deg, M=None, N=None,
                        r2=None, theta2_deg=None):
    # r2/theta2_deg default to r/theta_deg (the symmetric shape the old
    # free-placement construction used). They are given explicitly since
    # 2026-08-14: a ghost is now the rotated image of a REAL reference
    # atom, so each side's radius/axial offset is that atom's own and the
    # two sides generally differ. See emit_ghost_ring.
    # The 1D Fourier-space model used to CALIBRATE k/delta_base (see
    # calibrate_ring) -- the comb-filter sum over an MxN grid (defaults to
    # M=N=n when not given, the original construction), evaluated over an
    # abstract phi grid. NOT the real 3D placement (see emit_ghost_ring for
    # that, and its header comment for the sign convention linking the two
    # -- verified by direct 3D simulation against a real axis before being
    # trusted here, see doc/BondedForceFieldSprings.md).
    #
    # M,N need not equal n individually -- only LCM(M,N)=n matters (found
    # 2026-08-08/09: the M=N=n choice used everywhere until now is not a
    # mathematical requirement, just the simplest uniform convention. The
    # general comb-filter derivation -- S(phi) = sum_i sum_j f(phi-delta_base
    # -beta_i+gamma_j) with beta_i=i*360/M, gamma_j=j*360/N -- shows the
    # surviving harmonics are exactly the common multiples of M and N, i.e.
    # multiples of LCM(M,N), regardless of M,N individually; verified
    # numerically that an asymmetric M=n,N=1 ring reproduces a real target
    # to amplitude ratio 0.9998 / phase diff 0.000 deg, identical to M=N=n,
    # while using fewer ghosts/springs -- see doc/BondedForceFieldSprings.md).
    if M is None:
        M = n
    if N is None:
        N = n
    if r2 is None:
        r2 = r
    if theta2_deg is None:
        theta2_deg = theta_deg
    e = np.zeros_like(phi_grid_deg)
    for i in range(M):
        beta = i * 360.0 / M
        for j in range(N):
            gamma = j * 360.0 / N
            delta_ij = delta_base + beta - gamma
            d2v = closed_form_d2(L_axis, r, theta_deg, r2, theta2_deg, delta_ij, phi_grid_deg)
            d = np.sqrt(np.maximum(d2v, 0.0))
            e += 0.5 * k * (d - d0) ** 2
    return e

def geom_via_far_anchor(L_axis, bond_to_other, angle_at_other_deg):
    """(r, theta) of a reference atom relative to anchor B, when that atom is
    bonded to the OTHER anchor C instead of to B.

    Needed by the PLANARITY impropers: AMBER lists the hub third and bonds
    all three peripherals to it, so the quad's first atom hangs off C, not
    off B, and its distance/angle to B are 1-3 quantities. Plain triangle
    B-C-R with |BC|=L_axis, |CR|=bond_to_other and the angle at C between
    them: the law of cosines gives |BR|, then the angle at B. Both are still
    pure AMBER table values -- no structure is measured.

    (Proper torsions never take this path: their reference atom really is
    bonded to its own anchor.)"""
    g = np.radians(angle_at_other_deg)
    r = float(np.sqrt(L_axis ** 2 + bond_to_other ** 2 - 2.0 * L_axis * bond_to_other * np.cos(g)))
    if r < 1e-9:
        return 0.0, 0.0
    cos_theta = (L_axis ** 2 + r ** 2 - bond_to_other ** 2) / (2.0 * L_axis * r)
    return r, float(np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0))))

def choose_d0(L_axis, n, r1, theta1_deg, r2, theta2_deg, M, N):
    """Ghost-ghost equilibrium length: the RMS ghost separation sqrt(C1),
    which is the middle of the range the pair can actually reach.

    d0 does NOT trade off against harmonic purity, and a search over it is
    degenerate -- provably so. On an M=n,N=1 ring the sum of d^2 over the
    comb is exactly n*C1 (the comb kills the fundamental for n>=2), so it
    contributes only a constant; every harmonic comes from the -2*d0*sum(d)
    term alone and is therefore proportional to d0. All harmonic RATIOS are
    d0-independent -- verified numerically too: leakage stays at 0.1108% for
    every mu from 0.05 to 1.6 on a real chi1 geometry. d0 sets amplitude
    (hence k) and DC, nothing else.

    An earlier version of this function minimised leakage over d0 and, with
    nothing to discriminate, wandered to d0=4.29 A for a pair that never
    gets past 3.05 A -- a permanently stretched spring, which showed up as
    a NEGATIVE total dihedral energy and a 10x kinetic-energy blow-up on
    GKinase. sqrt(C1) is what the old RING_SHAPES table used as well (mu=1.0
    for n=2,3), for the same reason.

    n=1 keeps d0=0 exactly -- an algebraic identity, not a compromise: with
    d0=0 the energy is 0.5*k*d^2 = 0.5*k*(C1 - C2*cos(phi-delta)), a pure
    first harmonic with EXACTLY zero power anywhere else."""
    if n == 1:
        return 0.0
    c1 = (r1 * np.sin(np.radians(theta1_deg))) ** 2 + (r2 * np.sin(np.radians(theta2_deg))) ** 2 + (
        r1 * np.cos(np.radians(theta1_deg)) - L_axis + r2 * np.cos(np.radians(theta2_deg))) ** 2
    return float(np.sqrt(c1))

def calibrate_ring(L_axis, n, r, theta_deg, d0, target_complex, M=None, N=None, r2=None, theta2_deg=None):
    # Energy is exactly linear in k at fixed geometry (see
    # doc/BondedForceFieldSprings.md, step 3 of the closed-form derivation)
    # -- one evaluation at k=1 gives the ring's own harmonic-n Fourier
    # coefficient (magnitude AND phase), from which k (magnitude match) and
    # delta_base (phase match) solve directly, no iteration. M,N (default
    # n,n): see ring_curve_abstract's own docstring -- the phase/k formulas
    # below depend only on n (=LCM(M,N)), verified algebraically that
    # rotating delta_base by D shifts the ring's n-th harmonic phase by
    # -nD regardless of M,N individually, so nothing below needs to change
    # when M != N, only the curve construction they're evaluated against.
    e0 = ring_curve_abstract(L_axis, n, r, theta_deg, d0, 1.0, 0.0, PHI_GRID_DEG, M, N, r2, theta2_deg)
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
    # The ring's own n=0 (DC/mean) Fourier component at the CALIBRATED k --
    # exactly k times its value at k=1 (energy is linear in k, see above).
    # This is a pure construction artifact (a sum of squares can only ever
    # ADD a positive baseline, see RING_SHAPES's own header comment): it has
    # nothing to do with the real AMBER torsion's own physical mean
    # (target[0] in combined_target_for_axis, handled separately by the
    # caller) and must be tracked so it can be subtracted back out when
    # reporting an absolute energy (never affects forces -- a constant has
    # zero gradient).
    dc_ring = k * float(np.real(complex_fourier_coeff(e0, 0)))
    return k, delta_base_deg, dc_ring

n_dihedral_ok = 0
n_dihedral_skip = 0
n_ghost_particles = 0
dihedral_lines = []
ghostparticle_lines = []
# Pure logging side-channel (no effect on anything emitted): one entry per
# emit_ghost_ring call, i.e. one per (resname, axis, harmonic, group) ring
# actually deployed -- lets an external validation script reconstruct
# EXACTLY the same AMBER target / calibrated ring curves this file itself
# used, with zero risk of drifting from the real numbers (see
# doc/BondedForceFieldSprings.md's energy-curve validation section).
AXIS_RING_LOG = []

# (resname, anchors, r, theta, azimuth) -> already-emitted ghost name. Lets
# per-pair ring groups share the ghosts they have in common instead of
# emitting duplicates -- see emit_ghost_ring's own dedup comment.
_ghost_registry = {}

def emit_ghost_ring(resname, axis_label, family, n, L_axis, target_complex, atom_B, atom_C, atom_ref_b, atom_ref_c,
                    axis_dc_target=0.0, group_tag="", ref_geom_b=None, ref_geom_c=None):
    """Emits one full ring group (M+N ghost particles, M*N ghost-ghost
    springs -- M=n,N=1 for n>=2, the ghost/spring-minimal construction,
    see the M,N note below) reproducing one target Fourier harmonic
    exactly (up to the small, bounded leakage quantified in RING_SHAPES's
    own derivation). atom_B/atom_C are the
    real dihedral axis atoms (already +/- resolved); atom_ref_b/atom_ref_c
    are each side's own real reference atom (fixes what "azimuth=0" means
    for that side -- must be the SAME atoms used to measure the target's
    own phase convention, i.e. b_neighbors[0]/c_neighbors[0], see
    combined_target_for_axis).

    axis_dc_target: the real AMBER torsion's own physical mean for this
    axis (target[0] in the caller, sum of each matched term's own k --
    AMBER's PeriodicTorsionForce term k*(1+cos(...)) has mean exactly k).
    Attributed to exactly ONE ring per axis by the caller (every other
    ring for the same axis passes 0.0) so it's subtracted exactly once
    per axis, not once per harmonic. Combined with this ring's own
    (non-physical, construction-artifact) DC into a single per-spring
    dc_offset column: BondedForceFieldReader accumulates this across
    every DIHEDRAL entry it actually applies, giving an exact correction
    (subtracted only when reporting energy, never affecting forces) so
    BioSpring's absolute dihedral energy is directly comparable to
    AMBER's, instead of being inflated by the ring construction's own
    unavoidable positive baseline (see calibrate_ring's own comment).

    group_tag: disambiguates ghost/rule names when the caller emits more
    than one ring group for the SAME axis anchored on different real atoms
    (see generate_backbone_axis) -- must be unique per group sharing an
    axis_label, empty ("", the default) is fine whenever an axis only ever
    has one group (every sidechain/omega axis today).

    Placement/sign convention (verified by a full 3D rotating simulation
    against a real axis before being trusted here -- see
    doc/BondedForceFieldSprings.md): a b-side ghost i sits at azimuth
    i*360/M (measured from atom_B, relative to atom_ref_b, via
    spn::GhostParticle::computePosition(atom_B, atom_C, atom_ref_b, ...)).
    A c-side ghost j sits at azimuth delta_base - j*360/N (note the
    MINUS -- placed via computePosition(atom_C, atom_B, atom_ref_c, ...),
    i.e. B/C swapped since the c-side ghost's own axis direction points
    the other way). This combination reproduces exactly the SAME abstract
    delta_ij = delta_base + beta_i - gamma_j used to calibrate k/delta_base
    in Fourier space above -- get either the swap or the minus sign wrong
    and the real 3D energy has the wrong phase relative to the molecule's
    true dihedral angle, silently, without any local check catching it.

    M,N (found and switched 2026-08-09): only LCM(M,N)=n matters for the
    comb filter (see ring_curve_abstract's own docstring), not M=N=n --
    M=n,N=1 (or the mirror M=1,N=n) reproduces the exact same curve with
    fewer ghosts (M+N=n+1 instead of 2n) and fewer springs (M*N=n
    instead of n*n), verified against a real target to amplitude ratio
    0.9998/phase diff 0.000 deg. Applied system-wide: 984->737 ghost
    particles (-25%), 1166->492 springs (-58%) across the 245 rings this
    file emits. n=1 is already M=N=1 (the minimal case, LCM(1,1)=1, no
    asymmetric variant exists or is needed)."""
    global n_dihedral_ok, n_ghost_particles

    AXIS_RING_LOG.append({"resname": resname, "axis_label": axis_label, "family": family, "n": n,
                          "L_axis": L_axis, "target": target_complex, "axis_dc_target": axis_dc_target,
                          "group_tag": group_tag})

    # Since 2026-08-14 a ghost is the rotated image of its side's real
    # reference atom (spn::GhostParticle::computePositionByRotation), so
    # its radius and axial offset are that atom's own, taken from the AMBER
    # bond/angle tables like everything else here -- never chosen. Only d0
    # and k are still calibrated. ref_geom_b/ref_geom_c are (r, theta) in
    # the same convention as combined_target_for_axis's b_geom/c_geom.
    r, theta_deg = ref_geom_b
    r_c, theta_c_deg = ref_geom_c
    M, N = (n, 1) if n >= 2 else (1, 1)
    d0 = choose_d0(L_axis, n, r, theta_deg, r_c, theta_c_deg, M, N)
    k, delta_base, dc_ring = calibrate_ring(L_axis, n, r, theta_deg, d0, target_complex, M, N, r_c, theta_c_deg)
    # Split evenly across this ring's M*N springs -- BondedForceFieldReader
    # sums dc_offset per spring it actually applies, so this naturally
    # handles a partially-skipped ring too (e.g. Met's phi ghosts missing
    # their "-C" anchor at the true N-terminus): the offset for the springs
    # that never get created is correctly never added either.
    dc_offset_per_spring = (dc_ring - axis_dc_target) / (M * N)

    # group_tag disambiguates ghost/rule names when the SAME axis emits
    # more than one ring group anchored on different real atoms (see
    # generate_backbone_axis's per-owning-pair anchoring) -- without it,
    # two groups sharing the same non-grouped side's reference (e.g. both
    # of psi's O-group and +N-group use the same b-side CB reference) would
    # emit ghosts with IDENTICAL names ("GHpsin3B0" twice), and
    # BondedForceFieldReader's name-based resolution would silently pick
    # whichever one happens to match first -- a real, found-in-review bug,
    # not hypothetical (caught by tracing the actual .bi.ff output before
    # this parameter existed). Sanitized to strip "+"/"-": a ghost's own
    # name is never itself given a cross-residue prefix when referenced
    # later (ghosts are always resolved within their own creating residue).
    tag = re.sub(r"[+-]", "", group_tag)

    # Deduplication: a per-pair axis (see generate_sidechain_axis's per_pair)
    # emits one ring group per real substituent pair, and two groups sharing a
    # side's substituent ask for the SAME ghost -- same anchors, same
    # (r, theta, delta). Omega's n=2 is the clearest case: its 4 groups
    # (CA/CA, CA/H, O/CA, O/H) want only 6 distinct ghosts between them, not
    # 12, because groups CA/CA and CA/H share both their CA-anchored b-side
    # ghosts, and groups CA/CA and O/CA share their +CA-anchored c-side one.
    # A ghost is a massless virtual site fully determined by those 7 fields,
    # so two identical ones would sit at the same point and behave
    # identically -- reusing the first is a pure saving, with the springs
    # simply naming the shared ghost. Measured: 22% of all GHOSTPARTICLE
    # lines were exact duplicates before this.
    # gr/gtheta are the reference ATOM's own radius/angle for THAT side --
    # the two sides generally differ now that a ghost is that atom's rotated
    # image, so they are passed in rather than closed over.
    def ghost(name, res, a1, a2, ref, azim, gr, gtheta):
        global n_ghost_particles
        key = (res, a1, a2, ref, round(gr, 4), round(gtheta, 2), round(azim, 3))
        if key in _ghost_registry:
            return _ghost_registry[key]
        ghostparticle_lines.append(
            f"GHOSTPARTICLE {name:16s} {res:7s} {a1:5s} {a2:6s} {ref:6s} "
            f"{gr:7.4f} {gtheta:7.2f} {azim:9.3f}")
        _ghost_registry[key] = name
        n_ghost_particles += 1
        return name

    b_names, c_names = [], []
    for i in range(M):
        beta = i * 360.0 / M
        b_names.append(ghost(f"GH{axis_label}{tag}n{n}B{i}", resname, atom_B, atom_C, atom_ref_b, beta,
                             r, theta_deg))
    for j in range(N):
        gamma = delta_base - j * 360.0 / N
        c_names.append(ghost(f"GH{axis_label}{tag}n{n}C{j}", resname, atom_C, atom_B, atom_ref_c, gamma,
                             r_c, theta_c_deg))

    for bn in b_names:
        for cn in c_names:
            rule_name = f"{resname}_{axis_label}{tag}_n{n}_{bn[-1]}{cn[-1]}"
            dihedral_lines.append(
                f"DIHEDRAL {rule_name:28s} {resname:7s} {family:9s} {bn:16s} {cn:17s} "
                f"{d0:7.4f} {k:10.4f} {dc_offset_per_spring:10.5f}")
            n_dihedral_ok += 1

def emit_ghost_rings_for_axis(resname, axis_label, family, L_axis, target, dc_by_harmonic, atom_B, atom_C,
                              atom_ref_b, atom_ref_c, ref_geom_b, ref_geom_c):
    """Emits one ring (see emit_ghost_ring) per non-zero harmonic in
    `target` (n -> complex Fourier coefficient, see combined_target_for_axis),
    each calibrated against ITS OWN physical DC share (dc_by_harmonic[n] --
    NOT the axis's combined total attributed to a single ring: found and
    fixed 2026-08-05, see combined_target_for_axis's own comment for why
    that broke each ring's own curve even though it left the axis grand
    total correct) -- shared by generate_sidechain_axis, its GKinase
    variant, and generate_backbone_axis instead of repeating this
    bookkeeping 3 times."""
    for n, zt in target.items():
        if n == 0 or abs(zt) < 1e-6:
            continue
        emit_ghost_ring(resname, axis_label, family, n, L_axis, zt, atom_B, atom_C, atom_ref_b, atom_ref_c,
                       dc_by_harmonic.get(n, 0.0), ref_geom_b=ref_geom_b, ref_geom_c=ref_geom_c)

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
    nothing, not a gap.

    Returns (target, dc_by_harmonic): target[0] is still the axis's own
    total real DC (every matched term's own k, summed) for logging/sanity;
    dc_by_harmonic[n] is the sub-total contributed by matched terms whose
    OWN periodicity is exactly n -- i.e. the true physical DC that
    harmonic n's own ring should be calibrated against (see emit_ghost_ring
    -- every AMBER PeriodicTorsionForce term k*(1+cos(n*phi-phase)) has
    DC=k=its own amplitude, a 1:1 ratio for every n; a ring's own DC/
    amplitude ratio need not be 1:1 by construction -- e.g. the n=1 ring's
    exact-purity shape (d0=0) has a fixed, target-independent ratio of
    2:1 -- so subtracting the WRONG ring's own share, or the whole axis's
    combined DC, corrects the grand total but not each ring's own
    curve, found and fixed 2026-08-05: see doc/BondedForceFieldSprings.md)."""
    target = {0: 0j}
    dc_by_harmonic = {}
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
                dc_by_harmonic[n] = dc_by_harmonic.get(n, 0.0) + k
    return (target, dc_by_harmonic) if any_term else (None, None)

def generate_sidechain_axis(resname, b_name, c_name, axis_label, per_pair=False):
    """per_pair: emit ONE ring per real (b-substituent, c-substituent) pair,
    each anchored on that pair's own two real atoms, instead of one combined
    ring per harmonic on a shared generic anchor.

    Use it only for an axis whose pairs all reinforce IN PHASE under the ideal
    geometry (|target[n]| == sum of every pair's |k|, i.e. maximal coherence).
    There, a combined ring bakes that maximal coherence in and comes out
    systematically too stiff on real structures, whose substituents are never
    exactly at their ideal azimuths -- measured on omega: torque ratio
    |ring|/|AMBER| 1.51 median, ~14 kJ/mol/rad error. Per-pair removes the
    assumption instead of correcting for it: each ring reproduces one real
    AMBER term on ITS OWN real dihedral, so the sum is AMBER's sum by
    construction, whatever the real geometry is (verified on omega: energy
    error exactly 0.000 kJ/mol over 120 real instances, torque 14.25 -> 0.29).

    Do NOT use it where the pairs partly cancel -- there the combined ring is
    already unbiased (measured: psi 1.00, phi 0.85) and splitting only
    multiplies topology, which is exactly why it was reverted in 2026-08-07
    before the coherence criterion was understood. See
    doc/BondedForceFieldSprings.md."""
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
        # Formula-derived delta (found and fixed 2026-08-08, same
        # principle as generate_omega_axis/generate_backbone_axis -- no
        # PDB structure of any kind, see formula_derived_deltas/
        # ca_tetrahedral_delta's own docstrings): chi1's b-side is always
        # the backbone Cα stereocenter (reuses the already CIP-resolved
        # frame directly -- its 3 substituents, N/C/HA, are all distinct
        # classes, not the generic solvable case below); every other
        # vertex tries the generic trigonal/symmetric-tetrahedral solver.
        # Falls back to the old real-structure measurement only for a
        # genuine, still-unsolved 3-distinct-class stereocenter (Thr's
        # C-beta is the only one among every chi axis in this file).
        vertex_name = b_name if is_b_side else c_name
        axis_partner_name = c_name if is_b_side else b_name
        formula_deltas = None
        if vertex_name == "CA":
            try:
                formula_deltas = {n: ca_tetrahedral_delta(ref_name, vertex_name, axis_partner_name, n)
                                 for n in neighbors}
            except KeyError:
                formula_deltas = None
        if formula_deltas is None:
            atoms_with_classes = [(n, resolved_class(resname, n, False, False)) for n in neighbors]
            formula_deltas = formula_derived_deltas(vertex_class, other_class, atoms_with_classes, ref_name)

        # SIGN-CORRECTED on the b-side (found and fixed 2026-08-04, see
        # generate_omega_axis's own docstring for the full derivation):
        # delta must be ref_phi MINUS the measured value there, not the
        # other way round -- `phi_meas - ref_phi` on the b-side gives the
        # negative of the wanted delta. The c-side formula was already
        # correct (it's the b-side's asymmetric roles in
        # real_dihedral_deg's own argument order that flips the sign).
        geom = {}
        for n in neighbors:
            nc = resolved_class(resname, n, False, False)
            r = bond_len_A(nc, vertex_class)
            theta = valence_deg(vertex_class, nc, other_class)
            if formula_deltas is not None:
                delta = formula_deltas[n]
            elif n == ref_name:
                delta = 0.0
            elif is_b_side:
                delta = ref_phi - real_dihedral_deg(resname, n, b_name, c_name, c_ref)
            else:
                delta = real_dihedral_deg(resname, b_ref, b_name, c_name, n) - ref_phi
            geom[n] = (r, theta, delta)
        return geom

    b_geom = side_geom(b_neighbors, b_class, c_class, b_ref, True)
    c_geom = side_geom(c_neighbors, c_class, b_class, c_ref, False)

    class_of = lambda n: resolved_class(resname, n, False, False)
    target, dc_by_harmonic = combined_target_for_axis(b_geom, c_geom, b_class, c_class, class_of)
    if target is None:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}): no real substituent pair matches a "
             f"specific AMBER entry and no generic wildcard either")
        n_dihedral_skip += 1
        return

    if per_pair:
        emitted = 0
        for bn in b_neighbors:
            for cn in c_neighbors:
                t_pair, dc_pair = combined_target_for_axis({bn: (0, 0, 0.0)}, {cn: (0, 0, 0.0)},
                                                           b_class, c_class, class_of)
                if t_pair is None:
                    continue
                for n, tn in sorted(t_pair.items()):
                    if n == 0 or abs(tn) <= 1e-6:
                        continue
                    emit_ghost_ring(resname, axis_label, "SIDECHAIN", n, L_axis, tn, b_name, c_name,
                                    bn, cn, axis_dc_target=dc_pair.get(n, 0.0),
                                    group_tag=f"{bn}{cn}",
                                    ref_geom_b=b_geom[bn][:2], ref_geom_c=c_geom[cn][:2])
                    emitted += 1
        if emitted:
            n_dihedral_ok += 1
            return
        # Nothing resolved per pair: fall through to the combined construction
        # rather than silently emitting nothing for this axis.

    emit_ghost_rings_for_axis(resname, axis_label, "SIDECHAIN", L_axis, target, dc_by_harmonic, b_name, c_name,
                              b_ref, c_ref, b_geom[b_ref][:2], c_geom[c_ref][:2])

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

# Arg's guanidinium C-N rotations, past chi4's CD-NE: NOT chi axes (the
# textbook side-chain torsions stop at chi4) and not free rotamers either
# -- resonance keeps the group near-planar, the same spirit as omega
# keeping the peptide bond planar. But amber99sb.xml still parametrizes
# each of the three C-N bonds with a real generic X-CA-N2-X term (n=2,
# phase=180deg, k=10.0416 kJ/mol -- 12 real quadruplets per Arg once every
# substituent pair is counted), contributing energy that was entirely
# missed before (found 2026-08-09, same independent-OpenMM validation that
# surfaced the methyl rotors above). Every vertex involved is trigonal with
# exactly 2 non-axis substituents (NE: CD/HE; CZ: the 2 other N's; NH1/NH2:
# their own 2 H's), so all three axes route through the long-standing
# idealized_azimuth_deg path unchanged -- no new geometry needed, unlike
# the methyl case above. Requires the matching _CZ/_NH1/_NH2 split in
# ProteinAtomRigidGroups.rbody (see that file's own comment): with the old
# single rigid _GUAN clique, CZ-NH1/CZ-NH2 were frozen and these rings
# would have had nothing to act on.
# DELIBERATELY NOT GENERATED (decision 2026-08-10, after implementing and
# then measuring them). The three C-N bonds are chemically equivalent -- the
# protonated guanidinium is a fully delocalized resonance hybrid, and
# amber99sb.xml confirms it by giving all three identical parameters -- so
# there is no principled way to model two and skip one: it is all three or
# none. None was chosen, for two reasons.
#
# (1) What actually holds this group planar is AMBER's IMPROPER terms (hubs
#     NE/CZ/NH1/NH2, k=4.184/43.932/4.184/4.184 kJ/mol), and we model no
#     impropers at all -- that is the deprioritized PLANARITY work. The
#     single rigid _GUAN clique in ProteinAtomRigidGroups.rbody already
#     enforces the same planarity geometrically, exactly the argument that
#     deprioritized aromatic-ring planarity.
# (2) Measured against an independent OpenMM reference, these three axes were
#     the worst in the whole file -- CZ-NH1 off by x3.4 and CZ-NH2 outright
#     sign-flipped -- because a single ring per axis assumes the substituents
#     sit at their ideal relative azimuths and a real guanidinium is twisted
#     ~18 deg out of that. Freeing a rotation whose planarity is enforced by
#     terms we do not reproduce was the actual mistake.
#
# An earlier attempt split _GUAN into per-sp2-centre groups (_CZ/_NH1/_NH2)
# to free these rotations; that split was reverted with this. Note its own
# planarity argument was sound (each of AMBER's 4 impropers is fully
# contained in one of those groups, so every centre stayed planar) -- what it
# freed was the relative twist BETWEEN the three centres, which for a
# conjugated group is better left frozen than modelled badly.

# chi5 (NE-CZ) IS generated, and is a different case from the two CZ-NH bonds
# above -- checked directly in ProteinAtomRigidGroups.rbody rather than
# assumed: R_NE {CD,NE,CZ,HE} and R_GUAN {NE,CZ,NH1,NH2,HH*} overlap on
# exactly {NE,CZ}, which is this file's own "2 anchors, 1 free rotation"
# hinge pattern. So NE-CZ is ALREADY a free rotation under the unmodified
# rigid mesh -- unlike CZ-NH1/CZ-NH2, which sit wholly inside R_GUAN and are
# frozen, which is why rings there would have had nothing to act on. It
# carries a real 40.17 kJ/mol of AMBER torsion (4 pairs x k=10.0416, n=2,
# phase=180deg) that was governed only by --rigidbody's uniform stiffness
# until now. This is what classical nomenclature calls arginine's chi5; chi4
# (CD-NE) is a separate, also-free hinge whose AMBER k is a genuine zero.
#
# per_pair=True because all 4 of its substituent pairs reinforce in phase
# under ideal geometry -- exactly omega's signature, the condition under
# which a single combined ring assumes maximal coherence and comes out
# systematically too stiff (see generate_sidechain_axis's own docstring).
CHI5_AXIS = {"ARG": ("NE", "CZ")}
for _resname, (_b_name, _c_name) in CHI5_AXIS.items():
    generate_sidechain_axis(_resname, _b_name, _c_name, "chi5", per_pair=True)

# CZ-NH1 / CZ-NH2: the guanidinium's own two terminal-amine rotations. These
# are NOT free under the unmodified rigid mesh -- they sit wholly inside the
# single R_GUAN clique -- so they are only meaningful together with the
# _CZ/_NH1/_NH2 split in ProteinAtomRigidGroups.rbody (see that file's own
# comment). Kept in step with it deliberately: generating rings here without
# the split would put a soft torsion in competition with much stiffer rigid
# springs on an already over-constrained degree of freedom.
#
# per_pair=True is not optional here, it is the whole reason these are back:
# with a single combined ring per axis they were the worst axes in the file
# (x3.4 on CZ-NH1, outright sign flip on CZ-NH2), because all 4 substituent
# pairs reinforce in phase and the combined ring therefore bakes in maximal
# coherence while a real guanidinium is ~18 deg out of plane.
#
# Modelling choice, worth stating: a real guanidinium NH2 is neither a free
# rotor nor infinitely rigid -- it has a high but finite barrier. Governing it
# with AMBER's own 40.17 kJ/mol torsion is more faithful than freezing it into
# a uniform --stiffness clique, which is why the split is preferred now that
# per-pair anchoring makes the torsion accurate.
GUANIDINIUM_NH_AXIS = {"ARG": [("CZ", "NH1"), ("CZ", "NH2")]}
for _resname, _axes in GUANIDINIUM_NH_AXIS.items():
    for _b_name, _c_name in _axes:
        generate_sidechain_axis(_resname, _b_name, _c_name, f"guan_{_c_name}", per_pair=True)

# Asn/Gln side-chain amide rotation (CG-ND2 / CD-NE2): structurally the SAME
# torsion as the backbone peptide bond itself -- verified on a real built
# system: 4 substituent pairs at n=2, k=10.46, phase=180 deg, plus n=1 k=8.37
# carried by the (O, H) pairs, i.e. omega's exact parameter set. Frozen until
# now inside the single N_CG/Q_CD rigid cliques; requires the matching
# N_ND2/Q_NE2 split in ProteinAtomRigidGroups.rbody (see that file's own
# comment -- both AMBER impropers stay wholly inside one split group each, so
# sp2 planarity is untouched, only the C-N twist is freed). per_pair=True for
# the same reason as omega/guanidinium: all pairs reinforce in phase (the
# planar amide), so a single combined ring would bake in maximal coherence
# and come out systematically too stiff. Found via the all-20-AA GKinase
# validation: ~19 kJ/mol/instance of real AMBER torsion energy, the largest
# non-ring uncovered family (Gln 154.5 + Asn 39.0 kJ/mol/frame there).
AMIDE_AXIS = {"GLN": ("CD", "NE2"), "ASN": ("CG", "ND2")}
for _resname, (_b_name, _c_name) in AMIDE_AXIS.items():
    generate_sidechain_axis(_resname, _b_name, _c_name, f"amide_{_c_name}", per_pair=True)

# Hydroxyl/thiol/ammonium rotors (Ser OG, Thr OG1, Tyr OH, Cys SG, Lys NZ):
# the polar-hydrogen analogs of the terminal-methyl rotors above. Their
# hinges are ALREADY free in ProteinAtomRigidGroups.rbody (the S_OG/T_OG1/
# Y_OH/C_SG/K_NZ pivot groups, present from the start) -- only the energy was
# never generated, so these rotations were governed by nothing at all.
# amber99sb.xml gives each a real term (X-CT-OH-X k=0.70 n=3, X-CT-SH-X
# k=1.05 n=3, X-CT-N3-X n=3 for the ammonium; Tyr's aromatic C-OH matches
# its own specific entry, n=2 -- whatever matches, combined_target_for_axis
# resolves it). Small barriers (about kT or below), so the default combined
# ring is used, same policy as the methyls -- per_pair is reserved for the
# high-k in-phase families (omega, guanidinium, amides). Vertex geometry all
# routes through existing cases: 1-branch+2-identical (Ser/Cys/Thr/Lys CB or
# CE), all-same-class (Lys's 3 HZ), trigonal (Tyr CZ), single substituent
# (the hydroxyl H itself: it IS the side's reference atom, delta=0 by
# definition, nothing to derive).
ROTOR_AXIS = {"SER": ("CB", "OG"), "THR": ("CB", "OG1"), "TYR": ("CZ", "OH"),
              "LYS": ("CE", "NZ")}
for _resname, (_b_name, _c_name) in ROTOR_AXIS.items():
    generate_sidechain_axis(_resname, _b_name, _c_name, f"rotor_{_c_name}")
# Cys is absent from ubiquitin.pdb entirely -- sourced from the GKinase
# structure, same as its chi1 (see generate_sidechain_axis_gkinase above).

# Proline pyrrolidine ring torsions. The ring is deliberately NOT one rigid
# clique anymore (see ProteinAtomRigidGroups.rbody's own comment): real
# prolines pucker (Cgamma-endo/exo, few-kJ/mol barriers, coupled to phi),
# the real molecule has no improper there, and freezing it suppressed real
# conformational dynamics -- the largest single uncovered energy family on
# the all-20-AA validation (342 kJ/mol/frame on GKinase, 10 Pro). With the
# ring split into per-vertex hinge groups, its 4 side torsions become free
# and need their real AMBER terms (the 5th ring bond, N-CA, is phi's own
# axis and was already covered by PRO's phi rings). That the 4 dihedrals
# are coupled by ring closure is not a problem: AMBER itself treats each
# ring dihedral as an ordinary independent term, and so do these springs --
# closure and the bond/angle mesh do the rest, exactly as in a standard MD.
# CA-CB and CG-CD are pure generic n=3 aliphatic axes; CB-CG carries the
# usual 3-term CT-CT-CT-CT structure (same as any chi2) -- all verified by
# direct inspection of a real built system, and all standard combined-ring
# cases.
PRO_RING_AXIS = [("CA", "CB", "ring_CB"), ("CB", "CG", "ring_CG"), ("CG", "CD", "ring_CD")]
for _b_name, _c_name, _label in PRO_RING_AXIS:
    generate_sidechain_axis("PRO", _b_name, _c_name, _label)

def generate_pro_ring_ncd():
    """Pro's N-CD ring-closure torsion: owned by the SINGLE real pair
    (-C, CG) -- the proline-specific C-N-CT-CT AMBER entries (n=1 k=8.368,
    n=2 k=8.368, n=3 k=1.674, verified on a real built system); every
    other substituent pair around this bond is a genuine AMBER zero. So
    one ring per harmonic, anchored directly on the owning atoms, exactly
    like omega's own n=1 (O, +H) ring -- the original sole-owner rule, no
    per-pair machinery needed. Needs this dedicated emitter only because
    the -C anchor is cross-residue, which generate_sidechain_axis's
    intra-residue neighbor scan cannot see."""
    global n_dihedral_ok, n_dihedral_skip
    res = _residue_instance("PRO")
    if res.index == 0:
        print("SKIP DIHEDRAL (PRO ring_NCD): source instance is the chain start, no real -C")
        n_dihedral_skip += 1
        return
    n_atom = next(a for a in res.atoms if a.name == "N")
    cd_atom = next(a for a in res.atoms if a.name == "CD")
    cg_atom = next(a for a in res.atoms if a.name == "CG")
    prev_c = next(a for a in residues[res.index - 1].atoms if a.name == "C")
    b_class, c_class = bb_atom_class(n_atom), bb_atom_class(cd_atom)
    L_axis = bond_len_A(b_class, c_class)
    target, dc_by_harmonic = combined_target_for_axis({prev_c: (0, 0, 0.0)}, {cg_atom: (0, 0, 0.0)},
                                                      b_class, c_class, bb_atom_class)
    if target is None:
        print("SKIP DIHEDRAL (PRO ring_NCD): no AMBER entry for the (-C, CG) pair")
        n_dihedral_skip += 1
        return
    for n in sorted(target):
        if n == 0 or abs(target[n]) <= 1e-6:
            continue
        emit_ghost_ring("PRO", "ring_NCD", "SIDECHAIN", n, L_axis, target[n], "N", "CD", "-C", "CG",
                        axis_dc_target=dc_by_harmonic.get(n, 0.0),
                        ref_geom_b=(bond_len_A(bb_atom_class(prev_c), b_class),
                                    valence_deg(b_class, bb_atom_class(prev_c), c_class)),
                        ref_geom_c=(bond_len_A(bb_atom_class(cg_atom), c_class),
                                    valence_deg(c_class, bb_atom_class(cg_atom), b_class)))
    n_dihedral_ok += 1
# (called after the backbone section below -- bb_atom_class is defined there)

# "Methyl rotor" axes: NOT a chi (no rotamer identity -- see
# formula_derived_deltas' all-same-class case, and the CHI1_RESIDUES
# comment above excluding ALA for exactly this reason), but AMBER's own
# generic X-CT-CT-X wildcard (pure n=3, k=0.6508 kJ/mol, confirmed by
# direct XML inspection to have no n=1/n=2 companion term) still applies
# to every terminal methyl's own internal rotation, and contributes real,
# always-missed energy if ungenerated. Found 2026-08-09 via independent
# OpenMM/amber99sb.xml validation on real Fs-peptide trajectories (see
# doc/BondedForceFieldSprings.md and the fs-peptide validation report) --
# every one of these axes was simply never on any residue's list before,
# not a bug in an existing one. One entry per terminal CH3 in the
# standard 20 residues: Ala's only side-chain bond (CA-CB itself, unlike
# every other residue here where it's one step past an existing chi
# axis); Val's two CB branches; Thr's CG2 (alongside its already-covered
# genuine OG1 stereocenter); Ile's CG2 and its CD1 (one step past chi2);
# Leu's CD1/CD2 (one step past chi2); Met's CE (one step past chi3,
# through the CG-SD-CE thioether). Deliberately NOT extended to ACE/NME
# cap methyls (a different code path, smaller aggregate contribution,
# deferred). Each axis's "trunk" side reuses either the CA stereocenter
# or the already-solved 1-branch+2-identical case unchanged -- only the
# methyl's own 3 H's need the new all-same-class geometry, so no other
# function needed changing.
METHYL_AXIS = {"ALA": [("CA", "CB")],
               "VAL": [("CB", "CG1"), ("CB", "CG2")],
               "THR": [("CB", "CG2")],
               "ILE": [("CB", "CG2"), ("CG1", "CD1")],
               "LEU": [("CG", "CD1"), ("CG", "CD2")],
               "MET": [("SD", "CE")]}
for _resname, _axes in METHYL_AXIS.items():
    for _b_name, _c_name in _axes:
        generate_sidechain_axis(_resname, _b_name, _c_name, f"methyl_{_c_name}")

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

def generate_sidechain_axis_gkinase(resname, b_name, c_name, axis_label, per_pair=False):
    """per_pair: same semantics as generate_sidechain_axis's own flag (one
    ring per real substituent pair, for axes whose pairs all reinforce in
    phase) -- see that function's docstring."""
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
        # Formula-derived delta -- see generate_sidechain_axis's own
        # side_geom for the identical fix and its full derivation.
        vertex_name = b_name if is_b_side else c_name
        axis_partner_name = c_name if is_b_side else b_name
        formula_deltas = None
        if vertex_name == "CA":
            try:
                formula_deltas = {n: ca_tetrahedral_delta(ref_name, vertex_name, axis_partner_name, n)
                                 for n in neighbors}
            except KeyError:
                formula_deltas = None
        if formula_deltas is None:
            atoms_with_classes = [(n, resolved_class(resname, n, False, False)) for n in neighbors]
            formula_deltas = formula_derived_deltas(vertex_class, other_class, atoms_with_classes, ref_name)

        # SIGN-CORRECTED on the b-side -- see generate_sidechain_axis's own
        # side_geom for the identical fix and its derivation.
        geom = {}
        for n in neighbors:
            nc = resolved_class(resname, n, False, False)
            r = bond_len_A(nc, vertex_class)
            theta = valence_deg(vertex_class, nc, other_class)
            if formula_deltas is not None:
                delta = formula_deltas[n]
            elif n == ref_name:
                delta = 0.0
            elif is_b_side:
                delta = ref_phi - _gk_real_dihedral_deg(resname, n, b_name, c_name, c_ref)
            else:
                delta = _gk_real_dihedral_deg(resname, b_ref, b_name, c_name, n) - ref_phi
            geom[n] = (r, theta, delta)
        return geom

    b_geom = side_geom(b_neighbors, b_class, c_class, b_ref, True)
    c_geom = side_geom(c_neighbors, c_class, b_class, c_ref, False)

    class_of = lambda n: resolved_class(resname, n, False, False)
    target, dc_by_harmonic = combined_target_for_axis(b_geom, c_geom, b_class, c_class, class_of)
    if target is None:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}, GKinase): no real substituent pair "
             f"matches a specific AMBER entry and no generic wildcard either -- confirmed genuine "
             f"AMBER zero-barrier torsion, not a gap")
        n_dihedral_skip += 1
        return

    if per_pair:
        emitted = 0
        for bn in b_neighbors:
            for cn in c_neighbors:
                t_pair, dc_pair = combined_target_for_axis({bn: (0, 0, 0.0)}, {cn: (0, 0, 0.0)},
                                                           b_class, c_class, class_of)
                if t_pair is None:
                    continue
                for n, tn in sorted(t_pair.items()):
                    if n == 0 or abs(tn) <= 1e-6:
                        continue
                    emit_ghost_ring(resname, axis_label, "SIDECHAIN", n, L_axis, tn, b_name, c_name,
                                    bn, cn, axis_dc_target=dc_pair.get(n, 0.0),
                                    group_tag=f"{bn}{cn}",
                                    ref_geom_b=b_geom[bn][:2], ref_geom_c=c_geom[cn][:2])
                    emitted += 1
        if emitted:
            n_dihedral_ok += 1
            return

    emit_ghost_rings_for_axis(resname, axis_label, "SIDECHAIN", L_axis, target, dc_by_harmonic, b_name, c_name,
                              b_ref, c_ref, b_geom[b_ref][:2], c_geom[c_ref][:2])

generate_sidechain_axis_gkinase("CYS", "CA", "CB", "chi1")
generate_sidechain_axis_gkinase("CYS", "CB", "SG", "rotor_SG")  # thiol rotor, see ROTOR_AXIS above.
# Expected to SKIP on this source structure: GKinase's two CYS are 10 atoms
# each (no HG) -- disulfide-bonded cystines, no free thiol to rotate. A
# reduced-Cys rotor rule would need a source structure with a real HG;
# accepted residual (~2-3 kJ/mol per reduced Cys, X-CT-SH-X k=1.046 n=3).
generate_sidechain_axis_gkinase("TRP", "CA", "CB", "chi1")
generate_sidechain_axis_gkinase("TRP", "CB", "CG", "chi2")  # expected to correctly SKIP (aromatic
                                                            # ring, same zero-barrier family as
                                                            # Phe/Tyr/His chi2)

# ===========================================================================
# Aromatic ring torsions + PLANARITY impropers (Phe/Tyr/His/Trp). The rigid
# *_RING cliques were split into per-vertex groups (see
# ProteinAtomRigidGroups.rbody's own comment): unlike Pro's 5-ring, a 6-ring
# can fold ("boat" -- its para pairs are 1-4, invisible to bonds+angles) and
# an sp2 substituent can leave the plane at first order, so freeing these
# rings is only valid together with BOTH the ring proper torsions AND the
# improper family -- exactly the terms AMBER itself holds them with.
#
# Ring propers: every aromatic ring bond carries 4 substituent pairs, all on
# the same X-CA-CA-X-style entry (n=2, k=15.167 -- the largest k in this
# file), all reinforcing in phase in the planar ring: per_pair=True by the
# established criterion (see generate_sidechain_axis). Every ring vertex is
# trigonal with 2 non-axis substituents, so all deltas route through
# idealized_azimuth_deg unchanged. His is sourced from its real instance
# (whichever protonation it has); rings anchored on the absent protonation's
# H simply skip at deployment, same as every other missing-atom case.
AROMATIC_RING_AXIS = {
    "PHE": [("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ")],
    "TYR": [("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ")],
    "HIS": [("CG", "ND1"), ("CG", "CD2"), ("ND1", "CE1"), ("CD2", "NE2"), ("CE1", "NE2")],
}
for _resname, _axes in AROMATIC_RING_AXIS.items():
    for _b_name, _c_name in _axes:
        generate_sidechain_axis(_resname, _b_name, _c_name, f"aring_{_b_name}_{_c_name}", per_pair=True)
TRP_RING_AXIS = [("CG", "CD1"), ("CG", "CD2"), ("CD1", "NE1"), ("NE1", "CE2"), ("CD2", "CE2"),
                 ("CD2", "CE3"), ("CE2", "CZ2"), ("CE3", "CZ3"), ("CZ3", "CH2"), ("CH2", "CZ2")]
for _b_name, _c_name in TRP_RING_AXIS:
    generate_sidechain_axis_gkinase("TRP", _b_name, _c_name, f"aring_{_b_name}_{_c_name}", per_pair=True)

# PLANARITY impropers. An AMBER improper IS a 4-atom dihedral (i,j,k,l) with
# the hub listed THIRD (k bonded to i, j, l), measured about the real j-k
# bond -- so emit_ghost_ring models it as-is: axis=(j,k), refs=(i,l), one
# single-pair ring per improper (n=2, phase=180 deg). No reimplementation of
# amber99sb.xml's improper MATCHING rules (wildcards, trefoil ordering --
# the risky part task #36 flagged): the exact matched quadruplets and
# parameters are read from a real built OpenMM System per source structure,
# which uses AMBER's own assignment engine. This is a PARAMETER lookup, not
# a geometry calibration, so the formula-only rule for deltas is untouched
# (a single-pair ring has no free delta at all -- its refs ARE the anchors).
# Only side-chain hubs are emitted: the backbone C/N impropers of every
# residue stay geometrically enforced inside the _PHI/_PSI cliques, exactly
# like before. Trp's junction carbons CD2/CE2 genuinely carry no improper
# in AMBER (verified) -- nothing is fabricated for them.
def generate_planarity_impropers(resname, md_top, md_bonds):
    global n_dihedral_ok, n_dihedral_skip
    import openmm as _mm
    import openmm.app as _app
    import openmm.unit as _unit
    key = id(md_top)
    if key not in _improper_system_cache:
        omm_top = md_top.to_openmm()
        # Same free-N-terminus patch as every validation script: mdtraj's
        # standard-bond guesser misses MET1's N-H1 bond in ubiquitin.pdb.
        # Applied to the OPENMM copy only -- the generator's own md-level
        # bond graph (module-level `bonds`) must not gain a bond mid-run.
        res0 = next(iter(omm_top.residues()))
        r0_atoms = {a.name: a for a in res0.atoms()}
        if "H1" in r0_atoms and "N" in r0_atoms:
            r0_bonds = [(b[0], b[1]) for b in omm_top.bonds()
                        if b[0].residue is res0 and b[1].residue is res0]
            if not any({b[0].name, b[1].name} == {"N", "H1"} for b in r0_bonds):
                omm_top.addBond(r0_atoms["N"], r0_atoms["H1"])
        system = _app.ForceField("amber99sb.xml").createSystem(omm_top, nonbondedMethod=_app.NoCutoff,
                                                               constraints=None)
        tors = next(system.getForce(i) for i in range(system.getNumForces())
                    if system.getForce(i).__class__.__name__ == "PeriodicTorsionForce")
        bonded = set()
        for a1, a2 in md_bonds:
            bonded.add((a1.index, a2.index)); bonded.add((a2.index, a1.index))
        atoms = list(md_top.atoms)
        quads = []
        for t in range(tors.getNumTorsions()):
            i, j, k, l, per, ph, kv = tors.getTorsionParameters(t)
            a4 = [i, j, k, l]
            if not any(all((c, o) in bonded for o in a4 if o != c) for c in a4):
                continue
            quads.append((tuple(atoms[x] for x in a4), per,
                          ph.value_in_unit(_unit.radian), kv.value_in_unit(_unit.kilojoule_per_mole)))
        _improper_system_cache[key] = quads

    emitted = 0
    for (ai, aj, ak, al), per, phase, kv in _improper_system_cache[key]:
        if ak.residue.name != resname or ak.name in ("N", "C"):
            continue                       # side-chain hubs only, see above
        if any(x.residue.index != ak.residue.index for x in (ai, aj, ak, al)):
            continue                       # cross-residue = backbone territory
        # Dedup is GLOBAL (module-level, shared across calls) and keyed on the
        # HUB, not the exact quad: this function is deliberately called twice
        # for HIS (both source structures, hoping to union both protonations),
        # and a per-call set re-emitted the same 4 quads a second time --
        # duplicate DIHEDRAL rules under identical names, i.e. silently
        # DOUBLED improper energy (the identical ghosts deduplicated away,
        # which is what made it easy to miss). Hub-keyed also guards against
        # the same improper appearing with reordered peripherals in the other
        # source (same physics, different quad tuple): one hub = one improper.
        if (resname, ak.name) in _improper_emitted or kv <= 1e-6:
            continue
        _improper_emitted.add((resname, ak.name))
        names = (ai.name, aj.name, ak.name, al.name)
        L_axis = bond_len_A(resolved_class(resname, aj.name, False, False),
                            resolved_class(resname, ak.name, False, False))
        target = kv * np.exp(-1j * phase)
        cls = lambda nm: resolved_class(resname, nm, False, False)
        # names[0] (the quad's first peripheral) hangs off the HUB names[2],
        # not off the b-side anchor names[1] -- see geom_via_far_anchor.
        emit_ghost_ring(resname, f"imp_{ak.name}", "PLANARITY", per, L_axis, target,
                        names[1], names[2], names[0], names[3], axis_dc_target=kv,
                        ref_geom_b=geom_via_far_anchor(
                            L_axis, bond_len_A(cls(names[0]), cls(names[2])),
                            valence_deg(cls(names[2]), cls(names[0]), cls(names[1]))),
                        ref_geom_c=(bond_len_A(cls(names[3]), cls(names[2])),
                                    valence_deg(cls(names[2]), cls(names[3]), cls(names[1]))))
        emitted += 1
    if emitted:
        n_dihedral_ok += 1
    return emitted

_improper_system_cache = {}
_improper_emitted = set()
for _resname in ("PHE", "TYR", "HIS"):
    generate_planarity_impropers(_resname, top, bonds)
generate_planarity_impropers("TRP", _gk_top, _gk_bonds)
# His's OTHER protonation (its H sits on the other ring nitrogen): union the
# GKinase instances in too, so both HD1- and HE2-side impropers exist in the
# rule; whichever H is absent on a given deployed instance just skips.
generate_planarity_impropers("HIS", _gk_top, _gk_bonds)

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
    instance at all, so this is correct regardless of source_resname).

    Per-owning-pair anchoring (added 2026-08-04, same principle validated
    on generate_omega_axis's n=1 ring): phi's b-side neighbours are {-C, H}
    and psi's c-side neighbours are {O, +N} -- in both cases, ONE of the
    two neighbours (-C for phi, +N for psi) is itself an atom of the
    ADJACENT residue's own peptide plane, sitting right next to the (now
    separately modelled) omega axis. A single shared reference for that
    side means the deployed ghost assumes a FIXED offset between that
    cross-residue neighbour and whichever neighbour got picked as
    reference -- real (if very ordinary, bond-angle-scale) thermal noise
    then shows up as a visible energy residual simply because backbone
    torsion k's are large enough to make it visible (verified on ALA psi:
    up to ~90% error with a single shared c-side reference, even with the
    sign fix above). Splitting the CROSS-RESIDUE side into one group per
    atom, each with its own local reference and its own independently
    emitted ring(s), lets the ghost for the group anchored on that atom
    track its real position directly instead of assuming it's rigid
    relative to the other neighbour -- the SAME atom-vs-generic-anchor
    trade-off already made for omega's n=1, generalized to however many
    distinct atoms sit on the risky side. The non-risky side (phi's c-side
    {C, CB, HA}, psi's b-side {CB, HA, N} -- all directly bonded, ordinary
    substituents of the SAME vertex, no adjacent independent torsion) keeps
    a single shared reference, exactly as before: there is no equivalent
    risk there to correct."""
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

    def group_target(b_atoms, c_atoms):
        """Target for a LOCAL (b_atoms x c_atoms) submatrix, using each
        list's own first entry as this group's own local reference (a
        single-atom group is therefore anchored exactly on that atom --
        the omega-n=1-style exact case)."""
        ref_b, ref_c = b_atoms[0], c_atoms[0]
        ref_phi = real_dihedral_deg_atoms(ref_b, b_atom, c_atom, ref_c)

        # The CA stereocenter side (phi's c-side, psi's b-side -- the
        # non-cross-residue, tetrahedral side, never split into per-atom
        # groups) can be derived from AMBER's own angle tables + the
        # universal L-amino-acid chirality convention (ca_tetrahedral_
        # delta) instead of measured from any real PDB structure: found
        # and fixed 2026-08-08, same principle as generate_omega_axis's
        # idealized_azimuth_deg, extended past a planar vertex to a
        # tetrahedral one via the CIP/CORN rule (see
        # _build_ca_stereocenter_frame's own docstring for the full
        # derivation and its validation against real ubiquitin.pdb).
        # Restricted to the standard {N,C,CB,HA} substituent set -- excludes
        # Gly (no CB, HA2/HA3 instead: not a stereocenter at all, no
        # chirality to get right or wrong, falls back to the real-structure
        # measurement below unaffected by this change) and any future
        # non-standard case, rather than silently mismatching names.
        ca_side = "b" if b_atom.name == "CA" else "c" if c_atom.name == "CA" else None
        ca_atoms = b_atoms if ca_side == "b" else c_atoms if ca_side == "c" else None
        use_ca_formula = ca_atoms is not None and {a.name for a in ca_atoms} <= {"N", "C", "CB", "HA"}

        def delta_of(atom, is_b_side):
            if atom is (ref_b if is_b_side else ref_c):
                return 0.0
            if use_ca_formula and ((is_b_side and ca_side == "b") or (not is_b_side and ca_side == "c")):
                ref_name = (ref_b if is_b_side else ref_c).name
                return ca_tetrahedral_delta(ref_name, b_name, c_name, atom.name)
            if is_b_side:
                # SIGN-CORRECTED (found and fixed 2026-08-04, see
                # generate_omega_axis's own docstring for the derivation):
                # theta_atom - theta_ref_b, i.e. ref_phi MINUS the measured
                # value -- `measured - ref_phi` gives the negative of the
                # wanted delta, confirmed by direct algebraic derivation
                # and by a real-trajectory reconstruction test.
                return ref_phi - real_dihedral_deg_atoms(atom, b_atom, c_atom, ref_c)
            return real_dihedral_deg_atoms(ref_b, b_atom, c_atom, atom) - ref_phi

        b_geom = {a: (bond_len_A(bb_atom_class(a), b_class), valence_deg(b_class, bb_atom_class(a), c_class),
                     delta_of(a, True)) for a in b_atoms}
        c_geom = {a: (bond_len_A(bb_atom_class(a), c_class), valence_deg(c_class, bb_atom_class(a), b_class),
                     delta_of(a, False)) for a in c_atoms}
        target, dc_by_harmonic = combined_target_for_axis(b_geom, c_geom, b_class, c_class, bb_atom_class)
        # A ghost is the rotated image of its side's reference ATOM, so the
        # ring's radius/axial offset are that atom's own (see emit_ghost_ring).
        return target, dc_by_harmonic, ref_b, ref_c, b_geom[ref_b][:2], c_geom[ref_c][:2]

    b_is_cross = any(a.residue.index != res.index for a in b_neighbors)
    c_is_cross = any(a.residue.index != res.index for a in c_neighbors)
    if c_is_cross and not b_is_cross:
        groups = [(b_neighbors, [c]) for c in c_neighbors]
    elif b_is_cross and not c_is_cross:
        groups = [([b], c_neighbors) for b in b_neighbors]
    else:
        groups = [(b_neighbors, c_neighbors)]

    results = []
    for b_atoms, c_atoms in groups:
        target, dc_by_harmonic, ref_b, ref_c, rg_b, rg_c = group_target(b_atoms, c_atoms)
        if target is None:
            continue
        results.append((target, dc_by_harmonic, ref_b, ref_c, rg_b, rg_c))

    if not results:
        print(f"SKIP DIHEDRAL BACKBONE ({resname} {axis_label}): no real substituent pair matches "
             f"a specific AMBER Proper entry (all fall on the null wildcard)")
        n_dihedral_skip += 1
        return

    multi_group = len(results) > 1
    for target, dc_by_harmonic, ref_b, ref_c, rg_b, rg_c in results:
        _, rule_b_ref = rel_name(b_atom, ref_b)
        _, rule_c_ref = rel_name(b_atom, ref_c)
        # One atom is shared across every group (the non-risky side); the
        # OTHER is what actually distinguishes this group from the others
        # -- use whichever one varies as the disambiguating tag.
        group_tag = (ref_b.name if b_is_cross else ref_c.name) if multi_group else ""
        for n, zt in target.items():
            if n == 0 or abs(zt) < 1e-6:
                continue
            # axis_label is always exactly "phi" or "psi" here (see this
            # function's two call sites) -- .upper() gives the DihedralFamily
            # value directly (PHI/PSI each got their own runtime .msp toggle,
            # see BondedForceFieldReader.h's DihedralFamily comment).
            # dc_by_harmonic[n]: THIS group's own physical DC share for
            # harmonic n (see combined_target_for_axis's own comment) --
            # not the whole axis's combined DC dumped onto a single ring.
            emit_ghost_ring(resname, axis_label, axis_label.upper(), n, L_axis, zt, b_name, c_name, rule_b_ref,
                           rule_c_ref, axis_dc_target=dc_by_harmonic.get(n, 0.0), group_tag=group_tag,
                           ref_geom_b=rg_b, ref_geom_c=rg_c)

def generate_omega_axis(resname, source_resname=None):
    """omega: the peptide-bond torsion around C(i)-+N(i+1) -- unlike every
    other axis here, the AXIS ITSELF is cross-residue (BondedForceFieldReader
    already resolves +/- generically for atom_B/atom_C, same as it already
    does for atom_ref/atom_rotant, so no new C++ mechanism is needed).

    AMBER's real target for this axis is class-based (C/N/CT/H/O), not
    residue-specific, and has two harmonics: n=2 (k=10.46, phase=pi --
    trans/cis minima at 180/0 deg), carried EQUALLY by all 4 real
    substituent pairs (CA(i)/O(i) x CA(i+1)/H(i+1) -- same k/phase whether a
    pair matches the specific O-H entry or the generic wildcard, verified),
    and a much smaller n=1 (k=8.368, phase=0), carried by exactly ONE real
    pair, (O(i), +H(i+1)), alone -- the other 3 pairs' own match has zero k
    at n=1, a genuine AMBER asymmetry, not a simplification on our part.

    One ring per harmonic, never split further (found and fixed 2026-08-07,
    after a 2026-08-05 detour): each ring's target is the FULL
    combined_target_for_axis result (every real pair on both sides, exactly
    like every sidechain/backbone axis) -- the coherent sum of all
    contributing pairs, reproduced as a single sinusoid. A per-owning-pair
    split (one ring per real pair instead of one ring per harmonic) was
    tried for n=2 on 2026-08-05 and reverted here: it targets the same
    curve (by linearity, N rings each exactly reproducing their own pair's
    term sum to the same E(phi) as one ring reproducing the coherent sum),
    so it buys nothing physically, while multiplying the topology for no
    reason -- and since this file's rules are calibrated ONCE per resname
    from a single reference instance's real geometry and then replayed
    identically on every instance of that residue, 4 independently-scaled
    rings amplify residue-to-residue geometry variation instead of
    averaging over it the way one combined ring does: regression confirmed
    on ARG72 (gap got worse) while ARG42 improved, not a uniform fix. See
    doc/BondedForceFieldSprings.md.

    n=1's ring IS anchored directly on its own real owning pair (O, +H)
    rather than the generic CA/+CA reference used for n=2 -- this is not a
    per-pair split of one harmonic (it is still exactly one ring for n=1),
    it reflects that n=1 genuinely, physically belongs to that one real
    quartet alone (see above). n=2 has no single owner (all 4 pairs
    contribute the same term) and anchors on the generic CA/+CA reference
    instead -- the textbook definition of the omega dihedral itself.

    delta_o/delta_h (found and fixed 2026-08-08): derived from
    idealized_azimuth_deg, i.e. from AMBER's own real HarmonicAngleForce
    angles at C(i)/N(i+1) alone -- NOT measured from any real PDB
    structure. Both C(i) and N(i+1) are trigonal (exactly 2 non-axis
    substituents each: O/CA(i) at C, H/CA(i+1) at N), so this is exact and
    unambiguous (a planar 3-substituent vertex's relative azimuth is
    provably +-180 deg, no chirality/handedness ambiguity possible --
    unlike a tetrahedral vertex such as CA, see generate_backbone_axis).
    Why this matters, not just cosmetic: calibrating from ONE specific real
    instance's measured dihedral (crystal or MD frame, doesn't matter
    which) ties the generic rule to whichever structure happened to be
    used, instead of to the AMBER potential itself -- verified concretely:
    real delta_h differs by ~24 deg between ubiquitin.pdb's crystal ARG42
    and a relaxed MD frame of that SAME residue, which alone explains a
    ~22 deg phase / 10% amplitude mismatch in target[2] against the MD
    frame. Recalibrating a ring against the MD frame's own measured target
    (using the unchanged calibrate_ring/computePosition code) reproduces
    that frame's true curve to amplitude ratio 0.9998 / phase diff 0.000
    deg -- ruling out any bug in the ring construction itself; the earlier
    mismatch was 100% a calibration-source artifact, not a construction
    defect. Deriving delta from the formula instead of any one structure
    removes this dependency entirely and gives the SAME target for every
    amino acid's omega (barring Pro's missing H), matching AMBER's own
    class-based, residue-agnostic parameterization of this axis."""
    global n_dihedral_ok, n_dihedral_skip

    geom_resname = source_resname or resname
    res = _residue_instance(geom_resname)
    if res.index >= n_res - 1:
        print(f"SKIP DIHEDRAL OMEGA ({resname}): only real instance of {geom_resname} has no real "
             f"next residue to source from")
        n_dihedral_skip += 1
        return
    next_res = residues[res.index + 1]

    b_atom = next(a for a in res.atoms if a.name == "C")
    c_atom = next(a for a in next_res.atoms if a.name == "N")
    b_class = bb_atom_class(b_atom)
    c_class = bb_atom_class(c_atom)
    L_axis = bond_len_A(b_class, c_class)

    b_neighbors = real_bond_neighbors_atoms(b_atom, c_atom)   # CA(i), O(i)
    c_neighbors = real_bond_neighbors_atoms(c_atom, b_atom)   # CA(i+1), H(i+1) -- H missing if next is Pro
    if not b_neighbors or not c_neighbors:
        print(f"SKIP DIHEDRAL OMEGA ({resname}): axis C-+N has no real substituent on one side")
        n_dihedral_skip += 1
        return

    ca_b = next((a for a in b_neighbors if a.name == "CA"), None)
    o_b = next((a for a in b_neighbors if a.name == "O"), None)
    ca_c = next((a for a in c_neighbors if a.name == "CA"), None)
    h_c = next((a for a in c_neighbors if a.name == "H"), None)
    if ca_b is None or ca_c is None:
        print(f"SKIP DIHEDRAL OMEGA ({resname}): missing a real CA on one side")
        n_dihedral_skip += 1
        return

    # delta_o/delta_other_c: derived PURELY from AMBER's own real
    # HarmonicAngleForce angle values at C(i) and N(i+1) (idealized_azimuth_
    # deg), NOT measured from any real PDB structure (found and fixed
    # 2026-08-08: calibrating from a specific real instance's measured
    # dihedral -- crystal or MD frame, doesn't matter which -- ties the
    # generic rule to whichever structure happened to be used, instead of
    # to the AMBER potential itself; verified this is not a cosmetic
    # difference: real delta_h differs by ~24 deg between ubiquitin.pdb's
    # crystal instance and a relaxed MD frame of the SAME residue, which
    # alone explains a real ~22 deg phase / 10% amplitude mismatch in
    # target[2] -- confirmed by recalibrating a ring against the MD frame's
    # own measured target and finding a near-exact match (amplitude ratio
    # 0.9998, phase diff 0.000 deg), which also rules out any bug in
    # calibrate_ring/computePosition's own placement convention). CA is the
    # natural azimuth-zero reference on both sides (delta=0 trivially).
    #
    # other_c (generalizes h_c, found while wiring this in): the c-side
    # vertex (N(i+1)) always has exactly one other real substituent besides
    # CA -- H normally, but CD when the NEXT residue is Proline (its N
    # closes the pyrrolidine ring there instead of bearing an amide H, no
    # exception needed: idealized_azimuth_deg only needs THAT atom's own
    # AMBER class, whichever it is -- amber99sb.xml already has a real
    # CT-N-CT angle entry generic enough to cover CA/CD both being class
    # CT). Using h_c specifically here (as an earlier version of this fix
    # did) would silently pass None into idealized_azimuth_deg for this
    # Pro-next-residue case -- never exercised by ubiquitin.pdb's own
    # sequence, so it built and validated fine, but a latent bug for any
    # other input structure.
    other_b = next((a for a in b_neighbors if a is not ca_b), None)
    other_c = next((a for a in c_neighbors if a is not ca_c), None)
    delta_o = idealized_azimuth_deg(b_class, c_class, bb_atom_class(ca_b), bb_atom_class(other_b)) if other_b is not None else None
    delta_other_c = idealized_azimuth_deg(c_class, b_class, bb_atom_class(ca_c), bb_atom_class(other_c)) if other_c is not None else None

    def delta_of(atom, is_b_side):
        if is_b_side:
            return 0.0 if atom is ca_b else delta_o
        return 0.0 if atom is ca_c else delta_other_c

    b_geom = {a: (bond_len_A(bb_atom_class(a), b_class), valence_deg(b_class, bb_atom_class(a), c_class),
                 delta_of(a, True)) for a in b_neighbors}
    c_geom = {a: (bond_len_A(bb_atom_class(a), c_class), valence_deg(c_class, bb_atom_class(a), b_class),
                 delta_of(a, False)) for a in c_neighbors}

    target, dc_by_harmonic = combined_target_for_axis(b_geom, c_geom, b_class, c_class, bb_atom_class)
    if target is None:
        print(f"SKIP DIHEDRAL OMEGA ({resname}): no real substituent pair matches a specific/wildcard "
             f"AMBER entry")
        n_dihedral_skip += 1
        return

    # NOTE: every atom name emitted below must be relative to `res` (the
    # residue that owns this rule, i.e. `resname`'s own index) -- NOT
    # relative to whichever residue the atom itself happens to sit in.
    # Every other axis in this script has c_atom inside `res` too, so
    # rel_name(c_atom, ...) and rel_name(b_atom, ...) coincide there; omega
    # is the only axis where c_atom is itself in `next_res`, so the C-side
    # neighbours (ca_c, h_c) must be named relative to b_atom (in `res`),
    # never relative to c_atom -- using c_atom here previously produced an
    # unprefixed "CA"/"H" that silently resolved to the WRONG (this
    # residue's own) atom instead of failing loudly, since every residue
    # has one; only caught via GLY's own real unresolved-anchor warning.
    _, rule_ca_b = rel_name(b_atom, ca_b)
    _, rule_ca_c = rel_name(b_atom, ca_c)

    # Each ring below is calibrated against its OWN physical DC share
    # (dc_by_harmonic[n]: the real AMBER k's whose own periodicity is
    # exactly n -- see combined_target_for_axis's own comment), never the
    # axis's combined total dumped onto whichever ring happens to be
    # emitted first (found and fixed 2026-08-05: that bookkeeping left the
    # AXIS grand total correct but broke each ring's own curve -- e.g. the
    # n=1 ring's exact-purity shape has a fixed 2:1 DC/amplitude ratio,
    # unlike AMBER's universal 1:1 ratio for every real torsion term, so
    # subtracting the wrong constant left a systematic -amplitude bias).

    # n=2: ONE RING PER REAL SUBSTITUENT PAIR, each anchored on the two real
    # atoms its own AMBER term is written on (re-adopted 2026-08-10 after
    # being measured -- see this function's own docstring for the full
    # history, including why the same split was reverted in 2026-08-07).
    #
    # A single combined ring reproduces one function of ONE angle, which
    # equals AMBER's sum over the 4 pairs only if those pairs sit at their
    # ideal relative azimuths. Measured on real MD frames, they do not: the
    # amide nitrogen's own substituents (H vs CA) sit at 171.1 deg apart on
    # average instead of 180, which destroys part of the coherence AMBER
    # actually has. Since omega's 4 pairs all reinforce IN PHASE under the
    # ideal assumption (|target[2]| = 4k, the maximum possible), the
    # combined ring assumes maximal coherence and comes out systematically
    # too stiff: measured |ring|/|AMBER| torque ratio 1.51 median, ring
    # stiffer in 69% of instances, ~14 kJ/mol/rad mean error -- about kT
    # over a 10 deg rotation. No other axis shows this bias (psi 1.00, phi
    # 0.85), because no other axis has all its pairs reinforcing.
    #
    # Anchoring one ring per pair removes the assumption entirely rather
    # than correcting for it: each ring reproduces exactly one real AMBER
    # term evaluated on ITS OWN real dihedral, so the sum is AMBER's sum by
    # construction (to each ring's own ~0.5% leakage), whatever the real
    # geometry happens to be. This is why it is preferred over measuring the
    # real azimuths at build time and rescaling k/delta_base per instance:
    # same result, no per-instance calibration ported into C++, and it stays
    # correct if the geometry moves instead of being frozen at input.
    #
    # The 2026-08-07 revert was right on the evidence then: under IDEAL
    # geometry the two constructions are exactly equivalent by linearity, so
    # the split looked like pure added topology. What was missing was any
    # measurement of how far real geometry departs from ideal.
    # A ghost is the rotated image of its side's reference atom, so its
    # radius/axial offset are that atom's own (see emit_ghost_ring).
    rg_b = lambda a: (bond_len_A(bb_atom_class(a), b_class), valence_deg(b_class, bb_atom_class(a), c_class))
    rg_c = lambda a: (bond_len_A(bb_atom_class(a), c_class), valence_deg(c_class, bb_atom_class(a), b_class))

    n2_pairs = 0
    for b_sub in b_neighbors:
        for c_sub in c_neighbors:
            t_pair, dc_pair = combined_target_for_axis({b_sub: (0, 0, 0.0)}, {c_sub: (0, 0, 0.0)},
                                                       b_class, c_class, bb_atom_class)
            if t_pair is None or abs(t_pair.get(2, 0j)) <= 1e-6:
                continue
            _, rule_b_sub = rel_name(b_atom, b_sub)
            _, rule_c_sub = rel_name(b_atom, c_sub)
            emit_ghost_ring(resname, "omega", "OMEGA", 2, L_axis, t_pair[2], "C", "+N",
                            rule_b_sub, rule_c_sub, axis_dc_target=dc_pair.get(2, 0.0),
                            group_tag=f"{b_sub.name}{c_sub.name}",
                            ref_geom_b=rg_b(b_sub), ref_geom_c=rg_c(c_sub))
            n2_pairs += 1
    if n2_pairs == 0 and abs(target.get(2, 0j)) > 1e-6:
        # No pair resolved on its own (should not happen for a real peptide
        # bond) -- fall back to the old single combined ring rather than
        # silently emitting nothing for this harmonic.
        emit_ghost_ring(resname, "omega", "OMEGA", 2, L_axis, target[2], "C", "+N", rule_ca_b, rule_ca_c,
                       axis_dc_target=dc_by_harmonic.get(2, 0.0),
                       ref_geom_b=rg_b(ca_b), ref_geom_c=rg_c(ca_c))

    emitted_n1 = False
    if o_b is not None and h_c is not None:
        target_oh, dc_by_harmonic_oh = combined_target_for_axis({o_b: (0, 0, 0.0)}, {h_c: (0, 0, 0.0)}, b_class,
                                                                c_class, bb_atom_class)
        if target_oh is not None and abs(target_oh.get(1, 0j)) > 1e-6:
            _, rule_o_b = rel_name(b_atom, o_b)
            _, rule_h_c = rel_name(b_atom, h_c)
            emit_ghost_ring(resname, "omega", "OMEGA", 1, L_axis, target_oh[1], "C", "+N", rule_o_b, rule_h_c,
                           axis_dc_target=dc_by_harmonic_oh.get(1, 0.0),
                           ref_geom_b=rg_b(o_b), ref_geom_c=rg_c(h_c))
            emitted_n1 = True
    if not emitted_n1 and abs(target.get(1, 0j)) > 1e-6:
        # Fallback (generic CA/+CA anchor) for the rare case with no real
        # (O, +H) pair -- e.g. next residue is Pro (no amide H).
        emit_ghost_ring(resname, "omega", "OMEGA", 1, L_axis, target[1], "C", "+N", rule_ca_b, rule_ca_c,
                       axis_dc_target=dc_by_harmonic.get(1, 0.0),
                       ref_geom_b=rg_b(ca_b), ref_geom_c=rg_c(ca_c))

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
    generate_omega_axis(_resname, source_resname=_source)

# Pro's N-CD ring-closure torsion (see generate_pro_ring_ncd's own docstring
# for why it lives with the sidechain axes but is called here: its -C anchor
# needs bb_atom_class, defined just above for the backbone section).
generate_pro_ring_ncd()

lines.append("#")
lines.append("# GHOSTPARTICLE: massless virtual sites (see spn::GhostParticle), placed")
lines.append("# algebraically from 3 real anchor atoms. One M=N=n ring (2n particles)")
lines.append("# per real AMBER Fourier harmonic n of each covered axis -- see")
lines.append("# doc/BondedForceFieldSprings.md for the closed-form derivation.")
lines.append("#    type          name             resname atom_B atom_C atom_ref     r_A   theta_deg delta_deg")
lines.extend(ghostparticle_lines)
lines.append("#")
lines.append("# DIHEDRAL (proper): ghost-ghost springs, one ring group per real AMBER")
lines.append("# Fourier harmonic of each covered chi1-4 (family SIDECHAIN), phi (family PHI),")
lines.append("# psi (family PSI) and omega (family OMEGA) axis -- see")
lines.append("# doc/BondedForceFieldSprings.md. dc_offset")
lines.append("# (kJ.mol-1) is this spring's share of its axis's exact energy correction")
lines.append("# (ring construction artifact minus AMBER's own real DC, see")
lines.append("# calibrate_ring/emit_ghost_ring): BondedForceFieldReader sums it across")
lines.append("# every DIHEDRAL entry actually applied and subtracts the total only when")
lines.append("# *reporting* dihedral energy -- forces are unaffected (a constant has no")
lines.append("# gradient).")
lines.append("#    type     name                       resname family    atom_ref atom_rotant  d0_A       k        dc_offset")
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
