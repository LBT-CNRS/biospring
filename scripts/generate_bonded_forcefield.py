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
_arg42 = next(r for r in traj_ubq.topology.residues if r.name == "ARG" and r.resSeq == 42)
_atoms_by_name = {a.name: a for a in _arg42.atoms}

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

def bond_len_A(c1, c2):
    p = bond_params.get(frozenset((c1, c2)))
    return p[0] * 10.0 if p else None

def valence_deg(vertex_class, c1, c2):
    p = angle_params.get((vertex_class, frozenset((c1, c2))))
    return p[0] if p else None

def real_dihedral_deg(n1, n2, n3, n4):
    idx = [[_atoms_by_name[n].index for n in (n1, n2, n3, n4)]]
    return float(np.degrees(md.compute_dihedrals(traj_ubq, idx)[0, 0]))

def real_bond_neighbors(vertex_name, exclude_name):
    out = []
    for a1, a2 in bonds:
        if a1.residue.index != _arg42.index or a2.residue.index != _arg42.index:
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

def signed_projection(energy_over_phi_deg, n, phase_rad):
    # Real, signed component along cos(n*phi-phase) -- preserves sign/phase
    # alignment with AMBER's convention, so the linear solve below never
    # silently produces an unphysical negative stiffness without being caught.
    phi_rad = np.radians(PHI_GRID_DEG)
    return (2.0 / len(phi_rad)) * float(np.sum(energy_over_phi_deg * np.cos(n * phi_rad - phase_rad)))

def group_energy(pairs, k, L_axis, use_min_d0):
    total = np.zeros_like(PHI_GRID_DEG)
    for (bn, cn, r1, t1, r2, t2, delta) in pairs:
        d2v = closed_form_d2(L_axis, r1, t1, r2, t2, delta, PHI_GRID_DEG)
        anchor_phi = delta if use_min_d0 else 180.0 + delta
        d0 = np.sqrt(closed_form_d2(L_axis, r1, t1, r2, t2, delta, anchor_phi))
        total += 0.5 * k * (np.sqrt(np.maximum(d2v, 0.0)) - d0) ** 2
    return total

n_dihedral_ok = 0
n_dihedral_skip = 0
dihedral_lines = []

def generate_sidechain_axis(resname, b_name, c_name, axis_label):
    global n_dihedral_ok, n_dihedral_skip

    b_class = resolved_class(resname, b_name, False, False)
    c_class = resolved_class(resname, c_name, False, False)
    L_axis = bond_len_A(b_class, c_class)

    terms = lookup_torsion_wildcard(b_class, c_class)
    if terms is None:
        print(f"SKIP DIHEDRAL ({axis_label}): no generic X-{b_class}-{c_class}-X torsion entry")
        n_dihedral_skip += 1
        return

    b_neighbors = real_bond_neighbors(b_name, c_name)
    c_neighbors = real_bond_neighbors(c_name, b_name)
    b_ref, c_ref = b_neighbors[0], c_neighbors[0]
    ref_phi = real_dihedral_deg(b_ref, b_name, c_name, c_ref)

    def side_geom(neighbors, vertex_class, other_class, ref_name, is_b_side):
        geom = {}
        for n in neighbors:
            nc = resolved_class(resname, n, False, False)
            r = bond_len_A(nc, vertex_class)
            theta = valence_deg(vertex_class, nc, other_class)
            if n == ref_name:
                phi_meas = ref_phi
            elif is_b_side:
                phi_meas = real_dihedral_deg(n, b_name, c_name, c_ref)
            else:
                phi_meas = real_dihedral_deg(b_ref, b_name, c_name, n)
            geom[n] = (r, theta, phi_meas - ref_phi)
        return geom

    b_geom = side_geom(b_neighbors, b_class, c_class, b_ref, True)
    c_geom = side_geom(c_neighbors, c_class, b_class, c_ref, False)

    pairs = [(bn, cn, r1, t1, r2, t2, d2 - d1)
             for bn, (r1, t1, d1) in b_geom.items()
             for cn, (r2, t2, d2) in c_geom.items()]

    for (n, k_amber, phase_rad) in terms:
        e_max = group_energy(pairs, 1.0, L_axis, use_min_d0=False)
        proj_max = signed_projection(e_max, n, phase_rad)
        use_min_d0, proj = False, proj_max
        if proj_max <= 0:
            e_min = group_energy(pairs, 1.0, L_axis, use_min_d0=True)
            proj_min = signed_projection(e_min, n, phase_rad)
            if proj_min > 0:
                use_min_d0, proj = True, proj_min
            else:
                print(f"SKIP DIHEDRAL term ({axis_label} n={n}): no physical (positive-k) solution "
                     f"with either d0 convention -- projections {proj_max:.4g}/{proj_min:.4g}")
                n_dihedral_skip += 1
                continue

        k_solved = k_amber / proj
        e_final = group_energy(pairs, k_solved, L_axis, use_min_d0)
        reproduced = signed_projection(e_final, n, phase_rad)
        fft_final = np.fft.rfft(e_final) / len(e_final)
        residual = {m: round(2.0 * abs(fft_final[m]), 4)
                   for m in range(1, 4 * max(n, 1) + 1)
                   if m != n and 2.0 * abs(fft_final[m]) > 0.005 * abs(reproduced)}
        print(f"DIHEDRAL {axis_label} n={n} phase={np.degrees(phase_rad):.0f}deg "
             f"AMBER_k={k_amber:.4f} solved_k={k_solved:.4f} d0={'min' if use_min_d0 else 'max'} "
             f"reproduced={reproduced:.4f} residual(>0.5%)={residual}")

        for (bn, cn, r1, t1, r2, t2, delta) in pairs:
            anchor_phi = delta if use_min_d0 else 180.0 + delta
            d0 = np.sqrt(closed_form_d2(L_axis, r1, t1, r2, t2, delta, anchor_phi))
            rule_name = f"{resname}_{axis_label}_{bn}_{cn}_n{n}"
            dihedral_lines.append(
                f"DIHEDRAL {rule_name:28s} {resname:7s} SIDECHAIN {bn:5s} {cn:6s} "
                f"{d0:7.4f} {k_solved:10.4f}")
            n_dihedral_ok += 1

generate_sidechain_axis("ARG", "CA", "CB", "chi1")
generate_sidechain_axis("ARG", "CB", "CG", "chi2")

lines.append("#")
lines.append("# DIHEDRAL (proper, side chain): ghost springs reproducing Arg chi1/chi2's")
lines.append("# generic C-C rotation term, see doc/BondedForceFieldSprings.md Section 3.2.")
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

print(f"\nWrote {n_ok} unique STRETCH/BEND entries ({n_skip} skipped) and "
     f"{n_dihedral_ok} DIHEDRAL entries ({n_dihedral_skip} skipped) to {OUT}")
