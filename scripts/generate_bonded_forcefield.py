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
# class1/class2/[class3]/...), matched against the bond connectivity the
# XML residue templates state directly. STRETCH, BEND and the phi/psi/omega/
# chi rings are built on a SYNTHETIC residue chain assembled from those
# templates -- no PDB structure at all (see the _build_chain section below).
# Cys/Trp chi1, the Trp/His planarity impropers and the ACE/NME caps are the
# last three paths still reading a real structure; migrating them is
# in progress. BEND angles are found generically: any real
# atom with >=2 bonded neighbours is a valence vertex, its neighbour pairs
# are its real angles -- no residue-specific angle list needed. See
# doc/BondedForceFieldSprings.md for the full derivation and validation of
# the ghost-particle methodology (why ghost PARTICLES rather than springs
# directly between real atoms: a real-atom spring for BEND leaks bond-
# stretch noise, and for DIHEDRAL cannot always reach a real AMBER
# multi-term target without an unphysical negative stiffness).
#
# Requires `openmm` (pip install openmm) -- not a
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
# 0.5*k*(r-r0)^2 over the same bonds.
import os
import sys
import re
import warnings
import xml.etree.ElementTree as ET

import numpy as np
import openmm.app as openmm_app

# Ring construction and calibration live in a module shared with the
# nucleic-acid generator -- same code, not a copy, so a fix to the
# sign conventions or the DC bookkeeping applies to both.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bonded_rings import (PHI_GRID_DEG, closed_form_d2, complex_fourier_coeff,
                          ring_curve_abstract, choose_d0, calibrate_ring)

warnings.filterwarnings("ignore")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
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
    # A terminal residue is typed from its OWN template (NXXX/CXXX), not from
    # the interior one. This is not cosmetic: AMBER retypes atoms that exist
    # in both, and getting it wrong silently picks the neutral parameters for
    # a charged group. The C-terminal carbonyl O is class O2, not O (the
    # carboxylate is symmetric and delocalised, so O-C-OXT is really the
    # O2-C-O2 angle), and the N-terminal N is N3, not N (ammonium, not amide).
    # Reading the interior template first used to be a documented
    # approximation, forced by the terminal residues coming from whichever
    # ones a PDB happened to end with; now that they are built from AMBER's
    # own terminal templates, the approximation has no reason to stay.
    base = base_variant(resname)
    if is_nterm:
        c = atom_class("N" + base, atomname)
        if c is not None:
            return c
    if is_cterm:
        c = atom_class("C" + base, atomname)
        if c is not None:
            return c
    return atom_class(base, atomname)

# --- Synthetic residue chain built from the AMBER XML templates -----------
#
# Replaces loading real PDB structures. The generator only ever reads five
# attributes off these objects (.name/.index/.residue/.atoms/.resSeq --
# counted, nothing else), and it needs a chain so the +/- cross-residue
# neighbours resolve, so a synthetic chain carrying every residue type once
# is a complete stand-in.
#
# Why this replaced ubiquitin.pdb + GKinase/model.pdb: every parameter
# emitted already comes from the AMBER tables and the CIP rule (measured: the
# structure-reading fallback branches were taken 0 times over a full run), so
# the structures were only enumerating which atoms a residue has -- which the
# XML templates state exactly, and for EVERY residue rather than only those
# the two files happened to contain. The old arrangement silently dropped
# real parameters: Cys's thiol rotor was skipped because both cysteines in
# GKinase are disulfide-bonded and carry no HG, though AMBER defines both the
# template (with HG) and the torsion (X-CT-SH-X, k=1.046, n=3).

class _Atom:
    __slots__ = ("name", "index", "residue")
    def __init__(self, name, index, residue):
        self.name, self.index, self.residue = name, index, residue
    def __repr__(self):
        return "%s-%s%d" % (self.name, self.residue.name, self.residue.resSeq)

class _Residue:
    __slots__ = ("name", "index", "resSeq", "atoms", "is_nterm", "is_cterm")
    def __init__(self, name, index):
        self.name, self.index, self.resSeq, self.atoms = name, index, index + 1, []
        self.is_nterm = self.is_cterm = False
    def __repr__(self):
        return "%s%d" % (self.name, self.resSeq)

def _build_chain(residue_names, template_of, first_atom_index=0, first_residue_index=0):
    """One residue per name, linked C(i)-N(i+1) into a single chain.

    template_of(name) -> the XML <Residue> element to take atoms/bonds from
    (lets HIS be spelled HIS here while carrying HID's atoms, as before).
    """
    residues, bonds, counter = [], [], first_atom_index
    for i, name in enumerate(residue_names):
        tpl = template_of(name)
        res = _Residue(name, first_residue_index + i)
        names = [a.get("name") for a in tpl.findall("Atom")]
        for n in names:
            res.atoms.append(_Atom(n, counter, res))
            counter += 1
        by_name = {a.name: a for a in res.atoms}
        for b in tpl.findall("Bond"):
            f, t = b.get("from"), b.get("to")
            if f is not None:                      # amber99sb: indices
                bonds.append((res.atoms[int(f)], res.atoms[int(t)]))
            else:                                  # amber14: names
                bonds.append((by_name[b.get("atomName1")], by_name[b.get("atomName2")]))
        residues.append(res)
    # Peptide bond C(i)-N(i+1). Residues with no C or no N (caps) just skip.
    for a, b in zip(residues, residues[1:]):
        ca = next((x for x in a.atoms if x.name == "C"), None)
        nb = next((x for x in b.atoms if x.name == "N"), None)
        if ca is not None and nb is not None:
            bonds.append((ca, nb))
    return residues, bonds

# Every residue the generator asks about, in one chain. ACE/NME bracket the
# list so no residue of interest is terminal (the generator prefers interior
# instances), and the two caps are themselves generated from their own
# templates like any other residue.

class _SyntheticTopology:
    """Just enough of an mdtraj topology for generate_planarity_impropers:
    an ordered atom list and a to_openmm() that AMBER can parameterise."""
    def __init__(self, residues, bonds):
        self._residues = residues
        self._bonds = bonds
        self.atoms = [a for r in residues for a in r.atoms]

    def to_openmm(self):
        import openmm.app as _app
        top = _app.Topology()
        chain = top.addChain()
        self._omm = {}
        for r in self._residues:
            orr = top.addResidue(r.name, chain)
            for a in r.atoms:
                el = _app.element.get_by_symbol(_element_symbol(a.name))
                self._omm[a.index] = top.addAtom(a.name, el, orr)
        for a, b in self._bonds:
            top.addBond(self._omm[a.index], self._omm[b.index])
        return top


def _element_symbol(atom_name):
    # AMBER atom names start with their element except for hydrogens written
    # with a leading digit (1HB...), which no XML template uses.
    n = atom_name.lstrip("0123456789")
    return "Cl" if n.startswith("CL") else ("Na" if n.startswith("NA") else n[0])

# Each residue of interest sits between two ALA spacers, rather than in one
# long chain of all of them. That matters because several rules depend on the
# NEIGHBOURING residue: omega's n=1 harmonic is carried by the (O, +H) pair,
# so a residue followed by proline -- which has no amide H -- resolves
# differently and emits a differently-named rule. Chaining the 20 types back
# to back would just swap ubiquitin's arbitrary sequence for this list's
# arbitrary order; ALA spacers give every residue the same, neutral
# environment. _residue_instance takes the first interior instance, so the
# spacers themselves are never the ALA that ALA's own rules come from --
# hence the explicit ALA triplet first in the list.
_INTERIOR = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
             "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
             "TYR", "VAL"]
CHAIN_RESIDUES = (["ACE"]
                  + [r for name in _INTERIOR for r in ("ALA", name, "ALA")]
                  + ["NME"])

# A second, separate chain carrying the terminal variants: a free N-terminus
# (H1/H2/H3 ammonium) and a free C-terminus (OXT carboxylate). Both exist in
# AMBER for all 20 residues (NALA..NVAL, CALA..CVAL) but appear in no capped
# chain, so their bonds and angles would otherwise never be emitted -- they
# were previously picked up only because ubiquitin happens to have a free
# terminus at each end.
TERMINAL_RESIDUES = [t + r for t in ("N", "C") for r in _INTERIOR]

_xml_residues = {r.get("name"): r for r in root.find("Residues").findall("Residue")}
residues, bonds = _build_chain(CHAIN_RESIDUES, lambda n: _xml_residues[base_variant(n)])
_main_residues, _main_bonds = list(residues), list(bonds)

# One free N-terminal and one free C-terminal instance of every amino acid,
# each as its own two-residue chain so the terminal one really is at an end.
# Their templates are AMBER's own NXXX/CXXX (NALA carries H1/H2/H3, CALA
# carries OXT), so this is still tables-only -- it just stops the terminal
# rules from depending on which residues a particular PDB happened to end
# with.
for _base in _INTERIOR:
    for _prefix, _nterm in (("N", True), ("C", False)):
        _tpl = _prefix + base_variant(_base)
        if _tpl not in _xml_residues:
            continue
        _pair = [_base, "ALA"] if _nterm else ["ALA", _base]
        _extra, _extra_bonds = _build_chain(
            _pair, lambda n, t=_tpl, b=_base: _xml_residues[t if n == b else base_variant(n)],
            first_atom_index=sum(len(r.atoms) for r in residues),
            first_residue_index=len(residues))
        _term = _extra[0] if _nterm else _extra[-1]
        _term.is_nterm, _term.is_cterm = _nterm, not _nterm
        residues.extend(_extra)
        bonds.extend(_extra_bonds)
_hie_residues, _hie_bonds = _build_chain(
    ["ACE", "ALA", "HIS", "ALA", "NME"],   # capped: OpenMM rejects a bare chain
    lambda n: _xml_residues["HIE" if n == "HIS" else base_variant(n)],
    first_atom_index=sum(len(r.atoms) for r in residues),
    first_residue_index=len(residues))
_hie_top = _SyntheticTopology(_hie_residues, _hie_bonds)

n_res = len(residues)

def _is_nterm(res):
    """True for a FREE N-terminus (NH3+), i.e. one built from AMBER's NXXX
    template -- not for the ACE cap, which is not a free terminus."""
    return res.is_nterm


def _is_cterm(res):
    """True for a FREE C-terminus (COO-), built from AMBER's CXXX template."""
    return res.is_cterm

top = _SyntheticTopology(_main_residues, _main_bonds)

# No bond guessing to patch: the XML templates state every bond.

lines = [
    "# STRETCH (real 1-2 bonds), BEND (real 1-3 valence angles) and DIHEDRAL",
    "# ghost-ring springs: real AMBER ff99SB parameters (OpenMM",
    "# amber99sb.xml), for the 18 residue types, both free termini of each,",
    "# and the ACE/NME caps -- all enumerated from amber99sb.xml's own",
    "# residue templates, no structure.",
    "# Generated by scripts/generate_bonded_forcefield.py -- do not",
    "# hand-edit individual values, regenerate instead. BEND's theta0/k are",
    "# converted into a 1-3 distance spring's equilibrium/stiffness at build",
    "# time by BondedForceFieldReader, not by this script.",
    "#",
    "# STRETCH and BEND were removed for a while, because the ghost rings",
    "# applied 96 % of their force as radial and axial leak that deformed a",
    "# flexible mesh instead of twisting it. dihedral.tangentialonly fixes",
    "# that at the source, so soft bonds and angles are viable again -- see",
    "# doc/BondedForceFieldSprings.md.",
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

# Proline's ring stays rigid, so nothing retunes what is inside it: the
# .rbody's uniform springs are what hold those 5 atoms (the measurement
# behind this is with PRO_RING_AXIS, further down).
# A "+"/"-" prefixed atom belongs to another residue and is never
# ring-internal.
PRO_RING_ATOMS = frozenset(("N", "CA", "CB", "CG", "CD"))


def _pro_ring_internal(resname, *atoms):
    return resname == "PRO" and all(a in PRO_RING_ATOMS for a in atoms)

def emit_stretch(rule_name, resname, atom1, atom2_display, r0_A, k_biospring):
    global n_ok
    if _pro_ring_internal(resname, atom1, atom2_display):
        return
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
    if _pro_ring_internal(resname, atom1, atom2, atom3_display):
        return
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
    # Terminality is a property of the residue, not of its position in
    # this list: the chain below carries a free N-terminal and a free
    # C-terminal instance of EVERY amino acid, where a real PDB only
    # ever gave two (whichever residues happened to sit at its ends --
    # ubiquitin gave MET and GLY, so the other 18 got no terminal
    # rules at all).
    is_nterm = getattr(res, "is_nterm", ridx == 0)
    is_cterm = getattr(res, "is_cterm", ridx == n_res - 1)

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
    class_resolver=lambda a: resolved_class(a.residue.name, a.name,
                                            _is_nterm(a.residue), _is_cterm(a.residue)),
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

# The ghost-ring construction itself lives in bonded_axes.py, shared verbatim
# with the nucleic-acid generator: same comb/LCM theorem, same closed-form
# d(phi), same linear solve for k, same table-derived azimuths. Only the
# chemistry that genuinely differs stays here.
import bonded_axes as axes
from bonded_axes import (bond_len_A, valence_deg, idealized_azimuth_deg,        # noqa: F401
                         lookup_torsion_wildcard, lookup_torsion_specific,
                         formula_derived_deltas, geom_via_far_anchor,
                         combined_target_for_axis, emit_ghost_ring,
                         emit_ghost_rings_for_axis, _dihedral_from_points,
                         _build_two_vectors, _build_symmetric_tetrahedral_frame)

axes.configure(bond_params, angle_params, torsion_params)

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

# Thr's C-beta is the only genuinely chiral tetrahedral centre here that the
# angle tables cannot resolve on their own (the backbone C-alpha is handled
# by its own CIP frame at the call site, not through the generic solver).
# Set after the tables above because it is defined in terms of them.
axes.stereocenter_delta = thr_cb_tetrahedral_delta

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

    b_class = resolved_class(resname, b_name, False, False)
    c_class = resolved_class(resname, c_name, False, False)
    L_axis = bond_len_A(b_class, c_class)

    b_neighbors = real_bond_neighbors(resname, b_name, c_name)
    c_neighbors = real_bond_neighbors(resname, c_name, b_name)
    if not b_neighbors or not c_neighbors:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}): axis {b_name}-{c_name} has no real "
             f"substituent on one side")
        axes.bump_axis_skip()
        return
    b_ref, c_ref = b_neighbors[0], c_neighbors[0]

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
                # A vertex whose only non-axis substituent IS the reference
                # (Met's SD, for one) needs no solver: its azimuth is the
                # origin by definition.
                delta = 0.0
            else:
                # Every other vertex in this file is one of the solvable
                # cases; a new one must be DERIVED, never measured off a
                # structure -- that is what this whole generator is for.
                raise RuntimeError(
                    f"{resname} {axis_label}: no formula-derived azimuth for "
                    f"vertex {b_name if is_b_side else c_name} -- see "
                    f"formula_derived_deltas for the cases it solves")
            geom[n] = (r, theta, delta)
        return geom

    b_geom = side_geom(b_neighbors, b_class, c_class, b_ref, True)
    c_geom = side_geom(c_neighbors, c_class, b_class, c_ref, False)

    class_of = lambda n: resolved_class(resname, n, False, False)
    target, dc_by_harmonic = combined_target_for_axis(b_geom, c_geom, b_class, c_class, class_of)
    if target is None:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}): no real substituent pair matches a "
             f"specific AMBER entry and no generic wildcard either")
        axes.bump_axis_skip()
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
            axes.bump_axis_ok()
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

# Proline's pyrrolidine ring is left RIGID: no torsion about any of its
# bonds, and no STRETCH/BEND inside it either (see PRO_RING_ATOMS) -- the
# .rbody's uniform springs keep the 5 atoms at the shape the input
# structure had.
#
# This reverses the 2026-08-12 decision to free it, and the reason is
# measured rather than argued. A ghost ring's energy is not a function of
# its dihedral alone: the ghost is a REAL atom rotated about the axis, so
# deforming the ring's angles moves the ghosts, and there is a relaxation
# channel where every ghost-ghost distance reaches its d0 and the torsion
# energy collapses. One deformation relaxes every ring axis at once, so
# the gain outgrows what bonds and angles charge for it, and the ring
# flattens. Quenched, ubiquitin's three prolines went from 0.130 A out of
# plane to 0.021 -- flat -- and did so identically with every version of
# this file since the ring was freed. Freeing it therefore did not restore
# the pucker it was meant to restore; it produced a flat ring instead.
#
# Aromatic rings are NOT frozen and keep their torsions and impropers:
# for them the idealised geometry IS the real one, so the same channel
# leads to the correct structure (measured: 0.000 A out of plane after the
# same quench). Only SATURATED rings, whose truth is puckered while their
# ideal is flat, are affected.
PRO_RING_AXIS = [("CA", "CB", "ring_CB"), ("CB", "CG", "ring_CG"), ("CG", "CD", "ring_CD")]

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
    res = _residue_instance("PRO")
    if _is_nterm(res):
        print("SKIP DIHEDRAL (PRO ring_NCD): source instance is the chain start, no real -C")
        axes.bump_axis_skip()
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
        axes.bump_axis_skip()
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
    axes.bump_axis_ok()
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
# The two caps' methyls are hindered rotors in amber99sb exactly like the
# side-chain ones, but they sit on the caps' own axes: ACE's methyl turns
# about C-CH3, NME's about N-C. Without a rule they rotate FREELY -- the
# rigid mesh leaves them one degree of freedom by design (ACE_C and ACE_CH3
# overlap on {C, CH3}, a hinge) and nothing then supplies the barrier.
METHYL_AXIS = {"ACE": [("C", "CH3")],
               "NME": [("N", "CH3")],
               "ALA": [("CA", "CB")],
               "VAL": [("CB", "CG1"), ("CB", "CG2")],
               "THR": [("CB", "CG2")],
               "ILE": [("CB", "CG2"), ("CG1", "CD1")],
               "LEU": [("CG", "CD1"), ("CG", "CD2")],
               "MET": [("SD", "CE")]}
for _resname, _axes in METHYL_AXIS.items():
    for _b_name, _c_name in _axes:
        generate_sidechain_axis(_resname, _b_name, _c_name, f"methyl_{_c_name}")

generate_sidechain_axis("CYS", "CA", "CB", "chi1")
generate_sidechain_axis("CYS", "CB", "SG", "rotor_SG")  # thiol rotor, see ROTOR_AXIS above.
# The thiol rotor is real here: the synthetic chain's CYS carries HG, so
# X-CT-SH-X (k=1.046, n=3) is emitted. It used to be skipped because both
# cysteines in the source structure were disulfide-bonded cystines.
generate_sidechain_axis("TRP", "CA", "CB", "chi1")
generate_sidechain_axis("TRP", "CB", "CG", "chi2")  # expected to correctly SKIP (aromatic
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
    generate_sidechain_axis("TRP", _b_name, _c_name, f"aring_{_b_name}_{_c_name}", per_pair=True)

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
                                    valence_deg(cls(names[2]), cls(names[3]), cls(names[1]))),
                        dc_align="minimum")
        emitted += 1
    if emitted:
        axes.bump_axis_ok()
    return emitted

_improper_system_cache = {}
_improper_emitted = set()
for _resname in ("PHE", "TYR", "HIS"):
    generate_planarity_impropers(_resname, top, bonds)
generate_planarity_impropers("TRP", top, bonds)
# His's OTHER protonation (its H sits on the other ring nitrogen): union the
# HIE chain in too, so both HD1- and HE2-side impropers exist in the rule;
# whichever H is absent on a given deployed instance just skips.
generate_planarity_impropers("HIS", _hie_top, _hie_bonds)

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

def bb_atom_class(a):
    return resolved_class(a.residue.name, a.name, _is_nterm(a.residue), _is_cterm(a.residue))

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

    geom_resname = source_resname or resname
    res = _residue_instance(geom_resname)
    if _is_nterm(res) or _is_cterm(res):
        print(f"SKIP DIHEDRAL BACKBONE ({resname} {axis_label}): only real ubiquitin instance "
             f"of {geom_resname} is a chain terminus, no real neighbouring residue to source from")
        axes.bump_axis_skip()
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
        axes.bump_axis_skip()
        return

    def group_target(b_atoms, c_atoms):
        """Target for a LOCAL (b_atoms x c_atoms) submatrix, using each
        list's own first entry as this group's own local reference (a
        single-atom group is therefore anchored exactly on that atom --
        the omega-n=1-style exact case)."""
        ref_b, ref_c = b_atoms[0], c_atoms[0]

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

        # Glycine's CA is the one backbone vertex the CIP path above cannot
        # type: it carries HA2/HA3 where every other residue has CB/HA, so
        # use_ca_formula's name test rejects it. It is also not a
        # stereocenter at all -- its two hydrogens are chemically identical
        # -- which makes it exactly the symmetric-tetrahedral case
        # formula_derived_deltas already solves for every chi axis. Both
        # sides go through it here, same convention as the chi caller.
        def _side_formula(atoms, ref_atom, vertex_class, other_class, is_ca_side):
            if use_ca_formula and is_ca_side:
                return None
            return formula_derived_deltas(vertex_class, other_class,
                                          [(a, bb_atom_class(a)) for a in atoms], ref_atom)

        b_formula = _side_formula(b_atoms, ref_b, b_class, c_class, ca_side == "b")
        c_formula = _side_formula(c_atoms, ref_c, c_class, b_class, ca_side == "c")

        def delta_of(atom, is_b_side):
            if atom is (ref_b if is_b_side else ref_c):
                return 0.0
            if use_ca_formula and ((is_b_side and ca_side == "b") or (not is_b_side and ca_side == "c")):
                ref_name = (ref_b if is_b_side else ref_c).name
                # VERTEX first, then its axis partner -- the same convention
                # side_geom uses (see generate_sidechain_axis) and the one
                # combined_target_for_axis assumes: each side's azimuths are
                # measured about its OWN outward axis. Passing (b_name,
                # c_name) unconditionally measured the c-side about B->C
                # instead of C->B, i.e. negated, which is only ever the case
                # for phi (the sole backbone axis whose CA sits on the c
                # side). It cancelled exactly against the pair-combination
                # sign error in combined_target_for_axis, so both were
                # invisible until either was fixed alone.
                vertex_name = b_name if is_b_side else c_name
                partner_name = c_name if is_b_side else b_name
                return ca_tetrahedral_delta(ref_name, vertex_name, partner_name, atom.name)
            side_formula = b_formula if is_b_side else c_formula
            if side_formula is None:
                raise RuntimeError(
                    f"{resname} {axis_label}: no formula-derived azimuth for the "
                    f"{'b' if is_b_side else 'c'} side -- see formula_derived_deltas")
            return side_formula[atom]

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
        axes.bump_axis_skip()
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

    geom_resname = source_resname or resname
    res = _residue_instance(geom_resname)
    if res.index >= n_res - 1:
        print(f"SKIP DIHEDRAL OMEGA ({resname}): only real instance of {geom_resname} has no real "
             f"next residue to source from")
        axes.bump_axis_skip()
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
        axes.bump_axis_skip()
        return

    ca_b = next((a for a in b_neighbors if a.name == "CA"), None)
    o_b = next((a for a in b_neighbors if a.name == "O"), None)
    ca_c = next((a for a in c_neighbors if a.name == "CA"), None)
    h_c = next((a for a in c_neighbors if a.name == "H"), None)
    if ca_b is None or ca_c is None:
        print(f"SKIP DIHEDRAL OMEGA ({resname}): missing a real CA on one side")
        axes.bump_axis_skip()
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
        axes.bump_axis_skip()
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
# A backbone axis needs N, CA and C. The caps have none of that (ACE is
# CH3-C=O, NME is N-CH3), and never reached this loop before because no
# ubiquitin residue was a cap -- with the synthetic chain they do, so the
# precondition is now stated rather than left to the input's contents.
def _has_backbone(resname):
    names = {a.name for a in _residue_instance(resname).atoms}
    return {"N", "CA", "C"} <= names

BACKBONE_RESIDUES = sorted(r for r in {r.name for r in residues} if _has_backbone(r))
GENERIC_BACKBONE_SOURCE = {"MET": "LEU", "CYS": "LEU", "TRP": "LEU"}
for _resname in BACKBONE_RESIDUES + ["CYS", "TRP"]:
    _source = GENERIC_BACKBONE_SOURCE.get(_resname)
    generate_backbone_axis(_resname, "N", "CA", "phi", source_resname=_source)
    generate_backbone_axis(_resname, "CA", "C", "psi", source_resname=_source)
    generate_omega_axis(_resname, source_resname=_source)

# Pro's N-CD ring-closure torsion is frozen with the rest of the ring (see
# PRO_RING_ATOMS); generate_pro_ring_ncd is kept, unused, because it is the
# only place the proline-specific C-N-CT-CT entries are worked out and it
# is what would be called again if the ring were ever freed.

lines.append("#")
lines.append("# GHOSTPARTICLE: massless virtual sites (see spn::GhostParticle), placed")
lines.append("# algebraically from 3 real anchor atoms. One M=N=n ring (2n particles)")
lines.append("# per real AMBER Fourier harmonic n of each covered axis -- see")
lines.append("# doc/BondedForceFieldSprings.md for the closed-form derivation.")
lines.append("#    type          name             resname atom_B atom_C atom_ref     r_A   theta_deg delta_deg")
lines.extend(axes.ghostparticle_lines)
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
lines.extend(axes.dihedral_lines)

# ACE/NME under the PDB spelling of their methyls. The templates ACE/NME are
# already in the synthetic chain, so their bonds and angles come out of the
# main loop above like any other residue -- but amber99sb.xml names the cap
# methyls CH3/HH31/HH32/HH33, while PDB files (and every reader that loads
# one) write NME's methyl carbon C and all three cap hydrogens H1/H2/H3.
# Same atoms, different convention. A rule listing only one spelling matches
# NOTHING on a file using the other, with no error at all -- the same silent
# class as His's HD1/HE2 above. Emitting both is free: the reader skips
# whichever atom is absent.
#
# Replaces reading the Fs-peptide by absolute path from outside the repo,
# which is what previously supplied these names (and only these -- the
# chemistry was always the tables').
CAP_ATOM_ALIASES = {
    "ACE": {"HH31": "H1", "HH32": "H2", "HH33": "H3"},
    "NME": {"CH3": "C", "HH31": "H1", "HH32": "H2", "HH33": "H3"},
}

def emit_cap_name_aliases():
    def alias(resname, name):
        # "+X"/"-X" names a NEIGHBOURING residue's atom, so it is never
        # renamed here: NME's "-C" is the previous residue's carbonyl, not
        # this cap's own methyl, even though the alias maps CH3 -> C.
        prefix = name[0] if name[:1] in "+-" else ""
        bare = name[len(prefix):]
        return prefix + CAP_ATOM_ALIASES[resname].get(bare, bare)

    for line in list(lines):
        f = line.split()
        if not f or f[0] not in ("STRETCH", "BEND"):
            continue
        resname = f[2]
        if resname not in CAP_ATOM_ALIASES:
            continue
        atoms = f[3:-2]
        renamed = [alias(resname, a) for a in atoms]
        if renamed == atoms:
            continue
        rule = "%s_%s" % (resname, "_".join(a.lstrip("+-") for a in renamed))
        if f[0] == "STRETCH":
            emit_stretch(rule, resname, renamed[0], renamed[1], float(f[-2]), float(f[-1]))
        else:
            emit_bend(rule, resname, renamed[0], renamed[1], renamed[2], float(f[-2]), float(f[-1]))

emit_cap_name_aliases()

with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")

print(f"\nWrote {axes.n_ghost_particles} GHOSTPARTICLE entries and {axes.n_dihedral_ok} "
     f"DIHEDRAL entries ({axes.n_dihedral_skip} skipped) and {n_ok} STRETCH/BEND entries "
     f"({n_skip} skipped) to {OUT}")
