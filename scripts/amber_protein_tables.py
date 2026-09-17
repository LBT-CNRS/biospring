#!/usr/bin/env python3
"""AMBER ff99SB's protein tables, and a synthetic residue chain to enumerate on.

A library, not a generator: it parses OpenMM's bundled amber99sb.xml into the
bond, angle and torsion parameter tables, resolves each residue's atoms to
their AMBER classes, and builds a synthetic chain from the XML residue
templates so callers can walk real connectivity without reading any structure.

What used to live here -- the STRETCH/BEND emission and the ghost-particle
dihedral rings -- is gone with the ghost model itself. `amber_tabulate.py` is
the generator now, and it needs exactly four things from this module:
`residues`, `bonds`, `resolved_class` and `root`.

Requires `openmm` (pip install openmm), not a BioSpring build dependency.

Units: amber99sb.xml already expresses k in the convention BioSpring uses
(E = 0.5*k*(q-q0)^2, OpenMM's own), so every constant is used exactly as read.
Angles are converted rad -> deg for readability, nothing else.
"""

import os
import sys
import re
import warnings
import xml.etree.ElementTree as ET

import numpy as np
import openmm.app as openmm_app

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

warnings.filterwarnings("ignore")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XML = os.path.join(os.path.dirname(openmm_app.__file__), "data", "amber99sb.xml")

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


# lookup_torsion_specific/_wildcard read these; inject them once here so a
# caller only has to import this module.
import bonded_axes as axes  # noqa: E402
from bonded_axes import lookup_torsion_specific, lookup_torsion_wildcard  # noqa: F401,E402

axes.configure(bond_params, angle_params, torsion_params)
