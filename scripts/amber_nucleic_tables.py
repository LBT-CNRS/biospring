#!/usr/bin/env python3
"""AMBER's nucleic tables: DNA.OL15 and RNA.OL3, parsed straight from the XML.

A library, not a generator. ForceFieldTables reads one of OpenMM's bundled
nucleic force fields into its bond, angle and torsion parameter tables, keyed
the way that file keys them -- amber14's nucleic files key by TYPE where
amber99sb keys by CLASS, and both spellings are accepted so this reads either
generation of file.

What used to live here -- the STRETCH/BEND emission, the sugar and base
ghost-particle rings, the CIP machinery that oriented them -- is gone with the
ghost model. `amber_tabulate_nuc.py` is the generator now, and it needs
ForceFieldTables, DATA and GLYCOSIDIC_N.

Requires `openmm` (pip install openmm), not a BioSpring build dependency.
"""

import os
import sys
import xml.etree.ElementTree as ET

import numpy as np
import openmm.app as openmm_app

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# The ghost-ring construction, shared verbatim with the protein generator.
import bonded_axes as axes

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(os.path.dirname(openmm_app.__file__), "data")

# Alternative spellings for the same physical atom. Emitting a rule under
# every spelling is free (the .rbody/.bi.ff readers skip an atom that is
# absent) and is the only thing standing between a correct model and a
# silently empty one.
# Direction matters and is easy to get backwards: the amber14 files spell
# every atom in the PDB v3 convention (OP1/OP2, H5'/H5'', H2'/H2'', HO2'),
# so the key is what the XML says and the value is the OTHER spelling --
# amber99sb.xml's, which older files and the older BioSpring rules use.
ATOM_ALIASES = {"H5'": ["H5'1"], "H5''": ["H5'2"],
                "H2'": ["H2'1"], "H2''": ["H2'2"],
                "HO2'": ["HO'2"], "OP1": ["O1P"], "OP2": ["O2P"]}


class ForceFieldTables:
    """Bond/angle/torsion parameter tables of one AMBER XML."""

    def __init__(self, path):
        root = ET.parse(path).getroot()
        self.root = root
        # Kept so the same file can be handed to OpenMM's own ForceField --
        # see base_improper_quads, which needs AMBER's assignment engine
        # rather than a reimplementation of its improper matching rules.
        self.path = path
        self.type_to_class = {t.get("name"): t.get("class")
                              for t in root.find("AtomTypes").findall("Type")}
        self.residues = {r.get("name"): r for r in root.find("Residues").findall("Residue")}

        # WHICH IDENTIFIER THE PARAMETERS ARE KEYED ON DIFFERS BY FILE, and
        # reading the wrong one fails silently -- every lookup just misses
        # and the generator writes an empty file. amber99sb.xml keys its
        # bonds and angles by CLASS (class1/class2); the amber14 nucleic
        # files key them by TYPE (type1/type2, e.g. "DNA-C"). So the file
        # itself decides, and atom_key_of below returns the matching
        # identifier for an atom.
        first_bond = root.find("HarmonicBondForce").find("Bond")
        self.keyed_by_type = first_bond.get("type1") is not None
        n = (lambda e, i: e.get("type%d" % i)) if self.keyed_by_type else (lambda e, i: e.get("class%d" % i))

        self.bond_params = {}
        for b in root.find("HarmonicBondForce").findall("Bond"):
            self.bond_params[frozenset((n(b, 1), n(b, 2)))] = (
                float(b.get("length")), float(b.get("k")))

        # The middle identifier is the vertex by AMBER convention; the two
        # outer ones are interchangeable (same physical angle read either
        # way), the vertex is not.
        self.angle_params = {}
        for a in root.find("HarmonicAngleForce").findall("Angle"):
            key = (n(a, 2), frozenset((n(a, 1), n(a, 3))))
            self.angle_params[key] = (float(a.get("angle")) * 180.0 / np.pi, float(a.get("k")))

        # (id1, id2, id3, id4) -> [(periodicity, k, phase_rad)], "" = wildcard.
        # A k of exactly 0 is AMBER stating a deliberate zero barrier, not a
        # gap, and is dropped here so it can never be mistaken for one.
        self.torsion_params = {}
        for t in root.find("PeriodicTorsionForce").findall("Proper"):
            key = tuple(n(t, i) or "" for i in (1, 2, 3, 4))
            terms, i = [], 1
            while t.get("periodicity%d" % i) is not None:
                k = float(t.get("k%d" % i))
                if k != 0.0:
                    terms.append((int(t.get("periodicity%d" % i)), k,
                                  float(t.get("phase%d" % i))))
                i += 1
            if terms:
                self.torsion_params.setdefault(key, []).extend(terms)

    def use(self):
        """Point the shared core at this force field's tables."""
        axes.configure(self.bond_params, self.angle_params, self.torsion_params)

    def atoms_of(self, resname):
        """{atom name: the identifier this file's parameters are keyed on}."""
        res = self.residues.get(resname)
        if res is None:
            return {}
        return {a.get("name"): (a.get("type") if self.keyed_by_type
                                else self.type_to_class[a.get("type")])
                for a in res.findall("Atom")}

    def bonds_of(self, resname):
        """[(name1, name2)] internal bonds, plus which atoms bond externally."""
        res = self.residues.get(resname)
        if res is None:
            return [], []
        # amber14's nucleic files name the two atoms of a bond directly
        # (atomName1/atomName2); amber99sb.xml indexes them (from/to). Both
        # spellings are accepted so this reads either generation of file.
        names = [a.get("name") for a in res.findall("Atom")]
        internal = []
        for b in res.findall("Bond"):
            if b.get("atomName1") is not None:
                internal.append((b.get("atomName1"), b.get("atomName2")))
            else:
                internal.append((names[int(b.get("from"))], names[int(b.get("to"))]))
        external = [e.get("atomName") if e.get("atomName") is not None
                    else names[int(e.get("from"))]
                    for e in res.findall("ExternalBond")]
        return internal, external

    def neighbours_of(self, resname):
        adj = {n: set() for n in self.atoms_of(resname)}
        for a, b in self.bonds_of(resname)[0]:
            adj[a].add(b)
            adj[b].add(a)
        return adj


# ---------------------------------------------------------------------------
# Emission
# ---------------------------------------------------------------------------



GLYCOSIDIC_N = {"DA": "N9", "DG": "N9", "DC": "N1", "DT": "N1",
                "A": "N9", "G": "N9", "C": "N1", "U": "N1"}
