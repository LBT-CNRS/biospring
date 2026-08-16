#!/usr/bin/env python3
"""Generates DNA and RNA bonded-force-field rules (.bi.ff) from AMBER.

Companion to generate_bonded_forcefield.py, which does the same for
proteins; the ring construction they share lives in bonded_rings.py.

Two differences from the protein generator, both deliberate:

  * Parameters come from the NUCLEIC-ACID-SPECIFIC force fields, not
    amber99sb.xml. DNA uses OL15 (ff99 + parmbsc0 + the OL refinements of
    epsilon/zeta, chi and beta -- parmbsc0 is what fixed the alpha/gamma
    flips that used to wreck long DNA runs); RNA uses OL3 (ff99 + bsc0 +
    chiOL3, the glycosidic correction that made RNA simulable at all).
    amber99sb.xml does carry nucleic templates, but with the uncorrected
    ff99 torsions -- using them would silently reintroduce exactly the
    pathologies those force fields exist to remove.

  * Atoms and connectivity are read straight from the XML residue
    templates, with NO reference structure. The protein generator picks a
    real residue instance out of ubiquitin/GKinase to enumerate atoms; here
    the templates carry the full bond graph already (checked: DA is 32
    atoms / 34 internal bonds / 3 rings), so nothing needs measuring. That
    is strictly better and matches the rule this project settled on --
    derive from the tables, never from a structure.

Residue naming: OL15 spells DNA DA/DC/DG/DT and OL3 spells RNA A/C/G/U,
each with 5'/3'/nucleoside variants. Rules are emitted under BOTH the OL3
spelling and amber99sb's RA/RC/RG/RU, because real files use either and a
rule that knows only one matches nothing at all in the other -- measured, 0
springs out of 963 particles on a real rRNA fragment. Same reason the sugar
hydrogens are emitted under both H2'1/H2'2 (amber99sb) and H2'/H2'' (PDB
v3) spellings.
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

lines = []
counters = {"stretch": 0, "bend": 0, "skip": 0}
_seen = set()


def _aliases(name):
    return [name] + ATOM_ALIASES.get(name, [])


def emit_stretch(resname, a1, a2, r0_A, k):
    for n1 in _aliases(a1):
        for n2 in _aliases(a2):
            key = (resname, tuple(sorted((n1, n2.lstrip("+-")))), n2[:1] in "+-")
            if key in _seen:
                continue
            _seen.add(key)
            lines.append(f"STRETCH {resname}_{n1}_{n2.lstrip('+-'):<12s} {resname:7s} "
                         f"{n1:6s} {n2:7s} {r0_A:7.4f} {k:10.2f}")
            counters["stretch"] += 1


def emit_bend(resname, a1, vertex, a3, theta0, k):
    for n1 in _aliases(a1):
        for n3 in _aliases(a3):
            key = (resname, vertex, tuple(sorted((n1, n3))), "bend")
            if key in _seen:
                continue
            _seen.add(key)
            lines.append(f"BEND {resname}_{n1}_{vertex}_{n3.lstrip('+-'):<10s} {resname:7s} "
                         f"{n1:6s} {vertex:6s} {n3:7s} {theta0:7.3f} {k:10.4f}")
            counters["bend"] += 1


# ---------------------------------------------------------------------------
# Per-residue generation
# ---------------------------------------------------------------------------

# The 3'-to-5' link: O3' of residue i bonds to P of residue i+1. Written from
# the OWNING residue as "+P", exactly like the protein generator writes psi's
# "+N" -- the .bi.ff reader resolves +/- against the chain.
NEXT_P = "+P"
PREV_O3 = "-O3'"


def generate_residue(ff, resname, emit_as):
    """Emits STRETCH and BEND for one residue template.

    emit_as is the list of residue names the rules are written under (the
    same physical residue is spelled differently by different force fields
    and PDB generations -- see the module docstring).
    """
    atoms = ff.atoms_of(resname)
    if not atoms:
        print(f"SKIP residue {resname}: no template")
        return
    internal, external = ff.bonds_of(resname)
    adj = ff.neighbours_of(resname)

    # Internal bonds, plus the one cross-residue bond O3'(i) -> P(i+1).
    bonds = list(internal)
    if "O3'" in atoms and "P" in atoms:
        bonds.append(("O3'", NEXT_P))

    for a1, a2 in bonds:
        c1 = atoms[a1]
        c2 = atoms[a2.lstrip("+-")] if a2.lstrip("+-") in atoms else None
        if c2 is None:
            continue
        params = ff.bond_params.get(frozenset((c1, c2)))
        if params is None:
            print(f"SKIP STRETCH {resname} {a1}({c1})-{a2}({c2}): no HarmonicBondForce entry")
            counters["skip"] += 1
            continue
        length_nm, k_nm2 = params
        for name in emit_as:
            emit_stretch(name, a1, a2, length_nm * 10.0, k_nm2 / 100.0)

    # Angles: every pair of neighbours around each vertex. The cross-residue
    # ones (X-O3'-+P and O3'-+P-...) belong to the next residue's own rules,
    # so only intra-residue vertices are walked here -- same split the
    # protein generator uses.
    for vertex, nbrs in adj.items():
        nb = sorted(nbrs)
        for i in range(len(nb)):
            for j in range(i + 1, len(nb)):
                a1, a3 = nb[i], nb[j]
                key = (atoms[vertex], frozenset((atoms[a1], atoms[a3])))
                params = ff.angle_params.get(key)
                if params is None:
                    print(f"SKIP BEND {resname} {a1}-{vertex}-{a3}: no HarmonicAngleForce entry")
                    counters["skip"] += 1
                    continue
                theta0, k = params
                for name in emit_as:
                    emit_bend(name, a1, vertex, a3, theta0, k)



# ---------------------------------------------------------------------------
# Dihedrals: the same ghost-ring construction as the protein generator, from
# the same shared core (bonded_axes.py). Nothing here re-derives anything --
# only the chemistry differs.
# ---------------------------------------------------------------------------

# Standard nucleic torsions, named as the literature names them, each written
# as the BOND it turns about (the ghost-ring construction is per-axis, not
# per-4-atom-tuple):
#     alpha   O3'(i-1)-P-O5'-C5'        axis P-O5'
#     beta    P-O5'-C5'-C4'             axis O5'-C5'
#     gamma   O5'-C5'-C4'-C3'           axis C5'-C4'
#     delta   C5'-C4'-C3'-O3'           axis C4'-C3'   (a RING bond: this is
#                                       the pucker's readout, not a free hinge)
#     epsilon C4'-C3'-O3'-P(i+1)        axis C3'-O3'
#     zeta    C3'-O3'-P(i+1)-O5'(i+1)   axis O3'-P(i+1) -- the only one whose
#                                       axis crosses a residue boundary
BACKBONE_AXES = [("P", "O5'", "alpha"), ("O5'", "C5'", "beta"),
                 ("C5'", "C4'", "gamma"), ("C4'", "C3'", "delta"),
                 ("C3'", "O3'", "epsilon"), ("O3'", NEXT_P, "zeta")]

# The furanose's other four bonds. They matter for the same reason the
# protein's aromatic rings needed theirs: with the ring's angles fixed and
# nothing else, the pucker cannot move. These are what --bending acts on.
SUGAR_RING_AXES = [("C1'", "C2'", "ring_C1_C2"), ("C2'", "C3'", "ring_C2_C3"),
                   ("C4'", "O4'", "ring_C4_O4"), ("O4'", "C1'", "ring_O4_C1")]

# chi, the glycosidic torsion: O4'-C1'-N9-C4 on a purine, O4'-C1'-N1-C2 on a
# pyrimidine. Plus RNA's 2'-OH rotor.
GLYCOSIDIC_N = {"DA": "N9", "DG": "N9", "DC": "N1", "DT": "N1",
                "A": "N9", "G": "N9", "C": "N1", "U": "N1"}


def neighbours(ff, resname, vertex, exclude):
    """Real bonded substituents of `vertex` other than `exclude`, in template
    order -- the same stable, structure-free ordering the protein generator
    gets from its synthetic chain."""
    adj = ff.neighbours_of(resname)
    order = list(ff.atoms_of(resname))
    out = [n for n in adj.get(vertex, ()) if n != exclude.lstrip("+-")]
    return sorted(out, key=order.index)


def generate_axis(ff, resname, b_name, c_name, axis_label, family, emit_as):
    """One axis, one call -- mirrors generate_sidechain_axis exactly."""
    atoms = ff.atoms_of(resname)
    b_bare, c_bare = b_name.lstrip("+-"), c_name.lstrip("+-")
    if b_bare not in atoms or c_bare not in atoms:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}): axis atom missing")
        axes.bump_axis_skip()
        return
    b_class, c_class = atoms[b_bare], atoms[c_bare]
    L_axis = axes.bond_len_A(b_class, c_class)
    if L_axis is None:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}): no bond length for the axis")
        axes.bump_axis_skip()
        return

    b_neighbors = neighbours(ff, resname, b_bare, c_bare)
    c_neighbors = neighbours(ff, resname, c_bare, b_bare)
    if not b_neighbors or not c_neighbors:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}): no real substituent on one side")
        axes.bump_axis_skip()
        return
    b_ref, c_ref = b_neighbors[0], c_neighbors[0]

    def side_geom(nbrs, vertex_class, other_class, ref_name):
        with_classes = [(n, atoms[n]) for n in nbrs]
        deltas = axes.formula_derived_deltas(vertex_class, other_class, with_classes, ref_name)
        geom = {}
        for n in nbrs:
            nc = atoms[n]
            r = axes.bond_len_A(nc, vertex_class)
            theta = axes.valence_deg(vertex_class, nc, other_class)
            if r is None or theta is None:
                return None
            if deltas is not None:
                delta = deltas[n]
            elif n == ref_name:
                delta = 0.0
            else:
                print(f"SKIP DIHEDRAL ({resname} {axis_label}): vertex has no "
                      f"formula-derived azimuth (4 distinct substituent types)")
                return None
            geom[n] = (r, theta, delta)
        return geom

    b_geom = side_geom(b_neighbors, b_class, c_class, b_ref)
    c_geom = side_geom(c_neighbors, c_class, b_class, c_ref)
    if b_geom is None or c_geom is None:
        axes.bump_axis_skip()
        return

    target, dc_by_harmonic = axes.combined_target_for_axis(
        b_geom, c_geom, b_class, c_class, lambda n: atoms[n])
    if target is None:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}): AMBER gives this axis a "
              f"zero barrier (a real, deliberate null -- not a gap)")
        axes.bump_axis_skip()
        return

    for name in emit_as:
        axes.emit_ghost_rings_for_axis(name, axis_label, family, L_axis, target,
                                       dc_by_harmonic, b_name, c_name, b_ref, c_ref,
                                       b_geom[b_ref][:2], c_geom[c_ref][:2])


def main():
    dna = ForceFieldTables(os.path.join(DATA, "amber14", "DNA.OL15.xml"))
    rna = ForceFieldTables(os.path.join(DATA, "amber14", "RNA.OL3.xml"))

    header = [
        "#Bonded parameters for DNA and RNA, generated by",
        "#scripts/generate_nucleic_forcefield.py -- do not edit by hand.",
        "#",
        "#Sources: amber14/DNA.OL15.xml and amber14/RNA.OL3.xml, the",
        "#nucleic-acid-specific AMBER force fields (see the script's docstring",
        "#for why amber99sb.xml's own nucleic templates are the wrong choice).",
        "#",
        "#STRETCH <name> <resname> <atom1> <atom2> <r0_A> <k_kJ.mol-1.A-2>",
        "#BEND    <name> <resname> <atom1> <vertex> <atom3> <theta0_deg> <k>",
        "#",
        "#'+X' is atom X of the NEXT residue: the only cross-residue bond in a",
        "#nucleic acid is O3'(i)-P(i+1), the phosphodiester link.",
        "#",
        "#Every rule is emitted under each known spelling of its residue and",
        "#atoms (DA vs A, H2'1 vs H2', O1P vs OP1...). A rule that knows only",
        "#one spelling matches NOTHING on a file written in another, with no",
        "#error at all -- measured, 0 springs out of 963 particles.",
        "#",
        "#DIHEDRAL COVERAGE IS PARTIAL, and the gap is one single cause. Emitted:",
        "#alpha, beta, zeta (the phosphate torsions) and RNA's 2'-OH rotor.",
        "#Missing: gamma, delta, epsilon, chi and the four sugar-ring bonds --",
        "#every axis with a sugar carbon (C1'..C4') at one end. Each of those is",
        "#a tetrahedral vertex whose three non-axis substituents have three",
        "#DISTINCT parameter types, so its azimuths are a real chirality that no",
        "#angle table can resolve: it needs the ribose's own (beta-D) reference",
        "#frame, the exact counterpart of the L-amino-acid frame the protein",
        "#generator uses for C-alpha. Those axes are SKIPPED and reported, never",
        "#emitted from a guessed sign.",
        "#",
        "#Companion rigid-body files: DNAAtomRigidGroups.rbody and",
        "#RNAAtomRigidGroups.rbody. Those keep the sugar RIGID, so no",
        "#dihedral is emitted here yet for delta or the ring interior: with a",
        "#rigid furanose they would have no effect. Freeing the pucker is the",
        "#open question -- see the .rbody headers.",
        "",
    ]

    for base in ("DA", "DC", "DG", "DT"):
        generate_residue(dna, base, [base])
    # OL3 spells them A/C/G/U; amber99sb spells the same residues RA/RC/RG/RU.
    for base in ("A", "C", "G", "U"):
        generate_residue(rna, base, [base, "R" + base])

    # Dihedrals. The .bi.ff family column is what BioSpring gates at runtime,
    # and its five values (PHI/PSI/OMEGA/SIDECHAIN/PLANARITY) are named for
    # protein chemistry -- there is no alpha..zeta. Every nucleic proper
    # torsion therefore goes under SIDECHAIN, i.e. one gate
    # (dihedralchi.enable) for all of them, rather than borrowing
    # protein-backbone labels that would misdescribe them in the file.
    # Splitting them into their own families is the same C++ chain PHI/PSI/
    # OMEGA already went through (enum, parser, Topology collection, NetCDF,
    # .msp key, SpringNetwork gate) and is deliberately left as its own change.
    for ff, bases, spellings in ((dna, ("DA", "DC", "DG", "DT"), lambda b: [b]),
                                 (rna, ("A", "C", "G", "U"), lambda b: [b, "R" + b])):
        ff.use()
        for base in bases:
            emit_as = spellings(base)
            for b_name, c_name, label in BACKBONE_AXES:
                generate_axis(ff, base, b_name, c_name, label, "SIDECHAIN", emit_as)
            for b_name, c_name, label in SUGAR_RING_AXES:
                generate_axis(ff, base, b_name, c_name, label, "SIDECHAIN", emit_as)
            generate_axis(ff, base, "C1'", GLYCOSIDIC_N[base], "chi", "SIDECHAIN", emit_as)
            if "O2'" in ff.atoms_of(base):
                generate_axis(ff, base, "C2'", "O2'", "rotor_O2", "SIDECHAIN", emit_as)

    out = os.path.join(REPO_ROOT, "data/reducerules/NucleicAtomBonded.bi.ff")
    body = lines + [
        "#",
        "#GHOSTPARTICLE <name> <resname> <atom_B> <atom_C> <atom_ref> <r_A> <theta_deg> <delta_deg>",
        "#Massless virtual sites, placed algebraically from 3 real anchors.",
    ] + axes.ghostparticle_lines + [
        "#",
        "#DIHEDRAL <name> <resname> <family> <atom_ref> <atom_rotant> <d0_A> <k> <dc_offset>",
        "#One ring group per real AMBER Fourier harmonic of each covered axis.",
    ] + axes.dihedral_lines
    with open(out, "w") as f:
        f.write("\n".join(header + body) + "\n")
    print(f"\nWrote {counters['stretch']} STRETCH, {counters['bend']} BEND "
          f"({counters['skip']} skipped), {axes.n_ghost_particles} GHOSTPARTICLE and "
          f"{axes.n_dihedral_ok} DIHEDRAL entries ({axes.n_dihedral_skip} axes skipped) to {out}")


if __name__ == "__main__":
    main()
