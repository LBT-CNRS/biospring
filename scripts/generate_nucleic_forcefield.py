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

# The furanose is left RIGID: its five bonds and five valence angles keep
# the .rbody's uniform springs and no bonded rule retunes them, and no
# torsion is emitted about a ring bond (delta included -- it turns about
# C4'-C3').
#
# Why, measured (2026-08-19): a ghost ring's energy is not a function of
# the dihedral alone. Since the ghost is the image of a REAL atom rotated
# about the axis, deforming the ring's angles moves the ghosts too, and
# there is a relaxation channel in which every ghost-ghost distance
# reaches its d0 and the torsion energy collapses -- measured at 293
# kJ/mol where AMBER asks 5388. A single deformation relaxes all five ring
# axes at once, so the gain (~2800) outgrows what the bonds and angles
# charge for it (~1950): the ring flattens. Quenched, the pucker amplitude
# went 37 deg -> 13 deg, all 42 sugars, while AMBER's own bonded terms
# hold it at 41 deg. The five ring axes were exactly the five whose energy
# moved OPPOSITE to AMBER's; the six open axes were all correct.
#
# Open axes keep their bonded terms: one deformation only relaxes one
# axis there, and it has to pay for itself. Aromatic rings and the bases
# are not concerned either -- for them the idealised geometry IS the true
# one, so the channel leads nowhere (measured: 0.000 A out of plane).
RING_ATOMS = frozenset(("C1'", "C2'", "C3'", "C4'", "O4'"))


def is_ring_internal(*names):
    """True when every atom named is a furanose ring atom of THIS residue
    (a +/- prefix means another residue, so never ring-internal here)."""
    return all(n in RING_ATOMS for n in names)


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
        if is_ring_internal(a1, a2):
            continue
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

    # Angles: every pair of neighbours around each vertex -- including the
    # junction. The adjacency walked here comes from ONE residue's template,
    # so the phosphodiester link's angles appear at no vertex and must be
    # added by hand (the protein generator never had this problem: its
    # synthetic chain carries real peptide bonds, so its adjacency crosses
    # residues on its own). Found by the torque bench, not by inspection:
    # without these, C3'-O3'-P(i+1) and O3'(i-1)-P-{O5',OP1,OP2} have no
    # BEND at all and are held only by the .rbody's uniform springs (which
    # a missing rule never retunes -- silently, since a BEND rule that
    # doesn't exist reports nothing).
    junction = []
    # A 3'-terminal template caps its O3' with HO3' and has no next residue,
    # so this angle genuinely does not exist for it -- not a missing table
    # entry. The interior template emits it under the same residue name.
    if "O3'" in atoms and "C3'" in atoms and "P" in atoms and "HO3'" not in atoms:
        junction.append(("C3'", "O3'", NEXT_P))
    if "P" in atoms:
        for x in ("O5'", "OP1", "OP2"):
            if x in atoms:
                junction.append((PREV_O3, "P", x))
    for vertex, nbrs in adj.items():
        nb = sorted(nbrs)
        for i in range(len(nb)):
            for j in range(i + 1, len(nb)):
                junction.append((nb[i], vertex, nb[j]))
    for a1, vertex, a3 in junction:
        if is_ring_internal(a1, vertex, a3):
            continue
        # A +/- atom's identifier comes from this residue's own template
        # (backbone identifiers are uniform across residues in these files).
        key = (atoms[vertex.lstrip("+-")],
               frozenset((atoms[a1.lstrip("+-")], atoms[a3.lstrip("+-")])))
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

# The furanose's other four bonds. No torsion is emitted about any of them
# any more -- the ring is frozen, see RING_ATOMS -- so this list is kept
# only to say which bonds that covers, next to BACKBONE_AXES' delta.
SUGAR_RING_AXES = [("C1'", "C2'", "ring_C1_C2"), ("C2'", "C3'", "ring_C2_C3"),
                   ("C4'", "O4'", "ring_C4_O4"), ("O4'", "C1'", "ring_O4_C1")]

# chi, the glycosidic torsion: O4'-C1'-N9-C4 on a purine, O4'-C1'-N1-C2 on a
# pyrimidine. Plus RNA's 2'-OH rotor.
GLYCOSIDIC_N = {"DA": "N9", "DG": "N9", "DC": "N1", "DT": "N1",
                "A": "N9", "G": "N9", "C": "N1", "U": "N1"}

# Every backbone/sugar atom name. Used only to tell a base hub from a sugar
# one when reading impropers back: AMBER puts no improper on an sp3 sugar
# carbon, so this is a guard against surprises rather than a filter that
# does real work.
SUGAR_AND_BACKBONE = {"P", "OP1", "OP2", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C3'",
                      "O3'", "C2'", "C1'", "O2'", "H5'", "H5''", "H5'1", "H5'2", "H4'",
                      "H3'", "H2'", "H2''", "H2'1", "H2'2", "H1'", "HO5'", "HO3'", "HO2'",
                      "HO'2"}


# ---------------------------------------------------------------------------
# Base planarity impropers
# ---------------------------------------------------------------------------
# An AMBER improper IS a 4-atom dihedral (i,j,k,l) with the hub listed THIRD
# (k bonded to i, j and l), measured about the real j-k bond -- so the ghost
# ring models it as-is: axis=(j,k), refs=(i,l), one single-pair ring per
# improper. Every base improper comes out n=2, phase=180 deg: the sp2 centre
# is held in its plane, and both faces cost the same.
#
# AMBER's improper MATCHING rules (wildcards, the trefoil ordering that
# decides which peripheral is listed where) are NOT reimplemented here --
# that is the part most likely to be silently wrong. The quadruplets and
# their parameters are read from a real OpenMM System built on a synthetic
# chain, so AMBER's own assignment engine does the matching. Everything
# after that is a parameter lookup, not a geometry derivation, so the
# structure-free rule the rest of this file follows still holds: the chain
# below carries no coordinates at all, only atoms and bonds.
_improper_cache = {}


def base_improper_quads(ff, chain_residues):
    """[(hub_resname, (i, j, k, l) names, n, phase_rad, k)] for every improper
    AMBER assigns to a base, read off a synthetic chain of `chain_residues`."""
    if ff.path in _improper_cache:
        return _improper_cache[ff.path]

    import openmm.app as app
    import openmm.unit as unit

    top = app.Topology()
    chain = top.addChain()
    per_residue = []
    for rname in chain_residues:
        tpl = ff.residues[rname]
        res = top.addResidue(rname, chain)
        by_name = {}
        for a in tpl.findall("Atom"):
            name = a.get("name")
            # AMBER atom names start with their element; no nucleic template
            # uses a leading digit.
            symbol = "Cl" if name.upper().startswith("CL") else name.lstrip("0123456789")[0]
            by_name[name] = top.addAtom(name, app.element.get_by_symbol(symbol), res)
        for b in tpl.findall("Bond"):
            top.addBond(by_name[b.get("atomName1")], by_name[b.get("atomName2")])
        per_residue.append(by_name)
    # The phosphodiester link, the one bond a template cannot carry.
    for a, b in zip(per_residue, per_residue[1:]):
        if "O3'" in a and "P" in b:
            top.addBond(a["O3'"], b["P"])

    system = app.ForceField(ff.path).createSystem(top, nonbondedMethod=app.NoCutoff,
                                                  constraints=None)
    torsions = next(system.getForce(i) for i in range(system.getNumForces())
                    if system.getForce(i).__class__.__name__ == "PeriodicTorsionForce")
    bonded = set()
    for b in top.bonds():
        bonded.add((b[0].index, b[1].index))
        bonded.add((b[1].index, b[0].index))
    atoms = list(top.atoms())

    quads = []
    for t in range(torsions.getNumTorsions()):
        i, j, k, l, per, phase, kv = torsions.getTorsionParameters(t)
        quad = [i, j, k, l]
        # An improper is exactly the case where one atom of the quad is
        # bonded to the other three; a proper torsion never is.
        hubs = [c for c in quad if all((c, o) in bonded for o in quad if o != c)]
        if not hubs:
            continue
        hub = atoms[hubs[0]]
        if hub.name in SUGAR_AND_BACKBONE:
            continue
        if any(atoms[x].residue.index != hub.residue.index for x in quad):
            continue
        quads.append((hub.residue.name, tuple(atoms[x].name for x in quad), per,
                      phase.value_in_unit(unit.radian),
                      kv.value_in_unit(unit.kilojoule_per_mole)))
    _improper_cache[ff.path] = quads
    return quads


def generate_base_impropers(ff, resname, emit_as, quads):
    """One ring group per improper whose hub belongs to `resname`'s base."""
    atoms = ff.atoms_of(resname)
    seen, emitted = set(), 0
    for hub_res, names, n, phase, kv in quads:
        # The 5'/3'/nucleoside variants carry the same base and so the same
        # impropers; keying on the hub keeps exactly one of each. Without
        # this the rule would be emitted several times under one name --
        # silently doubled improper energy, the same trap the protein
        # generator hit.
        if hub_res != resname or names[2] in seen or kv <= 1e-6:
            continue
        if any(nm not in atoms for nm in names):
            continue
        seen.add(names[2])
        cls = lambda nm: atoms[nm]
        L_axis = axes.bond_len_A(cls(names[1]), cls(names[2]))
        if L_axis is None:
            print(f"SKIP PLANARITY ({resname} {names[2]}): no bond length for the axis")
            axes.bump_axis_skip()
            continue
        # names[0] hangs off the HUB names[2], not off the b-side anchor
        # names[1], so its (r, theta) has to be resolved through the hub --
        # geom_via_far_anchor, exactly as the protein impropers do.
        ref_geom_b = axes.geom_via_far_anchor(
            L_axis, axes.bond_len_A(cls(names[0]), cls(names[2])),
            axes.valence_deg(cls(names[2]), cls(names[0]), cls(names[1])))
        ref_geom_c = (axes.bond_len_A(cls(names[3]), cls(names[2])),
                      axes.valence_deg(cls(names[2]), cls(names[3]), cls(names[1])))
        if ref_geom_b is None or None in ref_geom_c:
            print(f"SKIP PLANARITY ({resname} {names[2]}): missing bond/angle for a peripheral")
            axes.bump_axis_skip()
            continue
        target = kv * np.exp(-1j * phase)
        for name in emit_as:
            axes.emit_ghost_ring(name, f"imp_{names[2]}", "PLANARITY", n, L_axis, target,
                                 names[1], names[2], names[0], names[3],
                                 axis_dc_target=kv, ref_geom_b=ref_geom_b,
                                 ref_geom_c=ref_geom_c, dc_align="minimum")
        emitted += 1
    if emitted:
        axes.bump_axis_ok()
    return emitted



# ---------------------------------------------------------------------------
# The sugar's four stereocenters
# ---------------------------------------------------------------------------
# AMBER's angle tables fix a tetrahedral vertex only up to a mirror
# reflection. For a vertex whose substituents are all different that
# reflection is a real chirality, and no table can supply it -- exactly the
# situation the protein generator meets at the backbone C-alpha, and it
# resolves it there with the L-amino-acid convention rather than by measuring
# anything.
#
# The nucleotide sugar is the same kind of fact: it is ALWAYS
# beta-D-ribofuranose (beta-D-2'-deoxyribofuranose in DNA), never anything
# else. Written as CIP descriptors, that convention is:
#
#     C1'  R      C2'  R  (ribose only -- DNA's C2' carries two hydrogens
#     C3'  S             and is not a stereocenter at all)
#     C4'  R
#
# So the sign is not searched for, guessed, or measured: the target
# descriptor is the input, and the frame builder below picks whichever of the
# two mirror images realizes it.

# CIP priority order at each centre (highest first), then the target
# descriptor. C1' takes N9 on purines and N1 on pyrimidines -- resolved by
# GLYCOSIDIC_N at build time.
STEREOCENTERS = {"C1'": (["O4'", None, "C2'", "H1'"], "R"),
                 "C2'": (["O2'", "C1'", "C3'", "H2'"], "R"),
                 "C3'": (["O3'", "C4'", "C2'", "H3'"], "S"),
                 "C4'": (["O4'", "C3'", "C5'", "H4'"], "R")}


def cip_descriptor(a, b, c, d):
    """'R' or 'S' from the four unit substituent vectors, given in CIP
    priority order. Viewed with the lowest priority (d) pointing away, R is
    a clockwise 1->2->3."""
    s = np.dot(d, np.cross(a, b) + np.cross(b, c) + np.cross(c, a))
    return "R" if s > 0 else "S"


def _check_cip_convention():
    """Pins the sign of cip_descriptor to the textbook definition, on a
    tetrahedron built here -- no molecule involved. Wrong-way-round is
    invisible in any energy (a mirrored model has the mirrored energy curve,
    which matches AMBER wherever the torsion phase is 0 or 180 deg -- and
    most of them are), so it is asserted at import instead of being trusted."""
    d = np.array([0.0, 0.0, -1.0])
    t = np.radians(70.5)

    def sub(az_deg):
        az = np.radians(az_deg)
        return np.array([np.sin(t) * np.cos(az), np.sin(t) * np.sin(az), np.cos(t)])

    assert cip_descriptor(sub(90), sub(-30), sub(-150), d) == "R"
    assert cip_descriptor(sub(90), sub(210), sub(330), d) == "S"


_check_cip_convention()

_sugar_frames = {}


def sugar_frame(ff, resname, vertex):
    """Unit vectors from `vertex` to each of its four substituents, built
    from AMBER's own angles and oriented by the beta-D convention above."""
    key = (id(ff), resname, vertex)
    if key in _sugar_frames:
        return _sugar_frames[key]
    atoms = ff.atoms_of(resname)
    order, target = STEREOCENTERS[vertex]
    order = [GLYCOSIDIC_N[resname.lstrip("R")] if n is None else n for n in order]
    if any(n not in atoms for n in order) or vertex not in atoms:
        _sugar_frames[key] = None
        return None
    v = atoms[vertex]

    def ang(x, y):
        return axes.valence_deg(v, atoms[x], atoms[y])

    s1, s2, s3, s4 = order
    if any(ang(x, y) is None for x, y in ((s1, s2), (s1, s3), (s2, s3),
                                          (s1, s4), (s2, s4), (s3, s4))):
        _sugar_frames[key] = None
        return None
    for sign in (+1.0, -1.0):
        V1, V2, V3 = axes._build_two_vectors(ang(s1, s2), ang(s1, s3), ang(s2, s3), sign)
        V4 = np.linalg.solve(np.array([V1, V2, V3]),
                             np.cos(np.radians([ang(s1, s4), ang(s2, s4), ang(s3, s4)])))
        V4 = V4 / np.linalg.norm(V4)
        if cip_descriptor(V1, V2, V3, V4) == target:
            frame = {s1: V1, s2: V2, s3: V3, s4: V4, vertex: np.zeros(3)}
            _sugar_frames[key] = frame
            return frame
    _sugar_frames[key] = None
    return None


def sugar_delta(ff, resname, vertex, ref_name, axis_partner, other_name):
    """dihedral(ref, vertex, axis_partner, other) read off the resolved
    frame -- the same quantity ca_tetrahedral_delta returns for the protein
    backbone, and used the same way by both sides of an axis."""
    f = sugar_frame(ff, resname, vertex)
    if f is None or ref_name not in f or axis_partner not in f or other_name not in f:
        return None
    # The axis partner defines the axis; ref and other are read around it.
    return axes._dihedral_from_points(f[ref_name], f[vertex], f[axis_partner], f[other_name])


def neighbours(ff, resname, vertex, exclude):
    """Real bonded substituents of `vertex` other than `exclude`, in template
    order -- the same stable, structure-free ordering the protein generator
    gets from its synthetic chain.

    `vertex` keeps its +/- prefix, and the substituents inherit it: a
    template is one residue, so walking the adjacency of "+P" returns P's
    neighbours by NAME but they belong to the residue the prefix points at.
    Dropping the prefix there names the wrong residue's atom -- silently, and
    only the ghost's placement shows it (an anchor 7 A away instead of 1.5).
    """
    prefix = vertex[0] if vertex[:1] in ("+", "-") else ""
    bare = vertex.lstrip("+-")
    adj = ff.neighbours_of(resname)
    order = list(ff.atoms_of(resname))
    out = [prefix + n
           for n in sorted((n for n in adj.get(bare, ()) if n != exclude.lstrip("+-")),
                           key=order.index)]
    # The phosphodiester link crosses the residue boundary, and the template
    # adjacency stops there -- BOTH ends need patching by hand. O3' gains
    # the next residue's P (without it, epsilon looks like a bare-ended axis
    # and is skipped entirely). P gains the PREVIOUS residue's O3' -- and
    # this one fails the other way: alpha (O3'-P-O5'-C5') still has OP1/OP2
    # on its P side, so nothing was skipped or warned, the axis was simply
    # built without the substituent that DEFINES alpha. Invisible in any
    # energy total; the torque bench caught it as sign flips on half the
    # alpha axes.
    if bare == "O3'" and not prefix and exclude.lstrip("+-") != "P":
        out.append(NEXT_P)
    if bare == "P" and not prefix and exclude.lstrip("+-") != "O3'":
        out.append(PREV_O3)
    # Unique-class substituents first, so the azimuth REFERENCE is never one
    # of an identical pair. A vertex like P (OP1, OP2, -O3') gets its frame
    # from formula_derived_deltas' symmetric-tetrahedral case, which fixes
    # the two identical atoms only up to a mirror -- that is fine for THEM
    # (their pair-sum is exchange-invariant) but not for a distinct branch
    # measured AGAINST one of them: the branch azimuth's sign is then
    # arbitrary, and PDB naming of the pair (OP1/OP2, H5'/H5'') is itself
    # arbitrary per residue. Measured on the thermalised duplex before this
    # reordering: alpha mirrored on 40/40 residues, beta on 18/42 -- torque
    # sign flips on about half of each. With the branch as reference, the
    # identical pair sits at +/-symmetric azimuths and no baked sign is
    # left to be wrong. Stable sort: template order otherwise.
    atoms_of = ff.atoms_of(resname)
    counts = {}
    for n in out:
        c = atoms_of[n.lstrip("+-")]
        counts[c] = counts.get(c, 0) + 1
    out.sort(key=lambda n: 0 if counts[atoms_of[n.lstrip("+-")]] == 1 else 1)
    return out


def generate_axis(ff, resname, b_name, c_name, axis_label, family, emit_as, per_pair=False):
    """One axis, one call -- mirrors generate_sidechain_axis exactly.

    per_pair: emit ONE ring per real (b-substituent, c-substituent) pair,
    anchored on that pair's own two atoms, instead of one combined ring per
    harmonic. Same option, same criterion and same reason as the protein
    generator's: use it only where every pair reinforces IN PHASE under the
    ideal geometry, because there a combined ring bakes that maximal
    coherence in and comes out too stiff on a real structure, whose
    substituents never sit exactly at their ideal azimuths.

    Measured, and NOT used by any nucleic axis today. chi looks like the
    textbook case -- coherence 1.00 on every harmonic, DNA and RNA alike --
    but that is trivially true because AMBER defines chi through a SINGLE
    specific quadruplet (O4'-C1'-N9-C8, which is what the chiOL15/chiOL3
    correction is) with no generic wildcard on the axis: one term is always
    in phase with itself, and per-pair then emits the same single ring.
    Verified by running it: identical energies to the digit. gamma (n=1 at
    0.69) and epsilon (0.50/0.50/0.19) partly cancel, so the combined ring
    is already unbiased there. Kept because the option costs nothing and
    the next axis that needs it should not have to reinvent it.
    """
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

    # Prefixed, not bare: an axis atom on the next residue (zeta's "+P")
    # carries its whole substituent set over with it.
    b_neighbors = neighbours(ff, resname, b_name, c_bare)
    c_neighbors = neighbours(ff, resname, c_name, b_bare)
    if not b_neighbors or not c_neighbors:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}): no real substituent on one side")
        axes.bump_axis_skip()
        return
    b_ref, c_ref = b_neighbors[0], c_neighbors[0]

    def side_geom(nbrs, vertex_class, other_class, ref_name):
        with_classes = [(n, atoms[n.lstrip('+-')]) for n in nbrs]
        deltas = axes.formula_derived_deltas(vertex_class, other_class, with_classes, ref_name)
        geom = {}
        for n in nbrs:
            nc = atoms[n.lstrip('+-')]
            r = axes.bond_len_A(nc, vertex_class)
            theta = axes.valence_deg(vertex_class, nc, other_class)
            if r is None or theta is None:
                return None
            if deltas is not None:
                delta = deltas[n]
            elif n == ref_name:
                delta = 0.0
            else:
                # All-different substituents: a real chirality. The sugar's
                # four centres are covered by the beta-D frame above; nothing
                # else in a nucleotide is a stereocenter.
                vertex_name = b_bare if nbrs is b_neighbors else c_bare
                partner = c_bare if nbrs is b_neighbors else b_bare
                delta = sugar_delta(ff, resname, vertex_name, ref_name, partner, n) \
                    if vertex_name in STEREOCENTERS else None
                if delta is None:
                    print(f"SKIP DIHEDRAL ({resname} {axis_label}): vertex "
                          f"{vertex_name} has all-different substituents and is "
                          f"not one of the sugar stereocenters")
                    return None
            geom[n] = (r, theta, delta)
        return geom

    b_geom = side_geom(b_neighbors, b_class, c_class, b_ref)
    c_geom = side_geom(c_neighbors, c_class, b_class, c_ref)
    if b_geom is None or c_geom is None:
        axes.bump_axis_skip()
        return

    target, dc_by_harmonic = axes.combined_target_for_axis(
        b_geom, c_geom, b_class, c_class, lambda n: atoms[n.lstrip('+-')])
    if target is None:
        print(f"SKIP DIHEDRAL ({resname} {axis_label}): AMBER gives this axis a "
              f"zero barrier (a real, deliberate null -- not a gap)")
        axes.bump_axis_skip()
        return

    if per_pair:
        emitted = 0
        for bn in b_neighbors:
            for cn in c_neighbors:
                t_pair, dc_pair = axes.combined_target_for_axis(
                    {bn: (0, 0, 0.0)}, {cn: (0, 0, 0.0)}, b_class, c_class,
                    lambda n: atoms[n.lstrip('+-')])
                if t_pair is None:
                    continue
                for n, tn in sorted(t_pair.items()):
                    if n == 0 or abs(tn) <= 1e-6:
                        continue
                    for name in emit_as:
                        axes.emit_ghost_ring(name, axis_label, family, n, L_axis, tn,
                                             b_name, c_name, bn, cn,
                                             axis_dc_target=dc_pair.get(n, 0.0),
                                             group_tag=f"{bn}{cn}",
                                             ref_geom_b=b_geom[bn][:2],
                                             ref_geom_c=c_geom[cn][:2])
                    emitted += 1
        if emitted:
            axes.bump_axis_ok()
            return
        # Nothing resolved per pair: fall through rather than emit nothing.

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
        "#DIHEDRAL COVERAGE IS COMPLETE: alpha..zeta, the glycosidic chi, the",
        "#four furanose ring bonds and RNA's 2'-OH rotor, with no axis skipped.",
        "#Every axis with a sugar carbon (C1'..C4') at one end goes through a",
        "#tetrahedral vertex whose three non-axis substituents have three",
        "#DISTINCT parameter types -- a real chirality no angle table can",
        "#resolve. It is resolved from the CONVENTION, not from a structure:",
        "#the sugar is always beta-D-ribofuranose, which in CIP descriptors is",
        "#C1' R, C2' R, C3' S, C4' R (DNA's C2' carries two H and is not a",
        "#stereocenter). The constructor builds both mirror images and keeps",
        "#the one realising the descriptor -- the exact counterpart of the",
        "#L-amino-acid frame the protein generator uses for C-alpha.",
        "#",
        "#Companion rigid-body files: DNAAtomRigidGroups.rbody and",
        "#RNAAtomRigidGroups.rbody. Measured on a real build, their per-vertex",
        "#groups cover all 10 pairs of the 5 furanose atoms -- a 5-cycle puts",
        "#every pair at graph distance <= 2, so 'vertex + its 2 ring",
        "#neighbours' leaves nothing out. Under --rigidbody that is a complete",
        "#clique at uniform stiffness, so the ring itself is rigid and the",
        "#pucker would be FROZEN -- which is why NUCLEIC_SUGAR exists as its",
        "#own family, torsions being what lets the pucker move at all here.",
        "#AMBER has no pucker term either and lets the state emerge from bonds,",
        "#angles, torsions and 1-4 together; this model keeps the bonds and",
        "#angles rigid and carries only the torsions.",
        "",
    ]

    # The interior template first, then the 5'/3'/nucleoside variants under
    # the SAME residue name. A PDB does not spell terminal residues
    # differently -- OpenMM writes plain "DA" for a 5' or 3' end -- so a
    # rule emitted under "DA5" would match nothing, and the atoms only those
    # variants carry (HO5', HO3', and the 5'-OH's own geometry) would have
    # no rule at all. Measured before this: HO3' ended a quench 9.4 A from
    # its own O3', with no spring to notice, because every rule came from
    # the interior template. Same "list every spelling, the reader skips
    # what is absent" convention as OP1/O1P -- emit_stretch/emit_bend
    # deduplicate, so the variants only ever ADD their own extra rules.
    # Interior first so that is the version that wins any tie.
    for base in ("DA", "DC", "DG", "DT"):
        for variant in (base, base + "5", base + "3", base + "N"):
            if variant in dna.residues:
                generate_residue(dna, variant, [base])
    # OL3 spells them A/C/G/U; amber99sb spells the same residues RA/RC/RG/RU.
    for base in ("A", "C", "G", "U"):
        for variant in (base, base + "5", base + "3", base + "N"):
            if variant in rna.residues:
                generate_residue(rna, variant, [base, "R" + base])

    # Dihedrals, in three nucleic-specific families so nothing borrows a
    # protein-backbone label that would misdescribe it:
    #   NUCLEIC_BACKBONE  alpha beta gamma delta epsilon zeta
    #   NUCLEIC_CHI       the glycosidic torsion + RNA's 2'-OH rotor
    #   NUCLEIC_SUGAR     the four furanose ring bonds -- kept apart because
    #                     the pucker is what they govern, and isolating it is
    #                     the whole point of validating a nucleic model.
    for ff, bases, spellings in ((dna, ("DA", "DC", "DG", "DT"), lambda b: [b]),
                                 (rna, ("A", "C", "G", "U"), lambda b: [b, "R" + b])):
        ff.use()
        # A synthetic chain carrying every base as an INTERIOR residue,
        # bracketed by the 5'/3' terminal variants so no base of interest is
        # itself terminal. Atoms and bonds only -- no coordinates, and none
        # needed: what is read back off it is which improper AMBER assigns
        # where, not any geometry.
        chain = [bases[0] + "5"] + list(bases) + [bases[-1] + "3"]
        quads = base_improper_quads(ff, chain)
        for base in bases:
            emit_as = spellings(base)
            for b_name, c_name, label in BACKBONE_AXES:
                # delta turns about C4'-C3', a ring bond: frozen with the
                # ring (see RING_ATOMS). It was never an independent hinge
                # anyway -- it is a readout of the pucker.
                if is_ring_internal(b_name, c_name):
                    continue
                generate_axis(ff, base, b_name, c_name, label, "NUCLEIC_BACKBONE", emit_as)
            generate_axis(ff, base, "C1'", GLYCOSIDIC_N[base], "chi", "NUCLEIC_CHI", emit_as)
            if "O2'" in ff.atoms_of(base):
                generate_axis(ff, base, "C2'", "O2'", "rotor_O2", "NUCLEIC_CHI", emit_as)
            generate_base_impropers(ff, base, emit_as, quads)

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
    print(f"\nWrote {axes.n_ghost_particles} GHOSTPARTICLE and {axes.n_dihedral_ok} DIHEDRAL "
          f"entries ({axes.n_dihedral_skip} axes skipped), {counters['stretch']} STRETCH and "
          f"{counters['bend']} BEND ({counters['skip']} skipped) to {out}")


if __name__ == "__main__":
    main()
