#!/usr/bin/env python3
"""Generate BioSpring's Martini 3 reduce rules (.grp) and non-bonded force
field (.nbi.ff / .ff) from the published Martini 3.0.0 parameters.

PROVENANCE -- every number below was read out of one of these two sources,
not chosen here:

  * martini_v3.0.0.itp (Marrink lab, martini-forcefields, file dated
    2021-03-29): 844 bead types and 355746 [ nonbond_params ] pairs.
    Columns 4 and 5 of that section give sigma_ij and epsilon_ij.
    TYPE below holds the DIAGONAL (i == j) entries of the types a protein
    uses, plus the mass from [ atomtypes ].

  * vermouth / martinize2 force field "martini3001": aminoacids.ff gives
    each residue's bead types and charges, modifications.ff gives the
    terminal ones, and data/mappings/martini3001/*.amber.map gives which
    atoms make up each bead.

WHY THE MAPPING IS HEAVY-ATOM ONLY. A Martini bead is the MASS-weighted
average of the atoms mapped onto it: martinize2 runs vermouth's
DoAverageBead, which reads the force field variable center_weight -- "mass"
in martini3001's general.ff -- and computes
np.average(positions, weights=mapping_weight * atom_mass). A BioSpring grain
is measure::centroid, an UNWEIGHTED mean, and no .grp can weight anything.

So the two agree only when the mapped atoms have similar masses, which heavy
atoms do (12, 14, 16) and hydrogens emphatically do not: in a mass-weighted
average a hydrogen pulls with 1/13 of a carbon's weight, while an unweighted
centroid gives it the same weight. Measured against the true mass-weighted
bead over the 401 beads of example/001_GKinase/model.pdb:

  listing Martini's weight-1 atoms, hydrogens included   median 0.342 A,
                                                         p90 0.537, max 0.730
  listing the heavy ones only                            median 0.028 A,
                                                         p90 0.134, max 0.415

Dropping the hydrogens is therefore not an approximation of Martini but a
TWELVEFOLD improvement on including them, and it is what this file does. The
residual 0.415 A worst case is on the sulfur-bearing beads -- CYS SC1 is CB
at 12 Da and SG at 32, so the mass-weighted centre sits 0.41 A from the
midpoint that BioSpring computes. Nothing in a .grp can fix that: an atom
name cannot be repeated to fake a weight, since ReduceRule stores its atoms
in an unordered_set.

WHAT THIS FILE DOES NOT CARRY, and why -- read data/martini3/README.md for
the measured size of each gap. The short version: Martini 3 stores sigma
and epsilon PER PAIR, and a BioSpring .ff stores them PER TYPE, so the
pair matrix cannot be written here at all. See README.md.

TWO NORMALISATIONS were applied to the upstream mapping files, which
contradict each other and themselves on these points:

  1. Beta hydrogens carry weight 0. Eighteen of the twenty .amber.map
     files write "HB1 !SC1" (the "!" means null weight); phe.amber.map
     writes "HB1 SC1" and trp.amber.map writes "HB2 !SC2" for an atom
     bonded to an SC1 carbon. The majority spelling is used everywhere.
  2. A ring or polar hydrogen belongs to the bead holding the atom it is
     bonded to. his/hid/hie/hip.amber.map disagree about this: his omits
     HD1 and HE1 altogether, hid puts HD1 in SC1 while CG is in SC2, and
     hie lists HD1 on a tautomer that has no delta proton. The rule here
     is chemical and is applied uniformly to all four.

Both normalisations are visible in the generated files and testable: see
README.md for what honouring the hydrogens is worth in Angstrom.
"""

import collections
import functools
import math
import os

# 2^(1/6): the factor between a Lennard-Jones sigma and its rmin.
TWO_ONE_SIXTH = 2.0 ** (1.0 / 6.0)

# Bead type -> (mass in Da, sigma_ii in A, epsilon_ii in kJ.mol-1).
# Diagonal of martini_v3.0.0.itp's [ nonbond_params ]; mass from [ atomtypes ].
TYPE = {
    "P6":    (72.0, 4.700, 4.990),
    "P5":    (72.0, 4.700, 4.480),
    "P2":    (72.0, 4.700, 4.060),
    "P1":    (72.0, 4.700, 3.880),
    "SP5":   (54.0, 4.100, 3.820),
    "SP2":   (54.0, 4.100, 3.310),
    "SP2a":  (54.0, 4.100, 3.070),
    "SP1":   (54.0, 4.100, 3.070),
    "TP1":   (36.0, 3.400, 2.280),
    "TP1dq": (36.0, 3.400, 2.280),
    "Q5":    (72.0, 4.700, 6.450),
    "Q5n":   (72.0, 4.700, 5.700),
    "SQ5n":  (54.0, 4.100, 4.080),
    "SQ4p":  (54.0, 4.100, 3.690),
    "SQ3p":  (54.0, 4.100, 3.310),
    "C6":    (72.0, 4.700, 3.520),
    "SC4":   (54.0, 4.100, 2.350),
    "SC3":   (54.0, 4.100, 2.350),
    "SC2":   (54.0, 4.100, 2.350),
    "SC1":   (54.0, 4.100, 2.350),
    "TC6":   (36.0, 3.400, 1.770),
    "TC5":   (36.0, 3.400, 1.510),
    "TC4":   (36.0, 3.400, 1.510),
    "TC3":   (36.0, 3.400, 1.510),
    "TN6":   (36.0, 3.400, 2.280),
    "TN6d":  (36.0, 3.400, 2.020),
    "TN5a":  (36.0, 3.400, 1.770),
    "SN6d":  (54.0, 4.100, 2.840),
}

# Backbone bead type, from aminoacids.ff. Everything else uses the macro
# prot_default_bb_type, which that file sets to P2.
BB_TYPE = {"GLY": "SP1", "ALA": "SP2", "VAL": "SP2", "PRO": "SP2a", "HYP": "P1"}
BB_TYPE_DEFAULT = "P2"

# Terminal backbone, from modifications.ff: N-ter and C-ter both replace the
# backbone type with Q5 and set the charge to +1 and -1.
BB_TERMINAL_TYPE = "Q5"

# One-letter prefix for the grain names, so a grain is keyed by residue as
# well as by bead: GLY's BB and ALA's BB are different types (SP1 vs SP2).
CODE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASH": "DH",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLH": "EH", "GLY": "G",
    "HIS": "H", "HIE": "HE", "HIP": "HP", "ILE": "I", "LEU": "L",
    "LYS": "K", "LYN": "KN", "MET": "M", "PHE": "F", "PRO": "P",
    "HYP": "O", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
    "VAL": "V",
}

# Hydrogen spellings that differ between conventions. "amber" is what tleap
# and data/reducerules/amber.grp use; "pdbv2" is the older PDB spelling that
# puts the branch digit first, which example/001_GKinase/model.pdb carries.
# A bead whose hydrogens are spelled the same either way gets ONE rule.
BOTH = "both"

# resname -> list of side chain beads, each
#   (label, type, charge, heavy atoms, {convention: [hydrogens]})
# Heavy atoms and weight-1 hydrogens from mappings/martini3001/*.amber.map,
# with the two normalisations the module docstring describes.
SIDECHAIN = {
    "ALA": [("SC1", "TC3", 0.0, ["CB"], {})],
    "GLY": [],
    "CYS": [("SC1", "TC6", 0.0, ["CB", "SG"], {})],
    "VAL": [("SC1", "SC3", 0.0, ["CB", "CG1", "CG2"], {})],
    "LEU": [("SC1", "SC2", 0.0, ["CB", "CG", "CD1", "CD2"], {})],
    "ILE": [("SC1", "SC2", 0.0, ["CB", "CG1", "CG2", "CD1"], {})],
    "MET": [("SC1", "C6", 0.0, ["CB", "CG", "SD", "CE"], {})],
    "PRO": [("SC1", "SC3", 0.0, ["CB", "CG", "CD"], {})],
    "HYP": [("SC1", "P1", 0.0, ["CB", "CG", "CD", "OD1"], {})],
    "ASN": [("SC1", "SP5", 0.0, ["CB", "CG", "OD1", "ND2"],
             {"amber": ["HD21", "HD22"], "pdbv2": ["1HD2", "2HD2"]})],
    "GLN": [("SC1", "P5", 0.0, ["CB", "CG", "CD", "OE1", "NE2"],
             {"amber": ["HE21", "HE22"], "pdbv2": ["1HE2", "2HE2"]})],
    "ASP": [("SC1", "SQ5n", -1.0, ["CB", "CG", "OD1", "OD2"], {})],
    "ASH": [("SC1", "SP2", 0.0, ["CB", "CG", "OD1", "OD2"], {BOTH: ["HD2"]})],
    "GLU": [("SC1", "Q5n", -1.0, ["CB", "CG", "CD", "OE1", "OE2"], {})],
    "GLH": [("SC1", "P2", 0.0, ["CB", "CG", "CD", "OE1", "OE2"], {BOTH: ["HE2"]})],
    "THR": [("SC1", "SP1", 0.0, ["CB", "OG1", "CG2"], {BOTH: ["HG1"]})],
    "SER": [("SC1", "TP1", 0.0, ["CB", "OG"], {BOTH: ["HG"]})],
    "LYS": [("SC1", "SC3", 0.0, ["CB", "CG", "CD"], {}),
            ("SC2", "SQ4p", 1.0, ["CE", "NZ"],
             {"amber": ["HZ1", "HZ2", "HZ3"], "pdbv2": ["1HZ", "2HZ", "3HZ"]})],
    "LYN": [("SC1", "SC3", 0.0, ["CB", "CG", "CD"], {}),
            ("SC2", "SN6d", 0.0, ["CE", "NZ"],
             {"amber": ["HZ2", "HZ3"], "pdbv2": ["2HZ", "3HZ"]})],
    "ARG": [("SC1", "SC3", 0.0, ["CB", "CG", "CD"], {}),
            ("SC2", "SQ3p", 1.0, ["NE", "CZ", "NH1", "NH2"],
             {"amber": ["HE", "HH11", "HH12", "HH21", "HH22"],
              "pdbv2": ["HE", "1HH1", "2HH1", "1HH2", "2HH2"]})],
    # HIS is AMBER's delta-protonated histidine, which is what both local
    # authorities carry under that name: data/reducerules/amber.grp lists
    # HD1 and no HE2, and example/001_GKinase/model.pdb does too. The bead
    # types are martini3001's HSD/HID entry. HIE and HIP are separate.
    "HIS": [("SC1", "TC4", 0.0, ["CB", "CG"], {}),
            ("SC2", "TN5a", 0.0, ["CD2", "NE2"], {BOTH: ["HD2"]}),
            ("SC3", "TN6d", 0.0, ["ND1", "CE1"], {BOTH: ["HD1", "HE1"]})],
    "HIE": [("SC1", "TC4", 0.0, ["CB", "CG"], {}),
            ("SC2", "TN6d", 0.0, ["CD2", "NE2"], {BOTH: ["HD2", "HE2"]}),
            ("SC3", "TN5a", 0.0, ["ND1", "CE1"], {BOTH: ["HE1"]})],
    "HIP": [("SC1", "TC4", 0.0, ["CB", "CG"], {}),
            ("SC2", "TP1dq", 0.5, ["CD2", "NE2"], {BOTH: ["HD2", "HE2"]}),
            ("SC3", "TP1dq", 0.5, ["ND1", "CE1"], {BOTH: ["HD1", "HE1"]})],
    "PHE": [("SC1", "SC4", 0.0, ["CB", "CG", "CD1"], {BOTH: ["HD1"]}),
            ("SC2", "TC5", 0.0, ["CE1"], {BOTH: ["HE1"]}),
            ("SC3", "TC5", 0.0, ["CZ", "CD2", "CE2"], {BOTH: ["HZ", "HD2", "HE2"]})],
    "TYR": [("SC1", "TC4", 0.0, ["CB", "CG"], {}),
            ("SC2", "TC5", 0.0, ["CD1", "CE1"], {BOTH: ["HD1", "HE1"]}),
            ("SC3", "TC5", 0.0, ["CD2", "CE2"], {BOTH: ["HD2", "HE2"]}),
            ("SC4", "TN6", 0.0, ["CZ", "OH"], {BOTH: ["HH"]})],
    # TRP SC3 is a VIRTUAL SITE in Martini 3 (mass 0, placed by
    # [ virtual_sitesn ] from SC1/SC2/SC4/SC5 -- the only one in the whole
    # of martini3001). BioSpring has no virtual sites, so it is a real bead
    # here, mapped onto the two ring-fusion carbons it sits between, and it
    # carries its size class's mass. README.md states what that costs.
    "TRP": [("SC1", "TC4", 0.0, ["CB", "CG"], {}),
            ("SC2", "TN6d", 0.0, ["CD1", "NE1"], {BOTH: ["HD1", "HE1"]}),
            ("SC3", "TC5", 0.0, ["CE2", "CD2"], {}),
            ("SC4", "TC5", 0.0, ["CE3", "CZ3"], {BOTH: ["HE3", "HZ3"]}),
            ("SC5", "TC5", 0.0, ["CZ2", "CH2"], {BOTH: ["HZ2", "HH2"]})],
}

# Backbone composition. BB is the centroid of N, CA, C, O and the amide
# hydrogen ("HN BB", weight 1, in every .amber.map). Proline and
# hydroxyproline have no amide hydrogen. A free N-terminus carries three
# hydrogens and a free C-terminus the extra OXT, both from
# mappings/martini3001/modifications.amber.mapping.
BB_HEAVY = ["N", "CA", "C", "O"]
NO_AMIDE_H = {"PRO", "HYP"}

INTERIOR, NTER, CTER = "", "NTER", "CTER"


def backbone(resname, position):
    """(heavy atoms, {convention: [hydrogens]}) of a residue's BB bead."""
    heavy = list(BB_HEAVY)
    if position == CTER:
        heavy = heavy + ["OXT"]
    if position == NTER:
        if resname in NO_AMIDE_H:
            return heavy, {"amber": ["H2", "H3"], "pdbv2": ["2H", "3H"]}
        return heavy, {"amber": ["H1", "H2", "H3"], "pdbv2": ["1H", "2H", "3H"]}
    if resname in NO_AMIDE_H:
        return heavy, {}
    return heavy, {BOTH: ["H"]}


def labels_of(resname):
    """A residue's bead labels, in Martini's own order: BB, SC1, SC2..."""
    return ["BB"] + [label for label, _t, _c, _h, _hy in SIDECHAIN[resname]]


@functools.lru_cache(maxsize=None)
def _grain_names(resname):
    """Grain name per bead label: Martini's own BEAD TYPE.

    The grain name is the force-field key (ReduceRule::forcefield_key), so
    naming a grain after its Martini type makes the .ff Martini's own type
    table instead of a per-residue expansion of it: 90 keys rather than 207,
    and "name TC5" in a viewer selects every aromatic ring vertex.

    Two properties of the engine constrain the spelling, and both were
    measured rather than assumed, on the 401 beads and 982 springs of
    example/001_GKinase/model.pdb:

    * A name is an ADDRESS. RigidBodyBuilder::_resolve_atom returns the FIRST
      particle bearing it, so beads of one residue sharing a type would share
      one address. Six residues do share one -- PHE and TYR (TC5 twice), TRP
      (TC5 three times), HIP (TP1dq twice), GLH (P2 on BB and SC1), HYP (P1
      likewise) -- and naming them by bare type left 10 beads of 401 with
      ZERO springs: the ring vertices of every PHE, TRP and TYR. A repeat
      inside a residue therefore takes a _2, _3 suffix.

    * A name must NOT change with chain position. Martini's terminal
      backbone bead is Q5 where the interior one is P2/SP1/SP2/SP2a/P1, but
      the .rbody reverse lookup (RigidBodyBuilder::_translate) goes through
      data/reducerules/martini3.names.grp, and the engine fetches that table by
      residue name alone -- ReduceRuleContainer::get_rules_for_residue
      ignores the position and _translate returns the first match. Renaming
      the terminal bead Q5 dropped 6 of the 982 springs, silently, all at the
      two chain ends. So the name here is always the INTERIOR type and the
      terminal change rides on the force-field KEY: "P2:NTER" carries Q5's
      +1 e. That is what ReduceRule::forcefield_key is for, and what
      amber.grp already does with "AN".

    The cost is that a grain name no longer says which residue it belongs to,
    and that at a free chain end it names the interior type while carrying
    Martini's terminal parameters. Nothing in the engine reads a reduce rule
    by name alone outside the tests, and the resulting model is identical:
    same 401 beads, the same 982-spring SET, and zero difference in charge,
    radius, epsilon, mass or position against the residue-prefixed naming.
    """
    seen, names = collections.Counter(), {}
    types = [BB_TYPE.get(resname, BB_TYPE_DEFAULT)]
    types += [btype for _l, btype, _c, _h, _hy in SIDECHAIN[resname]]
    for label, btype in zip(labels_of(resname), types):
        seen[btype] += 1
        names[label] = btype if seen[btype] == 1 else "%s_%d" % (btype, seen[btype])
    return names


def grain_name(resname, label):
    """The grain name a reduce rule gives this residue's bead."""
    return _grain_names(resname)[label]


def beads(resname, position):
    """Every bead of a residue at a chain position, in Martini's order.

    Yields (grain name, bead type, charge, heavy atoms, hydrogens-by-convention).
    """
    heavy, hydrogens = backbone(resname, position)
    if position == INTERIOR:
        bb_type, bb_charge = BB_TYPE.get(resname, BB_TYPE_DEFAULT), 0.0
    else:
        bb_type = BB_TERMINAL_TYPE
        bb_charge = 1.0 if position == NTER else -1.0
    yield grain_name(resname, "BB"), bb_type, bb_charge, heavy, hydrogens
    for label, btype, charge, sc_heavy, sc_h in SIDECHAIN[resname]:
        yield grain_name(resname, label), btype, charge, sc_heavy, sc_h


def residue_field(resname, position):
    return resname if position == INTERIOR else resname + ":" + position


def forcefield_key(name, position):
    return name if position == INTERIOR else name + ":" + position


def wrap(text, width=78, prefix="# "):
    out, line = [], prefix
    for word in text.split():
        if len(line) + len(word) + 1 > width and line != prefix:
            out.append(line.rstrip())
            line = prefix
        line += word + " "
    out.append(line.rstrip())
    return "\n".join(out)


GRP_HEADER = """\
Reduce rules mapping an all-atom protein onto Martini 3.0.0 coarse-grain
beads. Generated by data/martini3/prepare_martini3_data.py -- do not edit by
hand. Pair it with data/forcefield/martini3.nbi.ff, which declares every grain
name this file produces, terminal qualifiers included.

GRAIN NAMES ARE MARTINI'S BEAD TYPES: "SP2 ALA N CA C O" says alanine's
backbone bead is of type SP2. A grain name is the force-field key
(ReduceRule::forcefield_key), so the companion .nbi.ff is Martini's own type
table -- 90 keys for these 207 rules -- and "name TC5" in a viewer selects
every aromatic ring vertex. Two spellings are forced by the engine and both
were measured: a type repeating inside one residue takes a _2 or _3 suffix,
and at a free chain end the name stays the INTERIOR type while the qualified
key carries Martini's terminal parameters. See _grain_names in the generator
for what each costs when ignored -- 10 beads with no spring at all, and 6
springs lost silently at the chain ends.

Unlike data/martini/reducerules/martini.grp, this is a REAL reduction: it
takes an atomic structure and builds the beads itself, so martinize is not
in the pipeline. Atom names are AMBER's.

HEAVY ATOMS ONLY, AND THAT IS THE ACCURATE CHOICE, not a shortcut. A Martini
bead is the MASS-weighted average of its atoms (martinize2 runs vermouth's
DoAverageBead with center_weight = "mass"); a BioSpring grain is
measure::centroid, an unweighted mean. Over the 401 beads of
example/001_GKinase/model.pdb, listing the heavy atoms only lands 0.028 A
from the true Martini bead in the median and 0.415 A at worst, where listing
Martini's weight-1 atoms with their hydrogens lands 0.342 A and 0.730 A --
twelve times further, because an unweighted centroid gives a hydrogen the
same pull as a carbon and a mass-weighted one gives it 1/13 of that. The
worst residual is on the sulfur beads: CYS SC1 is CB at 12 Da and SG at 32.

The input therefore does NOT need hydrogens. Any it has are claimed by no
rule and reported one per atom, which is correct and expected -- they are not
part of any bead here.

A free chain end gets its own ":NTER"/":CTER" rules, which REPLACE the
unqualified ones rather than adding to them, and its own force-field keys so
that it carries Martini's terminal charge. Their atom lists differ from the
interior one only by OXT at the C end, the N end having no heavy atom of its
own.\
"""

FF_HEADER = """\
#type	charge(e)	radius(A)	epsilon(kJ.mol-1)	mass(Da)	transferIMP(kJ.mol-1.A-2)	Hydrophobicity(sqrt(kJ.mol-1.A-1))\
"""

FF_PREAMBLE = """\
Per-type non-bonded parameters for Martini 3.0.0 proteins. Generated by
data/martini3/prepare_martini3_data.py -- do not edit by hand. Pair it with
data/reducerules/martini3.grp, which produces exactly the grain names declared
here.

THE KEYS ARE MARTINI'S BEAD TYPES, because a grain name is the force-field
key and the .grp names each grain after its type. So this file is Martini's
type table and not a per-residue expansion of it: 90 keys for 207 reduce
rules. Where one residue holds two beads of a single type the second is
suffixed -- PHE and TYR TC5_2, TRP TC5_2 and TC5_3, HIP TP1dq_2, GLH P2_2,
HYP P1_2 -- since a grain name is an address and two beads sharing one would
leave the second with no spring at all.

READ THIS BEFORE TRUSTING THE EPSILON COLUMN.

Martini 3 does not define sigma and epsilon per bead type. It defines them
PER PAIR: martini_v3.0.0.itp holds 355746 [ nonbond_params ] lines, 276 of
them over the 23 types a protein uses. A BioSpring .ff has one line per
type, and BioSpring combines two of them with a fixed rule -- the geometric
mean on epsilon (biospring_lorentz_berthelot_epsilon, hard-coded, no .msp
setting reaches it) and, for the 12-6 law, the SUM on the radius column.

So this file cannot carry Martini's matrix, and does not pretend to. What
it carries instead:

  RADIUS. BioSpring's 12-6 is eps*(R/r)^12 - 2*eps*(R/r)^6, whose minimum
  sits at r = R, so the radius column is rmin/2, not sigma/2. It is written
  here as 2^(1/6) * sigma_ii / 2, which makes a like-like pair exact. On the
  protein subset the sum rule then reproduces Martini's sigma exactly for
  every same-size pair (R-R 4.700, S-S 4.100, T-T 3.400 A) and overshoots
  every mixed-size pair by exactly 0.100 A, because Martini subtracts a flat
  0.100 A there (R-S 4.300 where the mean is 4.400, R-T 3.950 where it is
  4.050, S-T 3.650 where it is 3.750). 174 of the 253 cross pairs carry that
  0.100 A, which is 2.5 % of sigma.

  EPSILON. This column is Martini's eps_ii, the diagonal of its own matrix.
  Like-like pairs are therefore exact. Cross pairs are sqrt(eps_i * eps_j),
  which agrees with Martini on 9 of 253 (3.6 %) and misses by 0.626 kJ/mol
  in the median, 2.436 at the 90th percentile and 3.300 at worst -- P5-SQ5n,
  where Martini holds 7.58 kJ/mol and the geometric mean gives 4.28. The
  worst four are all polar-with-charged pairs, which is what a salt bridge
  is made of, and they are all UNDERestimated.

  Fitting the column instead of taking the diagonal was measured and
  rejected: the best rank-1 fit in log space halves the median to 0.342 and
  the p90 to 1.367, but raises the worst case to 4.284 and has to put SQ3p
  at 7.80 kJ/mol where Martini says 3.31, because Q-Q pairs are weak while
  Q-P pairs are strong and one number per type cannot be both. The diagonal
  at least never states a value Martini does not hold.

  TRANSFERIMP is 0.0 on every line: the IMPALA term needs an oil/water
  partition dataset per bead type, and no such dataset was published with
  Martini 3.0.0. data/martini/forcefield/martini_{octanol,hexadecane}.ff
  carry one for the 3.0.b.3.2 beta -- those numbers belong to those bead
  names, not to these. Leave impala.enable at 0 with this file.

  HYDROPHOBICITY is 0.0 on every line, as in every force field shipped with
  BioSpring, so hydrophobicity.enable does nothing here.

  TERMINI carry Martini's own charges, under the qualified keys
  "<grain>:NTER" and "<grain>:CTER" that ReduceRule::forcefield_key asks
  for: a free N-terminus is Q5 at +1 e and a free C-terminus Q5 at -1 e,
  from modifications.ff. Measured on example/001_GKinase: 22 beads at +1 e
  and 22 at -1 e.

Masses are Martini's nominal 72/54/36 Da for regular, small and tiny beads,
from [ atomtypes ]. TRP SC3 is a massless virtual site in Martini and a real
36 Da bead here, so a tryptophan weighs 36 Da more than Martini's.\
"""


def rules_for(position):
    """Every .grp line for a chain position, as (grain, residue field, atoms).

    HEAVY ATOMS ONLY, and that is not a simplification -- it is the closest a
    BioSpring .grp can get to Martini. See the module docstring.
    """
    lines = []
    for resname in sorted(SIDECHAIN, key=lambda r: CODE[r]):
        field = residue_field(resname, position)
        for name, _type, _charge, heavy, _hydrogens in beads(resname, position):
            lines.append((name, field, list(heavy)))
    return lines


def write_grp(path, header):
    with open(path, "w") as out:
        out.write(wrap(" ".join(header.split("\n\n")[0].split())) + "\n#\n")
        for paragraph in header.split("\n\n")[1:]:
            out.write(wrap(paragraph) + "\n#\n")
        for position, title in ((INTERIOR, "Interior residues"),
                                (NTER, "Free N-terminus"),
                                (CTER, "Free C-terminus")):
            out.write("#\n# --- %s ---\n" % title)
            for name, field, atoms in rules_for(position):
                out.write("%s %s %s\n" % (name, field, " ".join(atoms)))
    return path


def write_ff(path):
    with open(path, "w") as out:
        for paragraph in FF_PREAMBLE.split("\n\n"):
            out.write(wrap(paragraph) + "\n#\n")
        out.write(FF_HEADER + "\n")
        seen = set()
        for position in (INTERIOR, NTER, CTER):
            for resname in sorted(SIDECHAIN, key=lambda r: CODE[r]):
                for name, btype, charge, _h, _hy in beads(resname, position):
                    key = forcefield_key(name, position)
                    if key in seen:
                        continue
                    seen.add(key)
                    mass, sigma, epsilon = TYPE[btype]
                    radius = TWO_ONE_SIXTH * sigma / 2.0
                    out.write("%s\t%.2f\t%.5f\t%.3f\t%.4f\t%.4f\t%.1f\n"
                              % (key, charge, radius, epsilon, mass, 0.0, 0.0))
    return path



# =====================================================================
# BONDED TOPOLOGY -- which springs exist, and nothing about how stiff.
#
# A .rbody line creates every pairwise spring between the atoms it names
# and leaves the particles distinct, so a TWO-ATOM group is exactly one
# spring. Every group below has two atoms: the network is then Martini's
# bonded graph and not a distance cutoff's guess at it.
#
# No force constant appears in this file, and none belongs here. A spring's
# stiffness is a CONSTRUCTION input: pdb2spn's -stiffness sets it, the .nc
# stores it per spring, the .msp's spring group cannot change it, and the
# only option that reads per-pair constants from a .bi.ff
# (-bondstiffness) is off by default. Martini's own constants, where it
# states any, go in the companion .bi.ff.
#
# WHAT IS IN HERE, from martini3001's aminoacids.ff:
#   - every side-chain [bonds] and [constraints] pair, deduplicated: the
#     two sections declare the SAME pair twice for most side chains, once
#     under "ifndef FLEXIBLE" and once under "ifdef FLEXIBLE";
#   - the 1-3 pair of every side-chain [angles] entry, which is how
#     BioSpring holds an angle;
#   - from the 39 [link] blocks, the backbone bond BB-+BB and the 1-3
#     pairs of the backbone angles.
#
# WHAT CANNOT BE IN HERE, and why:
#   - "BB ++BB" and "BB +++BB", the elastic bonds that stabilise extended
#     regions and helices. RigidBodyBuilder strips ONE +/- prefix
#     (name.substr(1)) and moves exactly one residue, so a two- or
#     three-residue offset resolves to an atom named "+BB" in the next
#     residue, which does not exist. Four bond patterns and two of the five
#     dihedral patterns are lost this way.
#   - "SC1 >SC1", the disulfide bridge, which is not between neighbours at
#     all. BioSpring has -static-disulfide for that.
#   - the secondary-structure choice. Martini declares three alternative
#     BB-+BB bonds (0.350 / 0.3375 / 0.365 nm) and three constraints
#     (0.310 / 0.33 / 0.3375 nm), and martinize2 picks by DSSP. A .rbody
#     keys on residue name and chain position, never on secondary
#     structure. It costs nothing here: a BioSpring spring equilibrates at
#     the distance found in the INPUT structure, which in a helix already
#     IS 0.310 nm -- the same measurement Martini's own table came from.
#     The force constant is 4000 kJ.mol-1.nm-2 on all three bond variants
#     anyway, so only the bond-or-constraint distinction varies, and that
#     is the construction choice above.
#
# TRP SC3 GETS SPRINGS MARTINI DOES NOT GIVE IT. There, SC3 is a virtual
# site placed from SC1/SC2/SC4/SC5 and carries no bonded term at all. Here
# it is a real bead (BioSpring has no virtual sites), so without springs it
# would be held by nothing but the non-bonded terms and would drift out of
# the ring. It is tied to the four beads its virtual-site definition uses.

# resname -> the side-chain pairs, bonds/constraints and angle 1-3 pairs
# together, since in this file they are the same thing: a spring.
SIDECHAIN_PAIRS = {
    "GLY": [],
    "ALA": [("BB", "SC1")],
    "CYS": [("BB", "SC1")],
    "VAL": [("BB", "SC1")],
    "LEU": [("BB", "SC1")],
    "ILE": [("BB", "SC1")],
    "MET": [("BB", "SC1")],
    "PRO": [("BB", "SC1")],
    "HYP": [("BB", "SC1")],
    "ASN": [("BB", "SC1")],
    "GLN": [("BB", "SC1")],
    "ASP": [("BB", "SC1")],
    "ASH": [("BB", "SC1")],
    "GLU": [("BB", "SC1")],
    "GLH": [("BB", "SC1")],
    "THR": [("BB", "SC1")],
    "SER": [("BB", "SC1")],
    "LYS": [("BB", "SC1"), ("BB", "SC2"), ("SC1", "SC2")],
    "LYN": [("BB", "SC1"), ("BB", "SC2"), ("SC1", "SC2")],
    "ARG": [("BB", "SC1"), ("BB", "SC2"), ("SC1", "SC2")],
    "HIS": [("BB", "SC1"), ("BB", "SC2"), ("BB", "SC3"),
            ("SC1", "SC2"), ("SC1", "SC3"), ("SC2", "SC3")],
    "PHE": [("BB", "SC1"), ("BB", "SC2"), ("BB", "SC3"),
            ("SC1", "SC2"), ("SC1", "SC3"), ("SC2", "SC3")],
    "TYR": [("BB", "SC1"), ("BB", "SC2"), ("BB", "SC3"),
            ("SC1", "SC2"), ("SC1", "SC3"), ("SC2", "SC3"),
            ("SC2", "SC4"), ("SC3", "SC4")],
    # SC3 is Martini's virtual site; the four pairs tying it in are ours.
    "TRP": [("BB", "SC1"), ("BB", "SC2"), ("BB", "SC4"),
            ("SC1", "SC2"), ("SC1", "SC4"), ("SC2", "SC4"),
            ("SC2", "SC5"), ("SC4", "SC5"),
            ("SC1", "SC3"), ("SC2", "SC3"), ("SC4", "SC3"), ("SC5", "SC3")],
}
SIDECHAIN_PAIRS["HIE"] = SIDECHAIN_PAIRS["HIS"]
SIDECHAIN_PAIRS["HIP"] = SIDECHAIN_PAIRS["HIS"]

# The backbone pairs every residue carries, from the [link] blocks.
# "BB +BB" is the bond; the rest are the 1-3 pairs of the backbone angles
# -BB/BB/+BB, -BB/BB/SC1 and SC1/BB/+BB. "BB +SC1" is deliberately absent:
# it is the same pair as the next residue's "-BB SC1".
BACKBONE_PAIRS = [("BB", "+BB"), ("-BB", "+BB")]
BACKBONE_PAIRS_WITH_SC = [("-BB", "SC1"), ("SC1", "+BB")]


def _group_name(code, a, b):
    sign = {"+": "p", "-": "m"}
    def tag(x):
        return sign.get(x[0], "") + x.lstrip("+-")
    return "M3_%s_%s_%s" % (code, tag(a), tag(b))


def rbody_lines():
    """Every .rbody line, as (group name, residue, [atom1, atom2])."""
    out = []
    for resname in sorted(SIDECHAIN_PAIRS, key=lambda r: CODE[r]):
        code = CODE[resname]
        has_sc = bool(SIDECHAIN[resname])
        pairs = list(BACKBONE_PAIRS)
        if has_sc:
            pairs += BACKBONE_PAIRS_WITH_SC
        pairs += SIDECHAIN_PAIRS[resname]
        for a, b in pairs:
            out.append((_group_name(code, a, b), resname, [a, b]))
    return out


RBODY_HEADER = """\
Bonded topology for the Martini 3.0.0 coarse-grain protein model: which
springs exist. Generated by data/martini3/prepare_martini3_data.py -- do not
edit by hand. Use it with -rigidbody, which is exclusive with -cutoff:
a distance cutoff would invent a network, and the point of this file is that
Martini's own bonded graph is used instead.

Every group here has exactly TWO atoms, so every group is exactly one
spring. Pair it with data/reducerules/martini3.grp, whose bead names these are.

NO FORCE CONSTANT APPEARS IN THIS FILE and none belongs here: a spring's
stiffness is a construction input (pdb2spn's -stiffness, stored per spring
in the .nc, unreachable from the .msp). Martini's own constants, where it
states any, are in the companion .bi.ff and are read only if
-bondstiffness is given.\
"""


def write_rbody(path):
    with open(path, "w") as out:
        for paragraph in RBODY_HEADER.split("\n\n"):
            out.write(wrap(paragraph) + "\n#\n")
        for paragraph in _RBODY_NOTES:
            out.write(wrap(paragraph) + "\n#\n")
        current = None
        for name, resname, atoms in rbody_lines():
            if resname != current:
                out.write("#\n")
                current = resname
            out.write("%s %s %s\n" % (name, resname, " ".join(atoms)))
    return path


_RBODY_NOTES = [
    "WHAT IS HERE, from martini3001's aminoacids.ff: every side-chain "
    "[bonds] and [constraints] pair, deduplicated because the two sections "
    "declare the same pair twice for most side chains (once under 'ifndef "
    "FLEXIBLE' and once under 'ifdef FLEXIBLE'); the 1-3 pair of every "
    "side-chain [angles] entry, which is how BioSpring holds an angle; and "
    "from the 39 [link] blocks, the backbone bond BB-+BB plus the 1-3 pairs "
    "of the backbone angles.",

    "WHAT CANNOT BE HERE: 'BB ++BB' and 'BB +++BB', the elastic bonds that "
    "stabilise extended regions and helices, because RigidBodyBuilder strips "
    "ONE +/- prefix and moves exactly one residue -- a two-residue offset "
    "resolves to an atom named '+BB' in the next residue, which does not "
    "exist. Four of the backbone bond patterns and two of the five dihedral "
    "patterns are lost that way. Nor is 'SC1 >SC1', the disulfide bridge, "
    "which is not between neighbours at all: -static-disulfide covers it.",

    "THE SECONDARY-STRUCTURE CHOICE IS ALSO ABSENT, and costs nothing. "
    "Martini declares three alternative BB-+BB bonds (0.350 / 0.3375 / "
    "0.365 nm) and three constraints (0.310 / 0.33 / 0.3375 nm), picked by "
    "DSSP; a .rbody keys on residue name and chain position, never on "
    "secondary structure. But a BioSpring spring equilibrates at the "
    "distance in the INPUT structure, which in a helix already IS 0.310 nm "
    "-- the same measurement Martini's table came from. And the force "
    "constant is 4000 kJ.mol-1.nm-2 on all three bond variants, so only the "
    "bond-or-constraint distinction varies, which is the construction "
    "choice above.",

    "TRP SC3 GETS FOUR SPRINGS MARTINI DOES NOT GIVE IT. There it is a "
    "virtual site placed from SC1/SC2/SC4/SC5 and carries no bonded term at "
    "all. Here it is a real bead, so with no spring it would be held by "
    "nothing but the non-bonded terms and would drift out of the ring. It is "
    "tied to the four beads its virtual-site definition uses.",
]



# =====================================================================
# BONDED FORCE CONSTANTS -- Martini's own, where Martini states any.
#
# A STIFFNESS record RETUNES an existing spring and creates none, so every
# record below names a pair that data/reducerules/martini3.rbody already springs.
# It is read only when -bondstiffness is given; without it the whole mesh
# stays at the run's uniform -stiffness.
#
# All of Martini's values are converted out of GROMACS units here: a bond
# force constant in kJ.mol-1.nm-2 divided by 100 is kJ.mol-1.A-2, and both
# engines write the bond as 0.5 * k * (r - r0)^2 (GROMACS funct 1;
# biospring_spring_energy in forcefield/shared/spring_shared.h), so there is
# no factor of two to put in.
# =====================================================================

# The backbone, for a COIL residue. Martini's backbone terms are chosen by
# DSSP and this file carries the coil set, because BioSpring assigns no
# secondary structure. The helix and extended values are in the header.
BB_BB_LENGTH = 3.50          # A, from "BB +BB 1 0.350 4000"
BB_BB_STIFFNESS = 40.0       # kJ.mol-1.A-2, the same record's 4000 / 100
BBB_ANGLE = (10, 127.0, 20.0)   # funct, theta0 in deg, k in kJ.mol-1
BBS_ANGLE = (2, 100.0, 25.0)    # "-BB BB SC1", BBS angles regular martini
SBB_ANGLE = (2, 100.0, 25.0)    # "SC1 BB +BB", First SBB regular martini

SC_BOND = {
    ('ARG', 'BB', 'SC1'):        (3.300, 50.000),
    ('ARG', 'SC1', 'SC2'):       (3.800, 50.000),
    ('ASH', 'BB', 'SC1'):        (3.520, 75.000),
    ('ASN', 'BB', 'SC1'):        (3.520, 50.000),
    ('ASP', 'BB', 'SC1'):        (3.520, 75.000),
    ('CYS', 'BB', 'SC1'):        (3.410, 75.000),
    ('GLH', 'BB', 'SC1'):        (4.000, 50.000),
    ('GLN', 'BB', 'SC1'):        (4.000, 50.000),
    ('GLU', 'BB', 'SC1'):        (4.000, 50.000),
    ('HIE', 'BB', 'SC1'):        (3.360, 75.000),
    ('HIP', 'BB', 'SC1'):        (3.360, 75.000),
    ('HIS', 'BB', 'SC1'):        (3.360, 75.000),
    ('HYP', 'BB', 'SC1'):        (3.000, 75.000),
    ('LEU', 'BB', 'SC1'):        (3.630, 75.000),
    ('LYN', 'BB', 'SC1'):        (3.300, 50.000),
    ('LYN', 'SC1', 'SC2'):       (3.600, 50.000),
    ('LYS', 'BB', 'SC1'):        (3.300, 50.000),
    ('LYS', 'SC1', 'SC2'):       (3.600, 50.000),
    ('MET', 'BB', 'SC1'):        (4.000, 25.000),
    ('PHE', 'BB', 'SC1'):        (3.250, 75.000),
    ('PRO', 'BB', 'SC1'):        (3.300, 75.000),
    ('SER', 'BB', 'SC1'):        (2.870, 75.000),
    ('TRP', 'BB', 'SC1'):        (3.150, 50.000),
    ('TYR', 'BB', 'SC1'):        (3.250, 50.000),
}

SC_CONSTRAINT = {
    ('ALA', 'BB', 'SC1'):        2.700,
    ('HIE', 'SC1', 'SC2'):       3.200,
    ('HIE', 'SC1', 'SC3'):       3.000,
    ('HIE', 'SC2', 'SC3'):       2.700,
    ('HIP', 'SC1', 'SC2'):       3.200,
    ('HIP', 'SC1', 'SC3'):       3.200,
    ('HIP', 'SC2', 'SC3'):       2.700,
    ('HIS', 'SC1', 'SC2'):       3.200,
    ('HIS', 'SC1', 'SC3'):       3.000,
    ('HIS', 'SC2', 'SC3'):       2.700,
    ('ILE', 'BB', 'SC1'):        3.410,
    ('PHE', 'SC1', 'SC2'):       3.400,
    ('PHE', 'SC1', 'SC3'):       3.400,
    ('PHE', 'SC2', 'SC3'):       2.900,
    ('THR', 'BB', 'SC1'):        3.050,
    ('TRP', 'SC1', 'SC2'):       3.350,
    ('TRP', 'SC1', 'SC4'):       4.040,
    ('TRP', 'SC2', 'SC4'):       4.700,
    ('TRP', 'SC2', 'SC5'):       4.120,
    ('TRP', 'SC4', 'SC5'):       2.930,
    ('TYR', 'SC1', 'SC2'):       3.000,
    ('TYR', 'SC1', 'SC3'):       3.000,
    ('TYR', 'SC2', 'SC3'):       3.000,
    ('TYR', 'SC2', 'SC4'):       2.850,
    ('TYR', 'SC3', 'SC4'):       2.850,
    ('VAL', 'BB', 'SC1'):        2.920,
}

SC_ANGLE = {
    ('ARG', 'BB', 'SC1', 'SC2'): (2, 180.0, 25.0),
    ('HIE', 'BB', 'SC1', 'SC2'): (2, 120.0, 50.0),
    ('HIE', 'BB', 'SC1', 'SC3'): (2, 120.0, 50.0),
    ('HIP', 'BB', 'SC1', 'SC2'): (2, 120.0, 50.0),
    ('HIP', 'BB', 'SC1', 'SC3'): (2, 120.0, 50.0),
    ('HIS', 'BB', 'SC1', 'SC2'): (2, 120.0, 50.0),
    ('HIS', 'BB', 'SC1', 'SC3'): (2, 120.0, 50.0),
    ('LYN', 'BB', 'SC1', 'SC2'): (2, 180.0, 25.0),
    ('LYS', 'BB', 'SC1', 'SC2'): (2, 180.0, 25.0),
    ('PHE', 'BB', 'SC1', 'SC2'): (2, 120.0, 50.0),
    ('PHE', 'BB', 'SC1', 'SC3'): (2, 120.0, 50.0),
    ('TRP', 'BB', 'SC1', 'SC2'): (2, 120.0, 60.0),
    ('TRP', 'BB', 'SC1', 'SC4'): (2, 130.0, 60.0),
    ('TYR', 'BB', 'SC1', 'SC2'): (2, 120.0, 60.0),
    ('TYR', 'BB', 'SC1', 'SC3'): (2, 120.0, 60.0),
}

# GROMACS improper harmonic (dihedral funct 2): 0.5 * k * (xi - xi0)^2, with
# k in kJ.mol-1.rad-2. These two hold the aromatic rings flat and are the only
# dihedrals Martini 3 applies to a protein without -scfix or a DSSP string.
SC_IMPROPER = {
    "TRP": ("SC5", "SC4", "SC2", "SC1", 180.0, 100.0),
    "TYR": ("SC4", "SC2", "SC3", "SC1", 180.0,  50.0),
}

IMPROPER_BINS = 512          # as in data/reducerules/ProteinAtomTorsion.bi.ff
KT300 = 2.494                # kJ.mol-1, k_B T at 300 K


def sc_length(resname, a, b):
    """Martini's equilibrium distance for a side-chain pair, bond or constraint."""
    for key in ((resname, a, b), (resname, b, a)):
        if key in SC_BOND:
            return SC_BOND[key][0]
        if key in SC_CONSTRAINT:
            return SC_CONSTRAINT[key]
    return None


def angle_to_distance(funct, theta0, r12, r23):
    """Turn an angle force constant into a 1-3 distance one.

    Returns (factor, r13_at_equilibrium): multiply Martini's k by the factor.

    DERIVATION. The law of cosines gives r13^2 = r12^2 + r23^2
    - 2 r12 r23 cos(theta), so with the two bonds held at their own
    equilibrium, dr13/dtheta = r12 r23 sin(theta) / r13. Matching the second
    derivative at the minimum, k_dist = V''(theta0) / (dr13/dtheta)^2.

    V''(theta0) is k sin^2(theta0) for the cosine-harmonic form (GROMACS
    funct 2, V = k/2 (cos theta - cos theta0)^2) and exactly k for the
    restricted bending form (funct 10, the same over sin^2 theta), which is
    why the sine cancels in one case and not the other:

        funct  2:  k_dist = k * (r13 / (r12 r23))^2
        funct 10:  k_dist = k * (r13 / (r12 r23))^2 / sin^2(theta0)

    At theta0 = 180 deg -- arginine's and lysine's BB-SC1-SC2 -- both
    V''(theta0) and dr13/dtheta vanish and the ratio is a 0/0 whose limit the
    funct-2 formula already is: expanding in delta = pi - theta gives
    V = k delta^4 / 8 and r13_0 - r13 = r12 r23 delta^2 / (2 r13_0), so the
    quartic in the angle is EXACTLY quadratic in the distance, with that same
    constant. No special case is needed.
    """
    t = math.radians(theta0)
    r13 = math.sqrt(r12 * r12 + r23 * r23 - 2.0 * r12 * r23 * math.cos(t))
    factor = (r13 / (r12 * r23)) ** 2
    if funct == 10:
        factor /= math.sin(t) ** 2
    return factor, r13


def angle_energy(theta_deg, funct, theta0, k):
    """Martini's own angle energy in kJ/mol, exactly as GROMACS writes it."""
    c, c0 = math.cos(math.radians(theta_deg)), math.cos(math.radians(theta0))
    if funct == 2:
        return 0.5 * k * (c - c0) ** 2
    s = math.sin(math.radians(theta_deg))
    return 0.5 * k * (c - c0) ** 2 / (s * s)


def conversion_error(funct, theta0, k, r12, r23):
    """(k_dist, Boltzmann-weighted RMS energy error at 300 K in kJ/mol).

    The two curves cannot be made equal by any choice of constant -- a 1-3
    distance spring is not an angle spring -- so what matters is how far
    apart they are where the model actually sits. Sampled over the angle,
    weighted by exp(-V/kT) sin(theta) out to 8 kT.
    """
    factor, r0 = angle_to_distance(funct, theta0, r12, r23)
    kd = k * factor
    num = den = 0.0
    for i in range(1, 1800):
        theta = i * 0.1
        if funct == 10 and not 2.0 <= theta <= 178.0:
            continue
        energy = angle_energy(theta, funct, theta0, k)
        if energy > 8.0 * KT300:
            continue
        t = math.radians(theta)
        r = math.sqrt(r12 * r12 + r23 * r23 - 2.0 * r12 * r23 * math.cos(t))
        weight = math.exp(-energy / KT300) * math.sin(t)
        num += weight * (0.5 * kd * (r - r0) ** 2 - energy) ** 2
        den += weight
    return kd, math.sqrt(num / den)


def bi_ff_records():
    """Every STIFFNESS record, as (residue, atom1, atom2, k, source, rms)."""
    out = []
    for resname in sorted(SIDECHAIN, key=lambda r: CODE[r]):
        out.append((resname, "BB", "+BB", BB_BB_STIFFNESS, "bond 3.500 A", None))

        funct, theta0, k = BBB_ANGLE
        kd, rms = conversion_error(funct, theta0, k, BB_BB_LENGTH, BB_BB_LENGTH)
        out.append((resname, "-BB", "+BB", kd,
                    "angle -BB-BB-+BB funct %d %.1f deg %.1f kJ.mol-1" % BBB_ANGLE, rms))

        side = sc_length(resname, "BB", "SC1")
        if side is not None:
            for atoms, angle in ((("-BB", "SC1"), BBS_ANGLE), (("SC1", "+BB"), SBB_ANGLE)):
                funct, theta0, k = angle
                r12, r23 = (BB_BB_LENGTH, side) if atoms[0] == "-BB" else (side, BB_BB_LENGTH)
                kd, rms = conversion_error(funct, theta0, k, r12, r23)
                name = "-BB-BB-SC1" if atoms[0] == "-BB" else "SC1-BB-+BB"
                out.append((resname, atoms[0], atoms[1], kd,
                            "angle %s funct %d %.1f deg %.1f kJ.mol-1" % (name, funct, theta0, k), rms))

        for (res, a, b), (r0, k) in sorted(SC_BOND.items()):
            if res == resname:
                out.append((resname, a, b, k, "bond %.3f A" % r0, None))

        for (res, a, b, c), (funct, theta0, k) in sorted(SC_ANGLE.items()):
            if res != resname:
                continue
            r12, r23 = sc_length(res, a, b), sc_length(res, b, c)
            kd, rms = conversion_error(funct, theta0, k, r12, r23)
            out.append((resname, a, c, kd,
                        "angle %s-%s-%s funct %d %.1f deg %.1f kJ.mol-1" % (a, b, c, funct, theta0, k), rms))
    return out


def improper_table(xi0, k):
    """(energy, torque) over xi in [-180, 180], 0.5 k (xi - xi0)^2 wrapped."""
    energy, torque = [], []
    for i in range(IMPROPER_BINS + 1):
        xi = -180.0 + 360.0 * i / IMPROPER_BINS
        delta = math.radians((xi - xi0 + 180.0) % 360.0 - 180.0)
        energy.append(0.5 * k * delta * delta)
        torque.append(-k * delta)
    return energy, torque


BI_FF_PREAMBLE = """\
Martini 3.0.0's own bonded force constants for proteins. Generated by
data/martini3/prepare_martini3_data.py -- do not edit by hand. It is read
only when -bondstiffness is given, alongside
-rigidbody data/reducerules/martini3.rbody, whose springs it retunes.

STIFFNESS <res> <a1> <a2> <k>, k in kJ.mol-1.A-2. TORSIONTABLE <id> <bins>
<E_0> <T_0> ... over xi in [-180, 180], E in kJ/mol and T = -dV/dxi in
kJ.mol-1.rad-1. TORSION <res> PLANARITY <a1> <a2> <a3> <a4> <table>.

A STIFFNESS RECORD RETUNES AND CREATES NOTHING, so every record here names a
pair that martini3.rbody already springs: 141 records over 141 of its 171
groups, and not one record naming a pair with no spring. Of the 30 groups
left, 26 are Martini [constraints] -- rigid, carrying NO force constant --
and 4 are the pairs tying in tryptophan's SC3, which is Martini's massless
virtual site and ours alone. Martini's "#ifdef FLEXIBLE" fallback of 1000000
kJ.mol-1.nm-2 is a workaround for engines without constraints, not a
force-field value, so nothing is written for any of them: they keep the run's
uniform -stiffness.

Every number is converted out of GROMACS units here. A force constant in
kJ.mol-1.nm-2 divided by 100 is kJ.mol-1.A-2, and both engines write a bond
as 0.5 k (r - r0)^2 -- GROMACS funct 1, and biospring_spring_energy in
forcefield/shared/spring_shared.h -- so there is no factor of two to put in.

WHY THE ANGLES ARE CONVERTED HERE WHEN ProteinAtomTorsion.bi.ff LEAVES THEM
OUT. That file says, correctly, that a 1-3 distance spring is not an angle
spring, and for AMBER that settles it: the mesh already holds valence angles
far stiffer than anything a retune would give them. Martini is the opposite
case. Its angle constants are 20 to 60 kJ.mol-1 and convert to 4.2 to 18.5
kJ.mol-1.A-2, against the 650 kJ.mol-1.A-2 that -rigidbody defaults to -- a
factor of 35 to 154. Leaving them at the default does not approximate
Martini; it replaces a flexible chain with a rigid one. So the conversion is
not a refinement here, it is the only way this file says anything at all.

THE CONVERSION IS DERIVED, AND WHAT IT COSTS IS MEASURED. See
angle_to_distance in the generator for the derivation, including why the
sine cancels for GROMACS funct 2 and not for funct 10, and why theta0 = 180
deg -- arginine's and lysine's BB-SC1-SC2 -- needs no special case even
though both V''(theta0) and dr13/dtheta vanish there. The cost is listed
term by term below: the Boltzmann-weighted RMS energy difference between
Martini's own angle curve and the distance spring, at 300 K, over the range
the term explores out to 8 kT. It runs from 0.271 to 1.685 kJ.mol-1 over the
91 converted records, worst on the backbone BBB angle, which is 0.68 of kT.
Fitting the constant by least squares instead of matching the curvature was
measured and rejected: it brings the worst case to 0.668 kJ.mol-1 but bakes
300 K into a force-field file and gives up the exact curvature, which is
what sets the small-amplitude fluctuation.

THE BACKBONE IS THE COIL SET. Martini picks its backbone terms by DSSP and
BioSpring assigns no secondary structure, so this file carries the
unconditional link -- bond 0.350 nm at 4000 kJ.mol-1.nm-2, hence 40
kJ.mol-1.A-2 -- and the coil angle, funct 10 at 127 deg and 20 kJ.mol-1.
Martini's alternatives, for anyone who wants to edit the generator: helix is
a 0.310 nm constraint with a funct 2 angle at 96 deg and 700 kJ.mol-1;
extended is 0.365 nm with funct 2 at 119.2 deg and 150; sheet is funct 10 at
134 deg and 25, bend 130 and 20, turn 100 and 20.

NOT WRITTEN, each for its own reason. Martini's EQUILIBRIUM LENGTHS: a
BioSpring spring equilibrates where the input put it, so the lengths in the
generator serve only to convert the angles. The -scfix terms, angles
SC1-BB-+BB and BB-+BB-+SC1 at funct 10, 100 deg and 15 kJ.mol-1 plus the
SC-BB-BB-SC dihedrals at 75 kJ.mol-1: optional in martinize2 and off unless
-scfix is passed. The elastic bonds BB-++BB and BB-+++BB of extended regions
(0.640 and 0.970 nm at 2500 kJ.mol-1.nm-2): RigidBodyBuilder moves exactly
one residue, so those pairs have no spring to retune. The CYS SC1->SC1
disulfide constraint at 0.240 nm: use -static-disulfide. And everything
behind cgidr, with every DSSP-conditional backbone dihedral.\
"""


def write_bi_ff(path):
    records = bi_ff_records()
    with open(path, "w") as out:
        for paragraph in BI_FF_PREAMBLE.split("\n\n"):
            out.write(wrap(paragraph) + "\n#\n")

        # One line per Martini angle term, not per record: the same term gives a
        # different constant in every residue, because the constant carries that
        # residue's own bond lengths. The spread is the range over the residues.
        out.write("# MEASURED CONVERSION ERROR, Boltzmann-weighted RMS energy difference\n"
                  "# from Martini's own angle curve at 300 K, where kT = 2.494 kJ.mol-1.\n"
                  "# One line per angle term, over the residues that carry it:\n#\n")
        out.write("#   %-34s %4s %17s %17s\n"
                  % ("Martini angle term", "n", "k (kJ.mol-1.A-2)", "rms (kJ.mol-1)"))
        groups = {}
        for res, a, b, k, source, rms in records:
            if rms is None:
                continue
            shape = source[len("angle "):]
            key = shape.split(" funct")[0]
            label = "%s %s" % (key, shape.split(" ", 1)[1])
            entry = groups.setdefault(label, [0, k, k, rms, rms])
            entry[0] += 1
            entry[1], entry[2] = min(entry[1], k), max(entry[2], k)
            entry[3], entry[4] = min(entry[3], rms), max(entry[4], rms)
        for label, (n, klo, khi, rlo, rhi) in sorted(groups.items(), key=lambda kv: -kv[1][4]):
            krange = "%.3f" % klo if khi - klo < 5e-4 else "%.3f-%.3f" % (klo, khi)
            rrange = "%.3f" % rlo if rhi - rlo < 5e-4 else "%.3f-%.3f" % (rlo, rhi)
            out.write("#   %-34s %4d %17s %17s\n" % (label, n, krange, rrange))
        out.write("#\n")

        for index, resname in enumerate(sorted(SC_IMPROPER)):
            _a, _b, _c, _d, xi0, k = SC_IMPROPER[resname]
            energy, torque = improper_table(xi0, k)
            values = " ".join("%.6g %.6g" % pair for pair in zip(energy, torque))
            out.write("TORSIONTABLE %d %d %s\n" % (index, IMPROPER_BINS, values))

        for index, resname in enumerate(sorted(SC_IMPROPER)):
            a, b, c, d, _xi0, _k = SC_IMPROPER[resname]
            out.write("TORSION\t%s\tPLANARITY\t%s\t%s\t%s\t%s\t%d\n" % (resname, a, b, c, d, index))

        for res, a, b, k, _source, _rms in records:
            out.write("STIFFNESS\t%s\t%s\t%s\t%.3f\n" % (res, a, b, k))
    return path



NAMES_HEADER = """\
Naming-translation table for use with -rigidbody: it maps a bead's Martini
label ("BB", "SC1") to the grain name the reduction gave it, which is the
bead's Martini TYPE ("SP2", "TC3"). Generated by
data/martini3/prepare_martini3_data.py -- do not edit by hand.

WHY THIS FILE EXISTS. A .rbody rule has to say "+BB" to reach the NEXT
residue's backbone bead, and it cannot know what residue that is. It cannot
write the grain name either, because the grain name is the bead's type and
the backbone type differs by residue: GLY is SP1, ALA SP2, PRO SP2a, HYP P1,
every other P2.

RigidBodyBuilder resolves exactly this: when an atom is not found under the
name the .rbody wrote, it translates that name through a table, IN THE TARGET
RESIDUE'S CONTEXT (RigidBodyBuilder::_translate takes target_resname). So
"+BB" on an alanine followed by a glycine resolves to "SP1". The table it
accepts must have EXACTLY ONE ATOM PER RULE --
_check_translation_is_one_atom_per_rule kills pdb2spn otherwise -- which is
why data/reducerules/martini3.grp cannot serve as the table and this file has to
exist.

ENUMERATING THE NEIGHBOUR INSTEAD DOES NOT WORK, and that was measured, not
assumed: writing "+P1 +P2 +SP1 +SP2 +SP2a +Q5" in place of "+BB" built 1021
springs on example/001_GKinase/model.pdb where 982 are right. Four
side-chain beads carry a backbone type -- THR SC1 is SP1, ASH SC1 SP2, GLH
SC1 P2, HYP SC1 P1 -- so both the neighbour's backbone AND its first side
chain answer, and _apply_rule springs every pair it resolves.

INTERIOR RULES ONLY, on purpose: a grain name does not change at a chain end
(only its force-field key does), and get_rules_for_residue fetches by
residue name alone, ignoring the position.

It is therefore NOT a reduction. Passing it to a run that would reduce would
destroy the model: the rule "SP2 ALA BB" looks for an atom named BB, and an
all-atom structure has none. Use it only with an input that is already
coarse-grained, where pdb2spn skips the reduction -- see the two-step recipe
in README.md.\
"""


def write_names_grp(path):
    with open(path, "w") as out:
        for paragraph in NAMES_HEADER.split("\n\n"):
            out.write(wrap(paragraph) + "\n#\n")
        for resname in sorted(SIDECHAIN, key=lambda r: CODE[r]):
            for label in labels_of(resname):
                out.write("%s %s %s\n" % (grain_name(resname, label), resname, label))
    return path


if __name__ == "__main__":
    # THE SIX FILES GO WHERE EVERY OTHER FORCE FIELD'S DO, not into a tree of
    # their own. data/forcefield/ holds amber.ff, martini.ff, model.ff and the
    # CAonly* set; data/reducerules/ holds their .grp, .rbody and .bi.ff. Martini
    # 3 sat in data/martini3/forcefield and data/martini3/reducerules until
    # 2026-11-26, and the only thing that cost was discoverability -- somebody
    # looking for the model where the models are concluded it was not in the
    # branch at all. This script keeps living here, as data/martini/'s own
    # prepare_martini_data.py does; it is the GENERATOR, not the force field.
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)                      # data/
    grp = write_grp(os.path.join(root, "reducerules", "martini3.grp"), GRP_HEADER)
    rbody = write_rbody(os.path.join(root, "reducerules", "martini3.rbody"))
    names = write_names_grp(os.path.join(root, "reducerules", "martini3.names.grp"))
    biff = write_bi_ff(os.path.join(root, "reducerules", "martini3.bi.ff"))
    ff = write_ff(os.path.join(root, "forcefield", "martini3.nbi.ff"))
    # .ff and .nbi.ff are the same file under two names, as everywhere else
    # in data/forcefield. Keep them in sync: pdb2spn's --ff reads one and
    # --nbi.ff the other, and a type missing from either is dropped silently.
    plain = os.path.join(root, "forcefield", "martini3.ff")
    with open(ff) as src, open(plain, "w") as dst:
        dst.write(src.read())
    for path in (grp, rbody, names, biff, ff, plain):
        print("wrote %s" % path)
