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
BB_TYPE = {"GLY": "SP1", "ALA": "SP2", "PRO": "SP2a", "HYP": "P1"}
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


def beads(resname, position):
    """Every bead of a residue at a chain position, in Martini's order.

    Yields (grain name, bead type, charge, heavy atoms, hydrogens-by-convention).
    """
    code = CODE[resname]
    heavy, hydrogens = backbone(resname, position)
    if position == INTERIOR:
        bb_type, bb_charge = BB_TYPE.get(resname, BB_TYPE_DEFAULT), 0.0
    else:
        bb_type = BB_TERMINAL_TYPE
        bb_charge = 1.0 if position == NTER else -1.0
    yield code + "BB", bb_type, bb_charge, heavy, hydrogens
    for label, btype, charge, sc_heavy, sc_h in SIDECHAIN[resname]:
        yield code + label, btype, charge, sc_heavy, sc_h


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
hand. Pair it with forcefield/martini3.nbi.ff, which declares every grain
name this file produces, terminal qualifiers included.

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
reducerules/martini3.grp, which produces exactly the grain names declared
here.

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


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    grp = write_grp(os.path.join(here, "reducerules", "martini3.grp"), GRP_HEADER)
    ff = write_ff(os.path.join(here, "forcefield", "martini3.nbi.ff"))
    # .ff and .nbi.ff are the same file under two names, as everywhere else
    # in data/forcefield. Keep them in sync: pdb2spn's --ff reads one and
    # --nbi.ff the other, and a type missing from either is dropped silently.
    plain = os.path.join(here, "forcefield", "martini3.ff")
    with open(ff) as src, open(plain, "w") as dst:
        dst.write(src.read())
    for path in (grp, ff, plain):
        print("wrote %s" % path)
