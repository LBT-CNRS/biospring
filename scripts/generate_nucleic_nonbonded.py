#!/usr/bin/env python3
"""Generates the NON-BONDED nucleic parameters: amber.dna.grp/.nbi.ff and
amber.rna.grp/.nbi.ff, from AMBER's own DNA.OL15 and RNA.OL3.

Companion to generate_nucleic_forcefield.py, which does the BONDED side.
Without these, a nucleic system has no charge, no radius, no epsilon and --
the part that bites first -- no mass: pdb2spn falls back to 1 Da for every
atom, so a phosphorus weighs what a hydrogen does and the dynamics is not
the molecule's.

WHERE EACH NUMBER COMES FROM
  charge   per (residue, atom), from the <Residues> templates. It has to be
           per residue: only the sugar/phosphate atoms share a charge
           across bases, and even they differ between DNA and RNA (P is
           1.1659 in OL15 and 1.1662 in OL3).
  sigma,
  epsilon  per TYPE, from <NonbondedForce>.
  mass     per TYPE, from <AtomTypes>.

UNITS -- the .nbi.ff's own header says epsilon is in kJ/mol and that is
WRONG; the file has always held kcal/mol and ForceFieldReader passes the
value through untouched. Checked against amber.nbi.ff rather than assumed:
its AN/AH/AO/AC rows read 0.17 / 0.0157 / 0.21 / 0.086, which are AMBER's
N/H/O/C epsilons in kcal/mol exactly. Same check for the radius column: it
is R* (rmin/2) in angstrom, sigma_nm * 10 * 2^(1/6) / 2, which turns
carbon's 0.339967 nm into the 1.908 that file already carries.

TYPE NAMES are capped at 4 characters: a type becomes the particle's name,
and that is written into a PDB atom-name field. Hence the scheme, which is
what the charge measurements above allow:
  * sugar and phosphate -> ONE type each, prefix 'd' (DNA) or 'r' (RNA),
    apostrophes dropped: dP, dOP1, dC5, dH52 for H5'', ...
  * C1' and H1' -> per base, since the base they carry shifts their charge,
    prefixed by a DIGIT (1C1, 2H1, ...). Not by the base letter: guanine
    carries an H1 on its N1, and dropping the apostrophe would have made
    the sugar's H1' collide with it -- caught by the self-check below,
    which had DG at -1.1774 e instead of -1.0000.
  * base atoms -> per base, prefix a/c/g/t (DNA) or a/c/g/u (RNA)
Lowercase prefixes on purpose: amber.grp's protein types are all prefixed
by an UPPERCASE residue letter, so a protein-DNA complex can concatenate
the two files with no collision -- which is the case that matters, 1A74
being exactly that. DNA and RNA types DO collide with each other (both use
a/c/g); a system holding both needs one of the two regenerated under
different prefixes.

IMPALA transfer energy and hydrophobicity are left at 0: there is no
published nucleic parameterisation for either, and inventing one would be
worse than declaring it absent.
"""

import os
import sys
import xml.etree.ElementTree as ET

import openmm.app as openmm_app

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(os.path.dirname(openmm_app.__file__), "data")

KCAL = 4.184                      # kJ/mol per kcal/mol
RSTAR = 10.0 * 2.0 ** (1.0 / 6.0) / 2.0   # sigma(nm) -> R*(A)

# Atom-name spellings, same convention and same reason as the bonded
# generator: a rule that knows one spelling matches nothing on a file
# written in the other, silently.
ALIASES = {"H5'": ["H5'1"], "H5''": ["H5'2"], "H2'": ["H2'1"], "H2''": ["H2'2"],
           "HO2'": ["HO'2"], "OP1": ["O1P"], "OP2": ["O2P"]}

# Sugar + phosphate. Shared across bases EXCEPT C1'/H1', which sit next to
# the base and take its charge with them (measured, not assumed).
BACKBONE = ["P", "OP1", "OP2", "O5'", "C5'", "H5'", "H5''", "C4'", "H4'",
            "O4'", "C3'", "H3'", "C2'", "H2'", "H2''", "O3'", "O2'", "HO2'",
            "HO5'", "HO3'"]
PER_BASE_SUGAR = ["C1'", "H1'"]


# Hydrogen-bond donor/acceptor roles, as (donor, acceptor).
#
# Same source as the protein side (HBPLUS / McDonald & Thornton, and
# Reduce/Probe), cited in ProteinDonorAcceptor.hbond's own header. Using one
# classification for both chemistries is the point: example 076 loads the
# protein and DNA tables together, and they must not disagree about what a
# donor is.
#
# N7 and N3 are the Hoogsteen and minor-groove faces, not Watson-Crick.
# They are included because HBPLUS includes them and because leaving them
# out would silently forbid every non-canonical pair -- the wobble and
# Hoogsteen geometries a folded RNA actually uses. Measured on 074's duplex,
# they cost almost nothing: BioSpring's reciprocal-best-hit rule still
# recovers 46 of the 48 Watson-Crick bonds.
PURINE = {"N3": (0, 1), "N7": (0, 1)}
BASE_ROLES = {
    "A": {"N1": (0, 1), "N6": (1, 0), **PURINE},          # N1/N6 Watson-Crick
    "G": {"N1": (1, 0), "N2": (1, 0), "O6": (0, 1), **PURINE},
    "C": {"N3": (0, 1), "N4": (1, 0), "O2": (0, 1)},
    "T": {"N3": (1, 0), "O4": (0, 1), "O2": (0, 1)},
    "U": {"N3": (1, 0), "O4": (0, 1), "O2": (0, 1)},
}
# Phosphate oxygens are strong acceptors; the esters and the furanose ring
# oxygen are weak ones. O2' is RNA's alone and is both roles at once.
SUGAR_ROLES = {"OP1": (0, 1), "OP2": (0, 1),
               "O5'": (0, 1), "O3'": (0, 1), "O4'": (0, 1),
               "O2'": (1, 1)}


def roles(base, atom):
    """(donor, acceptor) for one atom, or None when it has no H-bond role."""
    if atom in BASE_ROLES[base]:
        return BASE_ROLES[base][atom]
    return SUGAR_ROLES.get(atom)


def hbond_lines(grp, base_of):
    """<resname> <type> <donor> <acceptor>, one line per (residue, type).

    Derived from the .grp just built, so the two files cannot drift apart:
    a type that the .grp stops naming stops being classified here too.
    """
    out, emitted = [], {}
    for line in grp:
        type_name, resname, atom = line.split()
        role = roles(base_of(resname), atom)
        if role is None:
            continue
        key = (resname, type_name)
        if key in emitted:
            # Aliases (OP1/O1P) share a type and must share a role; anything
            # else means two different atoms collapsed onto one type name.
            if emitted[key] != role:
                raise SystemExit("type %s in %s would be both %s and %s"
                                 % (type_name, resname, emitted[key], role))
            continue
        emitted[key] = role
        out.append("%-4s %-6s %d %d" % (resname, type_name, role[0], role[1]))
    return out


def short(prefix, atom):
    """<=4 characters, apostrophes dropped and a doubled prime turned into a
    trailing 2 (H5'' -> H52), which is what keeps every name inside a PDB
    atom-name field. Uniqueness is not assumed from this: the self-check in
    main() re-reads the tables and refuses to write if any residue's charge
    stops matching AMBER's."""
    name = atom.replace("''", "2").replace("'", "")
    out = prefix + name
    if len(out) > 4:
        raise SystemExit("type name too long: %s (%s + %s)" % (out, prefix, atom))
    return out


def generate(xml_name, kind, sugar_prefix, base_prefix, spellings):
    root = ET.parse(os.path.join(DATA, "amber14", xml_name)).getroot()
    mass = {t.get("name"): float(t.get("mass")) for t in root.find("AtomTypes").findall("Type")}
    lj = {a.get("type"): (float(a.get("sigma")), float(a.get("epsilon")))
          for a in root.find("NonbondedForce").findall("Atom")}
    residues = {r.get("name"): r for r in root.find("Residues").findall("Residue")}

    nbi, grp, seen = [], [], {}

    def emit(type_name, atom_type, charge, resnames, atom_names):
        if type_name not in seen:
            sigma, eps = lj[atom_type]
            nbi.append("%-6s\t%.4f\t%.4f\t%.4f\t%.3f\t0.0\t0.0"
                       % (type_name, charge, sigma * RSTAR, eps / KCAL, mass[atom_type]))
            seen[type_name] = charge
        for rn in resnames:
            for an in atom_names:
                grp.append("%-6s %-4s %s" % (type_name, rn, an))

    for digit, (base, prefix) in enumerate(base_prefix.items(), start=1):
        digit_prefix = str(digit)
        # Every variant of this base: interior, 5', 3' and nucleoside. They
        # share the base and sugar chemistry; the terminal ones only add or
        # drop a phosphate, and a PDB spells them all with the plain name.
        variants = [v for v in (base, base + "5", base + "3", base + "N") if v in residues]
        emitted_here = set()
        for variant in variants:
            for a in residues[variant].findall("Atom"):
                name, atype, q = a.get("name"), a.get("type"), float(a.get("charge"))
                if name in emitted_here:
                    continue
                emitted_here.add(name)
                if name in PER_BASE_SUGAR:
                    tn = short(digit_prefix, name)
                elif name in BACKBONE:
                    tn = short(sugar_prefix, name)
                else:
                    tn = short(prefix, name)
                names = [name] + ALIASES.get(name, [])
                emit(tn, atype, q, spellings(base), names)

    # SELF-CHECK. A type name is 4 characters and several transformations
    # collapse into it, so uniqueness cannot be assumed -- it has to be
    # demonstrated. Rebuild each residue's total charge from the tables just
    # written and compare it with the template's: a collision shows up
    # immediately as a residue that no longer sums to AMBER's value (guanine
    # did, at -1.1774 e instead of -1.0000, before C1'/H1' were moved onto a
    # digit prefix).
    q_of = dict(seen)
    type_of = {}
    for line in grp:
        t, rn, an = line.split()
        type_of[(rn, an)] = t
    terminal_dev = []
    for base in base_prefix:
        for variant in (base, base + "5", base + "3", base + "N"):
            if variant not in residues:
                continue
            atoms = [(a.get("name"), float(a.get("charge")))
                     for a in residues[variant].findall("Atom")]
            amber = sum(q for _, q in atoms)
            mine = sum(q_of[type_of[(base, n)]] for n, _ in atoms if (base, n) in type_of)
            if variant == base:
                if abs(mine - amber) > 1e-4:
                    raise SystemExit(
                        "charge mismatch on %s: %.4f e against AMBER's %.4f -- two atoms "
                        "share a type name" % (variant, mine, amber))
            else:
                terminal_dev.append((variant, mine - amber))
    print("  %s : charges interieures exactes ; variantes terminales a %+.3f e "
          "(la .grp ne peut pas les distinguer, un PDB les nomme pareil)"
          % (kind.upper(), max(d for _, d in terminal_dev) if terminal_dev else 0.0))
    return nbi, grp


def write(path, header, lines):
    with open(path, "w") as fh:
        fh.write("\n".join(header) + "\n" + "\n".join(lines) + "\n")
    print("  %-46s %d lignes" % (os.path.relpath(path, REPO_ROOT), len(lines)))


def main():
    for xml_name, kind, sp, bp, spell in (
            ("DNA.OL15.xml", "dna", "d",
             {"DA": "a", "DC": "c", "DG": "g", "DT": "t"}, lambda b: [b]),
            ("RNA.OL3.xml", "rna", "r",
             {"A": "a", "C": "c", "G": "g", "U": "u"}, lambda b: [b, "R" + b])):
        nbi, grp = generate(xml_name, kind, sp, bp, spell)
        up = kind.upper()
        write(os.path.join(REPO_ROOT, "data/forcefield/amber.%s.nbi.ff" % kind),
              ["#Non-bonded parameters for %s, generated by" % up,
               "#scripts/generate_nucleic_nonbonded.py -- do not edit by hand.",
               "#Source: amber14/%s." % xml_name,
               "#",
               "#type\tcharge(e)\tradius(A)\tepsilon(kcal.mol-1)\tmass(g.mol-1)"
               "\ttransferIMP\tHydrophobicity",
               "#",
               "#epsilon is in kcal/mol, like amber.nbi.ff and unlike that file's own",
               "#header: ForceFieldReader passes the column through untouched, so the",
               "#unit is whatever amber.nbi.ff has always used (checked: its N/H/O/C",
               "#rows are AMBER's kcal/mol values exactly). radius is R* = rmin/2.",
               "#",
               "#transferIMP and Hydrophobicity are 0: no published nucleic",
               "#parameterisation exists for either, and a made-up number would be",
               "#worse than a declared absence.",
               "#"], nbi)
        write(os.path.join(REPO_ROOT, "data/reducerules/amber.%s.grp" % kind),
              ["#Atom-to-type mapping for %s, generated by" % up,
               "#scripts/generate_nucleic_nonbonded.py -- do not edit by hand.",
               "#Companion to data/forcefield/amber.%s.nbi.ff." % kind,
               "#",
               "#All-atom identity mapping (one type per atom, nothing merged), the",
               "#same shape as amber.grp -- which is what --bondedinteraction requires",
               "#of a translation table.",
               "#",
               "#Types are lowercase-prefixed and amber.grp's are uppercase-prefixed,",
               "#so this file can be concatenated with amber.grp for a protein-nucleic",
               "#complex. amber.dna.grp and amber.rna.grp DO collide with each other.",
               "#"], grp)
        # DA/DC/DG/DT, A/C/G/U and RA/RC/RG/RU all end on the base letter.
        base_of = lambda rn: rn[-1]
        write(os.path.join(REPO_ROOT, "data/reducerules/amber.%s.hbond" % kind),
              ["#Donor/acceptor classification for %s, generated by" % up,
               "#scripts/generate_nucleic_nonbonded.py -- do not edit by hand.",
               "#Keyed by amber.%s.grp's type names, so use it with that file;" % kind,
               "#it is the nucleic counterpart of amber.hbond and concatenates",
               "#with it for a protein-nucleic complex.",
               "#",
               "#Same classification as the protein table (HBPLUS / McDonald &",
               "#Thornton, Reduce/Probe) -- see ProteinDonorAcceptor.hbond's header.",
               "#Watson-Crick atoms plus the Hoogsteen (N7) and minor-groove (N3)",
               "#faces, the phosphate and ester oxygens, and the furanose O4'.",
               "#",
               "#Format: <resname> <type> <donor> <acceptor>, each 0 or 1.",
               "#",
               "#TWO LIMITS TO KNOW, both in the mechanism rather than this table:",
               "#",
               "#1. The columns are booleans, but real capacity is a COUNT: an",
               "#   amino nitrogen (N6, N4, N2) carries two donatable H, and a",
               "#   carbonyl oxygen (O6, O4, O2) has two lone pairs. BioSpring",
               "#   holds one partner per particle, so each is capped at one here.",
               "#   Canonical Watson-Crick pairing does not need more (its six",
               "#   atoms are each used once), but base triples do.",
               "#",
               "#2. O2' is the one atom listed as donor AND acceptor. In a folded",
               "#   RNA it routinely serves both roles AT ONCE, which one partner",
               "#   per particle cannot represent. Expect it to under-bond.",
               "#",
               "#The attractive term is also distance-only, with no donor-H...",
               "#acceptor angle. On 074's duplex that costs 2 of 48 Watson-Crick",
               "#bonds: two stacked same-strand groups sit closer to each other",
               "#(2.75 A) than to their real partners (2.90 and 2.94 A).",
               "#"], hbond_lines(grp, base_of))


if __name__ == "__main__":
    main()
