#ifndef __TOPOLOGY_TERMINALPOSITION_HPP__
#define __TOPOLOGY_TERMINALPOSITION_HPP__

#include <set>
#include <string>
#include <utility>

#include "utils/string.hpp"

namespace biospring
{
namespace topology
{

// Where in its chain a residue has to sit for a record to apply to it.
//
// AMBER types a free chain end from its own residue template (NALA/CALA, not
// ALA), and the retyping is neither cosmetic nor local. It reaches EVERY atom
// of the residue: of NALA's 12 atoms, 9 of the 10 it shares with ALA carry a
// different charge (N goes from -0.4157 to +0.1414 e as it becomes an ammonium
// N3, but so does CB, -0.1825 to -0.0597, three bonds away), and the same holds
// at the other end. Over the 20 residues and both ends, 628 shared atoms are
// retyped and 79 exist only at a terminus. What that buys is the net charge: an
// AMBER terminal residue is exactly +1 e at an N-terminus and -1 e at a
// C-terminus, where the interior parameters give -1.2719 e and +1.0000 e.
//
// It also changes which bonded terms exist. amber99sb assigns NO torsion at all
// to N3-CT-C-N, HP-CT-C-O or H1-CT-C-O2 while assigning one to the identical
// interior quadruplets -- so an interior record invents a term the force field
// does not have (measured on 071: a psi well 17.8 kJ.mol-1 deep and
// 20.3 kJ.mol-1.rad-1 of torque on a residue whose psi AMBER leaves free), and
// the terms AMBER does add there, the ammonium rotor and the carboxylate
// improper, are missing.
//
// So a residue name may carry a terminal qualifier, "ALA:NTER" or "ALA:CTER",
// and a residue is parameterised from the records matching its own position:
// the interior set is NOT also applied at a terminus. Two files use it, with
// the same syntax and the same semantics, which is why both read it from here:
//
//   .bi.ff  TORSION <resname>[:NTER|:CTER] ...   (BondedForceFieldReader)
//   .grp    <grain> <resname>[:NTER|:CTER] ...   (ReduceRuleReader)
//
// The residue name itself stays AMBER's own in both -- "ALA", never "NALA" --
// like every other data file here.
enum TerminalPosition
{
    POSITION_INTERIOR = 0, // no qualifier: a residue with neighbours on both sides
    POSITION_NTER,         // ":NTER", a free N-terminus (NH3+), AMBER's NXXX
    POSITION_CTER          // ":CTER", a free C-terminus (COO-), AMBER's CXXX
};

// The qualifier that spells `position`, without its colon. Empty for the
// interior, which carries no qualifier at all.
inline const char * terminal_qualifier(TerminalPosition position)
{
    switch (position)
    {
    case POSITION_NTER:
        return "NTER";
    case POSITION_CTER:
        return "CTER";
    default:
        return "";
    }
}

// Splits "ALA" or "ALA:NTER" into the residue name and the position.
//
// Returns false, leaving both outputs untouched, if the field is not a residue
// name with at most one known qualifier. Callers report that themselves: they
// know their own file name and line number, and the message has to name the
// format the user actually wrote.
inline bool parse_residue_field(const std::string & field, std::string & resname, TerminalPosition & position)
{
    const auto parts = utils::string::split(field, ":");
    if (parts.empty() || parts.size() > 2 || parts[0].empty())
        return false;

    TerminalPosition parsed = POSITION_INTERIOR;
    if (parts.size() == 2)
    {
        if (parts[1] == "NTER")
            parsed = POSITION_NTER;
        else if (parts[1] == "CTER")
            parsed = POSITION_CTER;
        else
            return false;
    }

    resname = parts[0];
    position = parsed;
    return true;
}

// The atoms a residue has ONLY at a free chain end: the ammonium's hydrogens
// and the carboxylate's second oxygen, under their PDB names.
//
// Being first in a chain is not the same as being a free N-terminus, and the
// difference is what these answer. A structure carved out of a larger one -- a
// domain, a binding pocket, a single residue lifted out for a test -- has a
// first residue whose nitrogen is still an amide with its one hydrogen, and a
// last residue whose carbonyl is still a carbonyl. Typing those from AMBER's
// NXXX/CXXX template would put a charge on them that the structure does not
// have, and suppress backbone torsions that are really there. Measured on the
// examples here: 056's lone methionine (no H1/H2/H3, no OXT, a fragment at both
// ends) would come out at +0.60 e instead of 0, and 4 of 077's chain ends with
// it -- 077 being a ribosome pocket carved out of 1VY4, where every cut is one
// of these.
//
// So a free terminus is recognised the way AMBER and OpenMM recognise it, from
// the residue's own atom list, and a chain end that carries none of these atoms
// keeps the interior parameters. PRO is covered: NPRO has H2 and H3 and no H1.
// Only OXT is listed for the C end, not the older OT1/OT2 spellings -- amber.grp
// declares no rule for those at all, so a file using them already has its atoms
// reported as claimed by no grain.
inline bool is_nterminal_only_atom(const std::string & name)
{
    return name == "H1" || name == "H2" || name == "H3";
}

inline bool is_cterminal_only_atom(const std::string & name) { return name == "OXT"; }

// Which records a residue is entitled to, from where it sits in its chain, what
// its own atoms say it is, and what the file actually declares.
//
// `has_nterminal_atom`/`has_cterminal_atom` are what the residue's atom names
// say -- see is_nterminal_only_atom. Both readers compute them over the same
// residue, so the charges and the torsions can never disagree about where a
// terminus is.
//
// `declared` is the set of (resname, position) pairs the file carries, and a
// position the file says nothing about falls back to the interior records.
// Without that fallback a file generated before the qualifier existed -- the
// DNA and RNA .bi.ff still are, and so is every .grp but amber.grp -- would
// silently lose every record at every chain end instead of keeping today's
// behaviour.
//
// A residue that is both a free N-terminus and a free C-terminus (a genuine
// free amino acid, ammonium and carboxylate both present) is treated as an
// N-terminus: AMBER has no template for one, so there is nothing to generate
// the real combination from. For the torsions the ammonium rotor is then right,
// and so are both suppressions, since HP-CT-C-O and H1-CT-C-O2 are equally
// absent from the force field; only the carboxylate's own improper is missed.
// For the charges its OXT is then claimed by no rule, which is reported per
// atom rather than passed over (Reducer's _check_unknown_particles).
inline TerminalPosition wanted_position(const std::string & resname, bool is_nterm, bool is_cterm,
                                        bool has_nterminal_atom, bool has_cterminal_atom,
                                        const std::set<std::pair<std::string, TerminalPosition>> & declared)
{
    // N before C, so a free amino acid takes the N-terminal set (see above).
    if (is_nterm && has_nterminal_atom && declared.count({resname, POSITION_NTER}) != 0)
        return POSITION_NTER;
    if (is_cterm && has_cterminal_atom && declared.count({resname, POSITION_CTER}) != 0)
        return POSITION_CTER;
    return POSITION_INTERIOR;
}

} // namespace topology
} // namespace biospring

#endif // __TOPOLOGY_TERMINALPOSITION_HPP__
