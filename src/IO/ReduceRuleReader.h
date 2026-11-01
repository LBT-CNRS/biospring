#ifndef _REDUCEREADER_H_
#define _REDUCEREADER_H_

#include <string>
#include <vector>

#include "IO/ReaderBase.h"
#include "reduce/ReduceRule.h"
#include "reduce/ReduceRuleContainer.hpp"

namespace biospring
{
namespace reduce
{

// Parses a .grp file: which atoms of which residue make up each grain, and what
// that grain is called.
//
//   <grain> <resname>[:NTER|:CTER] <atom> [<atom> ...]
//
// One line per grain. A grain built from several atoms sits at their centroid
// (a real coarse-grain bead); a grain built from one is a renaming, which is
// what the all-atom files here do -- amber.grp turns alanine's CA into ACA so
// that the .ff can give it AMBER's own per-atom charge and mass.
//
// A grain name may be repeated, for two reasons. Different residues can reduce
// to the same physical bead (the nucleic backbone beads of CGZaccharias.grp),
// and one grain can accept several spellings of its atom, one rule each
// ("XH1 ACE H1" and "XH1 ACE HH31"); a rule naming both spellings would instead
// require both atoms to be present. Rules are looked up by residue name and
// chain position, never by name alone.
//
// The residue name may carry a terminal qualifier, ":NTER" or ":CTER", the same
// one the .bi.ff carries on the same column and with the same meaning: the
// rules matching a residue's own position replace the unqualified ones. See
// topology/TerminalPosition.hpp for the rule and ReduceRule::forcefield_key for
// what the position changes, which is the parameters and not the grain's name.
class ReduceRuleReader : public ReaderBase
{
  protected:
    ReduceRuleContainer _rules;

  public:
    ReduceRuleReader() : ReaderBase() {}
    ReduceRuleReader(const std::string & path) : ReaderBase(path) {}
    ReduceRuleReader(const char * const path) : ReaderBase(path) {}

    const ReduceRuleContainer & rules(void) const { return _rules; }

    void read();

  protected:
    ReduceRule parse_rule(const std::string & line, size_t line_id) const;
};

} // namespace reduce
} // namespace biospring

#endif
