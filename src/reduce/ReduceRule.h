#ifndef _REDUCERULE_H_
#define _REDUCERULE_H_

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "topology/TerminalPosition.hpp"

class Reduce;

namespace biospring
{
namespace reduce
{

class ReduceRule
{
  public:
    ReduceRule(const std::string & name, const std::string & resname,
               topology::TerminalPosition position = topology::POSITION_INTERIOR)
        : _name(name), _residue_name(resname), _position(position), _atomnames(), _isallatom(false)
    {
    }

    ReduceRule() : ReduceRule("", "") {}

    void setName(const std::string & name) { _name = name; }
    const std::string & getName() const { return _name; }
    const std::string & name() const { return _name; }

    void setResidueName(const std::string & resname) { _residue_name = resname; };
    const std::string & getResidueName() const { return _residue_name; }
    const std::string & residue_name() const { return _residue_name; }

    // Where in its chain the residue has to sit for this rule to apply, from
    // the ":NTER"/":CTER" qualifier on its residue name. A residue takes the
    // rules matching its own position INSTEAD of the unqualified ones -- see
    // topology/TerminalPosition.hpp and ReduceRuleContainer::get_rules_for_residue.
    void setPosition(topology::TerminalPosition position) { _position = position; }
    topology::TerminalPosition position() const { return _position; }

    // The force-field type this rule's grain is parameterised from.
    //
    // A positional rule keeps the grain NAME of the interior rule it replaces
    // -- a terminal alanine's nitrogen is still "AN" -- because that name is
    // the particle's identity for every file keyed on it: the .hbond donor and
    // acceptor tables, the .psite peptide sites, and the reverse lookup that
    // resolves a .rbody or .bi.ff atom name onto a renamed particle
    // (RigidBodyBuilder::_translate). Renaming the terminal grains would drop
    // every terminal residue out of all four, silently.
    //
    // What the position changes is the PARAMETERS, so it is the force-field key
    // that carries the qualifier: "AN:NTER" holds AMBER's ammonium charge
    // (+0.1414 e) where "AN" holds the amide's (-0.4157 e). Same spelling as
    // the .grp's own qualifier, so amber.nbi.ff reads as a variant table.
    std::string forcefield_key() const
    {
        const std::string qualifier = topology::terminal_qualifier(_position);
        return qualifier.empty() ? _name : _name + ":" + qualifier;
    }

    const std::unordered_set<std::string> & getAtomNames() const { return _atomnames; }

    // The atom name this rule was DECLARED with first, which is the canonical
    // spelling of what it renames. _atomnames is an unordered_set, so it can
    // say which spellings are accepted but not which came first -- and a type
    // may accept several: amber.grp has both "XH1 ACE H1" and "XH1 ACE HH31",
    // and "ZC NME C" alongside "ZC NME CH3". A parameter table is written
    // against one of them, so recording whichever spelling a particular
    // structure happened to use would match that table only by luck.
    const std::string & getCanonicalAtomName() const { return _canonicalatom; }

    void setAllAtom(bool allatom) { _isallatom = allatom; }
    bool isAllAtom() const { return _isallatom; }

    size_t getNumberOfAtoms(void) const { return _atomnames.size(); }
    size_t number_of_atoms(void) const { return _atomnames.size(); }

    void addAtom(const std::string & atomname)
    {
        if (_canonicalatom.empty())
            _canonicalatom = atomname;
        _atomnames.insert(atomname);
    };

    // Returns true if the rule contains an atom with the given name.
    bool hasAtomNamed(const std::string name) const { return _atomnames.find(name) != _atomnames.end(); }

    void print() const;

  private:
    std::string _name;
    std::string _canonicalatom;
    std::string _residue_name;
    topology::TerminalPosition _position;
    std::unordered_set<std::string> _atomnames;

    bool _isallatom;
};

} // namespace reduce
} // namespace biospring

#endif
