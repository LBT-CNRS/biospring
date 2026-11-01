#ifndef __REDUCE_RULE_CONTAINER_HPP__
#define __REDUCE_RULE_CONTAINER_HPP__

#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// ReduceRule used to arrive here through reduce/Reduce.h, which was legacy and
// is gone. The dependency is real, so it is declared.
#include "ReduceRule.h"

namespace biospring
{
namespace reduce
{

class ReduceRuleContainer
{
  public:
    using value_type = ReduceRule;
    using size_type = std::size_t;
    using reference = value_type &;
    using const_reference = const value_type &;
    using pointer = value_type *;
    using iterator = std::vector<value_type>::iterator;
    using const_iterator = std::vector<value_type>::const_iterator;

  private:
    std::vector<ReduceRule> _rules;
    std::unordered_map<std::string, size_t> _name_to_position; // maps rule name to position in `_rules`

  public:
    ReduceRuleContainer() {}
    ReduceRuleContainer(const std::vector<ReduceRule> & rules) : _rules(rules)
    {
        for (size_t i = 0; i < _rules.size(); i++)
            _name_to_position[_rules[i].getName()] = i;
    }

    // =============================================================================
    // Capacity.
    // =============================================================================

    size_type size() const { return _rules.size(); }
    bool empty() const { return _rules.empty(); }

    // =============================================================================
    // Modifiers.
    // =============================================================================

    // Adds a rule to the container.
    // Rule names are not required to be globally unique: the same grain name
    // can legitimately be reused across different residues (e.g. the nucleic
    // acid backbone beads in CGZaccharias.grp, chemically identical for
    // every base). Rules are selected by residue name (see
    // get_rules_for_residue), never by looking up a name across residues, so
    // this is safe.
    void append(const ReduceRule & rule)
    {
        _rules.push_back(rule);
        _name_to_position[rule.getName()] = _rules.size() - 1;
    }

    // =============================================================================
    // Access to elements.
    // =============================================================================

    reference at(size_type pos) { return _rules.at(pos); }
    const_reference at(size_type pos) const { return _rules.at(pos); }

    reference operator[](size_type pos) { return _rules[pos]; }
    const_reference operator[](size_type pos) const { return _rules[pos]; }

    reference operator[](const std::string & name) { return _rules[_name_to_position.at(name)]; }
    const_reference operator[](const std::string & name) const { return _rules[_name_to_position.at(name)]; }

    iterator begin() { return _rules.begin(); }
    const_iterator begin() const { return _rules.begin(); }
    const_iterator cbegin() const { return _rules.cbegin(); }

    iterator end() { return _rules.end(); }
    const_iterator end() const { return _rules.end(); }
    const_iterator cend() const { return _rules.cend(); }

    // =============================================================================
    // Lookup.
    // =============================================================================

    // Returns true if the rule name exists.
    bool contains(const std::string & name) const { return _name_to_position.find(name) != _name_to_position.end(); }

    // Returns the rules for a given residue name, whatever chain position they
    // are qualified for.
    // If the residue name is not found, returns an empty container.
    //
    // This is the right overload for the reverse lookup that resolves an
    // original atom name onto a renamed particle (RigidBodyBuilder::_translate,
    // BondedForceFieldReader::_translate): a positional rule keeps the grain
    // name of the interior rule it replaces, so the answer does not depend on
    // the position, and the caller does not have to know one.
    ReduceRuleContainer get_rules_for_residue(const std::string & residue_name) const
    {
        std::vector<ReduceRule> output;

        for (const ReduceRule & rule : _rules)
            if (rule.getResidueName() == residue_name)
                output.push_back(rule);

        return output;
    }

    // Returns the rules that parameterise a residue sitting at `position`.
    //
    // A residue takes the rules qualified for its own position INSTEAD of the
    // unqualified ones -- an N-terminal alanine is AMBER's NALA, not ALA with
    // three hydrogens bolted on, and every one of its atoms is retyped (see
    // topology/TerminalPosition.hpp). So pass POSITION_INTERIOR, or the result
    // of topology::wanted_position, which applies the fallback that keeps a
    // file declaring no qualifier at all behaving exactly as before.
    ReduceRuleContainer get_rules_for_residue(const std::string & residue_name,
                                              topology::TerminalPosition position) const
    {
        std::vector<ReduceRule> output;

        for (const ReduceRule & rule : _rules)
            if (rule.getResidueName() == residue_name && rule.position() == position)
                output.push_back(rule);

        return output;
    }

    // The (residue, position) pairs the file carries, which is what
    // topology::wanted_position needs to fall back on the interior rules for a
    // position the file says nothing about.
    std::set<std::pair<std::string, topology::TerminalPosition>> declared_positions() const
    {
        std::set<std::pair<std::string, topology::TerminalPosition>> declared;

        for (const ReduceRule & rule : _rules)
            declared.insert({rule.getResidueName(), rule.position()});

        return declared;
    }
};

} // namespace reduce
} // namespace biospring

#endif // __REDUCE_RULE_CONTAINER_HPP__