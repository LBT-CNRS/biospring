#ifndef __RIGIDBODYGROUP_RIGIDBODYRULECONTAINER_HPP__
#define __RIGIDBODYGROUP_RIGIDBODYRULECONTAINER_HPP__

#include <string>
#include <unordered_map>
#include <vector>

#include "RigidBodyRule.h"

namespace biospring
{
namespace rigidbodygroup
{

class RigidBodyRuleContainer
{
  public:
    using value_type = RigidBodyRule;
    using size_type = std::size_t;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator = std::vector<value_type>::iterator;
    using const_iterator = std::vector<value_type>::const_iterator;

  private:
    std::vector<RigidBodyRule> _rules;

  public:
    RigidBodyRuleContainer() {}
    RigidBodyRuleContainer(const std::vector<RigidBodyRule> & rules) : _rules(rules) {}

    size_type size() const { return _rules.size(); }
    bool empty() const { return _rules.empty(); }

    // Rule names are not required to be globally unique (the same pivot name,
    // e.g. "_CA", legitimately repeats across residue types). Rules are
    // always selected by residue name (see get_rules_for_residue), never by
    // name alone.
    void append(const RigidBodyRule & rule) { _rules.push_back(rule); }

    reference operator[](size_type pos) { return _rules[pos]; }
    const_reference operator[](size_type pos) const { return _rules[pos]; }

    iterator begin() { return _rules.begin(); }
    const_iterator begin() const { return _rules.begin(); }
    iterator end() { return _rules.end(); }
    const_iterator end() const { return _rules.end(); }

    // Returns the rules for a given residue name, or an empty container if none match.
    RigidBodyRuleContainer get_rules_for_residue(const std::string & residue_name) const
    {
        std::vector<RigidBodyRule> output;
        for (const RigidBodyRule & rule : _rules)
            if (rule.getResidueName() == residue_name)
                output.push_back(rule);
        return output;
    }
};

} // namespace rigidbodygroup
} // namespace biospring

#endif
