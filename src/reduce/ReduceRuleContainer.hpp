#ifndef __REDUCE_RULE_CONTAINER_HPP__
#define __REDUCE_RULE_CONTAINER_HPP__

#include <string>
#include <unordered_map>
#include <vector>

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
    // Dies if the rule name already exists.
    void append(const ReduceRule & rule)
    {
        if (contains(rule.getName()))
        {
            std::string msg = "ReduceRuleContainer::append: rule name already exists: " + rule.getName();
            throw std::runtime_error(msg);
        }
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

    // Returns the rules for a given residue name.
    // If the residue name is not found, returns an empty container.
    ReduceRuleContainer get_rules_for_residue(const std::string & residue_name) const
    {
        std::vector<ReduceRule> output;

        for (const ReduceRule & rule : _rules)
            if (rule.getResidueName() == residue_name)
                output.push_back(rule);

        return output;
    }
};

} // namespace reduce
} // namespace biospring

#endif // __REDUCE_RULE_CONTAINER_HPP__