#ifndef __RIGIDBODYGROUP_RIGIDBODYRULEREADER_H__
#define __RIGIDBODYGROUP_RIGIDBODYRULEREADER_H__

#include <string>

#include "IO/ReaderBase.h"
#include "rigidbodygroup/RigidBodyRule.h"
#include "rigidbodygroup/RigidBodyRuleContainer.hpp"

namespace biospring
{
namespace rigidbodygroup
{

// Parses a .rbody file: same columnar syntax as the coarse-grain reduce
// rules (name resname atom1 atom2 ...), but the resulting rules are consumed
// by RigidBodyBuilder to create springs, not to reduce particles into grains.
// Atom names may carry a "+"/"-" prefix (CHARMM-style: next/previous residue).
class RigidBodyRuleReader : public ReaderBase
{
  protected:
    RigidBodyRuleContainer _rules;

  public:
    RigidBodyRuleReader() : ReaderBase() {}
    RigidBodyRuleReader(const std::string & path) : ReaderBase(path) {}
    RigidBodyRuleReader(const char * const path) : ReaderBase(path) {}

    const RigidBodyRuleContainer & rules() const { return _rules; }

    void read();

  protected:
    RigidBodyRule parse_rule(const std::string & line, size_t line_id) const;
};

} // namespace rigidbodygroup
} // namespace biospring

#endif
