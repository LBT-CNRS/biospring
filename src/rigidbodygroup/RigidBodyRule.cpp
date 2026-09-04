
#include "RigidBodyRule.h"

#include <iostream>

namespace biospring
{
namespace rigidbodygroup
{

void RigidBodyRule::print() const
{
    std::cout << "rigidbody " << _name << " resname " << _residue_name;
    std::cout << " atom (";

    for (const auto & name : _atomnames)
        std::cout << name << " ";

    std::cout << ")" << std::endl;
}

} // namespace rigidbodygroup
} // namespace biospring
