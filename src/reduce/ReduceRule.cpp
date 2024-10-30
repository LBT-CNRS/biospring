
#include "ReduceRule.h"

#include <iostream>
#include <string>

namespace biospring
{
namespace reduce
{

void ReduceRule::print() const
{
    std::cout << "group " << _name << " resname " << _residue_name;
    std::cout << " atom (";

    for (std::unordered_set<std::string>::const_iterator it = _atomnames.begin(); it != _atomnames.end(); ++it)
    {
        std::cout << *it << " ";
    }

    std::cout << " )" << std::endl;
}

} // namespace reduce
} // namespace biospring