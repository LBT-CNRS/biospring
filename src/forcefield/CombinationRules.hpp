#ifndef __FORCEFIELD_COMBINATION_RULES_HPP__
#define __FORCEFIELD_COMBINATION_RULES_HPP__

#include <cmath>

namespace biospring
{
namespace forcefield
{

namespace combination_rules
{

// Lorentz-Berthelot Rules.
// https://en.wikipedia.org/wiki/Combining_rules#Lorentz-Berthelot_rules
namespace lorentz_berthelot
{
inline float epsilon(float epsilon_i, float epsilon_j) { return sqrt(epsilon_i * epsilon_j); }
} // namespace lorentz_berthelot

// Good-Hope Rules.
// https://en.wikipedia.org/wiki/Combining_rules#Good-Hope_rule
namespace good_hope
{
inline float radius(float radius_i, float radius_j) { return sqrt(radius_i * radius_j); }
} // namespace good_hope

namespace zacharias
{
inline float epsilon(float epsilon_i, float epsilon_j) { return epsilon_i * epsilon_j; }
inline float radius(float radius_i, float radius_j) { return radius_i * radius_j; }

} // namespace zacharias

} // namespace combination_rules
} // namespace forcefield
} // namespace biospring

#endif // __FORCEFIELD_COMBINATION_RULES_HPP__