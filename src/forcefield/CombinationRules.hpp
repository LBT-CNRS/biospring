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
// The arithmetic mean on sigma, which is a SUM when the parameter stored is
// AMBER's R* = rmin/2, as it is in data/forcefield/amber*.ff. This is the rule
// the Lennard-Jones laws use.
inline float radius(float radius_i, float radius_j) { return radius_i + radius_j; }
} // namespace lorentz_berthelot

// Good-Hope Rules.
// https://en.wikipedia.org/wiki/Combining_rules#Good-Hope_rule
// Good-Hope Rules. The geometric mean applies to sigma, not to R*, so no law in
// this code uses it; it is kept named for the record.
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