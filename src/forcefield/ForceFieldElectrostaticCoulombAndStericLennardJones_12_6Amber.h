#ifndef _FORCEFIELDELECTROSTATICCOULOMBANDSTERICLENNARDJONES_12_6_AMBER_H_
#define _FORCEFIELDELECTROSTATICCOULOMBANDSTERICLENNARDJONES_12_6_AMBER_H_

#include "ForceField.h"
#include "energy.hpp"

namespace biospring
{
namespace forcefield
{

// Real AMBER van der Waals form: a classic 12-6 Lennard-Jones potential
// (see steric_energy_amber/steric_force_module_amber in energy/steric.hpp,
// which do use exponents 12 and 6). Distinct from the genuinely-8-6
// ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Lewitt/Zacharias
// classes (separate files) -- this class used to be named "..._8_6Amber"
// (and the user-facing steric.mode value "lennard-jones-8-6Amber"), which
// was misleading about its own exponent.
class ForceFieldElectrostaticCoulombAndStericLennardJones_12_6Amber : public ForceField
{
  public:
    ForceFieldElectrostaticCoulombAndStericLennardJones_12_6Amber() : ForceField() {}
    virtual ~ForceFieldElectrostaticCoulombAndStericLennardJones_12_6Amber() {}

    // Assignement operator.
    using ForceField::operator=;

    virtual float computeStericEnergy(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                      float distance) const override
    {
        return _stericscale * steric_energy_amber(radius_i, radius_j, epsilon_i, epsilon_j, distance);
    }

    virtual float computeStericForceModule(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                           float distance) const override
    {
        return _stericscale * steric_force_module_amber(radius_i, radius_j, epsilon_i, epsilon_j, distance);
    }
};

} // namespace forcefield
} // namespace biospring

#endif
