#ifndef _FORCEFIELDELECTROSTATICCOULOMBANDSTERICLENNARDJONES_8_6_ZACHARIAS_H_
#define _FORCEFIELDELECTROSTATICCOULOMBANDSTERICLENNARDJONES_8_6_ZACHARIAS_H_

#include "ForceField.h"
#include "energy.hpp"

namespace biospring
{
namespace forcefield
{

class ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Zacharias : public ForceField
{
  public:
    ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Zacharias() : ForceField() {}
    virtual ~ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Zacharias() {}

    // Assignement operator.
    using ForceField::operator=;

    virtual float computeStericEnergy(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                      float distance) const override
    {
        return _stericscale * steric_energy_zacharias(radius_i, radius_j, epsilon_i, epsilon_j, distance);
    }

    virtual float computeStericForceModule(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                           float distance) const override
    {
        return _stericscale * steric_force_module_zacharias(radius_i, radius_j, epsilon_i, epsilon_j, distance);
    }
};

} // namespace forcefield
} // namespace biospring

#endif
