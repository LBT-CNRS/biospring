#ifndef _FORCEFIELDELECTROSTATICCOULOMBANDSTERICLENNARDJONES_8_6_LEWITT_H_
#define _FORCEFIELDELECTROSTATICCOULOMBANDSTERICLENNARDJONES_8_6_LEWITT_H_

#include "ForceField.h"
#include "energy.hpp"

namespace biospring
{
namespace forcefield
{

class ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Lewitt : public ForceField
{
  public:
    ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Lewitt() : ForceField() { _radiusrule = STERIC_RADIUS_RULE_LEWITT; }
    virtual ~ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Lewitt() {}

    // Assignement operator.
    using ForceField::operator=;

    virtual float computeStericEnergy(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                      float distance) const override
    {
        return _stericscale * steric_energy_lewitt(radius_i, radius_j, epsilon_i, epsilon_j, distance, _radiusrule);
    }

    virtual float computeStericForceModule(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                           float distance) const override
    {
        return _stericscale * steric_force_module_lewitt(radius_i, radius_j, epsilon_i, epsilon_j, distance, _radiusrule);
    }
};

} // namespace forcefield
} // namespace biospring

#endif
