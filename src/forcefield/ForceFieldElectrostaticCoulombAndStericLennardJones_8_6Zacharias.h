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
    // This law's radius rule is not separable from the law: it pairs a product
    // on the radius with a product on epsilon, and any other combination makes
    // the potential incoherent. _radiusrule is set so that getRadiusRule()
    // reports the truth, but the law does not read it and steric.radiusrule
    // cannot change it.
    ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Zacharias() : ForceField()
    {
        _radiusrule = STERIC_RADIUS_RULE_ZACHARIAS;
    }
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
