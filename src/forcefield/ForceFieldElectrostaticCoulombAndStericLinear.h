#ifndef _FORCEFIELDELECTROSTATICCOULOMBANDSTERICLINEAR_H_
#define _FORCEFIELDELECTROSTATICCOULOMBANDSTERICLINEAR_H_

#include "ForceField.h"
#include "energy.hpp"

namespace biospring
{
namespace forcefield
{

class ForceFieldElectrostaticCoulombAndStericLinear : public ForceField
{
  public:
    ForceFieldElectrostaticCoulombAndStericLinear() : ForceField() {}
    virtual ~ForceFieldElectrostaticCoulombAndStericLinear() {}

    // Assignement operator.
    using ForceField::operator=;

    virtual float computeStericEnergy(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                      float distance) const override
    {
        return _stericscale * steric_energy_linear(radius_i, radius_j, distance);
    }

    virtual float computeStericForceModule(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                           float distance) const override
    {
        return _stericscale * steric_force_module_linear(radius_i, radius_j, distance);
    }
};

} // namespace forcefield
} // namespace biospring

#endif
