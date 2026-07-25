#include "Spring.h"

#include "Vector3f.h"

namespace biospring
{
namespace spn
{

const float Spring::DEFAULT_EQUILIBRIUM = 1.0;
const float Spring::DEFAULT_STIFFNESS = 1.0;

void Spring::applyForceToParticle(const biospring::forcefield::ForceField & ff)
{
    const Vector3f force = computeForce(ff);
    _p1.addForce(force);
    _p2.addForce(-force);
}

Vector3f Spring::computeForce(const biospring::forcefield::ForceField & ff)
{
    _energy = 0.0f;

    if ((_p1.isRigid() && _p2.isRigid()) || (!_p1.isDynamic() && !_p2.isDynamic()))
        return {};

    const Vector3f displacement = _p2.getPosition() - _p1.getPosition();
    _length = displacement.norm();
    computeEnergy(ff);

    Vector3f direction = displacement;
    direction.normalize();
    return direction * ff.computeSpringForceModule(_length, _stiffness, _equilibrium);
}

void Spring::computeEnergy(const biospring::forcefield::ForceField & ff)
{
    _energy = ff.computeSpringEnergy(_length, _stiffness, _equilibrium);
}
void Spring::computeLength() { _length = Particle::distance(_p1, _p2); }

} // namespace spn
} // namespace biospring
