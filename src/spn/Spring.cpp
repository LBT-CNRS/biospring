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
    if (_p1.isRigid() && _p2.isRigid()) return;
    
    if (_p1.isDynamic() or _p2.isDynamic())
    {
        Vector3f diff = _p2.getPosition() - _p1.getPosition();
        diff.normalize();
        computeLength();
        computeEnergy(ff);

        float forcemodule = 0;
        Vector3f force = Vector3f();
        forcemodule = ff.computeSpringForceModule(_length, _stiffness, _equilibrium);
        force = diff * forcemodule;

        _p1.addForce(force);
        _p2.addForce(-force);
    }
}

void Spring::computeEnergy(const biospring::forcefield::ForceField & ff)
{
    _energy = ff.computeSpringEnergy(_length, _stiffness, _equilibrium);
}
void Spring::computeLength() { _length = Particle::distance(_p1, _p2); }

} // namespace spn
} // namespace biospring
