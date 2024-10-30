#include "ParticleGroup.h"
#include "Particle.h"
#include <iostream>

namespace biospring
{
namespace reduce
{
namespace legacy
{
bool ParticleGroup::contains(const std::string & name) const
{
    for (const biospring::spn::Particle * const p : _group)
        if (name == p->getName())
            return true;
    return false;
}

bool ParticleGroup::contains(const biospring::spn::Particle & p) const { return contains(p.getName()); }

void ParticleGroup::addParticle(const biospring::spn::Particle * p)
{
    _group.push_back(p);
    this->setChainName(p->getChainName());
}

void ParticleGroup::print()
{
    std::cout << "ParticleGroup " << getName() << " " << getX() << " " << getY() << " " << getZ() << " :";
    for (unsigned i = 0; i < _group.size(); i++)
    {
        std::cout << _group[i]->getExtid() << " ";
    }
    std::cout << std::endl;
}

void ParticleGroup::_updatePositionWithMassCenter()
{
    Vector3f bary = Vector3f();
    for (unsigned i = 0; i < _group.size(); i++)
    {
        bary = bary + _group[i]->getPosition();
    }

    if (_group.size() != 0)
        bary = bary * (1.0 / ((float)_group.size()));
    setPosition(bary);
}

void ParticleGroup::reduce() { _updatePositionWithMassCenter(); }

} // namespace legacy
} // namespace reduce
} // namespace biospring
