#ifndef _GROUP_H_
#define _GROUP_H_

#include <string>
#include <vector>

#include "Particle.h"

using namespace std;

namespace biospring
{
namespace reduce
{
namespace legacy
{

class ParticleGroup : public biospring::spn::Particle
{
  public:
    ParticleGroup() : _group() {}

    void addParticle(const biospring::spn::Particle * p);
    void print();
    void reduce();

    size_t size(void) const { return _group.size(); }

    const biospring::spn::Particle * getParticleOfGroupFromIndex(unsigned i) { return _group[i]; }

    bool contains(const std::string & name) const;
    bool contains(const biospring::spn::Particle & p) const;

  private:
    vector<const biospring::spn::Particle *> _group;

    void _updatePositionWithMassCenter();
};

} // namespace legacy
} // namespace reduce
} // namespace biospring
#endif
