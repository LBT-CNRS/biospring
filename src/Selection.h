#ifndef __SELECTION_H__
#define __SELECTION_H__

#include "Particle.h"
#include "Vector3f.h"
#include "measure.hpp"

#include <vector>

class Selection
{
  public:
    Selection(const std::string & name) : _name(name) {}

    std::string getName() const { return _name; }

    Vector3f getBarycentre() const { return biospring::measure::centroid(_particles); }

    void addParticle(biospring::spn::Particle * p) { _particles.push_back(p); }

    void addForce(const Vector3f & force)
    {
        Vector3f force_ = force / float(_particles.size());
        for (auto p : _particles)
            p->addForce(force_);
    }

  protected:
    std::string _name;
    std::vector<biospring::spn::Particle *> _particles;
};

#endif // __SELECTION_H__
