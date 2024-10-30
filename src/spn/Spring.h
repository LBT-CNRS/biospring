#ifndef __SPRING_H__
#define __SPRING_H__

#include "Particle.h"
#include "forcefield/ForceField.h"

namespace biospring
{
namespace spn
{

class Spring
{
  public:
    static const float DEFAULT_EQUILIBRIUM;
    static const float DEFAULT_STIFFNESS;

    Spring(Particle & p1, Particle & p2, float equilibrium, float stiffness)
        : _p1(p1), _p2(p2), _equilibrium(equilibrium), _stiffness(stiffness), _length(0.0), _energy(0.0)
    {
        computeLength();
    }

    Spring(Particle & p1, Particle & p2) : Spring(p1, p2, DEFAULT_EQUILIBRIUM, DEFAULT_STIFFNESS) {}

    Spring(Particle & p1, Particle & p2, float equilibrium) : Spring(p1, p2, equilibrium, DEFAULT_STIFFNESS) {}

    void setId(int id) { _id = id; }
    int getId() const { return _id; }

    const Particle & getParticle1() const { return _p1; }
    const Particle & getParticle2() const { return _p2; }
    Particle & getParticle1() { return _p1; }
    Particle & getParticle2() { return _p2; }

    void setStiffness(float stiffness) { _stiffness = stiffness; }
    float getStiffness() const { return _stiffness; }

    void setEquilibrium(float equilibrium) { _equilibrium = equilibrium; }
    float getEquilibrium() const { return _equilibrium; }

    float getLength() const { return _length; }

    float getEnergy() const { return _energy; }

    void computeEnergy(const biospring::forcefield::ForceField & ff);
    void computeLength();

    void applyForceToParticle(const biospring::forcefield::ForceField & ff);

  private:
    Particle & _p1;
    Particle & _p2;
    float _equilibrium;
    float _stiffness;
    float _length;
    float _energy;
    int _id;
};

} // namespace spn
} // namespace biospring

#endif // __SPRING_H__
