#ifndef __INSERTION_VECTOR_H__
#define __INSERTION_VECTOR_H__

#include "Particle.h"
#include "Vector3f.h"

namespace biospring
{
namespace spn
{
class SpringNetwork;
} // namespace spn
} // namespace biospring

class InsertionVector
{
  public:
    InsertionVector(const biospring::spn::SpringNetwork & sp, biospring::spn::Particle & p1,
                    biospring::spn::Particle & p2)
        : _sp(sp), _p1(p1), _p2(p2), _vector(Vector3f(0.0, 0.0, 0.0)), _angle(0.0)
    {
    }

    const biospring::spn::Particle & getParticle(unsigned i) const;  // changed order of const
    biospring::spn::Particle & getParticle(unsigned i); // Non const usefull to modify the particles!


    Vector3f getVector() const { return _vector; }
    float getAngle() const { return _angle; }
    float getRollAngle() const { return _rollAngle; }
    float getInsertionDepth() const;

    void computeVector();
    void computeAngle();
    void computeRollAngle();

  protected:
    const biospring::spn::SpringNetwork & _sp;
    biospring::spn::Particle & _p1;
    biospring::spn::Particle & _p2;
    Vector3f _vector;
    float _angle;
    float _rollAngle;
};

#endif // __INSERTION_VECTOR_H__
