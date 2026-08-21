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

    void setId(unsigned id) { _id = id; }
    unsigned getId() const { return _id; }

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

    // Only meaningful for a dihedral ghost-ghost spring: this spring's
    // share of its axis's exact dihedral-energy correction (ring
    // construction artifact minus AMBER's own real DC -- see
    // topology::Spring::_dc_offset and scripts/generate_bonded_forcefield.py's
    // calibrate_ring). Zero for every other spring. Never affects forces (a
    // constant has zero gradient) -- SpringNetwork::computeDihedralForces
    // subtracts it from this spring's own energy when accumulating the total.
    void setDcOffset(float dcOffset) { _dcOffset = dcOffset; }
    float getDcOffset() const { return _dcOffset; }

    void computeEnergy(const biospring::forcefield::ForceField & ff);
    void computeLength();

    // Computes and returns the force applied to particle 1. The opposite force
    // must be applied to particle 2. This split makes parallel force evaluation
    // possible without concurrently modifying particles.
    //
    // `ignoreDynamicState`: bypasses the "skip if both endpoints are
    // non-dynamic" early exit (see the .cpp). Needed for dihedral ghost-ghost
    // springs (SpringNetwork::computeDihedralForces): both endpoints are
    // always static/massless virtual sites (spn::GhostParticle) by design,
    // so that guard's original intent -- skip pointless work between two
    // real, frozen atoms -- does not apply: the force/energy still matters
    // here and gets redistributed onto the ghosts' real, dynamic anchors
    // afterward (SpringNetwork::redistributeGhostForces).
    Vector3f computeForce(const biospring::forcefield::ForceField & ff, bool ignoreDynamicState = false);

    void applyForceToParticle(const biospring::forcefield::ForceField & ff);

  private:
    Particle & _p1;
    Particle & _p2;
    float _equilibrium;
    float _stiffness;
    float _length;
    float _energy;
    float _dcOffset = 0.0f;
    unsigned _id;
};

} // namespace spn
} // namespace biospring

#endif // __SPRING_H__
