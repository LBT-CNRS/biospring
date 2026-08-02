#ifndef __GHOSTPARTICLE_H__
#define __GHOSTPARTICLE_H__

#include "Particle.h"

namespace biospring
{
namespace spn
{

// A massless virtual site used to build dihedral ghost-spring torsion
// terms (see doc/BondedForceFieldSprings.md). Its position is always a
// fixed algebraic function of 3 real anchor particles (B, C, Ref) and 3
// calibrated parameters (r, theta, delta) -- never integrated on its own:
// no independent mass/velocity, always isStatic(), recomputed every step
// by updatePositionFromAnchors(). Deliberately NOT a "real" dynamical
// particle attached via a stiff rigid spring: that design (considered
// earlier) would add a genuine fast vibrational mode (its own mass on a
// finite-stiffness spring), forcing a smaller global timestep and, on the
// OpenCL backend, picking up the uniform per-particle velocity damping
// like any other integrated particle -- an artificial dissipation channel
// with no counterpart in the AMBER model being reproduced. A massless
// algebraic virtual site (matching GROMACS/OpenMM/NAMD virtual
// sites/lonepairs) has none of these failure modes: no inertia to
// perturb, no new oscillator for the integrator, and no thermostat can
// couple to something that is never itself integrated.
//
// Any force accumulated on this particle by ordinary spring force
// computation (e.g. a dihedral ghost spring between two GhostParticles)
// is redistributed onto its 3 anchors by redistributeForceToAnchors(),
// via the transpose of the placement Jacobian (virtual work principle:
// dE = -F_ghost . dX_ghost = -F_ghost . (J_B dB + J_C dC + J_Ref dRef), so
// F_B = J_B^T F_ghost, etc. -- verified numerically against the true
// energy gradient, not just the Jacobian itself, before implementing:
// using J directly instead of J^T silently breaks force conservation).
// Its own force is then reset to zero: it never reaches
// SpringNetwork::updateParticlePositions's normal per-dynamic-particle
// resetForce(), since it is static and excluded from that loop.
class GhostParticle : public Particle
{
  public:
    GhostParticle(Particle * anchorB, Particle * anchorC, Particle * anchorRef, float r, float theta_deg,
                  float delta_deg);

    Particle * getAnchorB() const { return _anchorB; }
    Particle * getAnchorC() const { return _anchorC; }
    Particle * getAnchorRef() const { return _anchorRef; }

    float getR() const { return _r; }
    float getTheta() const { return _theta_deg; }
    float getDelta() const { return _delta_deg; }

    // Recomputes this particle's position from its 3 anchors' CURRENT
    // positions. Called once per step, after the anchors themselves have
    // been integrated (see SpringNetwork::computeStep).
    void updatePositionFromAnchors();

    // Redistributes this particle's currently accumulated force (from
    // spring computation) onto its 3 anchors via addForce(), then resets
    // its own force to zero. Called once per step, after spring forces
    // have been computed but before the anchors are integrated.
    void redistributeForceToAnchors();

    // Pure placement formula, exposed for testing: places a virtual site
    // at distance r from B, angle theta_deg from the B->C axis direction,
    // azimuthal angle delta_deg from the reference direction defined by
    // Ref's perpendicular component relative to that axis.
    static Vector3f computePosition(const Vector3f & B, const Vector3f & C, const Vector3f & Ref, float r,
                                     float theta_deg, float delta_deg);

    // Redistributes a force `f` acting on the virtual site placed by
    // computePosition (same B, C, Ref, r, theta_deg, delta_deg) onto its 3
    // anchors, via the transpose of the placement Jacobian.
    static void redistributeForce(const Vector3f & B, const Vector3f & C, const Vector3f & Ref, float r,
                                   float theta_deg, float delta_deg, const Vector3f & f, Vector3f & F_B,
                                   Vector3f & F_C, Vector3f & F_Ref);

  private:
    Particle * _anchorB;
    Particle * _anchorC;
    Particle * _anchorRef;
    float _r;
    float _theta_deg;
    float _delta_deg;
};

} // namespace spn
} // namespace biospring

#endif // __GHOSTPARTICLE_H__
