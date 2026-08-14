#ifndef __GHOSTPARTICLE_H__
#define __GHOSTPARTICLE_H__

#include "Vector3f.h"

namespace biospring
{
namespace spn
{

// Binds a ghost (massless virtual-site) particle -- already present in
// SpringNetwork's particle list as an ordinary, isStatic()=true, mass=0
// Particle (see doc/BondedForceFieldSprings.md for why a massless
// algebraic virtual site was chosen over a real dynamical particle on a
// rigid attachment spring: the latter would add a genuine fast
// vibrational mode and, on the OpenCL backend, pick up the uniform
// per-particle velocity damping like any other integrated particle) --
// to its 3 real anchor particles, by INDEX rather than by pointer or
// reference: SpringNetwork::_particles is a std::vector<Particle> that
// can reallocate as more particles are added, which would silently
// invalidate a cached pointer/reference but never invalidates an index.
// This is why GhostParticle does NOT derive from Particle and does not
// store anchor pointers: it is plain data (see GhostParticleBinding)
// plus the two static placement/force-redistribution formulas below,
// used by SpringNetwork::addGhostParticle/redistributeGhostForces/
// updateGhostPositions to drive the actual Particle entries by index.
// A ghost's position in its own local frame: X = B + a*u + b*v + c*w.
// These three depend only on (r, theta_deg, delta_deg), which never change
// once the ghost is built, so they are computed once instead of costing
// four trigonometric calls per ghost per step (13964 ghosts on
// example/072 alone, evaluated twice per step -- placement and force
// redistribution).
struct GhostLocalOffset
{
    float a;
    float b;
    float c;
};

struct GhostParticleBinding
{
    unsigned ownIndex;
    unsigned anchorBIndex;
    unsigned anchorCIndex;
    unsigned anchorRefIndex;
    float r;
    float theta_deg;
    float delta_deg;
    // Cached from the three above; kept alongside rather than replacing
    // them because r/theta/delta are what the .nc stores and what
    // NetCDFWriter reads back out of the binding.
    GhostLocalOffset offset;
};

// The pure geometry/algebra behind a ghost particle: a 3-point ("NeRF"
// style) virtual-site placement and its Jacobian-transpose force
// redistribution (virtual work principle: dE = -F_ghost . dX_ghost =
// -F_ghost . (J_B dB + J_C dC + J_Ref dRef), so F_B = J_B^T F_ghost, etc.
// -- verified numerically against the true energy gradient before
// implementing: using J directly instead of J^T silently breaks force
// conservation). Free of any Particle/SpringNetwork dependency so it can
// be tested in isolation (see test-GhostParticle.cpp).
class GhostParticle
{
  public:
    // Converts the (r, theta_deg, delta_deg) description of a ghost into
    // its local-frame coordinates. Pure trigonometry, no anchor involved.
    static GhostLocalOffset localOffset(float r, float theta_deg, float delta_deg);

    // Places a virtual site at distance r from B, angle theta_deg from
    // the B->C axis direction, azimuthal angle delta_deg from the
    // reference direction defined by Ref's perpendicular component
    // relative to that axis.
    static Vector3f computePosition(const Vector3f & B, const Vector3f & C, const Vector3f & Ref, float r,
                                     float theta_deg, float delta_deg);

    // Same placement, with the local-frame coordinates already known. This
    // is the form the simulation loop uses: placement needs only the frame
    // (u, v, w), never the placement Jacobian, so it must not pay for the
    // ~9 3x3 matrix products that building the Jacobian costs.
    static Vector3f computePosition(const Vector3f & B, const Vector3f & C, const Vector3f & Ref,
                                     const GhostLocalOffset & offset);

    // Redistributes a force `f` acting on the virtual site placed by
    // computePosition (same B, C, Ref, r, theta_deg, delta_deg) onto its 3
    // anchors, via the transpose of the placement Jacobian.
    static void redistributeForce(const Vector3f & B, const Vector3f & C, const Vector3f & Ref, float r,
                                   float theta_deg, float delta_deg, const Vector3f & f, Vector3f & F_B,
                                   Vector3f & F_C, Vector3f & F_Ref);

    // Same redistribution, with the local-frame coordinates already known.
    // Unlike placement this genuinely needs the Jacobian, so only the
    // trigonometry is saved here.
    static void redistributeForce(const Vector3f & B, const Vector3f & C, const Vector3f & Ref,
                                   const GhostLocalOffset & offset, const Vector3f & f, Vector3f & F_B,
                                   Vector3f & F_C, Vector3f & F_Ref);
};

} // namespace spn
} // namespace biospring

#endif // __GHOSTPARTICLE_H__
