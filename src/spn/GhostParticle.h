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
// How a ghost's position is derived from its anchors. The two bonded
// families need genuinely different constructions -- they are not sharing
// ghosts, they share this placement machinery -- so each ghost carries its
// own mode.
enum class GhostPlacement : unsigned
{
    // DIHEDRAL rings: the ghost is the image of a real reference atom under
    // a rotation of delta about the B->C axis, at that atom's own distance
    // from the axis and its own position along it.
    AxisRotation = 0,
};

struct GhostParticleBinding
{
    unsigned ownIndex;
    unsigned anchorBIndex;
    unsigned anchorCIndex;
    // The real atom this ghost is the rotated image of. It also remains
    // the azimuthal origin it always was, so delta_deg keeps its meaning:
    // a ring is simply this atom replicated at M equispaced azimuths
    // about the B->C axis.
    unsigned anchorRefIndex;
    float r;
    float theta_deg;
    float delta_deg;
    // Cached cos/sin of delta_deg, for the rotation placement.
    float cos_delta;
    float sin_delta;
    // Index into SpringNetwork's per-axis accumulators: every ghost of a
    // ring shares one axis, and the closed-form axis reaction is computed
    // once per axis rather than once per ghost (it is only valid over
    // complete spring pairs -- see redistributeAxisReaction).
    unsigned axisIndex;
    GhostPlacement placement;
};

// The pure geometry behind a ghost particle. Both constructions here are
// closed form: a ghost is either the image of a real atom under a rotation
// about the axis, or a point ON the axis. Neither needs a placement
// Jacobian, and the general 3-anchor placement that did was removed --
// measured at ~59% of a bonded-only step for its ~9 3x3 matrix products
// per ghost, which is more than the classical treatments it competes with
// ever cost. Free of any Particle/SpringNetwork dependency so it can be
// tested in isolation (see test-GhostParticle.cpp).
class GhostParticle
{
  public:
    // ---- Rotation placement ------------------------------------------
    //
    // Places the ghost as the image of a REAL atom under a rotation of
    // delta about the B->C axis. The ghost then sits at exactly that
    // atom's distance from the axis and at its axial position, differing
    // only in azimuth -- which is what makes the whole construction
    // collapse (see doc/BondedForceFieldSprings.md):
    //
    //  * dX_ghost/dX_atom is EXACTLY the rotation matrix R(axis, delta),
    //    so the force owed to the real atom is just R(axis, -delta) . f
    //    -- one Rodrigues application instead of building the placement
    //    Jacobian's ~9 3x3 matrix products;
    //  * the forces owed to the two axis atoms carry no torque about the
    //    axis at all (they do zero virtual work under the only motion
    //    --rigidbody leaves free), and are fully determined by global
    //    force/torque balance -- so they are reconstructed in closed form
    //    rather than differentiated. Verified over 200 random geometries
    //    and ring sizes: matches the true Jacobian to 3e-7, i.e. to the
    //    finite-difference noise floor of the check itself.
    //
    // cos_delta/sin_delta are passed precomputed: delta is fixed for the
    // lifetime of a ghost.
    static Vector3f computePositionByRotation(const Vector3f & B, const Vector3f & C, const Vector3f & atom,
                                              float cos_delta, float sin_delta);

    // Redistributes a force acting on a ghost placed by
    // computePositionByRotation onto the real atom it images and onto the
    // two axis atoms. `f_atom_total`/`torque_total` accumulate the pair's
    // two ends before the axis reconstruction, which is why the axis part
    // is computed by redistributeAxisReaction below rather than here: the
    // reconstruction is only valid once BOTH ends of a spring have been
    // accounted for (it enforces the balance of the whole interaction).
    static Vector3f rotateForceToAtom(const Vector3f & B, const Vector3f & C, const Vector3f & f, float cos_delta,
                                      float sin_delta);

    // Closed-form reaction on the two axis atoms. Redistribution must
    // TRANSFER what acted on the ghosts, not cancel it, so it solves
    //     F_B + F_C     = (sum of ghost forces)  - (sum of real-atom forces)
    //     (C-B) x F_C   = (sum of ghost torques) - (sum of real-atom torques)
    // both taken about B, whose minimum-norm solution is
    // F_C = (T x a)/|a|^2 with a = C-B.
    //
    // Always solvable: rotating about the axis preserves the axial
    // component of a torque, so ghost and real-atom axial torques cancel
    // ghost by ghost and T never has a component along a. And exact, not
    // approximate: the true Jacobian's answer has no axial component
    // either (measured over 200 random geometries: 1e-7), because a ghost
    // built by rotation is rigidly carried by its atom's group and adds no
    // freedom of its own.
    //
    // Note the sums are over a whole axis. On a complete ring the ghost
    // sums vanish (every spring is internal to it) and this reduces to a
    // pure reaction -- but the ghost terms must be carried anyway, or an
    // unpaired force would be annihilated instead of transferred.
    static void redistributeAxisReaction(const Vector3f & B, const Vector3f & C, const Vector3f & sum_ghost_forces,
                                          const Vector3f & sum_ghost_torques_about_B, const Vector3f & sum_atom_forces,
                                          const Vector3f & sum_atom_torques_about_B, Vector3f & F_B, Vector3f & F_C);
};

} // namespace spn
} // namespace biospring

#endif // __GHOSTPARTICLE_H__
