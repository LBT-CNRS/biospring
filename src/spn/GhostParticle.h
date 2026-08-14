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
    // BEND: the ghost sits ON the axis, at distance r from B along B->C.
    // No azimuth, no reference atom (r/theta/delta's azimuthal term is
    // multiplied by sin(theta=0)). Models a 1-3 distance without letting
    // the real 1-2 bond stretch leak into it -- see
    // BondedForceFieldReader's own comment.
    AxialOffset = 1,
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
    // Cached from the three above; kept alongside rather than replacing
    // them because r/theta/delta are what the .nc stores and what
    // NetCDFWriter reads back out of the binding.
    GhostLocalOffset offset;
    // Cached cos/sin of delta_deg, for the rotation placement.
    float cos_delta;
    float sin_delta;
    // Index into SpringNetwork's per-axis accumulators: every ghost of a
    // ring shares one axis, and the closed-form axis reaction is computed
    // once per axis rather than once per ghost (it is only valid over
    // complete spring pairs -- see redistributeAxisReaction). Unused for
    // AxialOffset ghosts, whose redistribution balances on its own.
    unsigned axisIndex;
    GhostPlacement placement;
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

    // ---- Rotation placement: same physics, no Jacobian ---------------
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

    // ---- Axial placement (BEND) --------------------------------------
    //
    // X = B + r * u, u = (C-B)/|C-B|. Its Jacobian is closed form and
    // needs no frame either: with P = I - u u^T the projector across the
    // axis, dX/dC = (r/L) P and dX/dB = I - (r/L) P, so
    //     F_C = (r/L) * (f - u (u.f)),   F_B = f - F_C
    // which transfers the force exactly and conserves the torque about B
    // by construction ((C-B) x F_C = r u x f, the ghost's own torque).
    // Nothing accumulates across ghosts here: each is already balanced on
    // its own, unlike the rotation placement.
    static Vector3f computePositionAxial(const Vector3f & B, const Vector3f & C, float r);

    static void redistributeForceAxial(const Vector3f & B, const Vector3f & C, float r, const Vector3f & f,
                                        Vector3f & F_B, Vector3f & F_C);
};

} // namespace spn
} // namespace biospring

#endif // __GHOSTPARTICLE_H__
