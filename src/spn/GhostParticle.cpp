#include "GhostParticle.h"

#include <cmath>

namespace biospring
{
namespace spn
{

// ---- Rotation placement ------------------------------------------------
//
// Rodrigues about the axis, applied to the atom's radial component only:
// the axial component is what makes the ghost share the atom's position
// along the axis, and rotating about that axis leaves it untouched.

namespace
{
inline Vector3f rodrigues(const Vector3f & v, const Vector3f & u, float c, float s)
{
    return v * c + (u ^ v) * s + u * (u.dot(v) * (1.0f - c));
}
} // namespace

Vector3f GhostParticle::computePositionByRotation(const Vector3f & B, const Vector3f & C, const Vector3f & atom,
                                                  float cos_delta, float sin_delta)
{
    Vector3f d = C - B;
    const Vector3f u = d / d.norm();

    const Vector3f w = atom - B;
    const Vector3f axial = u * u.dot(w);
    const Vector3f radial = w - axial;

    return B + axial + rodrigues(radial, u, cos_delta, sin_delta);
}

Vector3f GhostParticle::rotateForceToAtom(const Vector3f & B, const Vector3f & C, const Vector3f & f,
                                          float cos_delta, float sin_delta)
{
    Vector3f d = C - B;
    const Vector3f u = d / d.norm();
    // J^T with J = R(u, delta) is R(u, -delta): same cosine, opposite sine.
    return rodrigues(f, u, cos_delta, -sin_delta);
}

void GhostParticle::redistributeAxisReaction(const Vector3f & B, const Vector3f & C,
                                             const Vector3f & sum_ghost_forces,
                                             const Vector3f & sum_ghost_torques_about_B,
                                             const Vector3f & sum_atom_forces,
                                             const Vector3f & sum_atom_torques_about_B, Vector3f & F_B, Vector3f & F_C)
{
    const Vector3f a = C - B;
    // (C-B) x F_C = T  =>  minimum-norm F_C = (T x a)/|a|^2. The component
    // this drops is the one along the axis, which the true Jacobian does
    // not produce either.
    const Vector3f T = sum_ghost_torques_about_B - sum_atom_torques_about_B;
    F_C = (T ^ a) / a.dot(a);
    F_B = (sum_ghost_forces - sum_atom_forces) - F_C;
}

} // namespace spn
} // namespace biospring
