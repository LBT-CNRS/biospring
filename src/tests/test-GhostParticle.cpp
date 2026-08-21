//
// Tests for spn::GhostParticle -- the two closed-form virtual-site
// placements and their force redistribution (see GhostParticle.h and
// doc/BondedForceFieldSprings.md for the derivation).
//

#include <cmath>

#include <gtest/gtest.h>

#include "spn/GhostParticle.h"

using biospring::spn::GhostParticle;


// =====================================================================================
// Rotation placement: the ghost is the image of a real atom under a
// rotation about the axis, so it must sit at exactly that atom's distance
// from the axis and at its axial position.
TEST(TestGhostParticle, computePositionByRotation_preserves_radius_and_axial_position)
{
    Vector3f B(0.3f, -0.7f, 0.2f);
    Vector3f C(1.9f, 0.4f, -0.5f);
    Vector3f atom(-0.4f, 1.1f, 0.8f);
    const float delta = 137.0f * static_cast<float>(M_PI) / 180.0f;

    Vector3f G = GhostParticle::computePositionByRotation(B, C, atom, std::cos(delta), std::sin(delta));

    Vector3f d = C - B;
    Vector3f u = d / d.norm();
    auto axial = [&](const Vector3f & p) { return u.dot(p - B); };
    auto radius = [&](const Vector3f & p) { return ((p - B) - u * u.dot(p - B)).norm(); };

    EXPECT_NEAR(axial(G), axial(atom), 1e-5f);
    EXPECT_NEAR(radius(G), radius(atom), 1e-5f);
    // and it is genuinely elsewhere: a non-zero azimuth really moved it
    EXPECT_GT((G - atom).norm(), 1e-2f);
}

// =====================================================================================
// The point of the rotation placement: the force owed to the real atom is
// R(axis, -delta) applied to the ghost force, and that must equal what the
// true placement Jacobian gives. Checked against a numerical J^T.
TEST(TestGhostParticle, rotateForceToAtom_matches_numerical_jacobian_transpose)
{
    Vector3f B(0.1f, 0.2f, -0.3f);
    Vector3f C(1.7f, -0.4f, 0.6f);
    Vector3f atom(0.5f, 1.3f, -0.2f);
    const float delta = 71.0f * static_cast<float>(M_PI) / 180.0f;
    const float cd = std::cos(delta), sd = std::sin(delta);
    Vector3f force(2.3f, -1.1f, 0.7f);

    // numerical J = dX_ghost/dX_atom, column by column
    const float h = 1e-4f;
    Vector3f col[3];
    for (int k = 0; k < 3; ++k)
    {
        Vector3f e(k == 0 ? 1.0f : 0.0f, k == 1 ? 1.0f : 0.0f, k == 2 ? 1.0f : 0.0f);
        Vector3f plus = GhostParticle::computePositionByRotation(B, C, atom + e * h, cd, sd);
        Vector3f minus = GhostParticle::computePositionByRotation(B, C, atom - e * h, cd, sd);
        col[k] = (plus - minus) / (2.0f * h);
    }
    // (J^T . force)_k = column_k . force
    Vector3f expected(col[0].dot(force), col[1].dot(force), col[2].dot(force));

    Vector3f got = GhostParticle::rotateForceToAtom(B, C, force, cd, sd);
    EXPECT_NEAR(got.getX(), expected.getX(), 1e-3f);
    EXPECT_NEAR(got.getY(), expected.getY(), 1e-3f);
    EXPECT_NEAR(got.getZ(), expected.getZ(), 1e-3f);
}

// =====================================================================================
// The closed-form axis reaction must make the whole interaction balance:
// zero net force and zero net torque, which is exactly what it solves for.
TEST(TestGhostParticle, redistributeAxisReaction_balances_force_and_torque)
{
    Vector3f B(0.0f, 0.0f, 0.0f);
    Vector3f C(1.526f, 0.0f, 0.0f);
    Vector3f atomI(-0.36f, 1.03f, 0.0f);
    Vector3f atomL(1.89f, 0.52f, 0.89f);

    // A genuine spring between two genuine ghosts -- not fabricated forces.
    // Both preconditions of the reconstruction come from this being a real
    // pair: a central force carries no net torque, and rotating about the
    // axis preserves the axial component, so neither survives to the axis
    // atoms. See GhostParticle.h.
    const float di = 1.1f, dl = 2.4f;
    Vector3f G = GhostParticle::computePositionByRotation(B, C, atomI, std::cos(di), std::sin(di));
    // the far side hangs off the same axis taken the other way round
    Vector3f H = GhostParticle::computePositionByRotation(C, B, atomL, std::cos(dl), std::sin(dl));

    Vector3f sep = H - G;
    Vector3f f = sep / sep.norm() * 3.7f;   // central: along the line joining them

    Vector3f F_i = GhostParticle::rotateForceToAtom(B, C, f, std::cos(di), std::sin(di));
    Vector3f F_l = GhostParticle::rotateForceToAtom(C, B, -f, std::cos(dl), std::sin(dl));

    Vector3f sumAtomF = F_i + F_l;
    Vector3f sumAtomT = ((atomI - B) ^ F_i) + ((atomL - B) ^ F_l);
    // what actually acted on the ghosts: a complete spring, so these cancel
    Vector3f sumGhostF = f + (-f);
    Vector3f sumGhostT = ((G - B) ^ f) + ((H - B) ^ (-f));

    Vector3f F_B, F_C;
    GhostParticle::redistributeAxisReaction(B, C, sumGhostF, sumGhostT, sumAtomF, sumAtomT, F_B, F_C);

    // Redistribution TRANSFERS: what the anchors end up with must equal
    // what acted on the ghosts, in both force and torque.
    Vector3f netF = F_i + F_l + F_B + F_C - sumGhostF;
    EXPECT_NEAR(netF.norm(), 0.0f, 1e-4f);

    Vector3f netT = sumAtomT + ((B - B) ^ F_B) + ((C - B) ^ F_C) - sumGhostT;
    EXPECT_NEAR(netT.norm(), 0.0f, 1e-4f);

    // and the reaction never pulls along the bond -- see GhostParticle.h
    Vector3f u = (C - B) / (C - B).norm();
    EXPECT_NEAR(u.dot(F_B), 0.0f, 1e-4f);
    EXPECT_NEAR(u.dot(F_C), 0.0f, 1e-4f);
}
