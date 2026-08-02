//
// Tests for spn::GhostParticle -- the massless virtual-site placement
// formula and its force redistribution (see GhostParticle.h and
// doc/BondedForceFieldSprings.md for the derivation).
//

#include <cmath>

#include <gtest/gtest.h>

#include "spn/GhostParticle.h"

using biospring::spn::GhostParticle;

// =====================================================================================
// computePosition places the virtual site at the expected distance from B
// and at the expected angle from the B->C axis, for a simple axis-aligned
// case where the answer can be checked by hand.
TEST(TestGhostParticle, computePosition_distance_and_angle)
{
    Vector3f B(0.0f, 0.0f, 0.0f);
    Vector3f C(0.0f, 0.0f, 2.0f);  // axis along +z
    Vector3f Ref(1.0f, 0.0f, 0.0f); // reference direction along +x (perpendicular to axis already)

    float r = 1.5f;
    float theta_deg = 90.0f; // perpendicular to the axis -> lies in the B-plane
    float delta_deg = 0.0f;  // aligned with the reference direction (+x)

    Vector3f X = GhostParticle::computePosition(B, C, Ref, r, theta_deg, delta_deg);

    EXPECT_NEAR(B.distance(X), r, 1e-4);
    // theta=90deg, delta=0 -> exactly r along +x from B in this axis-aligned setup.
    EXPECT_NEAR(X.getX(), r, 1e-4);
    EXPECT_NEAR(X.getY(), 0.0f, 1e-4);
    EXPECT_NEAR(X.getZ(), 0.0f, 1e-4);
}

TEST(TestGhostParticle, computePosition_theta_zero_lands_on_axis)
{
    Vector3f B(0.0f, 0.0f, 0.0f);
    Vector3f C(0.0f, 0.0f, 2.0f);
    Vector3f Ref(1.0f, 0.3f, 0.0f);

    Vector3f X = GhostParticle::computePosition(B, C, Ref, 0.7f, 0.0f, 0.0f);

    // theta=0 -> along the B->C axis regardless of delta.
    EXPECT_NEAR(X.getX(), 0.0f, 1e-4);
    EXPECT_NEAR(X.getY(), 0.0f, 1e-4);
    EXPECT_NEAR(X.getZ(), 0.7f, 1e-4);
}

// =====================================================================================
// redistributeForce must match the true gradient of a simple test energy
// E(X) = 0.5 |X - X0|^2 with respect to each anchor, computed by central
// finite differences directly on computePosition -- i.e. it must use the
// TRANSPOSE of the placement Jacobian (verified analytically and
// numerically in Python before this was translated to C++; using the
// Jacobian directly instead of its transpose silently breaks this test).
namespace
{
Vector3f numericalForce(const Vector3f & B, const Vector3f & C, const Vector3f & Ref, float r, float theta_deg,
                        float delta_deg, const Vector3f & X0, char which, float eps = 1e-3f)
{
    float grad[3];
    for (int k = 0; k < 3; k++)
    {
        Vector3f Bp = B, Cp = C, Rp = Ref;
        Vector3f Bm = B, Cm = C, Rm = Ref;
        Vector3f * targetp = (which == 'B') ? &Bp : (which == 'C') ? &Cp : &Rp;
        Vector3f * targetm = (which == 'B') ? &Bm : (which == 'C') ? &Cm : &Rm;
        (*targetp)[k] += eps;
        (*targetm)[k] -= eps;

        Vector3f Xp = GhostParticle::computePosition(Bp, Cp, Rp, r, theta_deg, delta_deg);
        Vector3f Xm = GhostParticle::computePosition(Bm, Cm, Rm, r, theta_deg, delta_deg);
        float Ep = 0.5f * (Xp - X0).dot(Xp - X0);
        float Em = 0.5f * (Xm - X0).dot(Xm - X0);
        grad[k] = (Ep - Em) / (2.0f * eps);
    }
    // True force = -dE/dAnchor.
    return Vector3f(-grad[0], -grad[1], -grad[2]);
}
} // namespace

TEST(TestGhostParticle, redistributeForce_matches_numerical_gradient)
{
    Vector3f B(1.2f, -0.4f, 0.3f);
    Vector3f C(2.5f, 0.9f, -1.1f);
    Vector3f Ref(0.6f, 1.8f, 0.7f);
    float r = 1.3f, theta_deg = 72.0f, delta_deg = -35.0f;

    Vector3f X = GhostParticle::computePosition(B, C, Ref, r, theta_deg, delta_deg);
    Vector3f X0(0.5f, -0.2f, 1.0f); // arbitrary reference point for the test energy.
    Vector3f force = -(X - X0);     // force on the ghost = -dE/dX = -(X - X0).

    Vector3f F_B, F_C, F_Ref;
    GhostParticle::redistributeForce(B, C, Ref, r, theta_deg, delta_deg, force, F_B, F_C, F_Ref);

    Vector3f F_B_numeric = numericalForce(B, C, Ref, r, theta_deg, delta_deg, X0, 'B');
    Vector3f F_C_numeric = numericalForce(B, C, Ref, r, theta_deg, delta_deg, X0, 'C');
    Vector3f F_Ref_numeric = numericalForce(B, C, Ref, r, theta_deg, delta_deg, X0, 'R');

    EXPECT_NEAR(F_B.getX(), F_B_numeric.getX(), 1e-2);
    EXPECT_NEAR(F_B.getY(), F_B_numeric.getY(), 1e-2);
    EXPECT_NEAR(F_B.getZ(), F_B_numeric.getZ(), 1e-2);

    EXPECT_NEAR(F_C.getX(), F_C_numeric.getX(), 1e-2);
    EXPECT_NEAR(F_C.getY(), F_C_numeric.getY(), 1e-2);
    EXPECT_NEAR(F_C.getZ(), F_C_numeric.getZ(), 1e-2);

    EXPECT_NEAR(F_Ref.getX(), F_Ref_numeric.getX(), 1e-2);
    EXPECT_NEAR(F_Ref.getY(), F_Ref_numeric.getY(), 1e-2);
    EXPECT_NEAR(F_Ref.getZ(), F_Ref_numeric.getZ(), 1e-2);
}

// =====================================================================================
// Linear momentum conservation: since the placement formula only ever
// depends on position DIFFERENCES (C-B, Ref-B) plus an additive B, it is
// translation-invariant, so the 3 redistributed forces must sum back to
// exactly the original force on the virtual site (no force is created or
// lost by the redistribution).
TEST(TestGhostParticle, redistributeForce_conserves_total_force)
{
    Vector3f B(-0.3f, 2.1f, 0.9f);
    Vector3f C(1.1f, -0.5f, 1.6f);
    Vector3f Ref(0.2f, 0.4f, -1.3f);
    float r = 0.9f, theta_deg = 110.0f, delta_deg = 48.0f;

    Vector3f force(0.7f, -1.4f, 2.2f);
    Vector3f F_B, F_C, F_Ref;
    GhostParticle::redistributeForce(B, C, Ref, r, theta_deg, delta_deg, force, F_B, F_C, F_Ref);

    Vector3f total = F_B + F_C + F_Ref;
    EXPECT_NEAR(total.getX(), force.getX(), 1e-3);
    EXPECT_NEAR(total.getY(), force.getY(), 1e-3);
    EXPECT_NEAR(total.getZ(), force.getZ(), 1e-3);
}

// =====================================================================================
// Translation invariance: computePosition only ever depends on position
// DIFFERENCES (C-B, Ref-B) plus an additive B, so rigidly translating all
// 3 anchors must translate the placed virtual site by exactly the same
// amount (this is also what makes the total-force conservation test above
// hold, by Noether's theorem).
TEST(TestGhostParticle, computePosition_is_translation_invariant)
{
    Vector3f B(0.2f, 0.1f, -0.3f);
    Vector3f C(1.5f, 0.4f, 0.8f);
    Vector3f Ref(0.9f, -1.2f, 0.1f);
    float r = 1.1f, theta_deg = 65.0f, delta_deg = 20.0f;

    Vector3f X_before = GhostParticle::computePosition(B, C, Ref, r, theta_deg, delta_deg);

    Vector3f shift(2.0f, -1.0f, 0.5f);
    Vector3f X_after =
        GhostParticle::computePosition(B + shift, C + shift, Ref + shift, r, theta_deg, delta_deg);

    Vector3f expected = X_before + shift;
    EXPECT_NEAR(X_after.getX(), expected.getX(), 1e-4);
    EXPECT_NEAR(X_after.getY(), expected.getY(), 1e-4);
    EXPECT_NEAR(X_after.getZ(), expected.getZ(), 1e-4);
}
