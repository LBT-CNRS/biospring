#include "GhostParticle.h"

#include <cmath>

namespace biospring
{
namespace spn
{

namespace
{

// Minimal row-major 3x3 matrix helper, private to this translation unit.
// Only what the placement Jacobian derivation needs (see
// doc/BondedForceFieldSprings.md for the vector-calculus derivation and
// its numerical verification against finite differences).
struct Mat3
{
    float m[3][3];

    static Mat3 zero()
    {
        Mat3 r;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                r.m[i][j] = 0.0f;
        return r;
    }

    static Mat3 identity()
    {
        Mat3 r = zero();
        r.m[0][0] = r.m[1][1] = r.m[2][2] = 1.0f;
        return r;
    }

    // Outer product a * b^T.
    static Mat3 outer(const Vector3f & a, const Vector3f & b)
    {
        float av[3] = {a.getX(), a.getY(), a.getZ()};
        float bv[3] = {b.getX(), b.getY(), b.getZ()};
        Mat3 r;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                r.m[i][j] = av[i] * bv[j];
        return r;
    }

    // Skew-symmetric cross-product matrix: skew(a) * x == a x x.
    static Mat3 skew(const Vector3f & a)
    {
        Mat3 r = zero();
        r.m[0][1] = -a.getZ();
        r.m[0][2] = a.getY();
        r.m[1][0] = a.getZ();
        r.m[1][2] = -a.getX();
        r.m[2][0] = -a.getY();
        r.m[2][1] = a.getX();
        return r;
    }

    Mat3 operator+(const Mat3 & o) const
    {
        Mat3 r;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                r.m[i][j] = m[i][j] + o.m[i][j];
        return r;
    }

    Mat3 operator-(const Mat3 & o) const
    {
        Mat3 r;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                r.m[i][j] = m[i][j] - o.m[i][j];
        return r;
    }

    Mat3 operator*(float s) const
    {
        Mat3 r;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                r.m[i][j] = m[i][j] * s;
        return r;
    }

    Mat3 operator*(const Mat3 & o) const
    {
        Mat3 r = zero();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    r.m[i][j] += m[i][k] * o.m[k][j];
        return r;
    }

    Vector3f operator*(const Vector3f & v) const
    {
        float vv[3] = {v.getX(), v.getY(), v.getZ()};
        float rv[3] = {0.0f, 0.0f, 0.0f};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                rv[i] += m[i][j] * vv[j];
        return Vector3f(rv[0], rv[1], rv[2]);
    }

    Mat3 transpose() const
    {
        Mat3 r;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                r.m[i][j] = m[j][i];
        return r;
    }
};

// Shared geometry for both computePosition and redistributeForce: the
// local frame (u, v, w) built from the 3 anchors, plus the Jacobians of
// u, v (and, if needed, w) with respect to each anchor. Computing these
// once and reusing them avoids repeating the whole derivation twice.
struct LocalFrame
{
    Vector3f u, v, w;
    float L, Lv;
    Mat3 Ju_B, Ju_C;    // Ju_R is always zero (u does not depend on Ref).
    Mat3 Jv_B, Jv_C, Jv_R;
};

// Just the orthonormal frame (u, v, w), without any of the derivative
// machinery below. Placing a ghost needs nothing else: X = B + a*u + b*v
// + c*w. Splitting this out of computeLocalFrame is what keeps the
// placement pass from building -- and immediately discarding -- the ~9
// 3x3 matrix products of the Jacobian. Profiling example/072's
// bonded-only model put computeLocalFrame at ~59% of the step, called
// once per ghost from each of the two passes.
void computeFrameBasis(const Vector3f & B, const Vector3f & C, const Vector3f & Ref, Vector3f & u, Vector3f & v,
                       Vector3f & w)
{
    Vector3f d = C - B;
    u = d / d.norm();

    Vector3f t = Ref - B;
    Vector3f tperp = t - u * u.dot(t);
    v = tperp / tperp.norm();

    w = u ^ v;
}

LocalFrame computeLocalFrame(const Vector3f & B, const Vector3f & C, const Vector3f & Ref)
{
    LocalFrame f;

    Vector3f d = C - B;
    f.L = d.norm();
    f.u = d / f.L;

    Vector3f t = Ref - B;
    float s = f.u.dot(t);
    Vector3f tperp = t - f.u * s;
    f.Lv = tperp.norm();
    f.v = tperp / f.Lv;

    f.w = f.u ^ f.v;

    Mat3 I = Mat3::identity();
    Mat3 Pu = I - Mat3::outer(f.u, f.u);
    Mat3 Pv = I - Mat3::outer(f.v, f.v);

    f.Ju_B = Pu * (-1.0f / f.L);
    f.Ju_C = Pu * (1.0f / f.L);

    // M = u t^T + s I
    Mat3 M = Mat3::outer(f.u, t) + I * s;

    Mat3 Jtperp_B = Pu * (-1.0f) - M * f.Ju_B;
    Mat3 Jtperp_C = M * f.Ju_C * (-1.0f);
    Mat3 Jtperp_R = Pu;

    f.Jv_B = Pv * Jtperp_B * (1.0f / f.Lv);
    f.Jv_C = Pv * Jtperp_C * (1.0f / f.Lv);
    f.Jv_R = Pv * Jtperp_R * (1.0f / f.Lv);

    return f;
}

} // namespace

GhostLocalOffset GhostParticle::localOffset(float r, float theta_deg, float delta_deg)
{
    const float theta = theta_deg * (static_cast<float>(M_PI) / 180.0f);
    const float delta = delta_deg * (static_cast<float>(M_PI) / 180.0f);
    GhostLocalOffset o;
    o.a = r * std::cos(theta);
    o.b = r * std::sin(theta) * std::cos(delta);
    o.c = r * std::sin(theta) * std::sin(delta);
    return o;
}

Vector3f GhostParticle::computePosition(const Vector3f & B, const Vector3f & C, const Vector3f & Ref,
                                        const GhostLocalOffset & offset)
{
    Vector3f u, v, w;
    computeFrameBasis(B, C, Ref, u, v, w);
    return B + u * offset.a + v * offset.b + w * offset.c;
}

Vector3f GhostParticle::computePosition(const Vector3f & B, const Vector3f & C, const Vector3f & Ref, float r,
                                        float theta_deg, float delta_deg)
{
    return computePosition(B, C, Ref, localOffset(r, theta_deg, delta_deg));
}

void GhostParticle::redistributeForce(const Vector3f & B, const Vector3f & C, const Vector3f & Ref, float r,
                                      float theta_deg, float delta_deg, const Vector3f & force, Vector3f & F_B,
                                      Vector3f & F_C, Vector3f & F_Ref)
{
    redistributeForce(B, C, Ref, localOffset(r, theta_deg, delta_deg), force, F_B, F_C, F_Ref);
}

void GhostParticle::redistributeForce(const Vector3f & B, const Vector3f & C, const Vector3f & Ref,
                                      const GhostLocalOffset & offset, const Vector3f & force, Vector3f & F_B,
                                      Vector3f & F_C, Vector3f & F_Ref)
{
    // Unlike placement, this genuinely needs the derivatives: the force on
    // each anchor is J_anchor^T . F_ghost (see the class comment).
    LocalFrame f = computeLocalFrame(B, C, Ref);

    const float a = offset.a;
    const float b = offset.b;
    const float c = offset.c;

    Mat3 Sv = Mat3::skew(f.v);
    Mat3 Su = Mat3::skew(f.u);
    Mat3 Jw_B = Su * f.Jv_B - Sv * f.Ju_B;
    Mat3 Jw_C = Su * f.Jv_C - Sv * f.Ju_C;
    Mat3 Jw_R = Su * f.Jv_R; // Ju_R == 0

    Mat3 I = Mat3::identity();
    Mat3 JX_B = I + f.Ju_B * a + f.Jv_B * b + Jw_B * c;
    Mat3 JX_C = f.Ju_C * a + f.Jv_C * b + Jw_C * c;
    Mat3 JX_R = f.Jv_R * b + Jw_R * c;

    // F_anchor = J_anchor^T . force (virtual work principle -- see
    // GhostParticle.h and doc/BondedForceFieldSprings.md; verified
    // numerically that using J directly instead of J^T is wrong).
    F_B = JX_B.transpose() * force;
    F_C = JX_C.transpose() * force;
    F_Ref = JX_R.transpose() * force;
}

} // namespace spn
} // namespace biospring
