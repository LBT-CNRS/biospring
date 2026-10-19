#ifndef __BIOSPRING_HBOND_SHARED_H__
#define __BIOSPRING_HBOND_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU -- see spring_shared.h for the rules this
// file obeys and why they exist.
//
// A hydrogen bond here is a Morse well in the donor-acceptor heavy-atom
// distance, weighted at BOTH ends by how well the bond lines up with the
// direction each side offers: the donor's hydrogen, standing in as
// "away from the antecedent", and the acceptor's lone pairs, standing in the
// same way. The physics behind the parameters is documented once, in
// forcefield/energy/hydrogenbond.hpp, which now calls straight into this file
// so there is one text rather than two.

/// Morse well between a donor and an acceptor heavy atom, no explicit
/// hydrogen.
/// @param distance    Donor-acceptor distance, in A.
/// @param welldepth   De, in kJ.mol-1 (already molar).
/// @param equilibrium re, in A.
/// @param width       a, in A-1.
/// @return Energy in kJ.mol-1: 0 at infinity, -De at equilibrium.
inline float biospring_hbond_energy(float distance, float welldepth, float equilibrium, float width)
{
    float u = exp(-width * (distance - equilibrium));
    float one_minus_u = 1.0f - u;
    return welldepth * one_minus_u * one_minus_u - welldepth;
}

/// dV/dr of the well above, in the "positive = attractive, on the FIRST
/// particle" convention the rest of the code uses.
/// @param convert Folded in by the caller rather than read from a global:
///     GLOBAL_SPRING_FORCE_CONVERT, which turns kJ.mol-1.A-1 into the
///     integrator's Da.A.fs-2. The well depth is molar, like a spring
///     stiffness, so it needs the same conversion.
/// @return Force module along the donor->acceptor axis.
inline float biospring_hbond_force_module(float distance, float welldepth, float equilibrium, float width,
                                          float convert)
{
    float u = exp(-width * (distance - equilibrium));
    float force_module = 2.0f * width * welldepth * u * (1.0f - u);
    return force_module * convert;
}

// ---- Angular weight -----------------------------------------------------
//
// w = max(cos theta, 0)^n, n = 2 by default. Clamped at zero rather than raised
// through: an acceptor sitting BEHIND the donor is not a weak hydrogen bond, it
// is not one at all, and an even power alone would revive it with the sign lost.
//
// n is a setting because HiRE-RNA v3 uses SIX where this has always used two
// (its nu(alpha) = cos^6(alpha - alpha0)), which is a far narrower angular
// acceptance. The exponent is a small integer, so the power is a loop rather
// than a pow() call: at n = 2 that is one multiply, and the function is on the
// innermost path of the term.

/// @param cos_theta Cosine of the angle between the side's own direction and
///     the bond axis.
/// @return Dimensionless weight in [0, 1].
inline float biospring_hbond_angular_factor_n(float cos_theta, int exponent)
{
    if (cos_theta <= 0.0f)
        return 0.0f;
    float w = 1.0f;
    for (int k = 0; k < exponent; ++k)
        w *= cos_theta;
    return w;
}

/// d(weight)/d(cos theta), for the force.
inline float biospring_hbond_angular_derivative_n(float cos_theta, int exponent)
{
    if (cos_theta <= 0.0f || exponent < 1)
        return 0.0f;
    float w = (float) exponent;
    for (int k = 1; k < exponent; ++k)
        w *= cos_theta;
    return w;
}

inline float biospring_hbond_angular_factor(float cos_theta)
{
    return biospring_hbond_angular_factor_n(cos_theta, 2);
}

inline float biospring_hbond_angular_derivative(float cos_theta)
{
    return biospring_hbond_angular_derivative_n(cos_theta, 2);
}

// ---- Planarity factor --------------------------------------------------
//
// HiRE-RNA v3's equation (5), which MULTIPLIES the well. Two bases pair when
// they meet edge to edge IN A SHARED PLANE, and the factor says how well they
// do: each of one base's plane atoms contributes e^(-(d/delta)^2) for its
// distance d out of the OTHER base's plane, the three are averaged, and the two
// sides multiply.
//
// Averaging rather than summing, unlike the reference's bare sum, so that a
// perfectly coplanar pair scores exactly 1. That matters here and not there:
// this multiplies a well depth that is already calibrated, and a factor that
// could exceed 1 would silently recalibrate it.
//
// It is a REDUCTION, never a reward: in [0, 1], 1 only at perfect coplanarity.

/// One plane atom's contribution.
/// @param out_of_plane Signed distance from the partner base's plane, in A.
/// @param delta        Width in A; the contribution is e^(-(d/delta)^2).
inline float biospring_hbond_planarity_term(float out_of_plane, float delta)
{
    const float u = out_of_plane / delta;
    return exp(-u * u);
}

/// d(term)/d(out_of_plane), in A-1.
inline float biospring_hbond_planarity_derivative(float out_of_plane, float delta)
{
    const float u = out_of_plane / delta;
    return -2.0f * out_of_plane / (delta * delta) * exp(-u * u);
}

#endif // __BIOSPRING_HBOND_SHARED_H__
