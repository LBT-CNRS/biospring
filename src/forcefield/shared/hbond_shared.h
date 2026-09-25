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
// w = max(cos theta, 0)^2. Clamped at zero rather than squared through: an
// acceptor sitting BEHIND the donor is not a weak hydrogen bond, it is not one
// at all, and cos^2 alone would revive it with the sign lost.

/// @param cos_theta Cosine of the angle between the side's own direction and
///     the bond axis.
/// @return Dimensionless weight in [0, 1].
inline float biospring_hbond_angular_factor(float cos_theta)
{
    return cos_theta > 0.0f ? cos_theta * cos_theta : 0.0f;
}

/// d(weight)/d(cos theta), for the force.
inline float biospring_hbond_angular_derivative(float cos_theta)
{
    return cos_theta > 0.0f ? 2.0f * cos_theta : 0.0f;
}

#endif // __BIOSPRING_HBOND_SHARED_H__
