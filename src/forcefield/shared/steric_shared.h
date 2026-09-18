#ifndef __BIOSPRING_STERIC_SHARED_H__
#define __BIOSPRING_STERIC_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU, on the same terms as spring_shared.h and
// electrostatic_shared.h: compiled twice, once as C++ into biospring-core and
// once as OpenCL C prepended to biospring.cl. The rules that keep it compilable
// by both toolchains are spelled out in spring_shared.h.
//
// One departure from the C++ this replaces, and it is deliberate: every literal
// carries the f suffix. The originals wrote `2.0`, `12.0`, which are DOUBLE
// literals and promote the whole expression, and an OpenCL device is not
// required to support double at all -- on one that does not, the kernel would
// simply fail to compile. The arithmetic is therefore float throughout on both
// sides, which moves the CPU's answers in their last bits and keeps the two
// backends running the same text. That is the trade spring_shared.h already
// made.
//
// Only the force modules. The energies stay in ../energy/steric.hpp: the host
// recomputes them from the state the device returns, so the device never needs
// them.

// Which law. Mirrors the .msp's steric.mode, which the host resolves to one of
// these before handing it down -- a kernel cannot be given a string.
#define BIOSPRING_STERIC_LINEAR        0
#define BIOSPRING_STERIC_AMBER_12_6    1
#define BIOSPRING_STERIC_LEWITT_8_6    2
#define BIOSPRING_STERIC_ZACHARIAS_8_6 3

// Combination rules, as in ../CombinationRules.hpp. Repeated here rather than
// included because a shared header may not include anything -- the CPU side
// keeps using the namespaced originals, which are the same one-liners.
inline float biospring_lorentz_berthelot_epsilon(float epsilon_i, float epsilon_j)
{
    return sqrt(epsilon_i * epsilon_j);
}
inline float biospring_good_hope_radius(float radius_i, float radius_j)
{
    return sqrt(radius_i * radius_j);
}
inline float biospring_zacharias_epsilon(float epsilon_i, float epsilon_j) { return epsilon_i * epsilon_j; }
inline float biospring_zacharias_radius(float radius_i, float radius_j) { return radius_i * radius_j; }

/// Linear overlap penalty: nothing until the two touch, then a spring.
/// @param stiffness The overlap stiffness, in kJ.mol-1.A-2. Passed rather than
///     read from a constant so the two sides cannot hold different values --
///     the energy and the force here once did, 100 against 1.0.
/// @param convert   kJ.mol-1.A-1 -> Da.A.fs-2.
inline float biospring_steric_force_module_linear(float radius_i, float radius_j, float distance,
                                                  float stiffness, float convert)
{
    float equilibrium = radius_i + radius_j;
    float distancevar = (distance - equilibrium);

    if (distancevar > 0.0f)
        return 0.0f;

    float force_module = -stiffness * fabs(distancevar);
    return force_module * convert;
}

/// AMBER's 12-6, with Lorentz-Berthelot on epsilon and Good-Hope on the radius.
/// @param mindistance Below this the pair contributes nothing: r^-13 at r = 0
///     is an infinity that would propagate through the whole step.
inline float biospring_steric_force_module_amber(float radius_i, float radius_j, float epsilon_i,
                                                 float epsilon_j, float distance,
                                                 float mindistance, float convert)
{
    if (distance < mindistance)
        return 0.0f;

    float epsilon_ij = biospring_lorentz_berthelot_epsilon(epsilon_i, epsilon_j);
    float radius_ij = biospring_good_hope_radius(radius_i, radius_j);

    float repulsive = -epsilon_ij * 12.0f * (pow(radius_ij, 12.0f) / pow(distance, 13.0f));
    float attractive = epsilon_ij * 2.0f * 6.0f * (pow(radius_ij, 6.0f) / pow(distance, 7.0f));

    float force_module = repulsive + attractive;
    return force_module * convert;
}

/// Lewitt's 8-6, same combination rules as AMBER's.
inline float biospring_steric_force_module_lewitt(float radius_i, float radius_j, float epsilon_i,
                                                  float epsilon_j, float distance,
                                                  float mindistance, float convert)
{
    if (distance < mindistance)
        return 0.0f;

    float epsilon_ij = biospring_lorentz_berthelot_epsilon(epsilon_i, epsilon_j);
    float radius_ij = biospring_good_hope_radius(radius_i, radius_j);

    float repulsive = -epsilon_ij * 3.0f * 8.0f * (pow(radius_ij, 8.0f) / pow(distance, 9.0f));
    float attractive = epsilon_ij * 4.0f * 6.0f * (pow(radius_ij, 6.0f) / pow(distance, 7.0f));

    float force_module = repulsive + attractive;
    return force_module * convert;
}

/// Zacharias's 8-6, which combines both parameters by plain product.
inline float biospring_steric_force_module_zacharias(float radius_i, float radius_j, float epsilon_i,
                                                     float epsilon_j, float distance,
                                                     float mindistance, float convert)
{
    if (distance < mindistance)
        return 0.0f;

    float epsilon_ij = biospring_zacharias_epsilon(epsilon_i, epsilon_j);
    float radius_ij = biospring_zacharias_radius(radius_i, radius_j);

    float repulsive = -epsilon_ij * 8.0f * (pow(radius_ij, 8.0f) / pow(distance, 9.0f));
    float attractive = epsilon_ij * 6.0f * (pow(radius_ij, 6.0f) / pow(distance, 7.0f));

    float force_module = repulsive + attractive;
    return force_module * convert;
}

/// Picks the law. The host resolves steric.mode to one of the constants above
/// once, at setup; this is what lets one kernel serve all four without a
/// separate build per mode.
inline float biospring_steric_force_module(int mode, float radius_i, float radius_j, float epsilon_i,
                                           float epsilon_j, float distance, float linearstiffness,
                                           float mindistance, float convert)
{
    if (mode == BIOSPRING_STERIC_AMBER_12_6)
        return biospring_steric_force_module_amber(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                                   mindistance, convert);
    if (mode == BIOSPRING_STERIC_LEWITT_8_6)
        return biospring_steric_force_module_lewitt(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                                    mindistance, convert);
    if (mode == BIOSPRING_STERIC_ZACHARIAS_8_6)
        return biospring_steric_force_module_zacharias(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                                       mindistance, convert);
    return biospring_steric_force_module_linear(radius_i, radius_j, distance, linearstiffness, convert);
}

#endif // __BIOSPRING_STERIC_SHARED_H__
