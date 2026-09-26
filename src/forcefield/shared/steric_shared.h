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

// How two radii combine. THIS IS A PROPERTY OF THE FORCE FIELD, not of the law,
// because the `radius` column of a .ff file does not hold the same quantity in
// every file, and no law can tell which it was handed:
//
//   data/forcefield/     radius column holds            rule
//   amber*.ff            AMBER's R* = rmin/2, 0.6-2.0 A  SUM      (rmin_ij = Ri + Rj)
//   model.ff             a bead radius, 2.6-3.9 A        SUM      (two beads touch)
//   CAonlyBaaden.ff      idem                            SUM
//   CAallresBaaden.ff    idem                            SUM
//   CAonlyLewitt.ff      a CONTACT DISTANCE, 5.0-7.4 A   GEOMETRIC_MEAN
//   CGZaccharias.ff      its own parameterisation        PRODUCT (the law's own)
//   martini.ff           sigma, 4.1/4.7 A                moot: every epsilon is 0
//
// The distinction is not academic. The shortest non-bonded CA-CA distance in a
// real protein is 3.7 A and the first percentile is 4.7 A, so Levitt's radii
// MEANED give a 5.1-6.5 A minimum, which is physical, while SUMMED they give
// 10.4-12.9 A, which would push every pair in the structure apart. Conversely
// AMBER's R* MEANED gives 1.9 A for a carbon pair, half its true contact
// distance, and a stack of nucleobases 3.4 A apart then sits on the attractive
// tail with no minimum at all -- which is the defect this enum exists to fix.
//
// Each law therefore carries the default of its own canonical parameter set
// (see STERIC_RADIUS_RULE_* in ../energy/steric.hpp), and steric.radiusrule in
// the .msp overrides it when a file is paired with a law it did not come from.
#define BIOSPRING_RADIUS_SUM             0
#define BIOSPRING_RADIUS_GEOMETRIC_MEAN  1
#define BIOSPRING_RADIUS_ARITHMETIC_MEAN 2
#define BIOSPRING_RADIUS_PRODUCT         3

// Combination rules, as in ../CombinationRules.hpp. Repeated here rather than
// included because a shared header may not include anything -- the CPU side
// keeps using the namespaced originals, which are the same one-liners.
inline float biospring_lorentz_berthelot_epsilon(float epsilon_i, float epsilon_j)
{
    return sqrt(epsilon_i * epsilon_j);
}
// SUM is Lorentz-Berthelot when the stored parameter is a radius: the arithmetic
// mean on sigma is a sum on R* = rmin/2, and it is what AMBER, CHARMM and OPLS
// all use. GEOMETRIC_MEAN is Good-Hope, right for a file whose column already
// holds a full contact distance. ARITHMETIC_MEAN is the same for such a file and
// is offered because it is what Lorentz-Berthelot means there. PRODUCT belongs
// to Zacharias' law, which pairs it with a product on epsilon too.
inline float biospring_combine_radius(int rule, float radius_i, float radius_j)
{
    if (rule == BIOSPRING_RADIUS_GEOMETRIC_MEAN)
        return sqrt(radius_i * radius_j);
    if (rule == BIOSPRING_RADIUS_ARITHMETIC_MEAN)
        return 0.5f * (radius_i + radius_j);
    if (rule == BIOSPRING_RADIUS_PRODUCT)
        return radius_i * radius_j;
    return radius_i + radius_j;
}
inline float biospring_zacharias_epsilon(float epsilon_i, float epsilon_j) { return epsilon_i * epsilon_j; }

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

/// AMBER's 12-6, with Lorentz-Berthelot on both epsilon and the radius.
/// @param mindistance Below this the pair contributes nothing: r^-13 at r = 0
///     is an infinity that would propagate through the whole step.
inline float biospring_steric_force_module_amber(float radius_i, float radius_j, float epsilon_i,
                                                 float epsilon_j, float distance,
                                                 float mindistance, float convert, int radiusrule)
{
    if (distance < mindistance)
        return 0.0f;

    float epsilon_ij = biospring_lorentz_berthelot_epsilon(epsilon_i, epsilon_j);
    float radius_ij = biospring_combine_radius(radiusrule, radius_i, radius_j);

    float repulsive = -epsilon_ij * 12.0f * (pow(radius_ij, 12.0f) / pow(distance, 13.0f));
    float attractive = epsilon_ij * 2.0f * 6.0f * (pow(radius_ij, 6.0f) / pow(distance, 7.0f));

    float force_module = repulsive + attractive;
    return force_module * convert;
}

/// Lewitt's 8-6, same combination rules as AMBER's. Its radius_ij is a minimum
/// too: dE/dr vanishes at r = radius_ij, where E = -epsilon_ij.
inline float biospring_steric_force_module_lewitt(float radius_i, float radius_j, float epsilon_i,
                                                  float epsilon_j, float distance,
                                                  float mindistance, float convert, int radiusrule)
{
    if (distance < mindistance)
        return 0.0f;

    float epsilon_ij = biospring_lorentz_berthelot_epsilon(epsilon_i, epsilon_j);
    float radius_ij = biospring_combine_radius(radiusrule, radius_i, radius_j);

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
    float radius_ij = biospring_combine_radius(BIOSPRING_RADIUS_PRODUCT, radius_i, radius_j);

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
                                           float mindistance, float convert, int radiusrule)
{
    if (mode == BIOSPRING_STERIC_AMBER_12_6)
        return biospring_steric_force_module_amber(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                                   mindistance, convert, radiusrule);
    if (mode == BIOSPRING_STERIC_LEWITT_8_6)
        return biospring_steric_force_module_lewitt(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                                    mindistance, convert, radiusrule);
    if (mode == BIOSPRING_STERIC_ZACHARIAS_8_6)
        return biospring_steric_force_module_zacharias(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                                       mindistance, convert);
    return biospring_steric_force_module_linear(radius_i, radius_j, distance, linearstiffness, convert);
}

// ---------------------------------------------------------------------------
// The ENERGIES of the four modes. They used to live only in
// ../energy/steric.hpp, so the device -- which walks every pair anyway -- had no
// way to report a steric energy, and the GPU path reported none.
//
// Each one is the integral of the force module above it, in kJ.mol-1, and they
// are written here so that the two backends run the same text: a force and an
// energy that are transcriptions of the same paper formula in two places drift,
// and this term has drifted before.

/// @return Overlap penalty, in kJ.mol-1. Zero when the two do not overlap.
inline float biospring_steric_energy_linear(float radius_i, float radius_j, float distance,
                                           float linearstiffness)
{
    float distancevar = distance - (radius_i + radius_j);
    if (distancevar > 0.0f)
        return 0.0f;
    return 0.5f * linearstiffness * distancevar * distancevar;
}

/// AMBER's 12-6, the integral of biospring_steric_force_module_amber.
inline float biospring_steric_energy_amber(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                          float distance, float mindistance, int radiusrule)
{
    if (distance < mindistance)
        return 0.0f;
    float epsilon_ij = biospring_lorentz_berthelot_epsilon(epsilon_i, epsilon_j);
    float radius_ij = biospring_combine_radius(radiusrule, radius_i, radius_j);
    // Written as the two terms the CPU has always summed, rather than factored,
    // so that sharing the text changes no digit of what it used to report.
    float repulsive = epsilon_ij * pow(radius_ij / distance, 12.0f);
    float attractive = -epsilon_ij * 2.0f * pow(radius_ij / distance, 6.0f);
    return repulsive + attractive;
}

/// Levitt's 8-6.
inline float biospring_steric_energy_lewitt(float radius_i, float radius_j, float epsilon_i, float epsilon_j,
                                           float distance, float mindistance, int radiusrule)
{
    if (distance < mindistance)
        return 0.0f;
    float epsilon_ij = biospring_lorentz_berthelot_epsilon(epsilon_i, epsilon_j);
    float radius_ij = biospring_combine_radius(radiusrule, radius_i, radius_j);
    float repulsive = epsilon_ij * 3.0f * pow(radius_ij / distance, 8.0f);
    float attractive = -epsilon_ij * 4.0f * pow(radius_ij / distance, 6.0f);
    return repulsive + attractive;
}

/// Zacharias' 8-6, with its own combining rules.
inline float biospring_steric_energy_zacharias(float radius_i, float radius_j, float epsilon_i,
                                              float epsilon_j, float distance, float mindistance)
{
    if (distance < mindistance)
        return 0.0f;
    float epsilon_ij = biospring_zacharias_epsilon(epsilon_i, epsilon_j);
    float radius_ij = biospring_combine_radius(BIOSPRING_RADIUS_PRODUCT, radius_i, radius_j);
    float repulsive = epsilon_ij * pow(radius_ij / distance, 8.0f);
    float attractive = -epsilon_ij * pow(radius_ij / distance, 6.0f);
    return repulsive + attractive;
}

/// The mode switch, in the same order as biospring_steric_force_module so that
/// a mode cannot pick one law's force and another's energy.
inline float biospring_steric_energy(int mode, float radius_i, float radius_j, float epsilon_i,
                                    float epsilon_j, float distance, float linearstiffness,
                                    float mindistance, int radiusrule)
{
    if (mode == BIOSPRING_STERIC_AMBER_12_6)
        return biospring_steric_energy_amber(radius_i, radius_j, epsilon_i, epsilon_j, distance, mindistance,
                                            radiusrule);
    if (mode == BIOSPRING_STERIC_LEWITT_8_6)
        return biospring_steric_energy_lewitt(radius_i, radius_j, epsilon_i, epsilon_j, distance, mindistance,
                                             radiusrule);
    if (mode == BIOSPRING_STERIC_ZACHARIAS_8_6)
        return biospring_steric_energy_zacharias(radius_i, radius_j, epsilon_i, epsilon_j, distance,
                                                 mindistance);
    return biospring_steric_energy_linear(radius_i, radius_j, distance, linearstiffness);
}

#endif // __BIOSPRING_STERIC_SHARED_H__
