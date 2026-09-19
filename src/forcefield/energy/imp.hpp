#ifndef __IMP_ENERGY_HPP__
#define __IMP_ENERGY_HPP__

// Vector3f is used below and was not included: this header only compiled when
// something before it had already pulled it in.
#include "Vector3f.h"
#include "../constants.hpp"

#include <cmath>

namespace biospring
{
namespace forcefield
{

static const float ALPHA = 1.99;  // A^-1
static const float Z0 = 15.75;    // A
static const float ALIP = -0.018; // kJ.mol^-1

/// @brief Compute IMPALA energy (double membrane version).
/// @link https://doi.org/10.3390/membranes13030362
/// @callergraph
/// @param x x coordinate, in Angstrom (A)
/// @param y y coordinate, in Angstrom (A)
/// @param z z coordinate, in Angstrom (A)
/// @param surface Solvent accessible surface of the particle, in A^2
/// @param transfer Transfer energy of the particle, in kJ.mol-1.A-2 -- PER UNIT
///     OF ACCESSIBLE SURFACE, which is the unit the .ff files declare for
///     their transferIMP column. It is multiplied by `surface` below, so
///     kJ.mol-1 alone would not balance.
/// @param offset IMPALA double membrane offset in angstrom
/// @param uppermembtubecurv Tube curvature of the upper membrane, in A^-1
/// @param lowermembtubecurv Tube curvature of the lower membrane, in A^-1
/// @return IMPALA energy of the particle in kJ.mol-1
inline float imp_energy(float x, float y, float z, 
                        float surface, float transfer, 
                        float uppermemboffset=0.0,
                        float lowermemboffset=0.0,
                        float uppermembtubecurv=0.0,
                        float lowermembtubecurv=0.0)
{
    // membrane radius based on its curvature 
    float uppermemb_radius = uppermembtubecurv == 0.0 ? 1000000 : std::abs(1 / uppermembtubecurv);
    float lowermemb_radius = lowermembtubecurv == 0.0 ? 1000000 : std::abs(1 / lowermembtubecurv);

    // sign of curv (-1 or 1)
    int uppermemb_curv_sign = (uppermembtubecurv > 0.0) - (uppermembtubecurv < 0.0);
    int lowermemb_curv_sign = (lowermembtubecurv > 0.0) - (lowermembtubecurv < 0.0);

    // "Center" of the the tube curved membrane relative to the particle
    Vector3f v_uppermemb_center = Vector3f(0.0, y,  uppermemboffset - uppermemb_curv_sign * uppermemb_radius);
    Vector3f v_lowermemb_center = Vector3f(0.0, y, -lowermemboffset - lowermemb_curv_sign * lowermemb_radius);

    // Vector from particle to the "center"
    Vector3f v_upper = Vector3f(x, y, z) - v_uppermemb_center;
    Vector3f v_lower = Vector3f(x, y, z) - v_lowermemb_center;

    // new z (insertion in upper membrane) based on its curvature 
    float z_upper = uppermemb_curv_sign == 0.0 ? z :
        z > v_uppermemb_center.getZ() ? 
            uppermemb_curv_sign * v_upper.norm() + uppermemboffset - uppermemb_radius : // Particle inside the upper tube membrane zone
            -uppermemb_curv_sign * v_upper.norm() + uppermemboffset - uppermemb_radius; // Particle outside the upper tube membrane zone
    
    float z_lower = lowermemb_curv_sign == 0.0 ? z :
        z > v_lowermemb_center.getZ() ? 
            lowermemb_curv_sign * v_lower.norm() - lowermemboffset - lowermemb_radius : // Particle inside the lower tube membrane zone
            -lowermemb_curv_sign * v_lower.norm() - lowermemboffset - lowermemb_radius; // Particle outside the lower tube membrane zone

    double cz_upper = 0.5 - 1.0 / (1.0 + exp(ALPHA * (std::abs(z_upper - uppermemboffset) - Z0)));
    double cz_lower = 0.5 - 1.0 / (1.0 + exp(ALPHA * (std::abs(z_lower + lowermemboffset) - Z0)));

    double hydro_upper = -surface * transfer * cz_upper;
    double hydro_lower = -surface * transfer * cz_lower;

    double lipid_upper = ALIP * surface * cz_upper;
    double lipid_lower = ALIP * surface * cz_lower;

        // Which membranes are actually there.
    //
    // A membrane whose offset AND curvature are both zero is not a membrane
    // sitting at z = 0, it is an absent one. The previous test asked only
    // whether ALL FOUR parameters were zero, so the moment a client nudged one
    // of them the other membrane -- still at zero, still flat -- started
    // contributing as well, and it contributes the SAME profile: with
    // lowermembtubecurv = 0 the code sets z_lower = z, and with
    // lowermemboffset = 0 that makes cz_lower bit-identical to cz_upper.
    //
    // The result was a jump by a factor of exactly two on both the energy and
    // the force, the instant an offset left zero -- measured at 1e-30 A, long
    // before the two membranes are physically distinct.
    //
    // With all four at zero, which is every batch run, the upper one is the
    // single flat membrane and this reduces to what it always returned.
    const bool anyset = uppermemboffset != 0.0 || lowermemboffset != 0.0 ||
                        uppermembtubecurv != 0.0 || lowermembtubecurv != 0.0;
    const bool upper_present = !anyset || uppermemboffset != 0.0 || uppermembtubecurv != 0.0;
    const bool lower_present = lowermemboffset != 0.0 || lowermembtubecurv != 0.0;

    double total = 0.0;
    if (upper_present)
        total += hydro_upper + lipid_upper;
    if (lower_present)
        total += hydro_lower + lipid_lower;
    return total;
}

/// @brief Compute IMPALA force module (double membrane version)
/// @link https://doi.org/10.3390/membranes13030362
/// @callergraph
/// @param x x coordinate
/// @param y y coordinate
/// @param z z coordinate
/// @param surface Solvent accessible surface of the particle
/// @param transfer Transfer energy of the particle, in kJ.mol-1.A-2 (see
///     imp_energy above).
/// @param offset IMPALA double membrane offset in angstrom
/// @param uppermembtubecurv Tube curvature of the upper membrane
/// @param lowermembtubecurv Tube curvature of the lower membrane
/// @return IMPALA force vector of the particle in Da.A.fs-2
inline Vector3f imp_force_vector(float x, float y, float z, 
                              float surface, float transfer, 
                              float uppermemboffset=0.0,
                              float lowermemboffset=0.0,
                              float uppermembtubecurv=0.0,
                              float lowermembtubecurv=0.0)
{
    // membrane radius based on its curvature 
    float uppermemb_radius = uppermembtubecurv == 0.0 ? 1000000 : std::abs(1 / uppermembtubecurv);
    float lowermemb_radius = lowermembtubecurv == 0.0 ? 1000000 : std::abs(1 / lowermembtubecurv);

    // sign of curv
    int uppermemb_curv_sign = (uppermembtubecurv > 0.0) - (uppermembtubecurv < 0.0);
    int lowermemb_curv_sign = (lowermembtubecurv > 0.0) - (lowermembtubecurv < 0.0);

    // "Center" of the the tube curved membrane relative to the particle
    Vector3f v_uppermemb_center = Vector3f(0.0, y,  uppermemboffset - uppermemb_curv_sign * uppermemb_radius);
    Vector3f v_lowermemb_center = Vector3f(0.0, y, -lowermemboffset - lowermemb_curv_sign * lowermemb_radius);

    // Vector from particle to the "center"
    Vector3f v_upper = Vector3f(x, y, z) - v_uppermemb_center;
    Vector3f v_lower = Vector3f(x, y, z) - v_lowermemb_center;

    // new z (insertion in upper membrane) based on its curvature 
    float z_upper = uppermemb_curv_sign == 0.0 ? z :
        z > v_uppermemb_center.getZ() ? 
            uppermemb_curv_sign * v_upper.norm() + uppermemboffset - uppermemb_radius : // Particle inside the upper tube membrane zone
            -uppermemb_curv_sign * v_upper.norm() + uppermemboffset - uppermemb_radius; // Particle outside the upper tube membrane zone
    
    float z_lower = lowermemb_curv_sign == 0.0 ? z :
        z > v_lowermemb_center.getZ() ? 
            lowermemb_curv_sign * v_lower.norm() - lowermemboffset - lowermemb_radius : // Particle inside the lower tube membrane zone
            -lowermemb_curv_sign * v_lower.norm() - lowermemboffset - lowermemb_radius; // Particle outside the lower tube membrane zone

    auto expo_side = [](float z, float offset) { return exp(ALPHA * (std::abs(z + offset) - Z0)); };

    double dcz_upper = (ALPHA * (z_upper - uppermemboffset) * expo_side(z_upper, -uppermemboffset)) / 
                (pow(expo_side(z_upper, -uppermemboffset) + 1, 2.0) * std::abs(z_upper - uppermemboffset));
    

    double dcz_lower = (ALPHA * (z_lower + lowermemboffset) * expo_side(z_lower, lowermemboffset)) / 
                (pow(expo_side(z_lower, lowermemboffset) + 1, 2.0) * std::abs(z_lower + lowermemboffset));
    
    if (std::isnan(dcz_upper) || !std::isfinite(dcz_upper))
        dcz_upper = 0.0;
    if (std::isnan(dcz_lower) || !std::isfinite(dcz_lower))
        dcz_lower = 0.0;

    Vector3f v_upper_dir = uppermemb_curv_sign == 0.0 ? Vector3f(0, 0, 1.0) : v_upper;
    v_upper_dir.normalize();

    Vector3f v_lower_dir = lowermemb_curv_sign == 0.0 ? Vector3f(0, 0, 1.0) : v_lower;
    v_lower_dir.normalize();

    double hydro_upper = -surface * transfer * dcz_upper;
    double hydro_lower = -surface * transfer * dcz_lower;

    double lipid_upper = ALIP * surface * dcz_upper;
    double lipid_lower = ALIP * surface * dcz_lower;

    Vector3f force_module_upper = v_upper_dir * ((hydro_upper + lipid_upper) * GLOBAL_IMP_FORCE_CONVERT);
    Vector3f force_module_lower = v_lower_dir * ((hydro_lower + lipid_lower) * GLOBAL_IMP_FORCE_CONVERT);

    // Which membranes are actually there.
    //
    // A membrane whose offset AND curvature are both zero is not a membrane
    // sitting at z = 0, it is an absent one. The previous test asked only
    // whether ALL FOUR parameters were zero, so the moment a client nudged one
    // of them the other membrane -- still at zero, still flat -- started
    // contributing as well, and it contributes the SAME profile: with
    // lowermembtubecurv = 0 the code sets z_lower = z, and with
    // lowermemboffset = 0 that makes cz_lower bit-identical to cz_upper.
    //
    // The result was a jump by a factor of exactly two on both the energy and
    // the force, the instant an offset left zero -- measured at 1e-30 A, long
    // before the two membranes are physically distinct.
    //
    // With all four at zero, which is every batch run, the upper one is the
    // single flat membrane and this reduces to what it always returned.
    const bool anyset = uppermemboffset != 0.0 || lowermemboffset != 0.0 ||
                        uppermembtubecurv != 0.0 || lowermembtubecurv != 0.0;
    const bool upper_present = !anyset || uppermemboffset != 0.0 || uppermembtubecurv != 0.0;
    const bool lower_present = lowermemboffset != 0.0 || lowermembtubecurv != 0.0;

    Vector3f total = Vector3f(0.0f, 0.0f, 0.0f);
    if (upper_present)
        total = total + force_module_upper;
    if (lower_present)
        total = total + force_module_lower;
    return total;
}

} // namespace forcefield
} // namespace biospring

#endif // __IMP_ENERGY_HPP__