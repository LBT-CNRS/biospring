#ifndef __IMP_ENERGY_HPP__
#define __IMP_ENERGY_HPP__

#include "../constants.hpp"

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
/// @param x x coordinate
/// @param y y coordinate
/// @param z z coordinate
/// @param surface Solvent accessible surface of the particle
/// @param transfer Transfer energy of the particle
/// @param offset IMPALA double membrane offset in angstrom
/// @param uppermembtubecurv Tube curvature of the upper membrane
/// @param lowermembtubecurv Tube curvature of the lower membrane
/// @return IMPALA energy of the particle in kJ.mol-1
inline float imp_energy(float x, float y, float z, 
                        float surface, float transfer, 
                        float uppermemboffset=0.0,
                        float lowermemboffset=0.0,
                        float uppermembtubecurv=0.0,
                        float lowermembtubecurv=0.0)
{
    // Initialize final Cz
    double cz = 0.0;

    // membrane radius based on its curvature 
    float uppermemb_radius = uppermembtubecurv == 0.0 ? 1000000 : abs(1/uppermembtubecurv);
    float lowermemb_radius = lowermembtubecurv == 0.0 ? 1000000 : abs(1/lowermembtubecurv);

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

    double cz_upper = 0.5 - 1.0 / (1.0 + exp(ALPHA * (abs(z_upper - uppermemboffset) - Z0)));
    double cz_lower = 0.5 - 1.0 / (1.0 + exp(ALPHA * (abs(z_lower + lowermemboffset) - Z0)));

    double hydro_upper = -surface * transfer * cz_upper;
    double hydro_lower = -surface * transfer * cz_lower;

    double lipid_upper = ALIP * surface * cz_upper;
    double lipid_lower = ALIP * surface * cz_lower;

        // If simple flat membrane
    if (uppermemboffset == 0.0 && lowermemboffset == 0.0 && uppermembtubecurv == 0.0 && lowermembtubecurv == 0.0)
        return hydro_upper + lipid_upper;
    else // Double membrane
        return hydro_upper + lipid_upper + hydro_lower + lipid_lower;
}

/// @brief Compute IMPALA force module (double membrane version)
/// @link https://doi.org/10.3390/membranes13030362
/// @callergraph
/// @param x x coordinate
/// @param y y coordinate
/// @param z z coordinate
/// @param surface Solvent accessible surface of the particle
/// @param transfer Transfer energy of the particle
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
    float uppermemb_radius = uppermembtubecurv == 0.0 ? 1000000 : abs(1/uppermembtubecurv);
    float lowermemb_radius = lowermembtubecurv == 0.0 ? 1000000 : abs(1/lowermembtubecurv);

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

    auto expo_side = [](float z, float offset) { return exp(ALPHA * (abs(z + offset) - Z0)); };

    double dcz_upper = (ALPHA * (z_upper - uppermemboffset) * expo_side(z_upper, -uppermemboffset)) / 
                (pow(expo_side(z_upper, -uppermemboffset) + 1, 2.0) * abs(z_upper - uppermemboffset));
    

    double dcz_lower = (ALPHA * (z_lower + lowermemboffset) * expo_side(z_lower, lowermemboffset)) / 
                (pow(expo_side(z_lower, lowermemboffset) + 1, 2.0) * abs(z_lower + lowermemboffset));
    
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

    // If simple flat membrane
    if (uppermemboffset == 0.0 && lowermemboffset == 0.0 && uppermembtubecurv == 0.0 && lowermembtubecurv == 0.0)
        return force_module_upper;
    else // Double membrane
        return force_module_upper + force_module_lower;
}

} // namespace forcefield
} // namespace biospring

#endif // __IMP_ENERGY_HPP__