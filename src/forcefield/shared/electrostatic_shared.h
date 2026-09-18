#ifndef __BIOSPRING_ELECTROSTATIC_SHARED_H__
#define __BIOSPRING_ELECTROSTATIC_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU, on the same terms as spring_shared.h:
// this file is compiled twice, once as C++ into biospring-core and once as
// OpenCL C, prepended to biospring.cl. Both sides run the same text rather than
// two transcriptions of the same paper formula. The rules that keep it
// compilable by both toolchains are spelled out in spring_shared.h.
//
// Only the FORCE lives here. The energy stays in ../energy/electrostatic.hpp,
// in double precision and with Avogadro's number: an OpenCL device is not
// required to support double at all, and the device does not need the energy
// anyway -- the host recomputes it from the state the device returns, exactly
// as it does for the springs (see SpringNetworkOpenCL::_computeEnergiesFrom-
// DeviceState). Only what the kernels integrate has to be shared.

/// Coulomb force module between two point charges.
///
/// @param charge1, charge2 Particle charges, in elementary charge units (e).
/// @param distance         Distance between them, in A.
/// @param dielectric       Relative dielectric constant (dimensionless).
/// @param mindistance      Below this, the pair contributes nothing. It is not
///     a physical cutoff but a guard against 1/r^2 at r = 0, which two
///     particles that have landed on each other would otherwise turn into an
///     infinity that propagates through the whole step.
/// @param fourpi           4*PI, passed rather than written, so that there is
///     one definition of PI (constants.hpp) and no way for the two sides to
///     round it differently.
/// @param convert          The kJ.mol-1.A-1 -> Da.A.fs-2 conversion.
/// @return Force module along the axis between the two, in Da.A.fs-2 once
///     `convert` carries the conversion. NEGATIVE for like charges, which with
///     the caller's axis convention (neighbour minus self) makes them push
///     apart.
///
/// The caller applies the force field's coulomb scale afterwards rather than
/// folding it in here, because that is the order ForceField::compute-
/// ElectrostaticForceModule has always used and floating point is not
/// associative.
inline float biospring_electrostatic_force_module(float charge1, float charge2, float distance,
                                                  float dielectric, float mindistance,
                                                  float fourpi, float convert)
{
    if (distance < mindistance)
        return 0.0f;

    float force_module = -(charge1 * charge2) / (fourpi * dielectric * distance * distance);
    return force_module * convert;
}

#endif // __BIOSPRING_ELECTROSTATIC_SHARED_H__
