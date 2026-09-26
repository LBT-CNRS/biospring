#ifndef __BIOSPRING_ELECTROSTATIC_SHARED_H__
#define __BIOSPRING_ELECTROSTATIC_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU, on the same terms as spring_shared.h:
// this file is compiled twice, once as C++ into biospring-core and once as
// OpenCL C, prepended to biospring.cl. Both sides run the same text rather than
// two transcriptions of the same paper formula. The rules that keep it
// compilable by both toolchains are spelled out in spring_shared.h.
//
// The force AND the energy live here. The energy used to stay in
// ../energy/electrostatic.hpp, in double precision, on the argument that the
// host would recompute it from the state the device returns "exactly as it does
// for the springs". It never did: _computeEnergiesFromDeviceState reads a
// per-particle buffer that the spring, torsion and hydrogen bond kernels fill,
// and no pairwise kernel filled one -- so the GPU path reported no steric, no
// Coulomb and no hydrophobic energy at all.
//
// It is written below so that single precision is enough: the charges stay in e
// and the distance in A, and one constant carries every unit at once, exactly as
// the force module does. Nothing in it goes near a float's limits -- the old
// double-precision version multiplied 1e-28 by Avogadro's number to get back to
// a number of order 1, which is what needed the double.

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

/// Coulomb energy between two point charges.
///
/// @param charge1, charge2 Particle charges, in elementary charge units (e).
/// @param distance         Distance between them, in A.
/// @param dielectric       Relative dielectric constant (dimensionless).
/// @param mindistance      Below this the pair contributes nothing, the same
///     guard the force module applies, so the two agree about which pairs
///     exist.
/// @param convert          (q1 q2)/(4.pi.eps0.dielectric.d) -> kJ.mol-1, i.e.
///     GLOBAL_ELECTROSTATIC_ENERGY_CONVERT, the familiar 1389.35.
/// @return Energy in kJ.mol-1, POSITIVE for like charges.
inline float biospring_electrostatic_energy(float charge1, float charge2, float distance, float dielectric,
                                           float mindistance, float convert)
{
    if (distance < mindistance)
        return 0.0f;
    return convert * (charge1 * charge2) / (dielectric * distance);
}

#endif // __BIOSPRING_ELECTROSTATIC_SHARED_H__
