#ifndef __BIOSPRING_HYDROPHOBIC_SHARED_H__
#define __BIOSPRING_HYDROPHOBIC_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU, on the same terms as spring_shared.h:
// compiled twice, once as C++ into biospring-core and once as OpenCL C
// prepended to biospring.cl. The rules that keep it compilable by both
// toolchains are spelled out there.
//
// Only the force. The energy stays in ../energy/hydrophobic.hpp: it carries
// Avogadro's number in double precision, and the host recomputes it from the
// state the device returns anyway.

/// Hydrophobic attraction between two particles.
///
/// @param hydrophobicity1, hydrophobicity2 Per-particle hydrophobicity, from
///     the force field. Zero on a particle that has none, which is why no
///     separate "is this particle hydrophobic" test is needed here: the
///     product does it.
/// @param distance Distance between them, in A.
/// @param convert  kJ.mol-1.A-1 -> Da.A.fs-2.
/// @return Force module along the axis between the two, in Da.A.fs-2 once
///     `convert` carries the conversion.
///
/// The caller applies the force field's hydrophobicity scale afterwards rather
/// than folding it in here, because that is the order
/// ForceField::computeHydrophobicityForceModule has always used and floating
/// point is not associative.
inline float biospring_hydrophobic_force_module(float hydrophobicity1, float hydrophobicity2,
                                                float distance, float convert)
{
    float force_module = (hydrophobicity1 * hydrophobicity2) * exp(-distance);
    return force_module * convert;
}

#endif // __BIOSPRING_HYDROPHOBIC_SHARED_H__
