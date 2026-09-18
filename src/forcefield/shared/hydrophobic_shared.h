#ifndef __BIOSPRING_HYDROPHOBIC_SHARED_H__
#define __BIOSPRING_HYDROPHOBIC_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU, on the same terms as spring_shared.h:
// compiled twice, once as C++ into biospring-core and once as OpenCL C
// prepended to biospring.cl. The rules that keep it compilable by both
// toolchains are spelled out there.
//
// THE LAW. This is the hydrophobic force law of Israelachvili & Pashley
// (Nature 300:341, 1982), whose published form is
//
//     W(D) = -22 exp(-D/10) mJ.m-2,   D in A
//
// an interaction free energy per unit area between two hydrophobic surfaces
// in water, decaying exponentially with a length of 10 A. Here it is written
// per pair of particles, with h_i h_j standing for the amplitude: a geometric
// (Berthelot) combining rule, so h_i is the square root of a per-particle
// hydrophobic energy. Because the product must stay positive for the pair to
// attract, h_i cannot be a signed hydrophobicity scale: a rectified one is
// required, and a particle that is not hydrophobic gets exactly 0.
//
// UNLIKE THE OTHER PAIRWISE LAWS, the energy lives here too rather than in
// ../energy/hydrophobic.hpp. It used to be kept out because it carried
// Avogadro's number in double precision -- but that was the bug: the force
// divides by Avogadro (through GLOBAL_SPRING_FORCE_CONVERT) while the energy
// multiplied by it, so the two were not the same function and disagreed by
// N_A * 1e-3 = 6.022e20. Written as below they are one function,
// differentiated and integrated, in single precision, and both sides run the
// same text.
//
//     E(r)   = -A * L * exp(-r/L)
//     -dE/dr = -A * exp(-r/L)
//
// so the force module, with the caller's axis convention (neighbour minus
// self, positive means attraction), is +A exp(-r/L).

/// Hydrophobic attraction between two particles.
///
/// @param hydrophobicity1, hydrophobicity2 Per-particle hydrophobicity, from
///     the force field's seventh column. Zero on a particle that has none,
///     which is why no separate "is this particle hydrophobic" test is needed
///     here: the product does it. Their product is in kJ.mol-1.A-1.
/// @param distance Distance between them, in A.
/// @param decaylength The distance over which the attraction decays, in A
///     (msp: hydrophobicity.decaylength). Passed rather than written because
///     it is a property of the solvent the modeller picks, in the same way
///     the dielectric constant is for the electrostatic law.
/// @param convert  kJ.mol-1.A-1 -> Da.A.fs-2.
/// @return Force module along the axis between the two, in Da.A.fs-2 once
///     `convert` carries the conversion.
///
/// The caller applies the force field's hydrophobicity scale afterwards rather
/// than folding it in here, because that is the order
/// ForceField::computeHydrophobicityForceModule has always used and floating
/// point is not associative.
inline float biospring_hydrophobic_force_module(float hydrophobicity1, float hydrophobicity2,
                                                float distance, float decaylength, float convert)
{
    float force_module = (hydrophobicity1 * hydrophobicity2) * exp(-distance / decaylength);
    return force_module * convert;
}

/// The potential the module above is the gradient of, in kJ.mol-1.
///
/// Note the factor `decaylength`: it is what makes this the integral of the
/// force rather than a transcription of the same expression. Dropping it
/// would leave an "energy" in kJ.mol-1.A-1, which is not an energy.
inline float biospring_hydrophobic_energy(float hydrophobicity1, float hydrophobicity2, float distance,
                                          float decaylength)
{
    return -(hydrophobicity1 * hydrophobicity2) * decaylength * exp(-distance / decaylength);
}

#endif // __BIOSPRING_HYDROPHOBIC_SHARED_H__
