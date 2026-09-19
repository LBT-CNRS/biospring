#ifndef __BIOSPRING_IMP_SHARED_H__
#define __BIOSPRING_IMP_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU, on the same terms as spring_shared.h:
// compiled twice, once as C++ into biospring-core and once as OpenCL C
// prepended to biospring.cl. The rules that keep it compilable by both
// toolchains are spelled out there.
//
// ONLY THE FLAT SINGLE MEMBRANE. IMPALA's general form allows a second
// membrane and a tube curvature for each (see ../energy/imp.hpp), and those
// need vectors -- the force is no longer along z, it points away from the tube
// axis. What is here is the case every batch run takes: the four geometry
// parameters default to 0.0 and no .msp key can change them, only MDDriver can,
// interactively, mid-run. imp.hpp tests for exactly that case and returns early,
// and this is the arithmetic of that early return.
//
// The model is a sigmoid in z alone. cz goes from -0.5 in the membrane core to
// +0.5 in bulk water, and each particle's contribution is its ACCESSIBLE
// SURFACE times an energy per unit of that surface -- which is why IMPALA
// implies --sasa and why the pairwise hydrophobic term, which has no surface in
// it, is a different model with a different column in the .ff.

/// The membrane profile at height z. -0.5 in the core, +0.5 in bulk water.
/// @param alpha 1.99 A-1, @param z0 15.75 A -- passed rather than written so
///     there is one definition, in ../energy/imp.hpp.
inline float biospring_imp_cz(float z, float alpha, float z0)
{
    return 0.5f - 1.0f / (1.0f + exp(alpha * (fabs(z) - z0)));
}

/// Its derivative, which is what the force is built from.
///
/// Written as alpha * sign(z) / (exp(u) + 2 + exp(-u)) rather than the
/// algebraically equal alpha * z * exp(u) / ((exp(u)+1)^2 * |z|), because the
/// second one is not computable in single precision. imp.hpp uses it in double
/// and catches what is left with isnan/isfinite; here u = alpha*(|z| - z0)
/// passes 88 at |z| = 60 A, exp overflows, (e+1)^2 overflows too, and inf/inf
/// is a NaN that spreads to the whole trajectory. Example 051 is long enough
/// in z to reach it, and did: "Found non-finite position for particle 1".
///
/// In this form the overflow is harmless -- one term of the denominator goes to
/// infinity and the quotient goes to zero, which is the right answer far from
/// the membrane.
///
/// The |z| in the original is a sign, and z = 0 is the one point where it is
/// not defined: the profile is even, so the force there is zero.
inline float biospring_imp_dcz(float z, float alpha, float z0)
{
    if (z == 0.0f)
        return 0.0f;
    float u = alpha * (fabs(z) - z0);
    float d = exp(u) + 2.0f + exp(-u);
    if (!(d > 0.0f))
        return 0.0f;
    return (z > 0.0f ? alpha : -alpha) / d;
}

/// IMPALA energy of one particle, in kJ.mol-1, before the force field's scale.
///
/// @param surface  Solvent-accessible surface of the particle, in A2.
/// @param transfer Its transfer energy PER UNIT of that surface, in
///     kJ.mol-1.A-2 -- the .ff's sixth column. Not to be confused with the
///     seventh, which feeds the pairwise hydrophobic term.
/// @param alip     The lipid term's own energy per unit surface,
///     -0.018 kJ.mol-1.A-2.
inline float biospring_imp_energy(float z, float surface, float transfer, float alip, float alpha, float z0)
{
    return surface * biospring_imp_cz(z, alpha, z0) * (alip - transfer);
}

/// The z component of the force, in Da.A.fs-2 once `convert` carries
/// GLOBAL_IMP_FORCE_CONVERT. The other two components are zero: a flat membrane
/// has no lateral gradient.
///
/// Note the sign. ForceField::computeIMPForceVector returns the NEGATIVE of
/// imp_force_vector, so the minus belongs to the force and not to the caller.
inline float biospring_imp_force_z(float z, float surface, float transfer, float alip, float alpha, float z0,
                                   float convert)
{
    return -surface * biospring_imp_dcz(z, alpha, z0) * (alip - transfer) * convert;
}

#endif // __BIOSPRING_IMP_SHARED_H__
