#ifndef __BIOSPRING_SPRING_SHARED_H__
#define __BIOSPRING_SPRING_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU. This file is compiled twice: once as C++
// as part of biospring-core, and once as OpenCL C, prepended to biospring.cl
// when the kernel source is embedded. Both sides therefore run the same text,
// not two transcriptions of the same paper formula.
//
// That distinction is not academic. The OpenCL spring kernel used to carry its
// own copy of this expression and was passed only the unit conversion, never
// the force field's spring scale; with spring.scale = 25 the GPU pulled 25
// times too weakly, and nothing in the build could notice.
//
// The rules that keep this file compilable by both toolchains:
//   - no namespaces, no classes, no templates, no references;
//   - no #include, and nothing from the C or C++ standard library;
//   - scalars in, scalar out -- vector types and memory layout differ between
//     the two sides and stay on their own side of the line;
//   - float literals carry the f suffix, so nothing is silently promoted to
//     double (which an OpenCL device is not required to support at all);
//   - names carry a biospring_ prefix, since there is no namespace to put
//     them in.
//
// Anything that needs a physical constant takes it as a PARAMETER rather than
// reading a global. The host owns the constants and passes them down, so there
// is exactly one definition of each and no way for the two sides to disagree.

/// Spring energy, unscaled, for E = 0.5.k.dr^2.
/// @param distance    Distance between the two particles, in A.
/// @param stiffness   Spring constant, in kJ.mol-1.A-2.
/// @param equilibrium Rest length, in A.
/// @return Energy in kJ.mol-1, before the force field's spring scale.
inline float biospring_spring_energy(float distance, float stiffness, float equilibrium)
{
    float displacement = distance - equilibrium;
    return 0.5f * stiffness * displacement * displacement;
}

/// Spring force module, the derivative of the above.
/// @param scale Everything the caller wants folded in: the force field's
///     spring scale times the kJ.mol-1.A-1 -> Da.A.fs-2 conversion. Passing it
///     rather than reading a constant is what keeps the CPU and the GPU from
///     drifting apart.
/// @return Force module along the axis, in Da.A.fs-2 once scale carries the
///     conversion.
inline float biospring_spring_force_module(float distance, float stiffness, float equilibrium, float scale)
{
    // Kept in this order -- the raw module first, the scale last -- because
    // that is the order the CPU used before this file existed, and floating
    // point is not associative.
    float force_module = stiffness * (distance - equilibrium);
    return force_module * scale;
}

#endif // __BIOSPRING_SPRING_SHARED_H__
