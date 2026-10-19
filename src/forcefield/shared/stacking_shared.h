#ifndef __BIOSPRING_STACKING_SHARED_H__
#define __BIOSPRING_STACKING_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU -- see spring_shared.h for the rules this
// file obeys and why they exist.
//
// BASE STACKING, after HiRE-RNA v3 (Cragnolini, Laurin, Derreumaux & Pasquali,
// arXiv:1404.0568). Two bases stack when their planes are close to a given
// distance apart, PARALLEL to each other, and FACE TO FACE rather than edge to
// edge. Their equation (4) is
//
//     E_st = eps e^(-(r-r_st)^2/sigma) (ni.nj)^2 (1-|ni x rhat|^4) (1-|nj x rhat|^4)
//
// and the three factors are exactly those three conditions: the Gaussian for
// the distance, (ni.nj)^2 for parallel, and one alignment factor per base for
// face to face.
//
// WHY THIS IS NOT THE HYDROGEN BOND TERM, AND CANNOT DOUBLE COUNT IT. With ni
// and rhat both unit, |ni x rhat|^2 = 1 - (ni.rhat)^2, so writing u = ni.rhat
// the alignment factor collapses to
//
//     A(u) = 1 - (1 - u^2)^2 = u^2 (2 - u^2),
//
// which is 1 when rhat lies ALONG the normal (one base stacked above the other)
// and exactly 0 when rhat lies IN the base plane (two bases side by side). A
// Watson-Crick pair is the second case, so it collects no stacking energy at
// all, and a stack collects no pairing energy because the hydrogen-bond term's
// own angular weight points along a base EDGE. The two terms describe disjoint
// geometries rather than sharing one.
//
// Sign conventions. eps is a positive well DEPTH in kJ.mol-1, like every other
// well in this code, and the energy carries the minus: E = -eps G P Ai Aj, so a
// perfect stack is at -eps and anything else is shallower. Both (ni.nj)^2 and
// A(u) are EVEN in each normal, so flipping a normal changes nothing -- which
// is required, because a base plane has no preferred face for stacking and the
// order the three frame atoms are named in must not be able to matter.

/// Gaussian in the plane-to-plane distance.
/// @param distance    Centre-to-centre distance of the two base frames, in A.
/// @param equilibrium r_st, in A.
/// @param width       sigma, in A^2 -- NOT squared again here, following the
///     reference: the factor is e^(-(r-r_st)^2/sigma), so the half-depth sits
///     at r_st +/- sqrt(sigma ln 2).
/// @return Dimensionless, in (0, 1].
inline float biospring_stacking_distance_factor(float distance, float equilibrium, float width)
{
    const float x = distance - equilibrium;
    return exp(-x * x / width);
}

/// d(factor)/d(distance), in A-1.
inline float biospring_stacking_distance_derivative(float distance, float equilibrium, float width)
{
    const float x = distance - equilibrium;
    return -2.0f * x / width * exp(-x * x / width);
}

/// Face-to-face alignment of ONE base: A(u) = u^2 (2 - u^2), u = normal.rhat.
/// 1 for a base directly above its partner, 0 for two bases edge to edge.
/// @return Dimensionless, in [0, 1] for |u| <= 1.
inline float biospring_stacking_alignment_factor(float u)
{
    const float u2 = u * u;
    return u2 * (2.0f - u2);
}

/// dA/du = 4u(1 - u^2).
inline float biospring_stacking_alignment_derivative(float u)
{
    return 4.0f * u * (1.0f - u * u);
}

#endif // __BIOSPRING_STACKING_SHARED_H__
