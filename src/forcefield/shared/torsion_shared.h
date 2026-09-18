#ifndef __BIOSPRING_TORSION_SHARED_H__
#define __BIOSPRING_TORSION_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU, on the same terms as spring_shared.h:
// compiled twice, once as C++ into biospring-core and once as OpenCL C
// prepended to biospring.cl.
//
// What is here and what is not, and why. A torsion is a four-atom term, so
// unlike the pairwise laws it cannot reduce to one scalar module: its result is
// four vectors. Vector types are the one thing this file may not carry -- a
// Vector3f and a float3 are different objects with different arithmetic -- so
// what is shared is every SCALAR the decomposition needs, and each side does
// the four-line vector assembly itself, spelled identically:
//
//     F1 = n1 * k1
//     F4 = n2 * k4
//     F2 = F1 * (-(c1 + 1)) + F4 * c3
//     F3 = F1 * c1 - F4 * (c3 + 1)
//
// That assembly is where the sign convention lives and it is not guessable: it
// is written in terms of r_ij = r_i - r_j and r_kl = r_k - r_l, which are the
// NEGATIVES of the b1 and b3 the callers build. The other convention passes the
// sum-to-zero test while being wrong by 5 rad/A, so a parity test against the
// CPU is what actually holds this together -- see
// SpringNetworkOpenCL.TorsionsMatchTheCPU.

/// The bin a dihedral angle falls in, for a table of `bins` intervals over
/// [-pi, pi]. The table holds bins + 1 samples, so bin + 1 is always readable.
/// @param pi PI, passed rather than written so there is one definition.
inline int biospring_torsion_bin(float phi, float pi, int bins)
{
    float x = (phi + pi) / (2.0f * pi) * (float)bins;
    if (!(x > 0.0f))
        x = 0.0f;
    int b = (int)x;
    return b > bins - 1 ? bins - 1 : b;
}

/// Where within that bin, in [0, 1).
inline float biospring_torsion_fraction(float phi, float pi, int bins)
{
    float x = (phi + pi) / (2.0f * pi) * (float)bins;
    if (!(x > 0.0f))
        x = 0.0f;
    return x - (float)biospring_torsion_bin(phi, pi, bins);
}

/// The coefficient of n1 in F1.
/// @param scale The force field's spring scale times the kJ.mol-1.A-1 ->
///     Da.A.fs-2 conversion times the torque read from the table.
inline float biospring_torsion_k1(float scale, float b2len, float n1sq)
{
    return -scale * b2len / n1sq;
}

/// The coefficient of n2 in F4. Opposite sign to k1, which is what makes the
/// four forces sum to zero identically once F2 and F3 are built from them.
inline float biospring_torsion_k4(float scale, float b2len, float n2sq)
{
    return scale * b2len / n2sq;
}

/// The projection coefficients that carry F1 and F4 onto the two inner atoms.
/// @param bdotb2   b1.b2 for c1, b3.b2 for c3.
/// @param b2lensq  b2.b2.
inline float biospring_torsion_c(float bdotb2, float b2lensq) { return bdotb2 / b2lensq; }

#endif // __BIOSPRING_TORSION_SHARED_H__
