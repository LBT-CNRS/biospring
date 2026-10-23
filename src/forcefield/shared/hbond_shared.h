#ifndef __BIOSPRING_HBOND_SHARED_H__
#define __BIOSPRING_HBOND_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU -- see spring_shared.h for the rules this
// file obeys and why they exist.
//
// A hydrogen bond here is a Morse well in the donor-acceptor heavy-atom
// distance, weighted at BOTH ends by how well the bond lines up with the
// direction each side offers: the donor's hydrogen, standing in as
// "away from the antecedent", and the acceptor's lone pairs, standing in the
// same way. The physics behind the parameters is documented once, in
// forcefield/energy/hydrogenbond.hpp, which now calls straight into this file
// so there is one text rather than two.

/// Morse well between a donor and an acceptor heavy atom, no explicit
/// hydrogen.
/// @param distance    Donor-acceptor distance, in A.
/// @param welldepth   De, in kJ.mol-1 (already molar).
/// @param equilibrium re, in A.
/// @param width       a, in A-1.
/// @return Energy in kJ.mol-1: 0 at infinity, -De at equilibrium.
inline float biospring_hbond_energy(float distance, float welldepth, float equilibrium, float width)
{
    float u = exp(-width * (distance - equilibrium));
    float one_minus_u = 1.0f - u;
    return welldepth * one_minus_u * one_minus_u - welldepth;
}

/// dV/dr of the well above, in the "positive = attractive, on the FIRST
/// particle" convention the rest of the code uses.
/// @param convert Folded in by the caller rather than read from a global:
///     GLOBAL_SPRING_FORCE_CONVERT, which turns kJ.mol-1.A-1 into the
///     integrator's Da.A.fs-2. The well depth is molar, like a spring
///     stiffness, so it needs the same conversion.
/// @return Force module along the donor->acceptor axis.
inline float biospring_hbond_force_module(float distance, float welldepth, float equilibrium, float width,
                                          float convert)
{
    float u = exp(-width * (distance - equilibrium));
    float force_module = 2.0f * width * welldepth * u * (1.0f - u);
    return force_module * convert;
}

// ---- Angular weight -----------------------------------------------------
//
// w = max(cos theta, 0)^2. Clamped at zero rather than squared through: an
// acceptor sitting BEHIND the donor is not a weak hydrogen bond, it is not one
// at all, and cos^2 alone would revive it with the sign lost.

/// @param cos_theta Cosine of the angle between the side's own direction and
///     the bond axis.
/// @return Dimensionless weight in [0, 1].
inline float biospring_hbond_angular_factor(float cos_theta)
{
    return cos_theta > 0.0f ? cos_theta * cos_theta : 0.0f;
}

/// d(weight)/d(cos theta), for the force.
inline float biospring_hbond_angular_derivative(float cos_theta)
{
    return cos_theta > 0.0f ? 2.0f * cos_theta : 0.0f;
}


// ---- Isostericity -------------------------------------------------------
//
// A REDUCTION, in [0, 1], on a pair of nucleotides whose base pair does not sit
// in the geometry a Watson-Crick pair occupies with respect to its own
// backbone. Like the well's angular weights it multiplies a depth that is
// already calibrated, so it can only take energy away: a real pair pays almost
// nothing and nothing can score above 1.
//
// WHAT IT IS FOR. The pairing term scores a bond by a distance and two angles,
// and all three are EXACTLY preserved when the partner nucleotide is turned 180
// degrees about the donor-acceptor axis -- both bonded atoms lie on that axis,
// and rotating an antecedent about it preserves the angle it makes with it. The
// measured degeneracy is -0.00 +/- 0.01 kJ.mol-1 over 49 bonds. So the term
// cannot tell a Watson-Crick pair from its mirror image, and that is why no
// pairing forms de novo in the right register.
//
// WHY THIS LIFTS IT. The flip moves the partner's C1' to the other side, which
// the base atoms alone never reveal. Measured on 8PKL, over 1121 strict
// Watson-Crick pairs (A.U or G.C, not sequence neighbours, their two central
// ring nitrogens within 3.2 A) and the same 1121 pairs with the partner flipped:
//
//                      C1'-C1' (A)       lambda (deg)
//   native             10.46 +/- 0.59    54.9 +/- 5.6
//   mirror             10.79 +/- 0.78    36.6 +/- 14.0
//
// The mirror KEEPS the C1'-C1' distance and destroys lambda. So lambda carries
// the discrimination and the distance only guards the helix width; both are
// kept because a pair that is isosteric must satisfy the two together, which is
// what isostericity means.
//
// lambda is the angle at one nucleotide between the C1'->C1' vector and its own
// glycosidic bond, C1'->N9 for a purine and C1'->N1 for a pyrimidine.
//
// WHAT IS BORROWED AND WHAT IS NOT. Isostericity is Leontis and Westhof's: base
// pairs are grouped into families that can substitute for one another in a
// motif without deforming it, and the C1'-C1' distance and lambda are their
// descriptors of a pair's geometry relative to its backbone (Leontis, Stombaugh
// and Westhof, NAR 30:3497, 2002). Their published METRIC, the IsoDiscrepancy
// Index (Stombaugh, Zirbel, Westhof and Leontis, NAR 37:2294, 2009), is NOT
// what is computed here and does not use lambda: it is built from the two
// pairs' C1'-C1' difference, from inter-pair C1'-C1' distances taken after
// SUPERPOSING one base of each, and from rotation angles about the normal to
// the pair plane, with thresholds 2.0 A isosteric, 3.3 A non-isosteric.
//
// The IDI cannot be a force-field term: it compares a candidate pair to an
// EXEMPLAR pair, so each evaluation needs a superposition, and this function
// sits in the innermost loop of the non-bonded. What is implemented is a
// two-descriptor reduction that keeps what the mirror breaks. The agreement is
// a real one as far as it goes -- the IDI's own "Case 2", the one with the
// larger angular coefficient, is exactly the case where a 180 degree flip is
// needed -- but it is an agreement in spirit, not the published number.
//
// So: the two quantities are theirs, the Watson-Crick centres are theirs and
// mine agree with them (they give about 10.5 A for the C1'-C1' of a WC pair,
// measured here at 10.46), and the Gaussian form and its two widths are MINE,
// chosen by the scan below. Anyone wanting the published classification should
// use FR3D and the Basepair Catalog, not this.
//
// THE OTHER ROUTE, and why this one. HiRE-RNA v3 breaks the same degeneracy
// with a published mechanism: its nu(alpha) = cos^6(alpha - alpha0) carries a
// per-pair-type alpha0, and a torsion tau flips alpha's sign so the cosine's
// two minima at +alpha and -alpha stop being equivalent (Cragnolini, Laurin,
// Derreumaux and Pasquali, arXiv:1404.0568, equations 6 to 8). That is the
// literature term and it is not what is here. Porting it needs the 22 per-pair
// geometries, which that paper does NOT tabulate -- it cites Lemieux and Major,
// NAR 30:4250, 2002 for the repertoire and obtains its own numbers from a
// genetic algorithm against decoys.
//
// The paper's own words on how binding its list is: "the choice of 22
// interactions is rather arbitrary and can be extended to any number of
// interactions as long as they are sufficiently distinct in interaction
// centers".
//
// WHAT THIS REPRODUCES AND WHAT IT DOES NOT, stated so no one has to find out
// by being surprised. alpha0 and tau together distinguish a pair from its
// mirror, and so does this, measured at 18-fold. But alpha0 is per pair TYPE,
// so the table also says WHICH family a pair belongs to and can reward a
// cis-Watson-Crick A.C at its own geometry. This factor cannot: it knows only
// how far a pair sits from the Watson-Crick geometry, so it is a Watson-Crick
// PREFERENCE, not a family classifier, and it attenuates every non-canonical
// pair to about a third whatever its family. Where that is wrong, it is the
// table that is needed, not a wider sigma.
//
// THE PARAMETERS ARE MEASURED, and the whole table is here so the choice is
// auditable. Centres are the native means above. The widths were scanned
// against the native, mirror and non-canonical populations at once:
//
//   sigma_d  sigma_lambda   native   mirror   non-canonical   discrimination
//    1.2 A      22 deg       85.6%     9.9%       36.1%            9x
//    1.5 A      18 deg       84.3%     4.7%       33.9%           18x   <- used
//    1.8 A      18 deg       85.0%     5.1%       35.6%           17x
//    2.2 A      22 deg       88.5%    12.4%       43.7%            7x
//
// A real pair keeps 84 % of its well, which is the attenuation the planarity
// factor was accepted at (87 %), and a mirror keeps 4.7 %. A NON-CANONICAL pair
// keeps 34 %: attenuated and not forbidden, which is deliberate -- the U32.A38
// that closes every tRNA's anticodon loop is a real, multi-bond, non-canonical
// pair, and a factor that killed it would unmake the loop it holds.
//
// ONLY WHERE BOTH SIDES DECLARE A GLYCOSIDIC BOND. A protein site has no C1'
// and no glycosidic nitrogen, so the factor is not defined there and must leave
// such a bond untouched to the bit, not merely within a tolerance. The caller
// is responsible for that test; this function is only the shape.

/// Gaussian acceptance in one coordinate, and its derivative.
/// @param value  The measured quantity.
/// @param centre Its value in a Watson-Crick pair.
/// @param width  The standard deviation of the acceptance, same unit.
/// @return Dimensionless weight in (0, 1], exactly 1 at the centre.
inline float biospring_isostericity_gaussian(float value, float centre, float width)
{
    float z = (value - centre) / width;
    return exp(-z * z);
}

/// d(weight)/d(value) of the above, for the force.
inline float biospring_isostericity_gaussian_derivative(float value, float centre, float width)
{
    float z = (value - centre) / width;
    return -2.0f * z / width * exp(-z * z);
}

/// The factor itself: the product of the three acceptances.
/// @param distance  C1'-C1' separation of the two nucleotides, in A.
/// @param lambda_i  Angle at the first nucleotide, in DEGREES.
/// @param lambda_j  Angle at the second, in degrees.
/// @return Dimensionless reduction in (0, 1].
inline float biospring_hbond_isostericity_factor(float distance, float lambda_i, float lambda_j,
                                                 float distance_centre, float lambda_centre,
                                                 float distance_width, float lambda_width)
{
    return biospring_isostericity_gaussian(distance, distance_centre, distance_width)
         * biospring_isostericity_gaussian(lambda_i, lambda_centre, lambda_width)
         * biospring_isostericity_gaussian(lambda_j, lambda_centre, lambda_width);
}

#endif // __BIOSPRING_HBOND_SHARED_H__
