"""Ring-construction mathematics shared by the bonded-force-field generators.

Pure geometry and Fourier algebra: nothing here knows about amino acids,
nucleotides, a particular AMBER XML or a particular structure. It takes bond
lengths, valence angles and AMBER torsion terms, and returns the ghost-ring
parameters that reproduce them.

Extracted from generate_bonded_forcefield.py so the nucleic-acid generator
can build on exactly the same construction instead of a copy: this code is
subtle (sign conventions, DC bookkeeping, the comb-filter derivation), and a
duplicate would drift the first time either side is fixed. Verified at
extraction time that the protein .bi.ff comes out byte-identical.

See doc/BondedForceFieldSprings.md for the derivation.
"""

import numpy as np

PHI_GRID_DEG = np.linspace(0.0, 360.0, 3600, endpoint=False)

def closed_form_d2(L, r1, theta1_deg, r2, theta2_deg, delta_deg, phi_deg):
    t1, t2, d, p = (np.radians(x) for x in (theta1_deg, theta2_deg, delta_deg, phi_deg))
    C1 = r1**2 * np.sin(t1)**2 + r2**2 * np.sin(t2)**2 + (r1 * np.cos(t1) - L + r2 * np.cos(t2))**2
    C2 = 2.0 * r1 * r2 * np.sin(t1) * np.sin(t2)
    return C1 - C2 * np.cos(p - d)

def complex_fourier_coeff(energy_over_phi_deg, n):
    # n=0 (DC/mean) needs a plain average, not the 2x-oscillation-amplitude
    # normalization that's correct for n>=1 (cos^2 integrates to 1/2 over a
    # period, a constant's own self-overlap integrates to 1 -- using the
    # n>=1 formula for n=0 silently doubles the DC target/coefficient, an
    # easy-to-miss bug found and fixed during the ghost-particle work, see
    # doc/BondedForceFieldSprings.md).
    phi_rad = np.radians(PHI_GRID_DEG)
    if n == 0:
        return complex(np.mean(energy_over_phi_deg))
    return (2.0 / len(phi_rad)) * np.sum(energy_over_phi_deg * np.exp(-1j * n * phi_rad))

# Deterministic, closed-form ghost-PARTICLE construction (see
# doc/BondedForceFieldSprings.md, "particules fantômes"/ghost-particle
# sections for the full derivation): replaces the earlier ghost-SPRING
# method (fixed real-atom geometry, only k free, sometimes had no
# non-negative-k solution for a multi-term real target) with springs
# between fully-free virtual points, ANCHORED to real atoms only through
# a 3-point placement (see spn::GhostParticle) -- r/theta/delta become
# real degrees of freedom instead of a choice among a few fixed shapes.
#
# For a target harmonic n>=2, a ring of M=N=n ghost-ghost springs (same
# base shape, evenly spaced by 360/n around the axis on each side) cancels
# -- by an exact roots-of-unity argument, verified numerically -- every
# harmonic that isn't a multiple of n, leaving a pure (up to a small,
# quantified leakage into 2n, 3n...) contribution at exactly n. n=1 is a
# special, EXACT case: a single ghost-ghost spring with d0=0 (zero rest
# length) gives a provably pure n=1 harmonic with ZERO leakage into any
# other harmonic (energy = 0.5*k*d(phi)^2, and d(phi)^2 is exactly
# C1-C2*cos(phi-delta) by the closed_form_d2 identity above -- no sqrt
# nonlinearity survives at d0=0) -- this is why n=1 uses mu=0 below, not a
# leak-minimized (rho,mu) like n>=2.
#
# (rho, mu) for n=2,3 were found by maximizing the achievable
# |a_n|/a0 ratio (the harmonic-to-DC-offset efficiency of one ring's own
# base shape) subject to leakage into the first unwanted multiple (2n)
# staying <=0.5% -- a small, one-time 2D grid search over the ABSTRACT
# shape function sqrt(1-rho*cos(x)), independent of any specific axis or
# residue (verified: reused identically across all 70 axes below).
def ring_curve_abstract(L_axis, n, r, theta_deg, d0, k, delta_base, phi_grid_deg, M=None, N=None,
                        r2=None, theta2_deg=None):
    # r2/theta2_deg default to r/theta_deg (the symmetric shape the old
    # free-placement construction used). They are given explicitly since
    # 2026-08-14: a ghost is now the rotated image of a REAL reference
    # atom, so each side's radius/axial offset is that atom's own and the
    # two sides generally differ. See emit_ghost_ring.
    # The 1D Fourier-space model used to CALIBRATE k/delta_base (see
    # calibrate_ring) -- the comb-filter sum over an MxN grid (defaults to
    # M=N=n when not given, the original construction), evaluated over an
    # abstract phi grid. NOT the real 3D placement (see emit_ghost_ring for
    # that, and its header comment for the sign convention linking the two
    # -- verified by direct 3D simulation against a real axis before being
    # trusted here, see doc/BondedForceFieldSprings.md).
    #
    # M,N need not equal n individually -- only LCM(M,N)=n matters (found
    # 2026-08-08/09: the M=N=n choice used everywhere until now is not a
    # mathematical requirement, just the simplest uniform convention. The
    # general comb-filter derivation -- S(phi) = sum_i sum_j f(phi-delta_base
    # -beta_i+gamma_j) with beta_i=i*360/M, gamma_j=j*360/N -- shows the
    # surviving harmonics are exactly the common multiples of M and N, i.e.
    # multiples of LCM(M,N), regardless of M,N individually; verified
    # numerically that an asymmetric M=n,N=1 ring reproduces a real target
    # to amplitude ratio 0.9998 / phase diff 0.000 deg, identical to M=N=n,
    # while using fewer ghosts/springs -- see doc/BondedForceFieldSprings.md).
    if M is None:
        M = n
    if N is None:
        N = n
    if r2 is None:
        r2 = r
    if theta2_deg is None:
        theta2_deg = theta_deg
    e = np.zeros_like(phi_grid_deg)
    for i in range(M):
        beta = i * 360.0 / M
        for j in range(N):
            gamma = j * 360.0 / N
            delta_ij = delta_base + beta - gamma
            d2v = closed_form_d2(L_axis, r, theta_deg, r2, theta2_deg, delta_ij, phi_grid_deg)
            d = np.sqrt(np.maximum(d2v, 0.0))
            e += 0.5 * k * (d - d0) ** 2
    return e

def choose_d0(L_axis, n, r1, theta1_deg, r2, theta2_deg, M, N):
    """Ghost-ghost equilibrium length: the RMS ghost separation sqrt(C1),
    which is the middle of the range the pair can actually reach.

    d0 does NOT trade off against harmonic purity, and a search over it is
    degenerate -- provably so. On an M=n,N=1 ring the sum of d^2 over the
    comb is exactly n*C1 (the comb kills the fundamental for n>=2), so it
    contributes only a constant; every harmonic comes from the -2*d0*sum(d)
    term alone and is therefore proportional to d0. All harmonic RATIOS are
    d0-independent -- verified numerically too: leakage stays at 0.1108% for
    every mu from 0.05 to 1.6 on a real chi1 geometry. d0 sets amplitude
    (hence k) and DC, nothing else.

    An earlier version of this function minimised leakage over d0 and, with
    nothing to discriminate, wandered to d0=4.29 A for a pair that never
    gets past 3.05 A -- a permanently stretched spring, which showed up as
    a NEGATIVE total dihedral energy and a 10x kinetic-energy blow-up on
    GKinase. sqrt(C1) is what the old RING_SHAPES table used as well (mu=1.0
    for n=2,3), for the same reason.

    n=1 keeps d0=0 exactly -- an algebraic identity, not a compromise: with
    d0=0 the energy is 0.5*k*d^2 = 0.5*k*(C1 - C2*cos(phi-delta)), a pure
    first harmonic with EXACTLY zero power anywhere else."""
    if n == 1:
        return 0.0
    c1 = (r1 * np.sin(np.radians(theta1_deg))) ** 2 + (r2 * np.sin(np.radians(theta2_deg))) ** 2 + (
        r1 * np.cos(np.radians(theta1_deg)) - L_axis + r2 * np.cos(np.radians(theta2_deg))) ** 2
    return float(np.sqrt(c1))

def calibrate_ring(L_axis, n, r, theta_deg, d0, target_complex, M=None, N=None, r2=None, theta2_deg=None):
    # Energy is exactly linear in k at fixed geometry (see
    # doc/BondedForceFieldSprings.md, step 3 of the closed-form derivation)
    # -- one evaluation at k=1 gives the ring's own harmonic-n Fourier
    # coefficient (magnitude AND phase), from which k (magnitude match) and
    # delta_base (phase match) solve directly, no iteration. M,N (default
    # n,n): see ring_curve_abstract's own docstring -- the phase/k formulas
    # below depend only on n (=LCM(M,N)), verified algebraically that
    # rotating delta_base by D shifts the ring's n-th harmonic phase by
    # -nD regardless of M,N individually, so nothing below needs to change
    # when M != N, only the curve construction they're evaluated against.
    e0 = ring_curve_abstract(L_axis, n, r, theta_deg, d0, 1.0, 0.0, PHI_GRID_DEG, M, N, r2, theta2_deg)
    z0 = complex_fourier_coeff(e0, n)
    phase_needed = np.angle(target_complex) - np.angle(z0)
    # Rotating delta_base by D shifts the ring's own harmonic-n phase by
    # -n*D, not +n*D (verified numerically -- an easy sign error, since the
    # naive expectation is +n*D). Using the wrong sign gives the complex
    # CONJUGATE of the target: invisible for a purely-real target (chi1's
    # single-term case) but ~45-60% wrong for any axis with a genuinely
    # complex target (phi/psi, chi2's 3-term case).
    delta_base_deg = np.degrees(-phase_needed / n)
    k = abs(target_complex) / abs(z0)
    # The ring's own n=0 (DC/mean) Fourier component at the CALIBRATED k --
    # exactly k times its value at k=1 (energy is linear in k, see above).
    # This is a pure construction artifact (a sum of squares can only ever
    # ADD a positive baseline, see RING_SHAPES's own header comment): it has
    # nothing to do with the real AMBER torsion's own physical mean
    # (target[0] in combined_target_for_axis, handled separately by the
    # caller) and must be tracked so it can be subtracted back out when
    # reporting an absolute energy (never affects forces -- a constant has
    # zero gradient).
    dc_ring = k * float(np.real(complex_fourier_coeff(e0, 0)))
    # The ring's own MINIMUM at the calibrated k, same linearity. The caller
    # needs it whenever "the artifact is a positive baseline added on top of
    # AMBER's mean" stops holding -- see emit_ghost_ring's dc_align. delta_base
    # only rotates the curve, so the set of values (hence the minimum) is the
    # same as at delta_base=0.
    dc_min = k * float(np.min(e0))
    return k, delta_base_deg, dc_ring, dc_min

