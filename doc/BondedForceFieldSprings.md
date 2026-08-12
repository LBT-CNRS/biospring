A Spring-Only Mechanical Model for AMBER's Bonded Interactions
================================================================

Abstract
--------

BioSpring represents molecular mechanics entirely through springs: particles
connected by distance restraints of the form `E = 0.5 * k * (r - r0)^2`. This
document derives, stage by stage, how a model built purely from such springs
can reproduce the mechanical behaviour of AMBER's four bonded interaction
terms — bond stretching, angle bending, proper dihedral torsion, and improper
dihedral torsion — starting from the simplest possible baseline (a uniform
rigid-body mesh) and layering real, spring-only parameters on top of it.
Bending and every dihedral family are represented using massless virtual
sites (`spn::GhostParticle`) rather than springs between real atoms directly
— both a real 1-3 angle and a real dihedral axis, represented on real atoms
alone, were found to leak unrelated real motion (bond stretch for bending, a
few degrees of ordinary thermal noise for backbone dihedrals) that a virtual
site's free, calibration-time-only geometry avoids by construction. Every
derivation is closed-form or an exact algebraic identity, verified either
symbolically or numerically (to machine precision where applicable); no step
is a curve fit. A microbenchmark comparing per-evaluation computational cost
against conventional trigonometric bonded force routines is included, along
with the full validation strategy used to check the model against
independent recomputation and real AMBER parameters — including, critically,
comparison against real molecular-dynamics trajectories rather than only
self-consistency at generation time, which is what actually found the two
leakage problems above.

1. Introduction
----------------

A production molecular force field (AMBER, CHARMM, and others) separates
interatomic interactions into *bonded* terms — computed only along the real
chemical bond graph (1-2, 1-3, 1-4 neighbours) — and *non-bonded* terms —
electrostatics and van der Waals, computed between all sufficiently distant
pairs. BioSpring already handles the non-bonded terms through a separate,
pre-existing mechanism outside springs; this document is concerned
exclusively with the four bonded terms:

| # | Term | What it restrains | AMBER functional form |
|---|------|--------------------|--------------------------|
| 1 | Bond stretching | distance between two bonded atoms | `0.5 k (r-r0)^2` |
| 2 | Angle bending | a 1-2-3 valence angle | `0.5 K_theta (theta-theta0)^2` |
| 3 | Proper dihedral | rotation about a real bond | `k (1 + cos(n*phi - phase))` |
| 4 | Improper dihedral | planarity of a trigonal centre | `k (1 + cos(n*phi - phase))` |

A spring can only restrain a *distance*. Term 1 needs no transformation at
all. Term 2 needs a geometric substitute (a 1-3 distance can stand in for an
angle) plus an argument for *why* that substitute is quantitatively correct,
not just qualitatively similar. Terms 3 and 4 are harder still: a dihedral
angle is not the distance between any single pair of atoms, so no direct
geometric substitute exists at all — a more elaborate construction ("ghost
springs", introduced in Section 3.2) is required.

BioSpring's existing, independent mechanism for approximating molecular
rigidity — `--rigidbody`, a uniform all-pairwise spring mesh within named
atom groups — is the natural point of departure: it already gives the right
*topology* of degrees of freedom (which axes rotate freely, which groups are
rigid), just not the right *numbers*. Section 3 below introduces the four
bonded terms as increasingly complete refinements layered on top of this
baseline, each stage motivated by exactly what the previous stage still
lacks. `--rigidbody` remains the permanent substrate throughout — these
refinements retune and extend its springs, they do not aim to replace it.

2. Baseline: the rigid-body model
-------------------------------------

`--rigidbody` builds, for every named group of a `.rbody` rule file (e.g. a
side-chain's chi-pivot group, or a backbone peptide-plane group), a spring
between *every pair* of atoms in that group, all at one uniform, non-physical
stiffness, with equilibrium distances read directly from the input
structure.

This has one purely geometric but powerful consequence: a complete graph of
distance constraints among four or more non-coplanar points removes every
internal degree of freedom of that set of points, leaving only the six
rigid-body degrees of freedom (three translations, three rotations) of the
group as a whole. Two adjacent groups that share exactly two atoms (the
hinge axis) are therefore rigid individually, but free to rotate relative to
each other about that shared axis — nothing more is needed for this
"two anchors, one free rotation" behaviour to appear.

This baseline already gives BioSpring the right qualitative mechanical
picture — rigid units connected by hinges at exactly the real rotatable
bonds — using nothing but springs. What it does not give is:

* **Real stiffness or equilibrium values** — every spring, whether it stands
  in for a real 1-2 bond, a 1-3 angle, or an incidental "extra" rigidifying
  pair with no direct chemical meaning, carries the same uniform
  `--stiffness`.
* **Any preferred rotation angle** — a hinge is completely free; there is no
  notion of a rotamer well or a preferred dihedral value.
* **A behaviour distinction between "genuinely rigid" and "genuinely free"**
  — a ring, for instance, is made fully rigid (no pucker at all) purely as a
  side effect of the all-pairwise construction, whether or not the real
  molecule is actually that stiff there.

Sections 3.1–3.3 remove these three limitations one at a time, always by
adding or retuning springs — never by introducing any other kind of
mechanism.

3. Spring-network transformation of the four bonded terms
--------------------------------------------------------------

### 3.1 Stage 1 — real stretching and bending, layered on the rigid-body mesh

The first refinement retunes the rigid-body mesh's *existing* springs with
real AMBER stiffness and equilibrium values, without changing which pairs of
atoms have a spring at all. `--rigidbody` has already created a spring for
every real 1-2 bond and every real 1-3 angle pair (they are simply two of
the many pairs inside a rigid group) — Stage 1 finds those specific springs
and replaces their uniform parameters with real ones.

**Term 1 — bond stretching (no transformation needed).** AMBER's
`HarmonicBondForce`, `E = 0.5 k (r-r0)^2`, already has the exact functional
form of a spring, in the same units and the same `0.5`-prefactor convention
BioSpring uses throughout (verified empirically against OpenMM's XML
tables — see the unit-convention note at the end of this section). `r0` and
`k` are copied directly: a stretching spring does not approximate a bond, it
*is* the bond, expressed in BioSpring's native representation.

**Term 2 — angle bending.** The equilibrium/stiffness conversion below (law
of cosines + curvature matching) is exact at `theta0` regardless of which
two points the resulting spring connects — the derivation was originally
applied directly to the two *real* 1-3 atoms, but that specific choice
turned out to leak an unrelated error into the model, described after the
derivation.

For a real 1-2-3 valence angle with vertex atom 2, the distance between the
two *outer* atoms, `r13`, is a monotonic function of `theta` for fixed real
bond lengths `r12`, `r23` (already correct, from Term 1):

    r13(theta)^2 = r12^2 + r23^2 - 2 * r12 * r23 * cos(theta)          (law of cosines)

At `theta = theta0` this fixes the spring's equilibrium directly:
`r13_0^2 = r12^2 + r23^2 - 2 r12 r23 cos(theta0)`.

The stiffness is fixed by requiring the *curvature* of the spring's energy
with respect to `theta` to match the real angle term's curvature at
`theta0` — both energies already agree in value and slope there (both are
stationary at equilibrium), so matching the second derivative as well is the
natural next-order equivalence, not a fit. Since
`E_spring(theta) = 0.5 k13 (r13(theta) - r13_0)^2`:

    dE_spring/dtheta = k13 (r13(theta)-r13_0) dr13/dtheta        -> 0 at theta0 (as required)
    d^2E_spring/dtheta^2 |_theta0 = k13 (dr13/dtheta |_theta0)^2  (the (r13-r13_0) term vanishes at theta0)

Differentiating the law of cosines, `dr13/dtheta = r12 r23 sin(theta) / r13`,
so at `theta0`: `dr13/dtheta|_theta0 = r12 r23 sin(theta0) / r13_0`. Matching
this curvature to `K_theta` (the real angle term's own, constant, curvature)
and solving for `k13`:

    k13 = K_theta * r13_0^2 / (r12 * r23 * sin(theta0))^2

This is a direct algebraic solve — no iteration, no fit. It is exact at
`theta0` (value, slope, curvature all match) and only approximate away from
it, since `r13(theta)` is not linear in `theta`; the deviation grows with
distance from equilibrium, the same category of approximation used
throughout the rest of this document.

**Bond-stretch leakage (found and fixed): why the spring connects two
*ghost* particles, not the real 1-3 atoms.** Connecting this spring
directly between the two real outer atoms initially seemed like the
obvious choice — but `r13`, the real 1-3 *distance*, depends on BOTH the
real angle `theta` AND the two adjacent real bond lengths `r12`/`r23`,
which independently flex via their own Term-1 springs during real
dynamics. A distance spring on the real atoms cannot distinguish "the
angle changed" from "a bond stretched" — both look identical as a change
in `r13`. Quantified by direct comparison against an independent OpenMM
(real MD frame) reference: **~60% systematic energy excess**, of which
~95% (59 of the 60 percentage points) is bond-stretch leakage and only
~5% is the genuine curvature-matching linearization error derived above.
The sensitivity `d(r13)/d(r12) = (r12 - r23*cos(theta))/r13` is not small
for realistic (obtuse, ~110-120deg) protein bond angles — confirmed both
analytically and by finite difference on two independent real examples
(0.81 and 0.87) — so bond stretch transmits almost 1:1 into the measured
`r13`.

Fix: for each real angle (vertex B, neighbours A, C), the spring connects
two massless virtual sites (`spn::GhostParticle`, Section 3.2 introduces
the mechanism in full) instead: `ghost_A = B + r0(bond BA) * normalize(A -
B)`, `ghost_C = B + r0(bond BC) * normalize(C - B)` — each rescaled to the
AMBER-ideal bond length but pointing along the REAL, current bond
direction. This reproduces the bonds-fixed synthetic case exactly during
real dynamics (leaving only the ~5% linearization error above), because
the two ghosts' distance from B is now pinned to the ideal bond length by
construction, immune to the real bonds' own independent stretching.
Force redistribution reuses the ghost-particle mechanism's own transpose-
Jacobian construction unchanged. **Validated end-to-end**: real AMBER
reference (Fs-peptide, 20 real trajectory frames), stretch+bend combined,
**mean relative error 1.03%, max 2.11%** — matching the theoretically
expected ~1-5% pure-linearization residual once the leakage is removed.

The real 1-2/1-3 pair that Term 1/Term 2 used to connect directly is never
removed once superseded this way (needed for nonbonded-exclusion
bookkeeping) — its stiffness is zeroed instead, and the real restraint
lives entirely in the new ghost-anchored spring, in its own dedicated
spring collection and energy channel (see the note on energy channels
below).

*Ring/planar-group caveat.* A rigid-body group with every pairwise distance
constrained (Section 2) has no internal degrees of freedom left at all —
this is stiffer than reality even after Stage 1, since bending alone does
not restore ring pucker or out-of-plane flexibility. Real planarity is a
distinct AMBER term (Term 4, Section 3.2) — bending is not meant to provide
it.

**Energy channels.** STRETCH, BEND, and DIHEDRAL each report their own
separate energy (`Stretch energy`, `Bend energy`, `Dihedral energy`), in
their own dedicated spring collection, independent of `--rigidbody`'s own
mesh. `Spring energy` was redefined to mean *only* that untouched
rigid-body/ENM baseline — it no longer includes the retuned/ghost-anchored
contributions, since "spring energy" as a single lumped number only ever
meant something for a plain elastic-network model in the first place.

*Unit convention (applies to every stage below).* All figures in this
document are sourced from `amber99sb.xml`'s `HarmonicBondForce`,
`HarmonicAngleForce` and `PeriodicTorsionForce` tables (bundled with the
`openmm` Python package), which already use BioSpring's own
`0.5`-prefactor convention (verified: e.g. a `k = 259408.0` equals
`2 * 4.184 * 31000` once converted from kcal — i.e. already in `0.5`-prefactor
form). No extra factor of 2 is applied anywhere — that would only be needed
if sourcing a raw AMBER `parm*.dat` file directly, which is not the case
here.

### 3.2 Stage 2 — dihedral torsion (proper and improper), added as an intermediate layer

Stage 1 retunes/extends springs anchored on real, existing atoms. Dihedral
torsion needs something categorically different: a dihedral angle `phi` is
not the distance between any single pair of real atoms, and — as Section
3.1's bending fix already found the hard way — even a real 1-3 *angle* leaks
unrelated bond-stretch noise when represented by a spring between the real
outer atoms. Both problems share the same fix: place the spring between
massless **virtual sites** (`spn::GhostParticle`) instead of real atoms —
points whose position is computed algebraically from 3 real anchor atoms
each simulation step, never integrated themselves, carrying no mass, and
exerting force back onto their anchors through the placement's own transpose
Jacobian (the virtual-work principle: verified numerically against finite
differences before being trusted, and that using `J^T` rather than `J`
matters — using `J` breaks energy conservation silently). This is the same
"dummy atom"/virtual-site technique used for e.g. TIP4P's M-site in
mainstream MD packages, not a BioSpring-specific invention.

`spn::GhostParticle::computePosition(B, C, Ref, r, theta, delta)` places a
ghost at distance `r`, polar angle `theta` from the B→C axis, azimuth
`delta` from the `Ref` direction, using a local frame built fresh from B/C/Ref's
*current* real positions every step — so a ghost always follows its anchors
correctly as the molecule moves, with `r`/`theta`/`delta` themselves fixed,
calibration-time constants. Crucially, **`r`, `theta`, and `delta` are free
choices**, unlike a spring anchored on two real atoms, whose geometry is
whatever the real chemistry happens to give. This one difference is what
separates the current mechanism from an earlier, abandoned "ghost-spring on
real atoms directly" attempt: that method's real substituent geometry
sometimes had no way to reproduce a real multi-term AMBER target without a
negative (unphysical) stiffness, or reproduced only 1 of several needed
harmonics — a ghost particle's free geometry removes that constraint
entirely (Steps 1-4 below).

This intermediate configuration — rigid-body mesh, Stage 1's real
stretching/bending, and now real dihedral wells layered on top — is already a
complete, independently useful, and independently testable model: bond and
angle vibration have real values, and side-chain/backbone rotation now has
real energetic preferences, while the underlying topology (which pairs have
springs at all) still traces back to the rigid-body mesh of Section 2.

In practice (`pdb2spn`), Stage 1 and Stage 2 are independent opt-in flags
on top of `-rigidbody -bondedinteraction <file>`: `-stretching`, `-bending`
(which requires `-stretching` too, since its conversion needs real bond
lengths), and, for dihedral, `-dihedralbackbone`/`-dihedralsidechain`
(one per proper-dihedral family — planarity has no flag yet, see Known
Limitations) plus `-dihedral` as a pure convenience alias for both
together. None is implied by the others or by `-bondedinteraction` alone —
combining only a dihedral flag, without `-stretching`/`-bending`, is
exactly this section's intermediate model. Once built, dihedral ghost
springs are ordinary springs at simulation time, gated by the same
`spring.enable` master switch as everything else in the network (the same
is true of stretch/bend's own ghost-anchored springs — see the
energy-channel note in Section 3.1) — on top of that, each proper family
also has its own independent runtime debug toggle
(`dihedralphi.enable`/`dihedralpsi.enable`/`dihedralomega.enable`/
`dihedralchi.enable`, plus `bending.enable`), all defaulting to enabled so
an `.msp` written before these existed keeps the same behaviour. There is
deliberately no equivalent toggle for stretching: a STRETCH spring shares
its real atom pair with the existing (zeroed) `-rigidbody` spring, so an
independent runtime switch would need a live link back to that sibling
spring, which does not exist today.

Both dihedral families — **proper** (side-chain chi1-4, backbone phi/psi,
the cross-residue peptide-bond torsion omega) and **improper** (planarity
of aromatic rings, guanidinium — still not implemented, see below) — use
exactly the same construction below. AMBER itself treats an improper
geometrically as an ordinary four-atom dihedral (the "axis" is the bond
between the 2nd and 3rd atom of its atom list, the other two are the
substituents); only the source table (`PeriodicTorsionForce`'s `<Proper>`
vs `<Improper>` entries) and the atoms involved differ.

#### Notation

For a dihedral defined by four points in a chain **R-B-C-S** (B and C
connected by a real bond, the axis; R attached to B; S attached to C) —
R and S may each be either a real atom (the improper-planarity case, not
yet implemented) or a ghost particle placed off B/C plus a third real
reference atom (every proper axis, as implemented):

* `L` — the B-C bond length.
* `r1`, `theta1` — R's distance/polar-angle from B (free for a ghost).
* `r2`, `theta2` — S's distance/polar-angle from C (free for a ghost).
* `phi` — the dihedral angle R-B-C-S.

#### Step 1 — closed-form distance between R and S as a function of phi

**Claim.** With B at the origin and C on the z-axis:

    d^2(phi) = C1 - C2 * cos(phi - delta)
    C1 = r1^2 sin^2(theta1) + r2^2 sin^2(theta2) + (r1 cos(theta1) - L + r2 cos(theta2))^2
    C2 = 2 * r1 * r2 * sin(theta1) * sin(theta2)

where `delta` is a fixed azimuthal offset depending only on how `phi`'s zero
reference is defined. This identity holds regardless of whether R/S are real
atoms or ghost particles — it is pure geometry, not a statement about
chemistry.

**Derivation.** Place B at the origin, C at `(0,0,L)`, R in the xz-plane
(fixing `phi`'s zero at R's azimuth): `R = (r1 sin(theta1), 0, r1 cos(theta1))`.
S is attached to C, at distance `r2`, angle `theta2` from the C→B direction,
azimuth `phi` relative to R:
`S = (r2 sin(theta2) cos(phi), r2 sin(theta2) sin(phi), L - r2 cos(theta2))`.
Expanding `d^2 = |S-R|^2` and using `cos^2+sin^2=1` to collapse the cross
terms gives exactly the formula above.

**Verification.** Implemented both ways — this closed form, and a direct 3D
vector construction taking the literal distance — for a representative sp3
carbon-chain geometry (`L=1.526`, `r1=1.526`, `theta1=109.5deg`, `r2=1.09`,
`theta2=109.5deg`), sampled at 37 values of `phi` across a full turn:
**maximum absolute difference 3.55e-15 A^2** — floating-point noise, not
model error. Separately verified that this same abstract formula (used only
to *calibrate* a ring, see Step 3) matches `GhostParticle::computePosition`'s
own literal 3D placement, driven by real, moving anchor atoms across a real
trajectory, to the same numerical precision.

A ghost-ghost ring member's energy is then simply `E(phi) = 0.5 k
(d(phi)-d0)^2` — an ordinary BioSpring distance spring, just between two
virtual points instead of two real atoms.

#### Step 2 — which harmonics survive averaging over several copies (Dirac comb)

A single spring's energy is **not** a pure cosine of `phi` (it is a
nonlinear function of `cos(phi)` through the square root and the square), so
it has content at every harmonic, not just one. Averaging several rotated
copies cancels the unwanted ones exactly, for *any* periodic function.

**Claim.** Let `g(phi)` have period `2*pi`, Fourier series
`g(phi)=sum_k ĝ_k exp(ikphi)`. Build `M` copies of the reference point
evenly spaced around B (`2*pi*i/M` apart) and `N` copies of the rotating
point evenly spaced around C (`2*pi*j/N` apart). Summing every `(i,j)`
pair's energy:

    E_total(phi) = sum_i sum_j g(phi + 2*pi*i/M - 2*pi*j/N) = M*N * sum over k that are multiples of LCM(M,N) of ĝ_k exp(ikphi)

**Proof.** Expand and swap summation order:
`E_total = sum_k ĝ_k exp(ikphi) * [sum_i exp(ik2pi*i/M)] * [sum_j exp(-ik2pi*j/N)]`.
Each bracket is a sum of `M` (resp. `N`) evenly spaced roots of unity to the
k-th power — a standard identity, equal to `M` (resp. `N`) if `M` (resp. `N`)
divides `k`, else `0`. Only `k` divisible by both — i.e. by `LCM(M,N)` —
survives, scaled by `M*N`. This holds for any periodic `g` and for either
real or virtual points equally: it is a property of averaging evenly-spaced
copies, not of the spring's specific shape.

**Verification (numerical, FFT of the real nonlinear ring energy):**

| M | N | LCM(M,N) | Harmonics surviving above 1% of peak (k=1..19) |
|---|---|----------|--------------------------------------------------|
| 3 | 3 | 3        | {3} |
| 1 | 2 | 2        | {2} |
| 2 | 3 | 6        | {6} |

Exactly as predicted. **Practical use:** for a ghost-ghost ring targeting a
real AMBER harmonic `n`, any `(M,N)` with `LCM(M,N)=n` works identically (the
ring shape is precomputed once per `n`, independent of axis or residue —
`RING_SHAPES` in `scripts/generate_bonded_forcefield.py` — since a ring's own
harmonic-to-DC efficiency depends only on its shape, not on which axis it
calibrates).

**Optimization: `M=n, N=1` instead of `M=N=n`.** Since only `LCM(M,N)`
matters (not `M` and `N` individually), the cheapest `(M,N)` pair for a given
`n>=2` is `M=n, N=1` (or the symmetric `M=1, N=n`) — one full ring on one
side, a single reference point on the other, instead of `n` points on both
sides. `emit_ghost_ring` uses this (`n=1` stays `M=N=1`, already minimal).
**Verified to give byte-identical results, not just an approximation of the
same quality:** point-by-point diff of the reconstructed energy curve for
every real (resname, axis) combination in the `.bi.ff`, `M=N=n` vs `M=n,N=1`,
max difference 0.000000 kJ/mol across all 90 curves x 91 sampled angles.
Reduces the deployed `.bi.ff` from 984 to 737 `GHOSTPARTICLE` entries
(-25.1%) and 1166 to 492 `DIHEDRAL` springs (-57.8%).

**`n=1` is a special, exact case.** A single ghost-ghost spring with `d0=0`
gives a *provably pure* n=1 harmonic with **zero** leakage into any other
harmonic — algebraically, `d0=0` makes the sqrt nonlinearity in Step 1's
identity vanish entirely (`E=0.5*k*d(phi)^2`, and `d(phi)^2` is already
exactly `C1-C2*cos(phi-delta)`, a pure single cosine). No ring/averaging is
needed for `n=1` at all: `M=N=1`, `d0=0`.

#### Step 3 — solving for stiffness and phase exactly (not an iterative fit)

For fixed ring geometry, a ring's energy is exactly linear in its shared
stiffness `k` (Step 1's identity has no `k`-dependence in `d(phi)` itself).
One evaluation at `k=1` gives the ring's own harmonic-`n` Fourier
coefficient, magnitude *and phase*, from which both the real target's
magnitude (`k = |target| / |achieved_at_k=1|`) and phase (via a free
azimuthal offset `delta_base`, rotating the whole ring rigidly) solve
directly — no iteration. This `delta_base` freedom is exactly what a
ghost-particle ring has that a real-atom-anchored ring does not: real atoms'
azimuth is whatever the real chemistry gives, with no equivalent free
parameter to rotate away a phase mismatch (this is precisely why a
real-atom-only attempt at omega, Section 3.2's later subsection, failed:
a ~169deg unfixable phase gap between the real geometry's own natural
harmonic content and AMBER's target).

*A sign pitfall, found and fixed.* Rotating a ring's base azimuth by `D`
shifts its own harmonic-`n` phase by `-n*D`, not `+n*D` — verified
numerically; using the naive `+n*D` sign gives the complex *conjugate* of
the target. Invisible for a purely real target (most single-term axes), but
~45-60% wrong for any axis with a genuinely complex target (phi/psi/omega,
or chi2's multi-term case).

#### Step 4 — the ring's own DC offset (a construction artefact, corrected exactly)

A ring made entirely of squared terms (`0.5*k*(d-d0)^2 >= 0`) can only ever
ADD a positive baseline to its own mean (DC) — it cannot subtract from it.
This mean has nothing to do with AMBER's own real torsion mean (`target[0]`,
exactly the sum of each matched real term's own `k`, since
`k*(1+cos(...))` has mean `k`): it is a pure side effect of using a
sum-of-squares construction to isolate one harmonic. Since it's an exact,
closed-form-computable constant (`dc_ring`, see `calibrate_ring`), it can be
subtracted back out *exactly* when reporting energy, with zero effect on
forces (a constant has no gradient): stored as a per-spring `dc_offset`
column in the `.bi.ff` `DIHEDRAL` line, summed by `BondedForceFieldReader`
across every entry actually applied and subtracted only in
`SpringNetwork::getDihedralEnergy()`. Multiple rings/harmonics on the same
axis attribute the real target's own DC to exactly the first ring emitted,
so it is subtracted exactly once per axis, not once per ring.

#### Step 5 — decoupling the ring from the axis bond length (`RING_SCALE`)

A ghost ring's energy is **not** a pure function of the torsion angle, unlike
the AMBER term it reproduces. With ghost planes separated by `Lam` and ring
radius `rho_r`, `d^2 = Lam^2 + 2*rho_r^2*(1-cos psi)`, so `d(d^2)/dL = 2*Lam`:
it also depends on the real distance `L` between the two axis atoms. Expanding
`E = sum_i 0.5k(d_i-d0)^2` over a ring of `M>=2` springs, the `sum_i cos(psi_i)`
term vanishes and at the torsion minimum everything collapses to

    E_min = 0.5 * k * M * (sqrt(2*rho_r^2 + Lam^2) - d0)^2

a positive-definite quadratic in `L`. Since `mu=1.0` sets `d0` to exactly that
square root evaluated at AMBER's *ideal* bond length, `E_min` vanishes only when
the real bond is ideal — any real deviation adds spurious energy.

**Measured on real MD frames before the fix:** replaying every deployed ring on
a real Fs-peptide frame and recomputing each with its axis rescaled to the ideal
length gave a total artifact of **+35.46 kJ/mol against an observed dihedral gap
of +38.06** — ~93% of the discrepancy. Per family: guanidinium +32.40, omega
+18.88, phi -11.09, methyl -2.94, chi1 -1.79. Note these largely *cancel*
(net +3.06 before the guanidinium axes existed), which is exactly why the defect
stayed invisible through every earlier validation. Arg's guanidinium is the
pathological case because `E_min >= 0` always, and resonance pins every real Arg
at its planar torsion minimum — so all 3 axes contribute a strictly one-signed
artifact that never cancels.

**The fix is free.** The harmonic content depends only on
`eps = 2*rho_r^2/(Lam^2 + 2*rho_r^2)`, which is exactly `RING_SHAPES`' own
`rho`. Scaling `rho_r` and `Lam` *together* leaves `eps` — hence the entire
torsion curve — untouched while multiplying the artifact by `1/scale^2`. `Lam`
is decoupled from the real bond by tilting the ghosts off their anchor's
perpendicular plane (`theta != 90 deg`, already supported by the file format and
by `GhostParticle::computePosition`, and already used by `n=1`), so **no format
change and no C++ change**. `n=1` needs no scaling: its exact `d0=0`/`theta=60`
construction already gives `Lam = L - 2*r*cos(60 deg) = 0` at the ideal length.

| scale | r (A) | theta | k | amplitude | artifact @ dL=-0.06 A |
|-------|-------|-------|-----|-----------|----------------------|
| 1 | 0.540 | 90.0 | 2208 | 80.33 | 14.87 |
| 2 | 1.271 | 121.8 | 552 | 80.33 | 5.81 |
| 4 | 2.950 | 133.0 | 138 | 80.33 | 2.50 |
| 6 | 4.660 | 136.0 | 61 | 80.33 | 1.58 |

Verified: the torsion curve is identical across all scales to `1.7e-13 kJ/mol`,
and a full regeneration at `RING_SCALE = 4.0` reproduces all 90 per-axis
validation curves to **0.000000 kJ/mol**. `ctest` 278/278, and a 5000-step
ubiquitin run at the usual 1.5 fs still relaxes monotonically (kinetic
38.01 -> 0.46 kJ/mol) despite the longer force-redistribution lever arms.

#### Residual harmonic content (quantified, small, and now a solved-once precomputation)

A ring's own leakage into unwanted harmonics (mostly `2n`, `3n`, ...) is
controlled by its shape (a 2-parameter family `(rho, mu)` — how far the
ring's own points sit from the axis, and the ring's own `d0`) — found once,
per `n`, by a small grid search on the *abstract* shape function alone
(independent of any specific axis or residue, reused identically across
every axis), targeting <=0.5-1% leakage into the first unwanted multiple.
Measured shape-only RMS error against the full target curve, for every axis
family currently implemented (chi1-4, phi/psi, omega): **0.07-0.2%** of the
axis's own peak-to-peak amplitude — two to three orders of magnitude
smaller than the abandoned real-atom method's 30-45% (Section 3.2's
historical note below).

#### Per-owning-pair anchoring: which real atoms a ring's `Ref` points use

A ring's calibrated `(r, theta, d0, k, delta_base)` only fixes its *shape*;
which 3 real atoms actually drive its placement at runtime (`B`, `C`, and
each side's own `Ref`) is a separate, free choice, and it matters: a
deployed ghost's real 3D motion depends **only** on its own 3 anchors, never
on any other real atom used solely to build the abstract calibration
target. If an axis's real AMBER target is built by combining several real
substituent pairs (e.g. psi's target mixes an `(N, +N)` pair and an `(HA,
O)` pair), but only ONE generic reference atom is picked per side for
*every* ring regardless of which pairs actually contributed, the deployed
ghost implicitly assumes every other real substituent sits at a FIXED
offset from that one reference — which is only as good as that assumption.

For ordinary same-residue substituents (e.g. chi1's CB-side neighbours, all
directly bonded to the same vertex, no intervening independent torsion)
this assumption holds fine, matching the small residual already quoted.
Found and fixed (2026-08-04/05) for two real cases where it didn't:

* **Omega's n=1 harmonic** is owned entirely by one real pair, `(O, +H)` —
  anchoring that ring directly on O and +H (instead of a generic `CA/+CA`
  reference shared with n=2) lets the deployed ghost track that pair's
  actual real motion instead of assuming it sits at a fixed offset from CA
  — reduces this ring's own reconstruction error, measured against a real
  trajectory, from a small-but-nonzero residual to **exactly zero**.
* **Phi's b-side `{-C, H}` and psi's c-side `{O, +N}`**: one atom on each of
  these sides belongs to the *adjacent* residue's own peptide plane. Using
  a single shared reference for that side (as the original implementation
  did) showed up as a real, large energy residual on psi specifically (ALA
  psi: up to ~90% error on some real conformations) — root-caused NOT to
  some hidden extra rotational freedom (both O and +N are directly,
  rigidly bonded to C, a trigonal-planar centre with no torsional freedom
  of its own substituents), but to ordinary bond-angle-scale thermal noise
  (a few degrees, present on every axis, always has been) becoming visible
  specifically here because backbone torsion `k`'s are large (psi's own
  n=2 term alone: k=6.61, an order of magnitude above a typical side-chain
  k). Fixed generically: `generate_backbone_axis` detects whichever side
  has a cross-residue neighbour and splits *that side only* into one group
  per atom, each independently anchored and independently emitted — the
  non-risky side keeps one shared reference, unaffected. This generalizes
  automatically and safely: applied to phi, it found phi's own `H` group
  contributes nothing (only the null wildcard matches `-C`'s partner atoms
  other than the canonical pair) and produced byte-identical output to
  before — the mechanism, not a per-axis decision, determined phi needed
  no change.

**Validated end-to-end** (Fs-peptide, real trajectory, 20 frames): isolating
ALA's psi axis specifically, mean/max error fell from +49.5/+88 kJ/mol
(worst-case 38% relative) to **+8.7/+18.1 kJ/mol (3.9% average, 7.3% max)**
after this fix — roughly a 5-6x reduction. Full-system validation (all
proper torsions including omega, fair comparison against an independent
OpenMM/`amber99sb.xml` reference) improved from 3.87%/10.25% (mean/max
relative error, corr=-0.27, i.e. barely tracking real conformational
change at all) to **3.53%/5.74% (corr=+0.67)** — the clearest evidence that
the fix addresses the actual mechanism rather than cosmetically shrinking
an aggregate number: the model went from mostly noise to a real signal
that meaningfully tracks the true energy across frames.

*A related sign bug, found and fixed at the same time.* The per-substituent
azimuthal offset (`delta_of` in the generator) had its b-side sign
backwards in every function that computed it (`generate_backbone_axis`,
`generate_sidechain_axis`, `generate_sidechain_axis_gkinase`) — confirmed
by direct algebraic derivation and by a real-trajectory reconstruction
test. Fixed in all three places; harmless for a purely-real single-term
axis (chi1) but materially wrong for any axis with a genuinely complex,
multi-pair target (backbone axes) — folded into the numbers above, since
both fixes were needed together to reach them.

#### Applying this to real AMBER data

**Side-chain chi1-4 — implemented, all 20 residue types (30 valid axes).**
Unlike the abandoned real-atom method (one shared generic term per axis,
Section 3.2's historical note below), the current generator matches each
axis's real AMBER-specific substituent pairs directly (specific-entry-first,
generic-wildcard fallback, exactly mirroring AMBER's own matching order),
decomposing multi-term axes (e.g. `CT-CT-CT-CT`'s real 3-term chi2/chi3
structure) into one ring per real harmonic rather than approximating with a
single shared term. 6 axes (Asp/His/Phe/Tyr chi2, Glu chi3, Arg chi4) are
confirmed, not assumed, genuinely flat in `amber99sb.xml` (k=0 for every
periodicity on the relevant wildcard) — correctly emit no ghost springs at
all, not a coverage gap.

**Backbone phi/psi — implemented, all 18 residue types (36 axes).** Checked
explicitly for a CMAP correction (a 2D phi-psi coupling table used in some
force fields) before generating anything: `amber99sb.xml` contains no
`CMAPTorsionForce` element — ff99SB describes phi/psi with independent 1D
terms only, so the per-axis construction above is complete for this force
field (would need rechecking for a revision that does have CMAP).

**Omega — implemented, the first cross-residue AXIS.** Every other axis has
both its axis atoms (`B`, `C`) inside the same residue; omega's own axis is
the peptide bond itself, `C(i)-+N(i+1)` — `BondedForceFieldReader`'s atom
resolution needed no change at all to support this (it already resolved
`+`/`-` prefixes independently per atom slot for every other purpose), but
the *generator* needed two real, cross-residue-specific bugs fixed while
wiring it in: (1) naming a ring's `Ref` atom relative to the wrong
residue (`c_atom`'s own residue, correct for every other axis since
`c_atom` is always inside the owning residue — wrong for omega, where it
isn't) silently resolved to the wrong, same-named real atom instead of
failing, caught only via a real unresolved-anchor warning on Gly; (2) two
ring groups sharing the same non-risky-side reference (see per-owning-pair
anchoring above) reused identical ghost/rule names, silently colliding —
fixed by a `group_tag` disambiguating ghost names whenever an axis emits
more than one group. Before choosing ghost particles for omega, a plain
real-atom-anchored attempt (the historical method, one ring using the 4
real substituents `CA(i)/O(i)` x `CA(i+1)/H(i+1)` directly) was tried first
and confirmed insufficient: a ~169deg phase mismatch between the real
geometry's own natural n=2 content and AMBER's target, with no free
parameter to fix it (see Step 3) — even a full grid search over (k,d0)
only reached 138% mean relative error. Ghost particles (n=2 generic
`CA/+CA` anchor, n=1 exact `O/+H` anchor per the rule above) validated at
mean relative error 3.87%/10.25% in isolation, improving to 3.53%/5.74%
system-wide after the phi/psi fixes above (omega itself unaffected by
those, already using per-owning-pair anchoring from the start).

Making omega flexible at all required a companion `--rigidbody` change:
the base `ProteinAtomRigidGroups.rbody`'s `_PHI`/`_PSI` groups treat the
*entire* peptide plane (6 atoms spanning both sides of the omega bond) as
one rigid clique, freezing omega completely regardless of any dihedral
term layered on top. A new file, `ProteinAtomRigidGroupsOmegaFree.rbody`
(the original left untouched, kept as the "no cis/trans freedom" model),
splits each group into two halves hinged at exactly the omega bond — the
same "2 shared atoms = 1 free hinge" trick already used for every other
rotatable axis. Verified two ways: an exact spring-pair diff against the
original file (294 pairs removed, 0 added, every removed pair of exactly
the 4 predicted "upstream-non-axis to downstream-non-axis" types, none
touching the hinge atoms themselves) and a perturbed-structure dynamics
run (clean relaxation, no divergence, comparable behaviour to the
original file).

**Improper (planarity) — deprioritized, not implemented.** Worth stating
precisely what this term would and would not fix: `--rigidbody`'s
all-pairwise mesh already geometrically guarantees ring/guanidinium
planarity on its own (fixing every pairwise distance in a group removes
every internal degree of freedom, including out-of-plane pucker), as long
as it remains the permanent substrate (Section 2) — which it always does
in this model. So an improper ghost spring would not fix a real risk of
losing planarity; it would only replace `--rigidbody`'s arbitrary uniform
stiffness with the real, physically-calibrated AMBER restoring force (a
real ring has small thermal out-of-plane fluctuations, not perfect
rigidity), and correctly except Proline (see below). Deprioritized by the
user (2026-07-31) once this distinction was clear — a refinement to
revisit later, not a gap being left open. The rest of this subsection
records what was worked out before that decision, for whenever it is
picked back up.

`amber99sb.xml`'s `<Improper>` table is
overwhelmingly uniform: n=2, phase 180deg, most commonly k=4.6024 (e.g.
`CA-CA-CA-CT`, the aromatic-ring-to-CB vertex shared by Phe/Tyr-type side
chains), with a few stiffer outliers (guanidinium's `CA-*-N2-N2`,
k=43.932). A single (2,2) group reproduces the common case cleanly. Matching
convention differs from proper dihedrals (wildcard/empty inner classes,
atom-order handling specific to the improper table) — verified per case, not
assumed to mirror the proper table. Proline is a deliberate exception: its
ring genuinely puckers in real AMBER (no improper restrains it) — unlike its
current `--rigidbody` treatment (full pairwise rigidity), which is *less*
physically accurate on this specific point.

#### Historical note: the abandoned real-atom-only method

Before ghost particles, dihedral rings were built directly between real
substituent atoms (Steps 1-2 above still apply verbatim — they are generic
geometric/Fourier identities — but `r`/`theta`/`delta` were then fixed by
real chemistry, and `M`/`N` capped by how many real substituents actually
exist on each side). This worked for single-term, evenly-spaced axes
(chi1's generic `X-CT-CT-X` term, real substituents close enough to 120deg
apart) but failed for anything more demanding: chi2's real 3-term structure
needed more independently-tunable degrees of freedom than 3 fixed real
hydrogens provide (residual 30-45% of target, and splitting further to
chase individual terms made small subgroups fail far worse, up to ~30x
target, from being too few/unevenly-spaced to cancel unwanted harmonics);
omega's real geometry turned out to have its own natural harmonic phase
~169deg away from AMBER's target with no free parameter to correct it.
Ghost particles remove the constraint that geometry must equal real
chemistry, at the cost of an extra virtual site per ring point — this is
why every axis implemented today uses ghost particles, not real atoms
directly (bending's own real-atom-to-real-atom leakage problem, Section
3.1, is the same lesson applied to Term 2).

Stages 1–2 retune and extend the rigid-body mesh, but always keep it as the
permanent substrate: the mesh's real 1-2/1-3 pairs are what Stage 1 retunes,
and every *other* pairwise spring the mesh creates inside a rigid group
(e.g. non-adjacent atoms inside a ring or a methyl group) stays at its
uniform value, providing the same all-pairwise rigidity as today wherever
a more specific term hasn't been layered on top. A model that drops
`--rigidbody` entirely is explicitly out of scope: it would need a full
coverage audit (every real 1-2 pair generated, no silent gaps) that hasn't
been attempted, and `--rigidbody` already does a perfectly good job of
covering whatever a given `.bi.ff` doesn't.

4. Computational performance
--------------------------------

A natural question once a term is reformulated purely as distance springs:
is the resulting force evaluation cheaper or more expensive than the
conventional trigonometric formula it replaces? A microbenchmark was built
to answer this directly (C++, `-O3`, `std::chrono`; both force formulas
checked against a finite-difference gradient of their own energy before any
timing was trusted — see below).

**Correctness check (prerequisite to trusting any timing number):**

    angle_force (trig): |analytic - finite-difference| = 1.27e-09
    dihedral_force (trig): |analytic - finite-difference| = 9.29e-11

Both formulas are correct to numerical precision.

**A methodological correction, reported because it materially changed the
conclusion.** The first version of this benchmark measured the ghost-spring
side using a small per-iteration `std::vector` allocation to collect forces
from a variable-size group of springs. That allocation overhead — not the
underlying arithmetic — dominated the measurement, making the ghost-spring
path look 2-4x *slower* than the traditional trigonometric formula. Rerunning
the identical arithmetic with a fixed-size, allocation-free accumulation (the
way BioSpring's actual `SpringNetwork` iterates its spring arrays — no
per-force-evaluation heap traffic) reversed the conclusion entirely. This is
reported explicitly because it is exactly the kind of benchmarking mistake
that produces a confident, wrong answer if not checked — the fix mattered
more than any other single number here.

**Results (20,000,000 iterations per case, randomised inputs to defeat
constant folding):**

| Case | Traditional (trig) | Springs (allocation-free) | Result |
|------|----------------------|-------------------------------|-----------|
| Angle bending (1 angle vs 1 spring) | 10.5 ns | 3.2 ns | **springs 3.3x faster** |
| Dihedral, single AMBER term (9 ghost springs, M=N=3) | 26.4 ns | 6.8 ns | **springs 3.9x faster** |
| Dihedral, full Arg chi2 (12 ghost springs, 3 terms) | 29.5 ns | 9.6 ns | **springs 3.1x faster** |

The spring-only formulation is faster in every case tested, by a consistent
factor of roughly 3-4x, *despite* needing 9-12 individual spring evaluations
to reproduce a single dihedral axis (vs. one shared, but transcendental-
function-heavy, dihedral-angle computation in the traditional formula). The
reason is architectural: a distance-spring force needs one square root and a
handful of multiplications; a trigonometric bonded force needs `acos`/`atan2`
plus several normalized cross products — each of those transcendental calls
individually costs more than an entire spring evaluation. Reformulating
every bonded term as a sum of plain distance springs replaces a small number
of expensive, branchy, transcendental-function-heavy evaluations with a
larger number of uniform, transcendental-function-free ones — and, measured
here, the uniform path wins even after accounting for the multiplicity a
dihedral term requires.

*Caveat.* This measures raw per-evaluation scalar cost on one machine, in
isolation. It does not measure whole-simulation throughput, where other
factors (memory layout, SIMD vectorisation of a uniform spring kernel across
heterogeneous term types, cache behaviour, thread scheduling) could shift
the picture further — plausibly in the spring-only network's favour, since
every term becomes the *same* uniform kernel rather than several
structurally different ones, but this has not been measured and is not
claimed here. It also predates the move from ghost springs on real atoms to
ghost *particles* (Section 3.2): a ghost particle's own placement
(`computePosition`) and force redistribution (the transpose-Jacobian
`redistributeForce`) add real per-ghost overhead this benchmark does not
account for, on top of the per-spring cost measured here — not remeasured
since, and worth doing before citing this table as still current for the
present mechanism.

5. Validation strategy
-------------------------

1. **Self-consistency at generation time.** The generator already asserts
   real bond/angle connectivity (mdtraj-inferred, never hand-typed) before
   emitting anything. For every dihedral ghost-particle ring, it
   additionally re-derives the ring's own reproduced harmonic amplitude
   (the FFT quadrature used above) and checks it against the AMBER target
   before writing anything.
2. **Full build + `ctest`** (278 tests, no regression) — every dihedral
   axis change so far (BEND's ghost-particle fix, omega, the sign-bug and
   per-owning-pair-anchoring fixes) needed no C++ changes at all, since the
   ghost-particle/`.bi.ff` mechanism is already fully generic; this step
   still catches anything unrelated it might disturb.
3. **Independent energy cross-check against real AMBER (OpenMM,
   `amber99sb.xml`), not just self-consistency.** Extract real frames from
   a real MD trajectory (Fs-peptide), build a matching BioSpring `.nc` per
   frame, and compare BioSpring's own reported energy (`Stretch energy`,
   `Bend energy`, `Dihedral energy` — see Section 3.1's energy-channel
   note) against an independent per-frame OpenMM computation. This is what
   actually found bending's bond-stretch leakage (Section 3.1) and psi's
   per-owning-pair-anchoring gap (Section 3.2) — self-consistency checks
   alone (point 1) had missed both, since they only confirm a ring
   reproduces its OWN calibration target, not that the target itself
   tracks real dynamics correctly.
4. **Per-axis, per-residue isolation** when an aggregate comparison shows
   an unexplained residual: filter the `.bi.ff` down to one axis (or one
   real substituent pair) at a time, rebuild, and compare against a
   matching real-AMBER decomposition for exactly that axis. This is how
   ALA's psi was isolated as the dominant contributor to an aggregate
   ~14-19% dihedral error that other axes' own small residuals didn't
   explain.
5. **Dynamical, qualitative validation on real data.** For Arg chi1/chi2:
   run BioSpring dynamics long enough to observe rotamer-well transitions,
   compare occupancy/dwell-time statistics qualitatively against the real
   measured Fs-peptide trajectory (`rotamer_hinge_vs_md_report.html`).
6. **Planarity check**, for a ring/guanidinium improper group (not yet
   implemented): verify the real dihedral stays near-planar with realistic
   fluctuation amplitude (not perfectly rigid as with `--rigidbody` alone,
   not unrestrained either).
7. **Intermediate-model check**, with `--dihedral` alone (Stage 2 without
   Stage 1, i.e. `-bondedinteraction` without `-stretching`/`-bending`):
   bonds/angles stay at the rigid-body uniform value, dihedral
   energy/behaviour is identical to the fully-refined model — the two
   concerns are independent by construction, confirmed experimentally.
8. **Structural diff for a `--rigidbody` change** (omega's
   `ProteinAtomRigidGroupsOmegaFree.rbody`): an exact before/after diff of
   every spring `--rigidbody` creates, checking that only the intended
   pairs disappear (and none appear) — a more rigorous check than a
   dynamics run alone, which can't easily isolate one specific hinge's
   freedom from whole-molecule relaxation noise.

6. Known limitations
------------------------

* **Stretching:** none — a direct, exact representation.
* **Bending:** ghost-anchored (Section 3.1) to remove real bond-stretch
  leakage; the remaining residual is the curvature-matching
  linearization error alone, validated at mean 1.03%/max 2.11% against
  real AMBER (Fs-peptide, 20 frames, stretch+bend combined).
* **Dihedral, all proper families (chi1-4, phi/psi, omega):** each ring's
  own shape-only residual is small (0.07-0.2% of peak-to-peak, Section
  3.2) — the dominant remaining error source is real thermal noise on
  whichever real atom a ring is anchored to, not yet correctable further
  without tracking additional real atoms per ring. Full real-AMBER
  validation (Fs-peptide, 20 frames, all proper torsions): mean relative
  error 3.53%, max 5.74%, corr=+0.67 (a real, meaningfully-tracking
  signal, not just noise around the right average). 6 side-chain axes
  (Asp/His/Phe/Tyr chi2, Glu chi3, Arg chi4) are confirmed genuinely flat
  in AMBER (k=0), not a gap.
* **Backbone phi/psi/omega:** no CMAP gap for ff99SB specifically,
  verified; a different AMBER/CHARMM revision with a real CMAP term would
  need this re-checked. Per-owning-pair anchoring (Section 3.2) covers the
  two known cross-residue-neighbour cases (phi's `-C`, psi's `+N`,
  omega's `+H`) found so far by direct real-trajectory validation, not by
  an exhaustive audit of every axis — a future axis with a similar
  cross-residue neighbour should be checked the same way, not assumed
  fixed by the same generic mechanism without verification.
* **Improper (planarity):** not implemented, deprioritized (Section 3.2) —
  `--rigidbody`'s mesh already geometrically guarantees planarity as long
  as it remains the permanent substrate; matching convention would need
  per-case verification against `amber99sb.xml` if picked back up, and
  Proline intentionally excluded (real ring pucker, no AMBER improper term
  there).
* **Scope:** `--rigidbody` remains a permanent, required substrate for this
  model — dropping it entirely is out of scope, not merely deferred.

7. Conclusion
----------------

All four of AMBER's bonded interaction terms can be expressed as ordinary
distance springs, with every conversion either exact (stretching), a
provable second-order match free of bond-stretch leakage (bending, via
ghost particles anchored along the real, current bond directions), or an
exact closed-form/Fourier ghost-particle-ring construction (every proper
dihedral family implemented: chi1-4, phi/psi, and omega, the first axis
whose own two axis atoms span a residue boundary) with quantified,
now-small residual error — validated end-to-end against real,
independent AMBER (OpenMM) computations on real MD trajectory frames, not
merely self-consistent with its own calibration target. Layering these
onto BioSpring's existing rigid-body mesh — which remains the permanent
substrate, not a step to be eventually removed — gives a practical,
incrementally-testable path from a purely topological approximation to a
physically-parameterised, spring-only model of the four bonded terms —
and, measured directly, the resulting force evaluations are not merely
equivalent but consistently faster than the conventional trigonometric
formulas they replace. The improper (PLANARITY) family is implemented
too, where it is physically load-bearing: an AMBER improper IS a 4-atom
dihedral about a real central bond (hub listed third), so
`emit_ghost_ring` models it as-is — one single-pair n=2 ring per
improper, axis = the central bond, refs = the two peripheral atoms,
parameters read from a real built OpenMM System (AMBER's own assignment
engine — no reimplementation of the improper matching rules). It gates
behind pdb2spn's `-dihedralplanarity` flag, folded into `-dihedral`.
Deployed for the aromatic rings (Phe/Tyr/His/Trp), whose rigid cliques
were split into per-vertex hinge groups at the same time: unlike Pro's
5-ring (freed earlier — bonds+angles hold every 5-ring distance, the
pucker survives as the real soft mode, and the real molecule has no
improper there), a 6-ring can fold (its para pairs are 1-4) and an sp2
substituent can leave the plane at first order, so freeing those rings is
only correct together with BOTH the ring proper torsions (X-CA-CA-X n=2
k=15.167, per-pair rings — measured exact, 0.00 torque error over 1140
instances) and these impropers (measured exact, 0.00 over 990 instances;
in real dynamics the freed rings stay planar to 0.15-0.66 deg mean
out-of-plane with stable para distances). Impropers on centres still held
rigid by cliques (the backbone peptide C/N inside `_PHI`/`_PSI`, the
split amide and guanidinium groups) are deliberately NOT emitted: their
planarity is enforced geometrically, and adding springs there would only
double-book energy.
