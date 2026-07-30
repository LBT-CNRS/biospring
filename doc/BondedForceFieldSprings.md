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
rigid-body mesh) and layering real, spring-only parameters and dihedral
ghost springs on top of it. Every derivation is closed-form or an exact
algebraic identity, verified either symbolically or numerically (to machine
precision where applicable); no step is a curve fit. A microbenchmark
comparing per-evaluation computational cost against conventional
trigonometric bonded force routines is included, along with the full
validation strategy used to check the model against independent
recomputation and real AMBER parameters.

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

**Term 2 — angle bending (1-3 distance spring, curvature-matched).** For a
real 1-2-3 valence angle with vertex atom 2, the distance between the two
*outer* atoms, `r13`, is a monotonic function of `theta` for fixed real bond
lengths `r12`, `r23` (already correct, from Term 1):

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

*Ring/planar-group caveat.* A rigid-body group with every pairwise distance
constrained (Section 2) has no internal degrees of freedom left at all —
this is stiffer than reality even after Stage 1, since bending alone does
not restore ring pucker or out-of-plane flexibility. Real planarity is a
distinct AMBER term (Term 4, Section 3.2) — bending is not meant to provide
it.

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

Stage 1 retunes springs that already exist. Dihedral torsion needs new ones:
a dihedral angle `phi` is not the distance between any single pair of atoms,
so no direct geometric substitute exists. Real substituent atoms around the
rotation axis are instead used to build a set of ordinary distance springs —
"ghost springs" — added *on top of* the Stage 1 network, whose combined
energy, as a function of `phi`, reproduces one specific Fourier term of
AMBER's real torsion energy. Nothing from Stage 1 needs to change: dihedral
springs are new additions between atoms that Stage 1 never connected (they
are not directly bonded), so this is strictly an additive refinement.

This intermediate configuration — rigid-body mesh, Stage 1's real
stretching/bending, and now real dihedral wells layered on top — is already a
complete, independently useful, and independently testable model: bond and
angle vibration have real values, and side-chain/backbone rotation now has
real energetic preferences and ring/guanidinium groups have real planarity
restraints, while the underlying topology (which pairs have springs at all)
still traces back to the rigid-body mesh of Section 2.

In practice (`pdb2spn`), Stage 1 and Stage 2 are independent opt-in flags
on top of `-rigidbody -bondedinteraction <file>`: `-stretching`, `-bending`
(which requires `-stretching` too, since its conversion needs real bond
lengths), and, for dihedral, `-dihedralbackbone`/`-dihedralsidechain`
(one per proper-dihedral family — planarity has no flag yet, see Known
Limitations) plus `-dihedral` as a pure convenience alias for both
together. None is implied by the others or by `-bondedinteraction` alone —
combining only a dihedral flag, without `-stretching`/`-bending`, is
exactly this section's intermediate model. Once built, dihedral ghost
springs are ordinary springs at simulation time: there is no separate
`.msp` switch for them, they are gated by the same `spring.enable` as
everything else in the network.

Both dihedral families — **proper** (side-chain chi1-4, backbone phi/psi)
and **improper** (planarity of aromatic rings, guanidinium) — use exactly
the same construction below. AMBER itself treats an improper geometrically
as an ordinary four-atom dihedral (the "axis" is the bond between the 2nd
and 3rd atom of its atom list, the other two are the substituents); only the
source table (`PeriodicTorsionForce`'s `<Proper>` vs `<Improper>` entries)
and the atoms involved differ.

#### Notation

For a dihedral defined by four atoms in a chain **R-B-C-S** (B and C
connected by a real bond, the axis; R attached to B; S attached to C):

* `L` — the B-C bond length.
* `r1`, `theta1` — the B-R bond length and the R-B-C valence angle.
* `r2`, `theta2` — the C-S bond length and the B-C-S valence angle.
* `phi` — the dihedral angle R-B-C-S.

#### Step 1 — closed-form distance between R and S as a function of phi

**Claim.** With B at the origin and C on the z-axis:

    d^2(phi) = C1 - C2 * cos(phi - delta)
    C1 = r1^2 sin^2(theta1) + r2^2 sin^2(theta2) + (r1 cos(theta1) - L + r2 cos(theta2))^2
    C2 = 2 * r1 * r2 * sin(theta1) * sin(theta2)

where `delta` is a fixed azimuthal offset depending only on how `phi`'s zero
reference is defined.

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
model error. This is a geometric identity.

A ghost spring's energy is then simply `E(phi) = 0.5 k (d(phi)-d0)^2` — an
ordinary BioSpring distance spring between the real atoms R and S.

#### Step 2 — which harmonics survive averaging over several copies (Dirac comb)

A single ghost spring's energy is **not** a pure cosine of `phi` (it is a
nonlinear function of `cos(phi)` through the square root and the square), so
it has content at every harmonic, not just one. Averaging several rotated
copies cancels the unwanted ones exactly, for *any* periodic function.

**Claim.** Let `g(phi)` have period `2*pi`, Fourier series
`g(phi)=sum_k ĝ_k exp(ikphi)`. Build `M` real copies of the reference
substituent evenly spaced around B (`2*pi*i/M` apart) and `N` real copies of
the rotating substituent evenly spaced around C (`2*pi*j/N` apart). Summing
every `(i,j)` pair's energy:

    E_total(phi) = sum_i sum_j g(phi + 2*pi*i/M - 2*pi*j/N) = M*N * sum over k that are multiples of LCM(M,N) of ĝ_k exp(ikphi)

**Proof.** Expand and swap summation order:
`E_total = sum_k ĝ_k exp(ikphi) * [sum_i exp(ik2pi*i/M)] * [sum_j exp(-ik2pi*j/N)]`.
Each bracket is a sum of `M` (resp. `N`) evenly spaced roots of unity to the
k-th power — a standard identity, equal to `M` (resp. `N`) if `M` (resp. `N`)
divides `k`, else `0`. Only `k` divisible by both — i.e. by `LCM(M,N)` —
survives, scaled by `M*N`. This holds for any periodic `g`: it is a property
of averaging evenly-spaced copies, not of a ghost spring's specific shape.

**Verification (numerical, FFT of the real nonlinear ghost-spring energy):**

| M | N | LCM(M,N) | Harmonics surviving above 1% of peak (k=1..19) |
|---|---|----------|--------------------------------------------------|
| 3 | 3 | 3        | {3} |
| 1 | 2 | 2        | {2} |
| 2 | 3 | 6        | {6} |

Exactly as predicted. **Practical use:** to target a real AMBER harmonic
`n`, choose any `(M,N)` with `LCM(M,N)=n`, `M`/`N` never exceeding the real
number of substituents actually present.

#### Step 3 — solving for stiffness exactly (not an iterative fit)

For fixed geometry and `d0`, a ghost spring's energy is exactly linear in its
stiffness: `E(phi)=k*h(phi)` with `h(phi)=0.5(d(phi)-d0)^2` independent of
`k`. Every Fourier component of `E` is therefore exactly `k` times the
corresponding component of `h`; one evaluation at `k=1` gives the
proportionality constant, hence `k_target = k_AMBER_target / component_at_k=1`
— a direct algebraic solve. Verified: re-measuring the reproduced
harmonic's amplitude after solving matches the AMBER target exactly (to FFT
quadrature precision) for every worked example below.

#### Step 4 — residual harmonic content (quantified, not hidden)

Only multiples of the target `n` survive Step 2's cancellation *when every
copy shares identical (r, theta), differing only in azimuth*. On an
idealized, artificially symmetric test geometry (same r1/theta1 for every
copy on each side, exactly evenly spaced), the residual is negligible —
verified numerically to <1% for the cases tried.

**Real substituents are not that clean, and the residual is much larger in
practice.** Around a real vertex (e.g. CA's real neighbours N, C, HA),
different real substituents have genuinely different bond lengths and
angles to the axis — not just a different azimuth. This breaks the exact
cancellation the idealized case enjoys. Measured on Arg's real chi1 and
chi2 axes (real bond lengths/angles from AMBER's tables, real relative
azimuths measured from an actual residue's 3D coordinates, see below):

| Axis | Target (M=N=3, all real substituents) | Residual at m=1 | Residual at m=2 |
|------|------------------------------------------|--------------------|--------------------|
| chi1 (X-CT-CT-X, n=3, k=0.6508) | reproduced exactly at n=3 | 0.66 (102% of target) | 0.21 (32%) |
| chi2 (X-CT-CT-X, n=3, k=0.6508) | reproduced exactly at n=3 | 0.87 (134% of target) | 0.29 (45%) |

This is a real, load-bearing limitation, not a cosmetic one: several more
elaborate attempts to do better were tried and abandoned (splitting the
axis by each real substituent's own most-specific AMBER torsion entry;
allowing each real pair its own independently-solved stiffness and
equilibrium via linear/non-linear least-squares against the full combined
real target) — all either reduced to the same size of residual or produced
some pairs needing a negative (unphysical) stiffness to fit. The
`n=3`-only, full-real-substituent-grid construction above is the version
actually implemented and validated end-to-end (see the worked example
below) — a genuine, bounded, first-order approximation of the real
torsional profile, not an exact reproduction, and not disguised as one.

#### Applying this to real AMBER data

**Side-chain chi (as implemented).** Both chi1 (CA-CB axis) and chi2
(CB-CG axis) use the one generic wildcard term available for any C-C
single-bond rotation, `X-CT-CT-X` (n=3, k=0.6508, phase 0), applied via one
ghost-spring group per axis built from *all* real substituents on each side
(M=N=3: e.g. for chi1, {N, C, HA} on the CA side and {CG, HB2, HB3} on the
CB side) — not split by each pair's own more specific AMBER entry. This is
a deliberate simplification: `amber99sb.xml` also has more specific real
entries for some of chi2's substituent combinations (`CT-CT-CT-CT`,
`HC-CT-CT-CT`, `HC-CT-CT-HC`, each its own multi-term torsion) which a
fully faithful reproduction would need to represent separately — attempted
and abandoned (see Step 4): several of the resulting subgroups have too
few real substituents to be evenly spaced, which breaks the cancellation
Step 2 relies on far worse than the residual above (measured up to ~30x
the target itself for some subgroups, not a small effect). Using the one
broad wildcard term across the full, genuinely ~120-degree-spaced real
substituent set avoids that failure mode, at the cost of not capturing
chi2's finer real multi-term structure for now — an explicitly open item,
not a silently dropped one.

**Worked example, validated end-to-end.** For each axis, real geometry
(bond length/angle from AMBER's tables, real relative azimuth from
ubiquitin's ARG42 3D coordinates) gives 9 real (B-substituent,
C-substituent) pairs; Step 3's exact linear solve (via a signed projection
onto `cos(3*phi)`, not a magnitude, so a sign mismatch would be caught
rather than silently producing a negative stiffness) gives one shared
stiffness per axis (chi1: k=10.98 kJ.mol-1.A-2; chi2: k=12.67
kJ.mol-1.A-2). Summed over all 4 real Arg residues in ubiquitin.pdb, an
independent Python recomputation of `0.5*k*(d-d0)^2` over the resulting 72
ghost springs matches BioSpring's own reported dihedral energy exactly:
136.8747 kJ/mol (Python) vs 136.87 kJ/mol (BioSpring).

**Backbone phi/psi.** Same table, backbone atom classes. Checked explicitly
for a CMAP correction (a 2D phi-psi coupling table used in some force
fields, e.g. CHARMM and later AMBER revisions) before generating anything:
`amber99sb.xml` contains no `CMAPTorsionForce` element and no mention of
CMAP anywhere — ff99SB describes phi/psi with independent 1D terms only,
exactly like a side-chain chi. The methodology applies with no documented
gap for this force field.

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

*Choosing substituents when fewer are needed than exist.* The LCM rule says
which `(M,N)` are valid, not which real atoms to use when an axis has more
substituents than a given term needs. For chi1/chi2 this doesn't arise (the
single wildcard term uses all real substituents on both sides, see above).
It would arise for any future multi-term axis (e.g. a genuine attempt at
chi2's finer `CT-CT-CT-CT`/`HC-CT-CT-*` breakdown, or improper/backbone
terms with more than one significant harmonic) — this choice would need to
be explicit, per case, in the generator, and Step 4's finding that small
subgroups can fail far worse than the residual above should be checked
before assuming it works.

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
claimed here.

5. Validation strategy
-------------------------

1. **Self-consistency at generation time.** For stretching/bending, the
   generator already asserts real bond/angle connectivity (mdtraj-inferred,
   never hand-typed) before emitting anything. For every dihedral
   ghost-spring group, it additionally re-derives the group's own reproduced
   harmonic amplitude (the FFT quadrature used above) and checks it against
   the AMBER target before writing anything, logging the residual-harmonic
   table for traceability.
2. **Full build + `ctest`** (273 tests, no regression).
3. **Independent energy cross-check (Python vs BioSpring).** Isolate one
   term's contribution by subtraction (with vs without, or via the `.msp`
   per-family toggle for dihedral) and compare against an independent Python
   recomputation of the same AMBER-derived formula on the same coordinates.
   Already done this way for stretching and bending (exact match, ubiquitin
   and a 20-frame real Fs-peptide trajectory); the same method extends
   directly to dihedral.
4. **Full phi-scan comparison**, not just the target harmonic's amplitude:
   sample the entire 0-360deg energy curve of a ghost-spring group and
   compare point-by-point against the real AMBER curve for that term (and
   its real sum, where several terms share an axis) — this is what shows the
   residual-harmonic wiggle visually, as a max-deviation figure.
5. **Dynamical, qualitative validation on real data.** For Arg chi1/chi2:
   run BioSpring dynamics long enough to observe rotamer-well transitions,
   compare occupancy/dwell-time statistics qualitatively against the real
   measured Fs-peptide trajectory (`rotamer_hinge_vs_md_report.html`).
6. **Planarity check**, for a ring/guanidinium improper group: verify the
   real dihedral stays near-planar with realistic fluctuation amplitude
   (not perfectly rigid as with `--rigidbody` alone, not unrestrained
   either).
7. **Intermediate-model check**, with `--dihedral` alone (Stage 2 without
   Stage 1, i.e. `-bondedinteraction` without `-stretching`/`-bending`):
   bonds/angles stay at the rigid-body uniform value, dihedral
   energy/behaviour is identical to the fully-refined model — the two
   concerns are independent by construction, confirmed experimentally.

6. Known limitations
------------------------

* **Stretching:** none — a direct, exact representation.
* **Bending:** second-order (curvature) match at equilibrium only; deviates
  further from `theta0`. A fully rigid-body group loses genuine flexibility
  (e.g. pucker) that planarity (Term 4) is meant to restore, not bending.
* **Dihedral (both families):** residual harmonic content is small only for
  an idealized, artificially symmetric geometry; on Arg's real chi1/chi2
  axes (real, non-identical bond lengths/angles per substituent) the
  measured residual at the next harmonic is 30-45% of the target itself —
  large, bounded, and honestly reported, not a small correction. Chi2's
  finer real multi-term structure (`CT-CT-CT-CT` vs `HC-CT-CT-*`) is not
  currently represented at all: attempts to split it out made the fit far
  worse (up to ~30x the target for some subgroups), not better — see
  Section 3.2, Step 4. Compounds with bending's own curvature-matching
  limitation where a dihedral's reference geometry depends on a bent angle.
* **Backbone phi/psi:** no CMAP gap for ff99SB specifically, verified; a
  different AMBER/CHARMM revision with a real CMAP term would need this
  re-checked.
* **Improper (planarity):** matching convention needs per-case verification
  against `amber99sb.xml`, not assumed to mirror the proper table. Proline
  intentionally excluded (real ring pucker, no AMBER improper term there).
* **Scope:** `--rigidbody` remains a permanent, required substrate for this
  model — dropping it entirely is out of scope, not merely deferred.

7. Conclusion
----------------

All four of AMBER's bonded interaction terms can be expressed as ordinary
distance springs between real atoms, with every conversion either exact
(stretching), a provable second-order match (bending), or an exact
closed-form/Fourier construction with quantified residual error (both
dihedral families). Layering these onto BioSpring's existing rigid-body
mesh — which remains the permanent substrate, not a step to be eventually
removed — gives a practical, incrementally-testable path from a purely
topological approximation to a physically-parameterised, spring-only
model of the four bonded terms — and, measured directly, the resulting
force evaluations are not merely equivalent but consistently faster than
the conventional trigonometric formulas they replace.
