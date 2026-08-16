#!/usr/bin/env python3
"""Force-field-agnostic core of the ghost-ring dihedral construction.

Everything here consumes AMBER parameter TABLES and nothing else -- no
structure, no residue templates, no notion of what a residue even is. That
is what lets the protein and the nucleic-acid generators share it verbatim
instead of growing two drifting copies of the same derivation: the comb/LCM
theorem, the closed-form d(phi), the linear solve for k, and the azimuths
derived from HarmonicAngleForce know nothing about amino acids.

Call configure() once with the tables the caller has parsed, then use
combined_target_for_axis to read an axis's real AMBER Fourier content and
emit_ghost_rings_for_axis to build the rings that reproduce it. The emitted
lines and counters accumulate here; the caller splices them into its own
output file (see dihedral_lines / ghostparticle_lines / the counters).

The one caller-specific hook is `stereocenter_delta`: a tetrahedral vertex
whose four substituents all have DISTINCT parameter identifiers has a real
chirality that no angle table can resolve, so the caller supplies it. Every
other vertex -- trigonal, or tetrahedral with a repeated substituent -- is
solved here from the tables alone (see formula_derived_deltas).

Shares bonded_rings.py with both callers; that module holds the pure ring
math (closed_form_d2, calibrate_ring, choose_d0).
"""

import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bonded_rings import (PHI_GRID_DEG, calibrate_ring, choose_d0,      # noqa: F401
                          closed_form_d2, complex_fourier_coeff, ring_curve_abstract)

# --- injected by configure() ------------------------------------------------
bond_params = {}
angle_params = {}
torsion_params = {}
stereocenter_delta = None


def configure(bond_params_, angle_params_, torsion_params_, stereocenter_delta_=None):
    """Point the core at one force field's tables.

    bond_params:    {frozenset((id1, id2)): (length_nm, k)}
    angle_params:   {(vertex_id, frozenset((id1, id2))): (angle_deg, k)}
    torsion_params: {(id1, id2, id3, id4): [(n, k, phase_rad), ...]}, "" wild
    stereocenter_delta: optional (vertex_ref_name, other_name) -> deg, for the
        genuinely chiral tetrahedral centres the tables cannot resolve.
    """
    global bond_params, angle_params, torsion_params, stereocenter_delta
    bond_params = bond_params_
    angle_params = angle_params_
    torsion_params = torsion_params_
    stereocenter_delta = stereocenter_delta_


def bump_axis_ok():
    global n_dihedral_ok
    n_dihedral_ok += 1


def bump_axis_skip():
    global n_dihedral_skip
    n_dihedral_skip += 1


def lookup_torsion_wildcard(c2, c3):
    for key in (("", c2, c3, ""), ("", c3, c2, "")):
        if key in torsion_params:
            return torsion_params[key]
    return None

def lookup_torsion_specific(c1, c2, c3, c4):
    # Exact (non-wildcarded) match for one specific real substituent pair,
    # tried in both read directions (a dihedral read backwards names the
    # same physical angle). Used for backbone phi/psi, which -- unlike
    # chi1's single generic X-CT-CT-X term -- amber99sb.xml encodes as a
    # handful of pair-specific entries (verified directly against the XML,
    # not assumed): the axis's own generic wildcard entry exists but carries
    # k=0 for every periodicity (a real, deliberate null placeholder, not a
    # missing table), so it never survives the `if _k != 0.0` filter above
    # and correctly never matches here either.
    for key in ((c1, c2, c3, c4), (c4, c3, c2, c1)):
        if key in torsion_params:
            return torsion_params[key]
    return None

def bond_len_A(c1, c2):
    p = bond_params.get(frozenset((c1, c2)))
    return p[0] * 10.0 if p else None

def valence_deg(vertex_class, c1, c2):
    p = angle_params.get((vertex_class, frozenset((c1, c2))))
    return p[0] if p else None

def idealized_azimuth_deg(vertex_class, ref_class, x_class, y_class):
    """Azimuthal offset of substituent y relative to x, viewed down the
    vertex->ref axis, derived PURELY from AMBER's own HarmonicAngleForce
    equilibrium angles at `vertex_class` -- no PDB structure of any kind.
    Spherical law of cosines: for 3 unit vectors from a common vertex at
    known pairwise angles (ref-x, ref-y, x-y), the azimuth phi of y around
    the ref axis (x's own azimuth taken as zero) satisfies
    cos(a_rx)*cos(a_ry) + sin(a_rx)*sin(a_ry)*cos(phi) = cos(a_xy).
    Assumes the vertex is planar (a_rx+a_ry+a_xy ~= 360 deg, true for any
    real sp2 center to within AMBER's own small parameterization
    residual) -- verified algebraically that exact planarity forces
    cos(phi)=-1 (phi=180 deg exactly); AMBER's real angle values give
    176-177 deg for omega's C(i)/N(i+1) vertices (~3 deg short of ideal
    360, a genuine feature of the real parameterization, not an error).
    Returns a value in [0,180] (sign/handedness needs one more real-world
    convention to resolve, irrelevant here since 180 deg is its own
    mirror image -- omega's case only, not a general-purpose solver)."""
    a_rx = np.radians(valence_deg(vertex_class, ref_class, x_class))
    a_ry = np.radians(valence_deg(vertex_class, ref_class, y_class))
    a_xy = np.radians(valence_deg(vertex_class, x_class, y_class))
    cos_phi = (np.cos(a_xy) - np.cos(a_rx) * np.cos(a_ry)) / (np.sin(a_rx) * np.sin(a_ry))
    return np.degrees(np.arccos(np.clip(cos_phi, -1.0, 1.0)))

def _dihedral_from_points(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x, y = np.dot(v, w), np.dot(np.cross(b1n, v), w)
    return np.degrees(np.arctan2(y, x))

def _build_two_vectors(a_rx_deg, a_ry_deg, a_xy_deg, sign):
    """ref along +z; x at angle a_rx from ref (azimuth 0); y at angle a_ry
    from ref, azimuth +-arccos(...) from x (sign picks which of the 2
    mirror-image solutions -- chirality-blind on its own, see callers for
    how each resolves it)."""
    a_rx, a_ry, a_xy = np.radians(a_rx_deg), np.radians(a_ry_deg), np.radians(a_xy_deg)
    ref = np.array([0.0, 0.0, 1.0])
    x = np.array([np.sin(a_rx), 0.0, np.cos(a_rx)])
    cos_phi = np.clip((np.cos(a_xy) - np.cos(a_rx) * np.cos(a_ry)) / (np.sin(a_rx) * np.sin(a_ry)), -1.0, 1.0)
    phi = sign * np.arccos(cos_phi)
    y = np.array([np.sin(a_ry) * np.cos(phi), np.sin(a_ry) * np.sin(phi), np.cos(a_ry)])
    return ref, x, y

_symmetric_tetrahedral_cache = {}

def _build_symmetric_tetrahedral_frame(vertex_class, axis_class, branch_class, h_class):
    """A tetrahedral vertex with exactly 3 non-axis substituents, 2 of
    which SHARE an AMBER class (e.g. CB's CG + HB2 + HB3 -- one real
    branch, two chemically-identical H's) -- the common case for every
    aliphatic chi vertex except the backbone C-alpha and Thr's C-beta
    (see formula_derived_deltas). Unlike the backbone stereocenter, the
    mirror-image ambiguity that 3 pairwise angles alone cannot resolve is
    PROVEN irrelevant here, not fixed by a chirality rule: swapping which
    of the 2 identical substituents is placed via the free sign choice
    (h1) vs solved from the other 3 (h2) is exactly a relabelling of two
    atoms that combined_target_for_axis sums over identically (their
    contribution is 1 + 2*cos(n*X), an even function of the mirror
    parameter X -- verified algebraically and, for a concrete case,
    numerically against combined_target_for_axis itself). sign=+1 is used
    arbitrarily; a genuine 3-distinct-class stereocenter (only the
    backbone Cα and Thr's Cβ, among every axis in this file) is NOT this
    case and must not be routed through this function."""
    key = (vertex_class, axis_class, branch_class, h_class)
    if key in _symmetric_tetrahedral_cache:
        return _symmetric_tetrahedral_cache[key]
    a_axis_branch = valence_deg(vertex_class, axis_class, branch_class)
    a_axis_h = valence_deg(vertex_class, axis_class, h_class)
    a_branch_h = valence_deg(vertex_class, branch_class, h_class)
    a_hh = valence_deg(vertex_class, h_class, h_class)
    axis_v, branch_v, h1_v = _build_two_vectors(a_axis_branch, a_axis_h, a_branch_h, +1.0)
    A = np.array([axis_v, branch_v, h1_v])
    b = np.array([np.cos(np.radians(a_axis_h)), np.cos(np.radians(a_branch_h)), np.cos(np.radians(a_hh))])
    h2_v = np.linalg.solve(A, b)
    h2_v = h2_v / np.linalg.norm(h2_v)
    frame = {"axis": axis_v, "branch": branch_v, "h1": h1_v, "h2": h2_v}
    _symmetric_tetrahedral_cache[key] = frame
    return frame

def formula_derived_deltas(vertex_class, axis_class, atoms_with_classes, ref_atom):
    """Attempts to derive delta (deg) for EVERY atom in atoms_with_classes
    (a list of (atom, class) -- every real non-axis substituent at this
    vertex) purely from AMBER's own angle tables, relative to ref_atom
    (matching whatever the caller's own ref_b/ref_c convention already
    is -- no PDB structure of any kind. Returns None if this vertex isn't
    one of the solvable cases:
    - exactly 2 substituents (trigonal/planar vertex): idealized_azimuth_deg,
      forced unambiguous by planarity (see that function).
    - exactly 3 substituents, all 3 sharing one class (a methyl's 3 H's,
      no branch atom at all -- see the "all-same-class" block below):
      pure 3-fold symmetry, no chirality question, no frame to build.
    - exactly 3 substituents with >=2 sharing a class (tetrahedral, one
      real branch + 2 chemically-identical substituents):
      _build_symmetric_tetrahedral_frame, sign-independent (see its own
      docstring).
    A genuine 3-distinct-class tetrahedral stereocenter (only the
    backbone Cα -- handled separately by ca_tetrahedral_delta -- and
    Thr's C-beta, among every axis in this file) returns None here; the
    caller falls back to real-structure measurement for that one case."""
    if len(atoms_with_classes) == 2:
        (a1, c1), (a2, c2) = atoms_with_classes
        delta_a2 = idealized_azimuth_deg(vertex_class, axis_class, c1, c2)
        deltas = {a1: 0.0, a2: delta_a2}
        shift = deltas[ref_atom]
        return {a: d - shift for a, d in deltas.items()}
    if len(atoms_with_classes) == 3:
        classes = [c for _, c in atoms_with_classes]
        if len(set(classes)) == 1:
            # A methyl's 3 H's (or any vertex with 3 non-axis substituents
            # that all share one class): an even MORE symmetric case than
            # the 1-branch+2-identical one below, with no branch atom to
            # build a frame from at all. combined_target_for_axis sums
            # this vertex's contribution over all 3 real substituents, and
            # every one of them looks up the identical (n,k,phase) AMBER
            # term when paired with any given other-side atom (same
            # class) -- so the sum only depends on the SET of 3 azimuths,
            # never on which specific atom gets which value. Placing them
            # at {0, +120, -120} (any assignment) is therefore exact, not
            # an approximation: algebraically, summing
            # exp(i*n*{0,+120,-120}) gives 1+2*cos(n*120deg), which is
            # exactly 3 when n is a multiple of 3 and exactly 0 otherwise
            # -- the vertex's own 3-fold symmetry forces every non-3x
            # harmonic to cancel on its own, matching amber99sb.xml's own
            # choice to parametrize this term as pure n=3 with no n=1/n=2
            # companion (confirmed by direct XML inspection, not assumed).
            a0, a1, a2 = (a for a, _ in atoms_with_classes)
            deltas = {a0: 0.0, a1: 120.0, a2: -120.0}
            shift = deltas[ref_atom]
            return {a: d - shift for a, d in deltas.items()}
        if set(classes) == {"OH", "CT", "H1"}:
            # Threonine's C-beta -- the one genuine 3-distinct-class
            # tetrahedral stereocenter among every chi vertex in this file
            # besides the backbone Cα (see thr_cb_tetrahedral_delta's own
            # docstring for the CIP derivation). Class-matched, not
            # resname-matched, so it fires wherever this exact substituent
            # pattern occurs -- today that is only Thr's CB.
            by_class = {c: a for a, c in atoms_with_classes}
            deltas = {by_class["CT"]: stereocenter_delta("OG1", "CG2"),
                     by_class["OH"]: 0.0,
                     by_class["H1"]: stereocenter_delta("OG1", "HB")}
            shift = deltas[ref_atom]
            return {a: d - shift for a, d in deltas.items()}
        common_class = next((c for c in classes if classes.count(c) >= 2), None)
        h_atoms = [a for a, c in atoms_with_classes if c == common_class]
        branch_atoms = [a for a, c in atoms_with_classes if c != common_class]
        if common_class is None or len(h_atoms) != 2 or len(branch_atoms) != 1:
            return None
        branch_atom = branch_atoms[0]
        h1_atom, h2_atom = h_atoms
        branch_class = next(c for a, c in atoms_with_classes if a is branch_atom)
        frame = _build_symmetric_tetrahedral_frame(vertex_class, axis_class, branch_class, common_class)
        origin = np.zeros(3)
        pts = {"origin": origin, "axis": origin + frame["axis"], "branch": origin + frame["branch"],
              "h1": origin + frame["h1"], "h2": origin + frame["h2"]}
        deltas = {
            branch_atom: 0.0,
            h1_atom: _dihedral_from_points(pts["branch"], pts["origin"], pts["axis"], pts["h1"]),
            h2_atom: _dihedral_from_points(pts["branch"], pts["origin"], pts["axis"], pts["h2"]),
        }
        shift = deltas[ref_atom]
        return {a: d - shift for a, d in deltas.items()}
    return None

def geom_via_far_anchor(L_axis, bond_to_other, angle_at_other_deg):
    """(r, theta) of a reference atom relative to anchor B, when that atom is
    bonded to the OTHER anchor C instead of to B.

    Needed by the PLANARITY impropers: AMBER lists the hub third and bonds
    all three peripherals to it, so the quad's first atom hangs off C, not
    off B, and its distance/angle to B are 1-3 quantities. Plain triangle
    B-C-R with |BC|=L_axis, |CR|=bond_to_other and the angle at C between
    them: the law of cosines gives |BR|, then the angle at B. Both are still
    pure AMBER table values -- no structure is measured.

    (Proper torsions never take this path: their reference atom really is
    bonded to its own anchor.)"""
    g = np.radians(angle_at_other_deg)
    r = float(np.sqrt(L_axis ** 2 + bond_to_other ** 2 - 2.0 * L_axis * bond_to_other * np.cos(g)))
    if r < 1e-9:
        return 0.0, 0.0
    cos_theta = (L_axis ** 2 + r ** 2 - bond_to_other ** 2) / (2.0 * L_axis * r)
    return r, float(np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0))))

n_dihedral_ok = 0
n_dihedral_skip = 0
n_ghost_particles = 0
dihedral_lines = []
ghostparticle_lines = []
# Pure logging side-channel (no effect on anything emitted): one entry per
# emit_ghost_ring call, i.e. one per (resname, axis, harmonic, group) ring
# actually deployed -- lets an external validation script reconstruct
# EXACTLY the same AMBER target / calibrated ring curves this file itself
# used, with zero risk of drifting from the real numbers (see
# doc/BondedForceFieldSprings.md's energy-curve validation section).
AXIS_RING_LOG = []

# (resname, anchors, r, theta, azimuth) -> already-emitted ghost name. Lets
# per-pair ring groups share the ghosts they have in common instead of
# emitting duplicates -- see emit_ghost_ring's own dedup comment.
_ghost_registry = {}

def emit_ghost_ring(resname, axis_label, family, n, L_axis, target_complex, atom_B, atom_C, atom_ref_b, atom_ref_c,
                    axis_dc_target=0.0, group_tag="", ref_geom_b=None, ref_geom_c=None):
    """Emits one full ring group (M+N ghost particles, M*N ghost-ghost
    springs -- M=n,N=1 for n>=2, the ghost/spring-minimal construction,
    see the M,N note below) reproducing one target Fourier harmonic
    exactly (up to the small, bounded leakage quantified in RING_SHAPES's
    own derivation). atom_B/atom_C are the
    real dihedral axis atoms (already +/- resolved); atom_ref_b/atom_ref_c
    are each side's own real reference atom (fixes what "azimuth=0" means
    for that side -- must be the SAME atoms used to measure the target's
    own phase convention, i.e. b_neighbors[0]/c_neighbors[0], see
    combined_target_for_axis).

    axis_dc_target: the real AMBER torsion's own physical mean for this
    axis (target[0] in the caller, sum of each matched term's own k --
    AMBER's PeriodicTorsionForce term k*(1+cos(...)) has mean exactly k).
    Attributed to exactly ONE ring per axis by the caller (every other
    ring for the same axis passes 0.0) so it's subtracted exactly once
    per axis, not once per harmonic. Combined with this ring's own
    (non-physical, construction-artifact) DC into a single per-spring
    dc_offset column: BondedForceFieldReader accumulates this across
    every DIHEDRAL entry it actually applies, giving an exact correction
    (subtracted only when reporting energy, never affecting forces) so
    BioSpring's absolute dihedral energy is directly comparable to
    AMBER's, instead of being inflated by the ring construction's own
    unavoidable positive baseline (see calibrate_ring's own comment).

    group_tag: disambiguates ghost/rule names when the caller emits more
    than one ring group for the SAME axis anchored on different real atoms
    (see generate_backbone_axis) -- must be unique per group sharing an
    axis_label, empty ("", the default) is fine whenever an axis only ever
    has one group (every sidechain/omega axis today).

    Placement/sign convention (verified by a full 3D rotating simulation
    against a real axis before being trusted here -- see
    doc/BondedForceFieldSprings.md): a b-side ghost i sits at azimuth
    i*360/M (measured from atom_B, relative to atom_ref_b, via
    spn::GhostParticle::computePosition(atom_B, atom_C, atom_ref_b, ...)).
    A c-side ghost j sits at azimuth delta_base - j*360/N (note the
    MINUS -- placed via computePosition(atom_C, atom_B, atom_ref_c, ...),
    i.e. B/C swapped since the c-side ghost's own axis direction points
    the other way). This combination reproduces exactly the SAME abstract
    delta_ij = delta_base + beta_i - gamma_j used to calibrate k/delta_base
    in Fourier space above -- get either the swap or the minus sign wrong
    and the real 3D energy has the wrong phase relative to the molecule's
    true dihedral angle, silently, without any local check catching it.

    M,N (found and switched 2026-08-09): only LCM(M,N)=n matters for the
    comb filter (see ring_curve_abstract's own docstring), not M=N=n --
    M=n,N=1 (or the mirror M=1,N=n) reproduces the exact same curve with
    fewer ghosts (M+N=n+1 instead of 2n) and fewer springs (M*N=n
    instead of n*n), verified against a real target to amplitude ratio
    0.9998/phase diff 0.000 deg. Applied system-wide: 984->737 ghost
    particles (-25%), 1166->492 springs (-58%) across the 245 rings this
    file emits. n=1 is already M=N=1 (the minimal case, LCM(1,1)=1, no
    asymmetric variant exists or is needed)."""
    global n_dihedral_ok, n_ghost_particles

    AXIS_RING_LOG.append({"resname": resname, "axis_label": axis_label, "family": family, "n": n,
                          "L_axis": L_axis, "target": target_complex, "axis_dc_target": axis_dc_target,
                          "group_tag": group_tag})

    # Since 2026-08-14 a ghost is the rotated image of its side's real
    # reference atom (spn::GhostParticle::computePositionByRotation), so
    # its radius and axial offset are that atom's own, taken from the AMBER
    # bond/angle tables like everything else here -- never chosen. Only d0
    # and k are still calibrated. ref_geom_b/ref_geom_c are (r, theta) in
    # the same convention as combined_target_for_axis's b_geom/c_geom.
    r, theta_deg = ref_geom_b
    r_c, theta_c_deg = ref_geom_c
    M, N = (n, 1) if n >= 2 else (1, 1)
    d0 = choose_d0(L_axis, n, r, theta_deg, r_c, theta_c_deg, M, N)
    k, delta_base, dc_ring = calibrate_ring(L_axis, n, r, theta_deg, d0, target_complex, M, N, r_c, theta_c_deg)
    # Split evenly across this ring's M*N springs -- BondedForceFieldReader
    # sums dc_offset per spring it actually applies, so this naturally
    # handles a partially-skipped ring too (e.g. Met's phi ghosts missing
    # their "-C" anchor at the true N-terminus): the offset for the springs
    # that never get created is correctly never added either.
    dc_offset_per_spring = (dc_ring - axis_dc_target) / (M * N)

    # group_tag disambiguates ghost/rule names when the SAME axis emits
    # more than one ring group anchored on different real atoms (see
    # generate_backbone_axis's per-owning-pair anchoring) -- without it,
    # two groups sharing the same non-grouped side's reference (e.g. both
    # of psi's O-group and +N-group use the same b-side CB reference) would
    # emit ghosts with IDENTICAL names ("GHpsin3B0" twice), and
    # BondedForceFieldReader's name-based resolution would silently pick
    # whichever one happens to match first -- a real, found-in-review bug,
    # not hypothetical (caught by tracing the actual .bi.ff output before
    # this parameter existed). Sanitized to strip "+"/"-": a ghost's own
    # name is never itself given a cross-residue prefix when referenced
    # later (ghosts are always resolved within their own creating residue).
    tag = re.sub(r"[+-]", "", group_tag)

    # Deduplication: a per-pair axis (see generate_sidechain_axis's per_pair)
    # emits one ring group per real substituent pair, and two groups sharing a
    # side's substituent ask for the SAME ghost -- same anchors, same
    # (r, theta, delta). Omega's n=2 is the clearest case: its 4 groups
    # (CA/CA, CA/H, O/CA, O/H) want only 6 distinct ghosts between them, not
    # 12, because groups CA/CA and CA/H share both their CA-anchored b-side
    # ghosts, and groups CA/CA and O/CA share their +CA-anchored c-side one.
    # A ghost is a massless virtual site fully determined by those 7 fields,
    # so two identical ones would sit at the same point and behave
    # identically -- reusing the first is a pure saving, with the springs
    # simply naming the shared ghost. Measured: 22% of all GHOSTPARTICLE
    # lines were exact duplicates before this.
    # gr/gtheta are the reference ATOM's own radius/angle for THAT side --
    # the two sides generally differ now that a ghost is that atom's rotated
    # image, so they are passed in rather than closed over.
    def ghost(name, res, a1, a2, ref, azim, gr, gtheta):
        global n_ghost_particles
        key = (res, a1, a2, ref, round(gr, 4), round(gtheta, 2), round(azim, 3))
        if key in _ghost_registry:
            return _ghost_registry[key]
        ghostparticle_lines.append(
            f"GHOSTPARTICLE {name:16s} {res:7s} {a1:5s} {a2:6s} {ref:6s} "
            f"{gr:7.4f} {gtheta:7.2f} {azim:9.3f}")
        _ghost_registry[key] = name
        n_ghost_particles += 1
        return name

    b_names, c_names = [], []
    for i in range(M):
        beta = i * 360.0 / M
        b_names.append(ghost(f"GH{axis_label}{tag}n{n}B{i}", resname, atom_B, atom_C, atom_ref_b, beta,
                             r, theta_deg))
    for j in range(N):
        gamma = delta_base - j * 360.0 / N
        c_names.append(ghost(f"GH{axis_label}{tag}n{n}C{j}", resname, atom_C, atom_B, atom_ref_c, gamma,
                             r_c, theta_c_deg))

    for bn in b_names:
        for cn in c_names:
            rule_name = f"{resname}_{axis_label}{tag}_n{n}_{bn[-1]}{cn[-1]}"
            dihedral_lines.append(
                f"DIHEDRAL {rule_name:28s} {resname:7s} {family:9s} {bn:16s} {cn:17s} "
                f"{d0:7.4f} {k:10.4f} {dc_offset_per_spring:10.5f}")
            n_dihedral_ok += 1

def emit_ghost_rings_for_axis(resname, axis_label, family, L_axis, target, dc_by_harmonic, atom_B, atom_C,
                              atom_ref_b, atom_ref_c, ref_geom_b, ref_geom_c):
    """Emits one ring (see emit_ghost_ring) per non-zero harmonic in
    `target` (n -> complex Fourier coefficient, see combined_target_for_axis),
    each calibrated against ITS OWN physical DC share (dc_by_harmonic[n] --
    NOT the axis's combined total attributed to a single ring: found and
    fixed 2026-08-05, see combined_target_for_axis's own comment for why
    that broke each ring's own curve even though it left the axis grand
    total correct) -- shared by generate_sidechain_axis, its GKinase
    variant, and generate_backbone_axis instead of repeating this
    bookkeeping 3 times."""
    for n, zt in target.items():
        if n == 0 or abs(zt) < 1e-6:
            continue
        emit_ghost_ring(resname, axis_label, family, n, L_axis, zt, atom_B, atom_C, atom_ref_b, atom_ref_c,
                       dc_by_harmonic.get(n, 0.0), ref_geom_b=ref_geom_b, ref_geom_c=ref_geom_c)

def combined_target_for_axis(b_geom, c_geom, b_class, c_class, class_of):
    """b_geom/c_geom: {atom: (r, theta, delta_deg)} for each side's real
    substituents (delta already phase-corrected relative to that side's
    own reference atom, i.e. b_neighbors[0]/c_neighbors[0] -- ref_phi
    itself cancels out of the pairwise delta and is never needed, see
    doc/BondedForceFieldSprings.md). class_of(atom) resolves an atom to
    its AMBER class. Tries each real substituent pair's specific AMBER
    entry first, falls back to the axis's generic wildcard (matches
    amber99sb.xml's own matching order) -- a pair whose specific entry
    exists but carries k=0 for every periodicity (a genuine AMBER zero,
    e.g. chi2 for Asp/His/Phe/Tyr, chi3 for Glu, chi4 for Arg, chi2 for
    Trp -- all confirmed by direct XML inspection) correctly contributes
    nothing, not a gap.

    Returns (target, dc_by_harmonic): target[0] is still the axis's own
    total real DC (every matched term's own k, summed) for logging/sanity;
    dc_by_harmonic[n] is the sub-total contributed by matched terms whose
    OWN periodicity is exactly n -- i.e. the true physical DC that
    harmonic n's own ring should be calibrated against (see emit_ghost_ring
    -- every AMBER PeriodicTorsionForce term k*(1+cos(n*phi-phase)) has
    DC=k=its own amplitude, a 1:1 ratio for every n; a ring's own DC/
    amplitude ratio need not be 1:1 by construction -- e.g. the n=1 ring's
    exact-purity shape (d0=0) has a fixed, target-independent ratio of
    2:1 -- so subtracting the WRONG ring's own share, or the whole axis's
    combined DC, corrects the grand total but not each ring's own
    curve, found and fixed 2026-08-05: see doc/BondedForceFieldSprings.md)."""
    target = {0: 0j}
    dc_by_harmonic = {}
    any_term = False
    for bn, (r1, t1, d1) in b_geom.items():
        for cn, (r2, t2, d2) in c_geom.items():
            terms = lookup_torsion_specific(class_of(bn), b_class, c_class, class_of(cn))
            if terms is None:
                terms = lookup_torsion_wildcard(b_class, c_class)
            if terms is None:
                continue
            any_term = True
            delta_pair_rad = np.radians(d2 - d1)
            for (n, k, phase) in terms:
                target[n] = target.get(n, 0j) + k * np.exp(-1j * phase) * np.exp(1j * n * delta_pair_rad)
                target[0] += k
                dc_by_harmonic[n] = dc_by_harmonic.get(n, 0.0) + k
    return (target, dc_by_harmonic) if any_term else (None, None)

