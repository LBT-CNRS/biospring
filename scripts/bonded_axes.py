#!/usr/bin/env python3
"""Look a torsion up in a force field's own table, the way AMBER matches it.

All that is left of what used to be a much larger module: the geometry it
carried -- idealised azimuths, tetrahedral frames, CIP stereocentres, the comb
theorem and the ring construction itself -- existed to place ghost particles,
and went with them.

Call configure() once with the tables the caller has parsed, then use the two
lookups. Both are shared verbatim between the protein and nucleic sides, which
is the only reason this is a module of its own.
"""

# --- injected by configure() ------------------------------------------------
bond_params = {}
angle_params = {}
torsion_params = {}


def configure(bond_params_, angle_params_, torsion_params_):
    """Point the lookups at one force field's tables.

    bond_params:    {frozenset((id1, id2)): (length_nm, k)}
    angle_params:   {(vertex_id, frozenset((id1, id2))): (angle_deg, k)}
    torsion_params: {(id1, id2, id3, id4): [(n, k, phase_rad), ...]}, "" wild
    """
    global bond_params, angle_params, torsion_params
    bond_params = bond_params_
    angle_params = angle_params_
    torsion_params = torsion_params_


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
    # missing table), so it never survives the caller's `k != 0` filter and
    # correctly never matches here either.
    for key in ((c1, c2, c3, c4), (c4, c3, c2, c1)):
        if key in torsion_params:
            return torsion_params[key]
    return None
