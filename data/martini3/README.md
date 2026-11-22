# Martini 3 in BioSpring

Configuration for the **Martini 3.0.0** coarse-grain protein model. Unlike
`data/martini/`, which assumes `martinize` has already built the beads, this
set makes BioSpring do the coarse-graining itself from an all-atom structure.

Built and measured against the code at **`69c277c`**, "Give the springs
AMBER's own bond stiffness, creating none". That is the earliest commit on the
rigid-body line carrying both features this force field needs:
`parse_residue_field` in the reduce rule reader, which arrived at `68290b8`
and makes `:NTER`/`:CTER` rules possible, and the `STIFFNESS` record in the
bonded reader, which arrives here. It stops before the sixteen
hydrogen-bond-cone, desolvation and power-stroke commits that follow, none of
which this work needs.

`feature/rigidbody-dihedral-dynamicbond` (tip `0b2ef7b`) is the OLDEST point
of that same line, 39 commits back, and neither feature exists there: a
`:NTER` rule would be filed under the literal residue name `"VAL:NTER"` and
`STIFFNESS` would be rejected as an unknown entry type.

## Files

| file | what it is | status |
|---|---|---|
| `reducerules/martini3.grp` | atom → bead grouping rules, AMBER atom names, heavy atoms | 0.028 Å from Martini's beads in the median |
| `forcefield/martini3.nbi.ff` | per-type charge / radius / epsilon / mass | charge, radius, mass, termini exact; **epsilon cannot be** |
| `forcefield/martini3.ff` | byte-identical copy of the above | keep the two in sync |
| `prepare_martini3_data.py` | generator; every number's provenance is in its docstring | — |

Not written yet: the bonded set (`.rbody` + `.bi.ff`) and a `.msp`. What it
would take is in **The bonded terms** below.

## Provenance

* `martini_v3.0.0.itp` (Marrink lab `martini-forcefields`, file dated
  2021-03-29): **844 bead types, 355746 `[ nonbond_params ]` pairs**.
* vermouth / martinize2 force field `martini3001`: `aminoacids.ff` for bead
  types and charges, `modifications.ff` for the terminal ones,
  `mappings/martini3001/*.amber.map` for the atom → bead composition,
  `general.ff` for the mapping weight.

The repository already vendors Martini **3.0.b.3.2** under `data/martini/`.
That beta and 3.0.0 differ in the mapping on exactly two residues: in the
beta, alanine's CB sits in BB (one bead, no side chain), and lysine counts
its three ammonium hydrogens. Everything else is identical. 3.0.0 was chosen
because it is the published version and because vermouth ships AMBER-named
mappings only for it — the beta has CHARMM names only, and BioSpring's atom
names are AMBER's.

## The mapping: 0.028 Å, and why it is not 0

**A Martini bead is a MASS-weighted average of the atoms mapped onto it.**
`martinize2` runs vermouth's `DoAverageBead`, which reads the force field
variable `center_weight` — `"mass"` in `martini3001`'s `general.ff` — and
computes `np.average(positions, weights = mapping_weight * atom_mass)`. A
BioSpring grain is `measure::centroid`, an **unweighted** mean, and nothing in
a `.grp` can weight anything.

The two therefore agree exactly only when the mapped atoms have equal masses,
and approximately when they have *similar* masses. Heavy atoms do (12, 14,
16); hydrogens emphatically do not. In a mass-weighted average a hydrogen
pulls with 1/13 of a carbon's weight, while an unweighted centroid gives it
the same weight as a carbon. Measured against the true mass-weighted bead over
the 401 beads of `example/001_GKinase/model.pdb`:

| what the `.grp` lists | median | p90 | worst |
|---|---|---|---|
| **heavy atoms only** *(what this file does)* | **0.028 Å** | 0.102 | 0.415 |
| Martini's weight-1 atoms, hydrogens included | 0.342 Å | 0.537 | 0.730 |

So dropping the hydrogens is not an approximation of Martini — it is a
**twelvefold improvement** on including them, and it means the input does not
need to be protonated at all.

The residual worst case, 0.415 Å, is on the sulfur-bearing beads: `CYS SC1` is
CB at 12 Da and SG at 32 Da, so Martini's centre of mass sits 0.41 Å away from
the midpoint BioSpring computes. Nothing in a `.grp` can fix that — an atom
name cannot be repeated to fake a weight, because `ReduceRule` stores its
atoms in an `unordered_set`. Closing it would mean teaching the reducer a
mass-weighted centroid, which is a one-line change in `Reducer.h` and would
bring every bead to 0.

Verified on `example/001_GKinase/model.pdb` (183 residues, 2814 atoms):

* **401 beads built of 401 expected**, nothing dropped.
* Positions **0.028 Å** from the true Martini beads in the median, 0.415 Å at
  worst, as above.
* **No heavy atom claimed by no grain.** Every unclaimed atom is a hydrogen,
  which is correct here: hydrogens are not part of any bead.
* The same 401 beads whether the file spells its N-terminal hydrogens
  `H1/H2/H3` or `1H/2H/3H`.
* Radii 2.6378 / 2.3010 / 1.9082 Å on 158 / 175 / 68 beads — the regular,
  small and tiny size classes.
* On the unprotonated maltoporin (`1AF6.pdb`, 10050 atoms): 3069 beads and
  **zero unclaimed atoms**.

**Charges, masses and termini**, from `aminoacids.ff` and
`modifications.ff`: 22 beads at +1 e (17 ARG, 4 LYS, 1 N-terminus) and 22 at
−1 e (10 ASP, 11 GLU, 1 C-terminus), net **+0.000 e**. A free N-terminus is
`Q5` at +1 e and a free C-terminus `Q5` at −1 e, reached through the qualified
force-field keys `<grain>:NTER` and `<grain>:CTER` that
`ReduceRule::forcefield_key` asks for. Histidine's protonated form splits
+0.5 e over SC2 and SC3.

### One caveat, and it is the engine's, not the file's

`topology::is_nterminal_only_atom` matches `H1`, `H2`, `H3` and nothing else,
so a structure spelling its terminal ammonium `1H/2H/3H` is not RECOGNISED as
a free N-terminus and its `:NTER` rules never apply. The bead is still built —
the interior rule lists only heavy atoms, so it matches at any chain position
— but it keeps the interior type, so it carries 0 e where Martini has +1 e.
Measured on `example/001_GKinase/model.pdb`: net **−1.000 e** instead of
+0.000 e, the C-terminus being detected normally because `OXT` is spelled the
same everywhere. Accepting `1H/2H/3H` in that one function is a one-line
change that can affect no AMBER-named file.

## What cannot be reproduced, and by how much

### Epsilon: Martini stores it per PAIR, a `.ff` stores it per TYPE

This is the real limit, and no configuration file can lift it. Martini 3
gives every pair of bead types its own σ and ε — 355746 of them, 276 over
the 23 types a protein uses. BioSpring reads one line per type and combines
two of them with a **hard-coded geometric mean** on ε
(`biospring_lorentz_berthelot_epsilon`; no `.msp` setting reaches it).

Martini's matrix is not of the form √(ε_i·ε_j), so the gap is structural.
Measured over the 253 cross pairs of the protein subset:

| per-type ε | pairs exact | median | p90 | worst |
|---|---|---|---|---|
| ε_i = ε_ii, the diagonal *(what this file ships)* | 9 (3.6 %) | 0.626 | 2.436 | **3.300 kJ/mol** |
| ε_i fitted, best rank-1 in log space | 3 (1.2 %) | 0.342 | 1.367 | **4.284 kJ/mol** |

The four worst cases under the diagonal are all polar-with-charged pairs and
all **underestimated**: P5–SQ5n is 7.58 kJ/mol in Martini and 4.28 by the
geometric mean, SP2a–SQ4p 6.57 against 3.37, SP2a–SQ3p 6.38 against 3.19,
P2–SQ5n 7.16 against 4.07. That pairing is what a salt bridge is made of.

Fitting was measured and rejected: it halves the median but has to put SQ3p
at 7.80 kJ/mol where Martini says 3.31, because Q–Q pairs are weak while
Q–P pairs are strong and one number per type cannot be both. The diagonal at
least never states a value Martini does not hold.

**Closing this gap needs code, not configuration**: a per-pair table read
alongside the `.ff`, consulted by `biospring_steric_*` instead of the
combining rule, and mirrored in the OpenCL kernel. 276 lines would cover
proteins.

### Radius: exact on same-size pairs, 0.100 Å out on mixed ones

BioSpring's 12-6 is `ε(R/r)¹² − 2ε(R/r)⁶`, whose minimum is at `r = R`, so
the radius column is **rmin/2**, not σ/2 — it is written as
`2^(1/6)·σ_ii/2`. Two of them are summed, which gives the arithmetic mean of
the two σ. Martini's own σ, on the protein subset:

| size classes | sum rule | Martini | error |
|---|---|---|---|
| R–R | 4.700 | 4.700 | 0.000 |
| S–S | 4.100 | 4.100 | 0.000 |
| T–T | 3.400 | 3.400 | 0.000 |
| R–S | 4.400 | 4.300 | **+0.100** |
| R–T | 4.050 | 3.950 | **+0.100** |
| S–T | 3.750 | 3.650 | **+0.100** |

Martini subtracts a flat 0.100 Å whenever the two sizes differ, which a sum
of two per-type radii cannot do. 174 of the 253 cross pairs carry it: 2.5 %
of σ.

### TRP SC3 is a virtual site in Martini and a real bead here

The only `[ virtual_sitesn ]` entry in all of `martini3001`. BioSpring has no
virtual sites, so SC3 is mapped onto the two ring-fusion carbons it sits
between and carries its size class's 36 Da. A tryptophan therefore weighs
36 Da more than Martini's, and has three translational degrees of freedom
Martini does not give it.

## The bonded terms

Not written yet. Martini 3's bonded set for the 22 standard residues is
small — 45 `bonds`, 26 `constraints`, 14 `angles`, 2 `dihedrals`, 1 virtual
site and 21 `exclusions`, with bond constants from 2500 to 7500
kJ·mol⁻¹·nm⁻² (25 to 75 kJ·mol⁻¹·Å⁻²), angles from 25 to 60 kJ·mol⁻¹ and the
two dihedrals at 50 and 100. The backbone terms are separate, in the 39
`[ link ]` blocks of `aminoacids.ff`. On this base the route is:

* **Which springs exist** comes from `--cutoff` or from a `.rbody` group
  file. A `.rbody` line creates every pairwise spring between the atoms it
  lists and leaves the particles distinct, so a two-atom group is exactly one
  spring — which is the control a faithful bonded topology needs, rather than
  letting a distance cutoff invent the network. Cross-residue `+`/`−` atom
  prefixes work there too (`RigidBodyBuilder.cpp`), so BB(i)–BB(i+1) is
  addressable as `+BB`.
* **Force constants: expressible, per pair.**
  `STIFFNESS <resname>[:NTER|CTER] <atom1> <atom2> <k>` takes k in
  kJ·mol⁻¹·Å⁻² and accepts the same `+`/`−` prefixes. Martini's constants are
  in kJ·mol⁻¹·nm⁻²: divide by 100.
* **Bond equilibrium lengths: not expressible.** `STIFFNESS` retunes springs
  that already exist; a spring's equilibrium is the distance in the input
  structure. Martini's fixed lengths (0.270, 0.315, 0.350 nm …) would not be
  imposed; the input geometry's would stand. Nor is there anything to express
  a `[ constraints ]` entry, which Martini uses for most ring side chains and
  which is rigid rather than stiff — though Martini itself offers a
  `#ifdef FLEXIBLE` alternative for every one of them, at `$stiff_fc` =
  1000000 kJ·mol⁻¹·nm⁻², i.e. 10000 kJ·mol⁻¹·Å⁻². That is fifteen times the
  650 this project settled on for a rigid-body mesh, so adopting it is a
  timestep decision as much as a stiffness one.
* **Angles.** Martini uses the cosine-harmonic form (GROMACS angle type 2).
  BioSpring holds an angle as a 1-3 distance spring, so each one needs a
  stiffness conversion that depends on the equilibrium geometry. That
  conversion has to be calibrated and measured, not assumed.
* **Dihedrals: expressible** as `TORSION` / `TORSIONTABLE` records, the same
  way `data/reducerules/ProteinAtomTorsion.bi.ff` tabulates AMBER's.

## Two upstream inconsistencies, and why they stopped mattering

The `*.amber.map` files contradict each other and themselves about hydrogens:
eighteen of the twenty write `HB1 !SC1` (the `!` means null weight) while
`phe.amber.map` writes `HB1 SC1` and `trp.amber.map` writes `HB2 !SC2` for an
atom bonded to an SC1 carbon; and the four histidine maps disagree about where
a ring proton goes — `his` omits HD1 and HE1 altogether, `hid` puts HD1 in SC1
while CG is in SC2, and `hie` lists HD1 on a tautomer with no delta proton.

Since this mapping lists **no hydrogens at all**, none of that reaches the
generated rules. It is recorded here because it is the reason to distrust any
hydrogen-level reading of those files.

What does still matter is the heavy-atom skeleton of histidine. `HIS` is taken
as AMBER's **delta-protonated** form, with `martini3001`'s HSD/HID bead types,
because that is what both local authorities carry under the name:
`data/reducerules/amber.grp` lists HD1 and no HE2, and
`example/001_GKinase/model.pdb` does the same. `HIE` and `HIP` are declared
separately, and they differ from `HIS` in their bead TYPES, not in their atoms.

## Usage

```
pdb2spn -s protein.pdb -o model.nc \
        --grp data/martini3/reducerules/martini3.grp \
        --ff  data/martini3/forcefield/martini3.nbi.ff \
        --cutoff 6.0
```

`impala.enable` must stay 0: the `transferIMP` column is 0.0 on every line,
because no oil/water partition dataset was published with Martini 3.0.0. The
numbers in `data/martini/forcefield/martini_{octanol,hexadecane}.ff` belong
to the 3.0.b.3.2 bead names, not to these.
