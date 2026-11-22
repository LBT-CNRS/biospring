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
| `reducerules/martini3.rbody` | bonded topology: which springs exist, 171 groups of two atoms | Martini's graph exactly, minus the ±2/±3 terms |
| `reducerules/martini3.names.grp` | naming table for `-rigidbody`, one atom per rule | required, see below |
| `forcefield/martini3.nbi.ff` | per-type charge / radius / epsilon / mass | charge, radius, mass, termini exact; **epsilon cannot be** |
| `forcefield/martini3.ff` | byte-identical copy of the above | keep the two in sync |
| `prepare_martini3_data.py` | generator; every number's provenance is in its docstring | — |

Not written yet: the `.bi.ff` carrying Martini's own force constants, and a
`.msp`.

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

`reducerules/martini3.rbody` carries Martini's bonded **topology** — which
springs exist. It carries no force constant, and none belongs in it: a
spring's stiffness is a construction input (`pdb2spn -stiffness`, stored per
spring in the `.nc`, unreachable from the `.msp`), and the only option that
reads per-pair constants from a `.bi.ff`, `-bondstiffness`, is off by default.

A `.rbody` line creates every pairwise spring between the atoms it names and
leaves the particles distinct, so a **two-atom group is exactly one spring**.
All 171 groups have two atoms, which is what makes the network Martini's graph
rather than a distance cutoff's guess at it.

What is in it, from `martini3001`: every side-chain `[bonds]` and
`[constraints]` pair, deduplicated — the two sections declare the same pair
twice for most side chains, once under `ifndef FLEXIBLE` and once under
`ifdef FLEXIBLE`; the 1-3 pair of every side-chain `[angles]` entry, which is
how BioSpring holds an angle; and from the 39 `[ link ]` blocks, the backbone
bond BB–+BB plus the 1-3 pairs of the backbone angles.

### Verified

Two-step build on `example/001_GKinase/model.pdb`, 183 residues:

* **982 springs created, 982 expected** from an independent count over the
  sequence.
* The graph, read back from the output's CONECT records:

  | class | n | median length |
  |---|---|---|
  | BB(i)–BB(i+1) | 182 | 3.23 Å |
  | BB(i)–BB(i+2) | 181 | 5.09 Å |
  | BB–SC1 intra | 168 | 2.73 Å |
  | SC1(i)–BB(i+1) | 168 | 5.24 Å |
  | BB(i)–SC1(i+1) | 167 | 4.80 Å |

  182 is every consecutive pair of the 183 residues, 181 every 1-3 pair, and
  168 is 183 minus the 15 glycines. The ring side chains account for the rest.

* The BB–BB bonds run from **2.84 to 3.78 Å**. Martini declares three
  alternatives — 3.10 Å (helix), 3.375 and 3.50 Å — and picks between them by
  DSSP. The input's own distances already span that range, which is the
  measurement behind the next paragraph.

### What cannot be in it

* **The secondary-structure choice**, and it costs nothing. A `.rbody` keys on
  residue name and chain position, never on secondary structure. But a
  BioSpring spring equilibrates at the distance in the INPUT structure, which
  in a helix already *is* 3.10 Å — the same measurement Martini's table came
  from, as the 2.84–3.78 Å spread above shows. The force constant is 4000
  kJ·mol⁻¹·nm⁻² on all three bond variants anyway, so only the
  bond-or-constraint distinction varies, and that is the construction choice.
* **`BB ++BB` and `BB +++BB`**, the elastic bonds that stabilise extended
  regions and helices. `RigidBodyBuilder` strips ONE `+`/`−` prefix
  (`name.substr(1)`) and moves exactly one residue, so a two-residue offset
  resolves to an atom named `+BB` in the next residue, which does not exist.
  Four of the backbone bond patterns and two of the five dihedral patterns are
  lost this way.
* **`SC1 >SC1`**, the disulfide bridge, which is not between neighbours at
  all. `-static-disulfide` covers it, at 1389.1 kJ·mol⁻¹·Å⁻².
* **Martini's equilibrium lengths.** A spring equilibrates where the input
  puts it. For the backbone that is a feature, as above; for a side chain it
  means 0.270 / 0.315 / 0.341 nm are not imposed.

### TRP SC3 gets four springs Martini does not give it

There, SC3 is a virtual site placed from SC1/SC2/SC4/SC5 and carries no bonded
term at all. Here it is a real bead, so with no spring it would be held by
nothing but the non-bonded terms and would drift out of the ring. It is tied
to the four beads its virtual-site definition uses.

### Still to write: the `.bi.ff`

Martini states a force constant for only **19 of the 59 side-chain pairs** —
25, 50 or 75 kJ·mol⁻¹·Å⁻² — and 40 kJ·mol⁻¹·Å⁻² for the backbone bond. The
other 40 pairs are `[constraints]`, which are rigid and carry **no constant at
all**; Martini's `#ifdef FLEXIBLE` fallback of 10000 kJ·mol⁻¹·Å⁻² is a
workaround for engines without constraints, not a force-field value, so
nothing should be written for them — they take the run's `-stiffness`. The
angles need one more step: Martini's is the cosine-harmonic form (GROMACS
type 2 or the restricted type 10) and BioSpring's is a 1-3 distance spring, so
the conversion has to be derived and measured, not assumed.

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

**Non-bonded only**, springs from a distance cutoff:

```
pdb2spn -s protein.pdb -o model.nc \
        --grp data/martini3/reducerules/martini3.grp \
        --ff  data/martini3/forcefield/martini3.nbi.ff \
        --cutoff 6.0
```

**With Martini's bonded topology**, which takes two steps. `--rigidbody`
needs a naming table with exactly one atom per rule
(`_check_translation_is_one_atom_per_rule` kills pdb2spn otherwise), so
`martini3.grp` cannot serve as that table and `martini3.names.grp` exists for
it. Passing the names table to a run that would REDUCE destroys the model —
`ABB ALA BB` looks for an atom named `BB`, which no all-atom structure has —
so the reduction happens first, and the second run takes the `.nc` it wrote,
where pdb2spn skips reduction and uses the table only for naming:

```
pdb2spn -s protein.pdb -o beads.nc --grp data/martini3/reducerules/martini3.grp --ff data/martini3/forcefield/martini3.nbi.ff
```

```
pdb2spn -s beads.nc -o model.nc --rigidbody data/martini3/reducerules/martini3.rbody --grp data/martini3/reducerules/martini3.names.grp --ff data/martini3/forcefield/martini3.nbi.ff
```

`--ff` is mandatory whenever `--grp` is given, even in the second run where it
is ignored. `--rigidbody` raises the default `-stiffness` to 650
kJ·mol⁻¹·Å⁻²; pass `-stiffness` explicitly to choose otherwise.

`impala.enable` must stay 0: the `transferIMP` column is 0.0 on every line,
because no oil/water partition dataset was published with Martini 3.0.0. The
numbers in `data/martini/forcefield/martini_{octanol,hexadecane}.ff` belong
to the 3.0.b.3.2 bead names, not to these.
