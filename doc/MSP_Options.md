BioSpring Simulation Parameter File (.msp file)
===============================================

Introduction
------------

The BioSpring simulation parameter file describes the parameters and settings of
the simulation. Its extension is `.msp`.

The syntax is:

    [Setting].[Parameter]<space>=<space>value

Example:

	viscosity.enable = 1
    
The parameters can be written in any order and none are compulsory as each
parameter has a default value.

You can comment a line using the # character.

Two mistakes are fatal rather than ignored, each reported with its line number: a parameter name
that does not exist (`SafeConfigurationReader: line N: invalid parameter "..."`) and the same
parameter given twice. A typo therefore never silently runs a different simulation -- and a group
this manual marks as *not implemented* cannot be left in a file "harmlessly".

Example of .msp file for BioSpring
----------------------------------

The following example file parameterizes a simulation run with a *2 fs* timestep
for *1000* steps. The viscosity is set to *0.1 Da.fs-1* and the steric 
interaction is set to the linear mode with a 8Å cutoff.These parameters are 
describe in the next section of this manual.

    #001_GKinase - Molecular Simulation Parameter file

   	#fs
   	simulation.timestep = 2.0
   	simulation.nbsteps = 1000
   	
   	#Da.fs-1
	viscosity.enable = 1
	viscosity.value = 0.1

   	#kJ.mol-1.A-2
	spring.enable = 1
	spring.scale = 1

	#A
	steric.enable = 1
	steric.gridscale = 1
	steric.cutoff = 8.0
	steric.mode = linear

General Simulation Parameters Description
-----------------------------------------

Define the steps and IO settings. The values written here are the default values.

* **simulation.timestep = 0.01** *(fs, float)* Timestep used in this run. The default is
`defaultConfiguration()`'s own value and not a recommendation: it sits far below anything usable,
so an `.msp` that forgets this line crawls instead of exploding. Every real run sets it -- see
*Choosing `--stiffness` and `simulation.timestep`* below for the measured values, **3.0 fs** on a
protein and **2.0 fs** on a nucleic acid, both at `--stiffness 650` kJ.mol-1.A-2.
* **simulation.nbsteps = -1** *(integer)* Number of simulation steps. The run stops after this step.
-1 defines an infinite run and you will have to kill the process manually. 
* **simulation.samplerate = 100** *(integer)* Frequence at which energies are printed on the standard
output.
* **simulation.neighborskin = -1** *(A, float)* Extra margin added to the steric, electrostatic
and hydrophobic cutoffs when a term builds its list of neighbour pairs. The list holds every
pair within `cutoff + skin`, so it stays valid until a particle has moved half the margin, and
is reused until then instead of being rebuilt. A negative value, the default, lets the backend
choose: **0.5 A on the CPU** and **1.0 A on the GPU**. `0` is a real setting and turns reuse
off, so every step rebuilds; on the CPU it falls back to walking the cells.

Reuse matters more than the exact margin, and it matters on both backends. On the GPU, median
steps/s over three runs of 1000 steps, against the same run at `0`:

| skin (A) | 023.Nucleosome | 024.CoarseGrain | 034.VirusCA | 042.FepA |
|---|---:|---:|---:|---:|
| 0 | 158.98 | 301.20 | 512.82 | 523.56 |
| 0.5 | +27.1 % | +35.5 % | **+38.3 %** | +45.8 % |
| **1.0** | **+30.0 %** | **+43.1 %** | +37.3 % | +49.2 % |
| 2.0 | +22.4 % | +40.1 % | +36.4 % | **+52.8 %** |

The optimum moves with how fast the structure drifts -- 023 rebuilds 165 times in 1000 steps at
1 A where the capsid rebuilds 9 -- but 1.0 A is the best of the four on two of them and within
3.5 % of it on the others, which is why it is the default.

The CPU wants a skin too, and a narrower one, so the two defaults differ on purpose:

| skin (A) | 023.Nucleosome | 024.CoarseGrain | 034.VirusCA | 042.FepA |
|---|---:|---:|---:|---:|
| 0 | 46.47 | 338.98 | 231.48 | 3846.15 |
| **0.5** | **+2.1 %** | **+4.6 %** | +57.1 % | +4.0 % |
| 1.0 | +0.6 % | +2.8 % | +58.8 % | +4.0 % |
| 2.0 | -5.5 % | +1.0 % | **+60.6 %** | +4.0 % |

0.5 A is the outright best on two of them and within 1.7 % on a third; 042's framerate print
quantises at 4000 steps/s, so its three margins are indistinguishable. Note 034: on the CPU any
margin at all is worth 57 to 61 %, which is what makes a list worth building there.
---
* **thermostat.enable = 0** *(boolean)* Turns the Langevin thermostat on. Off, the default,
`viscosity.value` is a pure brake: energy only ever leaves the system, so a run is a viscous
relaxation and not a dynamics at any temperature. Measured on 072.GK, the kinetic energy falls
from 15.93 to 0.63 kJ.mol-1 over 4000 steps -- 300 K on that system is around 10500.

  Its friction IS `viscosity.value`, not a second coefficient, so `viscosity.enable` must be on
too and setup refuses the combination otherwise. Friction and the random kicks are one physical
effect -- the solvent this model does not simulate -- and the fluctuation-dissipation relation
fixes the size of one from the other. Raising the viscosity therefore no longer just calms a
run: it couples it harder to the bath, in both directions, and the temperature it settles at is
unchanged.

* **thermostat.temperature = 300** *(K, float)* The bath temperature, used only when the
thermostat is enabled.

  The velocity update is the exact solution of the Langevin equation over one step,
`v <- c.v + sigma.xi` with `c = exp(-gamma.dt/m)`, not a first-order approximation of it, so the
temperature reached does not drift with the timestep. Measured on 042.FepA, asking for 300 K:

  | timestep (fs) | 2.0 | 1.0 | 0.5 | 0.25 |
  |---|---:|---:|---:|---:|
  | kinetic energy (kJ.mol-1) | 834.9 | 815.4 | 824.1 | 809.6 |

  The noise comes from a counter-based generator keyed on (step, particle, component), so a run
is reproducible from its seed, an OpenMP loop needs no per-thread state, and the CPU and the
device draw the SAME sequence -- which is what keeps the two backends comparable to the digit.

  **A timestep that was stable without the thermostat may not be with it.** The damping was
bleeding off the energy that makes a structure move; put it back and the integrator has to cope
with it. 072.GK ships with 3.0 fs and diverges at 300 K around step 3000, while 1.5 fs holds at
312 K indefinitely. The reported `Temperature:` line is there to be watched.
---
* **simulation.cellsize = 0** *(A, float)* Width of one neighbour-grid cell, for every term at
once. `0`, the default, gives each term cells the size of its own search radius, so each walks
the 27 cells around its own. Set it smaller to trade more cells walked for fewer candidates
tested, which pays on a dense structure: 041.Alphagalactosidase runs at 778 steps/s on the GPU
with the default and 844 at half its Coulomb cutoff. Set it larger on a hollow one, where most
cells are empty and visiting them is the cost. There is no value that suits every structure,
which is why this is a knob and not a heuristic.
---
* **pdbtrajectory.enable = 0** *(boolean)* Enables trajectory writing in pdb format.
* **pdbtrajectory.frequency = 100** *(integer)* Frequence at which frames are written.
* **pdbtrajectory.path = ""** *(string)* Name of the pdb trajectory file.
---
* **csvsampling.enable = 0** *(boolean)* Enables energies logging in csv format.
* **csvsampling.frequency = 100** *(integer)* Frequence of energy logging.
* **csvsampling.path = ""** *(string)* Name of the csv energies log.
---
* **xtctrajectory.enable = 0** *(boolean)* Enables trajectory writing in xtc format.
* **xtctrajectory.frequency = 100** *(integer)* Frequence at which frames are written.
* **xtctrajectory.path = ""** *(string)* Name of the xtc trajectory file.

The three writers share one rule: a frame goes out when the step counter is an exact multiple of
`frequency`, and that counter starts at **0**. So step 0 always produces a frame, and a run of
exactly `frequency` steps produces only that one -- `simulation.nbsteps = 1000` with
`frequency = 1000` writes a single-frame trajectory of the input structure.

`rigidbody.enable` and `rigidbody.enablesampling` together override `csvsampling.frequency` to
1000000 at setup: the line is then written at the end of every rotation about the insertion
vector, not on a step count.

Spring Network Parameters Description
-------------------------------------

Define how the spring network behavior. More details are given in the examples.

* **spring.enable = 0** *(boolean)* Enable spring forces.
* **spring.scale = 1.0** *(dimensionless factor, float)* Multiplier applied to the per-spring
stiffness (itself in kJ.mol-1.A-2, set at topology creation time by pdb2spn/editspn/mergespn's
`--stiffness` option). See doc/User_Manual.md for the usual literature value (~0.6 kcal.mol-1.A-2,
i.e. ~2.5 kJ.mol-1.A-2).
* **spring.cutoff = 15.0** *(Angstroms, float)* Currently unused: springs are created once, ahead of
time, by pdb2spn/editspn/mergespn's own `--cutoff` option (see there), not rebuilt at runtime from
this value.
---
---
* **dihedralphi.enable = 1** *(boolean)* Runtime on/off for the phi (backbone) torsions,
independently of psi/omega/chi. Only meaningful if the topology was actually built with
`-dihedralbackbone` (or `-dihedral`) in the first place -- a family not built has no torsions to
enable/disable either way. Defaults to enabled so an `.msp` written before this setting existed
keeps the same behaviour.
* **dihedralpsi.enable = 1** *(boolean)* Same as `dihedralphi.enable`, for the psi axis.
* **dihedralomega.enable = 1** *(boolean)* Same as `dihedralphi.enable`, for the omega
(peptide-bond) axis.
* **dihedralchi.enable = 1** *(boolean)* Same as `dihedralphi.enable`, for every side-chain
chi1-4 dihedral (the SIDECHAIN family in the `.bi.ff`).
* **dihedralplanarity.enable = 1** *(boolean)* Same as `dihedralphi.enable`, for the PLANARITY
impropers.
* **dihedralnucleicbackbone.enable = 1** *(boolean)* Same, for a nucleotide's alpha..zeta,
and for the 5'-OH and 3'-OH rotors of a chain end: those turn about the beta and epsilon axes
themselves, so they belong here rather than with the 2'-OH, which turns about C2'-O2'.
* **dihedralnucleicchi.enable = 1** *(boolean)* Same, for the glycosidic torsion and RNA's
2'-OH rotor.
* **dihedralnucleicsugar.enable = 1** *(boolean)* Same, for the four furanose ring bonds. A
family of its own rather than part of the nucleic backbone, because those bonds carry the sugar
pucker -- the single lever choosing the A or B helical form -- so its energy has to be readable,
and switchable, on its own. At build time it is opted in by `-dihedralbackbone`, not by
`-dihedralsidechain`: the chain runs *through* the ring (C4'-C3').

  These are not a refinement on top of the mesh, and they are not only about aromatic rings.
  Displace an sp2 hub by `z` out of its three substituents' plane: the three 1-3 distances do
  not change at all, and each bond changes by `z^2/2r`. Every pairwise distance among those
  four atoms is therefore an **even** function of `z`, stationary at `z = 0`, so a spring
  network's energy is `O(z^4)` -- it has no second-order stiffness in that mode at all, at any
  `--stiffness`. An improper's energy is `O(z^2)`.

  The mesh only holds a plane where the group is over-determined by a neighbour off it: a
  closed ring, or two overlapping `.rbody` groups. Where it is not -- a carboxylate, whose two
  oxygens appear in no other group -- the impropers are the only term doing the work.
  Ubiquitin, deviation from planarity in degrees, median (max), against OpenMM/amber99sb:

  | hub | AMBER | on, k=500 | on, k=100 | on, k=50 | off, k=500 | off, k=100 | off, k=50 |
  |---|---|---|---|---|---|---|---|
  | Asp CG | 0.02 | 0.04 | 0.03 | 0.02 | 6.20 | 17.35 | 24.31 |
  | Glu CD | 0.05 | 0.05 | 0.03 | 0.07 | 4.98 | 8.16 | 12.23 |
  | backbone C | 0.58 | 0.49 | 0.50 | 0.40 | 2.37 | 2.30 | 1.91 |
  | aromatic, Arg CZ, Asn/Gln | 0.04 | 0.02 | -- | -- | 0.01 | -- | -- |

  With the impropers on, planarity is **independent of `--stiffness`**; without them it is the
  mesh's job and it degrades as the mesh softens. That is what makes lowering `--stiffness`
  safe (see the next section). Widening the `.rbody` group instead was measured and rejected:
  it restores the plane but freezes the terminal torsion, because nothing beyond the hub's
  substituents can be reached without crossing the rotatable bond.

* **dihedralphi.scale = 1.0** *(dimensionless factor, float)* A multiplier on that family's
torsion energy and torque, and so for each of the seven others: `dihedralpsi.scale`,
`dihedralomega.scale`, `dihedralchi.scale`, `dihedralplanarity.scale`,
`dihedralnucleicbackbone.scale`, `dihedralnucleicchi.scale` and `dihedralnucleicsugar.scale`, all
`1.0` by default.

  It is the only way to weigh torsions against the rest of the model. `spring.scale` does reach
them, but the rigid-body mesh holding every bond and valence angle reads that same number, so
turning the torsions down with it turns the bonds down too.
* **dihedralphi.cutoff = 0** *(Angstroms, float)* Parsed, printed in the configuration dump, and
read by nothing -- as are `dihedralpsi.cutoff`, `dihedralomega.cutoff`, `dihedralchi.cutoff`,
`dihedralplanarity.cutoff`, `dihedralnucleicbackbone.cutoff`, `dihedralnucleicchi.cutoff` and
`dihedralnucleicsugar.cutoff`. Every dihedral family is an `EnergySetting`, which carries the
whole `enable`/`scale`/`cutoff` triplet, and a torsion is four named atoms rather than a search
radius.

**All eight default to enabled** (`defaultConfiguration()` sets them so deliberately): they are
an opt-OUT knob for isolating one family's contribution, not an opt-in feature switch, so
`-rigidbody ... -dihedral` needs nothing in the `.msp` to apply every family.

**They isolate a family's contribution; they do not undo the model.** Turning every one of them
off does *not* reproduce a topology built without the corresponding `pdb2spn` flags: a family
not requested at build time never gets a torsion or a NetCDF entry, so nothing at runtime can
bring it back.

To compare a rigid-body model against a bonded one, **build one `.nc` per stage**
(`--rigidbody` alone, then `+ --dihedral`) rather than toggling one `.nc` at runtime. See
`073.BondedStages` in the Biospring-Example repository, which does exactly that.

Bonds and valence angles have no `.msp` switch and no `pdb2spn` flag of their own: they are
held by the `--rigidbody` mesh at `--stiffness`. See `073.BondedStages`' README in the
Biospring-Example repository for
why the model is built that way, and for the `--stiffness` value it needs.

### Choosing `--stiffness` and `simulation.timestep`

The mesh is what limits the timestep -- verified, not assumed: with every torsion family
disabled the ceiling is unchanged (10 fs either way at k = 50), and `dt` tracks `1/sqrt(k)`
across the whole range. So the two settings are one choice, and the only question is how much
well fidelity a softer mesh costs.

**Use 650 kJ.mol-1.A-2, with `simulation.timestep = 3.0` on a protein and `2.0` on a nucleic
acid.** The rest of this section is how that was arrived at, and when to depart from it.

Measured by perturbing every chi1 (protein) or glycosidic chi (nucleic) by +40 deg off AMBER's
own bonded minimum, quenching 40 ps at `viscosity.value = 0.1`, and comparing every torsion
against OpenMM. `dt` is the largest value that survived 40 ps; cost is relative to the
recommended row at equal simulated time.

**What is being scored, and what is not.** A rotamer turns about a single bond and the mesh
must simply be rigid enough not to deform in its place. A furanose pucker is not a rotamer:
`nu1`, `nu2` and `delta` are ring coordinates, set by the ring's own bonds and angles, so the
mesh is their only support and they keep rewarding stiffness indefinitely -- see
`generate_ring_bonded.py`, which is the real answer to them. They are therefore reported apart:
the rotamers are the target, and the ring only has to stay a chemically sound furanose.

Ubiquitin, 1228 particles, against amber99sb. Its rings (Phe, Tyr, Trp, His, Pro) carry no
measured torsion -- chi2 turns about a single bond -- so there is nothing to separate:

| `--stiffness` (kJ.mol-1.A-2) | `timestep` (fs) | cost | median (deg) | p95 | within 10 deg | bond drift (A) |
|---|---|---|---|---|---|---|
| 250 | 4 | 0.75 | 1.15 | 9.98 | 95.0 % | 0.0022 |
| 500 | 3 | 1.00 | 0.85 | 7.27 | 96.2 % | 0.0013 |
| **650** | **3** | **1.00** | **0.75** | 6.99 | **96.6 %** | 0.0010 |
| 800 | 3 | 1.00 | 0.66 | 6.77 | 96.6 % | 0.0008 |
| 1200 | 2 | 1.50 | 0.82 | 7.13 | 96.6 % | 0.0006 |
| 1600 | 2 | 1.50 | 0.80 | 7.01 | 96.6 % | 0.0005 |

The protein has an interior optimum at 800: it is the last value that still runs at 3 fs, and
above it the timestep halves while the wells stop improving. 650 gives up 0.09 deg of median
for nothing else.

B-DNA duplex, 1270 particles, against amber14/DNA.OL15 -- rotamers (alpha, beta, gamma,
epsilon, zeta, chi; 252 angles) apart from the ring:

| `--stiffness` | `timestep` | cost | rotamers, median | within 10 deg | ring bonds (A) | ring angles (deg) |
|---|---|---|---|---|---|---|
| 250 | 3 | 0.67 | 2.00 | 94.0 % | 0.0266 | 7.90 |
| 500 | 2 | 1.00 | 1.40 | 94.4 % | 0.0163 | 7.51 |
| **650** | **2** | **1.00** | **1.30** | **93.7 %** | **0.0128** | 7.46 |
| 1000 | 1.5 | 1.33 | 1.17 | 94.4 % | 0.0087 | 7.33 |
| 2000 | 1 | 2.00 | 0.94 | 98.4 % | 0.0057 | 6.85 |
| 8000 | 0.5 | 4.00 | 0.67 | 100.0 % | 0.0028 | 5.80 |
| *AMBER itself* | | | | | *0.0013* | *5.12* |

RNA hairpin, 389 particles, against amber14/RNA.OL3:

| `--stiffness` | `timestep` | cost | rotamers, median | within 10 deg | ring bonds (A) | ring angles (deg) |
|---|---|---|---|---|---|---|
| 500 | 2 | 1.00 | 1.64 | 100.0 % | 0.0092 | 7.19 |
| **650** | **2** | **1.00** | **1.52** | **100.0 %** | **0.0079** | 6.48 |
| 1000 | 1.5 | 1.33 | 1.28 | 100.0 % | 0.0065 | 5.65 |
| 2000 | 1 | 2.00 | 1.06 | 100.0 % | 0.0045 | 3.88 |
| *AMBER itself* | | | | | *0.0004* | |

The nucleic rotamers are flat over the whole usable range -- 93.7 to 94.4 % on DNA from k = 100
to 500, and a clean 100 % on RNA from 500 to 2000 -- so there is no interior optimum there,
only a trade: everything that still improves with stiffness is the ring. RNA's ring is the
tighter of the two at equal k (0.0079 A against 0.0128) because the 2'-OH gives C2' one more
substituent, hence more 1-3 springs holding the ring closed.

So 650 is optimal nowhere and measurably worse nowhere, and it removes the need for three
values: it keeps each system's best timestep (3 fs and 2 fs), costs the protein 0.09 deg
against its own optimum, and holds the DNA ring 25 % tighter than 500 does.

**Depart from it** when the sugar pucker itself is the object of study: the ring wants 2000
(DNA 98.4 % of its ring coordinates, RNA 100 %) at half the speed. The right fix is the ring's
real bonded terms rather than a stiffer mesh -- with them, k = 250 beats k = 2000 at three
times the timestep.

Two things do *not* limit how far `--stiffness` can drop. The torsions do not: the timestep
ceiling is the same with them disabled. Planarity does not either, as long as the PLANARITY
impropers are on -- see `dihedralplanarity.enable` above, where the deviation is flat at
0.02-0.07 deg from k = 500 down to k = 50. What degrades is the mesh's grip on the torsion
wells themselves, and nothing else compensates for that.

`pdb2spn` applies this by itself: with `-rigidbody` and no explicit `-stiffness`, the default
is 650 rather than the generic 1.0, which suits a soft elastic network built from a distance
cutoff and does not hold a bond here, let alone a valence angle. An explicit `-stiffness` is
honoured whatever it says, so nothing that already passes one changes behaviour, and a
`-cutoff` network keeps its own 1.0.

Finally, `--rigidbody` takes every mesh spring's rest length from the **input structure's own
distances** (`RigidBodyBuilder` passes -1.0, meaning "use the current one"). Whatever geometry
the input has is thereby declared to be equilibrium. A torsion about a single bond does not
care, but a ring coordinate does: on an unminimised B-DNA the furanose pucker was frozen at
the input's value, putting nu1 31.6 deg off AMBER while raising `--stiffness` made it *worse*
(42.5 deg at k = 8000). Relax the structure under AMBER's bonded terms before building the
`.nc`, or the numbers above do not apply.
* **viscosity.enable = 0** *(boolean)* Enables a damping factor on the particles.
* **viscosity.value = 1.0** *(Da.fs-1, float)* Damping factor.

Force Field Parameters Description
----------------------------------

In addition to springs interaction, you can define steric (short-range)
and electrostatic (long-range) interactions. Steric interactions are simulated
by a linear repulsive force or a Lennard-Jones potential. Electrostatic interactions
are simulated by Coulomb equation.

A pre-computed electrostatic potential grid from APBS[1] can be used in some cases, as well as an
implicit membrane potential [2]. Both are shown by the examples that use them, in the
Biospring-Example repository: `011.DNAseI`, `012.MgSite_GK` and `021.Pax5` for the APBS grid,
`053.IMPALA_spheniscine`, `054.OmpA` and `057.TestIMP` for the implicit membrane.

* **steric.enable = 0** *(boolean)* Enable steric interaction.
* **steric.mode = linear** *(string)* Type of steric interaction Can be *linear, lennard-jones-8-6Lewitt, lennard-jones-8-6Zacharias, lennard-jones-12-6Amber*.
  `lennard-jones-8-6Amber` is the former spelling of `lennard-jones-12-6Amber`, still
  accepted: it is translated to the current name with a warning, so existing `.msp`
  files keep running the force field they ask for. Prefer the current name in new files.
* **steric.gridscale = 1** *(dimensionless factor, float)* Multiplier applied to steric forces
(steric stiffness in the *linear* mode is a fixed kJ.mol-1.A-2 constant in the code; the
Lennard-Jones modes use each particle's `epsilon`, in kJ.mol-1, from the .ff file).
* **steric.cutoff = 1** *(Angstroms, float)* Cutoff distance for steric calculation.
* **steric.radiusrule = default** *(string)* How a pair's two radii combine, for the
Lennard-Jones modes: *default, sum, geometric-mean, arithmetic-mean*. `default` means the law's
own rule, which is the rule of the parameter set that law was published with
(`STERIC_RADIUS_RULE_*` in `forcefield/energy/steric.hpp`). Set it explicitly only when pairing a
`.ff` with a law it did not come from, because **the `radius` column does not hold the same
quantity in every file**:

  | shipped `.ff` | what column 3 holds | rule its radii want |
  |---|---|---|
  | `amber*.ff` | AMBER's R\* = rmin/2, 0.6-2.0 A | sum |
  | `model.ff`, `CAonlyBaaden.ff`, `CAallresBaaden.ff` | bead radius, 2.6-3.9 A | sum |
  | `CAonlyLewitt.ff` | CA-CA **contact** distance, 5.0-7.4 A | geometric mean |
  | `CGZaccharias.ff` | its own parameterisation | product, inside the law |

  Each shipped file's header states its own. Getting it wrong is not subtle: summing Levitt's
contact distances puts a residue pair's minimum at 10-13 A where a real protein's closest
non-bonded CA pair is 3.71 A, and taking the geometric mean of AMBER's half-distances put two
carbons at 1.91 A instead of 3.82, which left nucleobases stacked 3.4 A apart with no minimum at
all -- the RNA hairpin's steric energy was -41.17 kJ.mol-1 where it should be -571.37. Zacharias'
law ignores this setting: its product on the radius is inseparable from its product on epsilon.
A multi-resolution network (a ribosome) needs both conventions at once, which is what made this
selectable rather than fixed.

  **Note that `steric.mode` defaults to `linear`, which is no Lennard-Jones at all**, so an
`.msp` that enables sterics without naming a mode never reads a `radius` column written as
rmin/2. The linear wall is soft: on Martini 3's large beads (rmin up to 5.28 A), with Coulomb on
at an 11 A cutoff and 300 K asked of the thermostat, glucokinase ran at 1377 K and two neutral,
unbonded beads were pushed from 7.69 A to 0.91 A -- where the same run under
`lennard-jones-12-6Amber` held 301.3 K.
---
* **coulomb.enable = 0** *(boolean)* Enables Coulomb interaction.
* **coulomb.scale = 1.0** *(dimensionless factor, float)* Multiplier applied to electrostatic
forces (charges, in elementary charge units *e*, come from the .ff file).
* **coulomb.cutoff = 16.0** *(Angstroms, float)* Cutoff distance for Coulomb.
calculation
* **coulomb.dielectric = 1.0** *(dimensionless, float)* Relative dielectric constant used in
Coulomb's equation.
* **coulomb.distancedependent = 0** *(boolean)* Makes the dielectric grow with distance,
`epsilon(r) = coulomb.dielectric * r` with r in Angstrom (Warshel & Levitt 1976), instead of
holding it constant.

  A constant 78 screens a pair at contact exactly as hard as one 30 A away, and that is wrong in
the direction that matters: no water fits between two atoms 3 A apart, so the short-range
interaction is barely screened in reality. Measured on a native alpha-helix hydrogen bond of
ubiquitin, the amide/carbonyl attraction is **-13.61 kJ.mol-1** unscreened and **-0.17 kJ.mol-1**
at a constant 78.

  Off by default: it changes every electrostatic energy in the model, so it is a modelling choice
and not a correction.
---
* **scale14.enable = 0** *(boolean)* Applies AMBER's 1-4 scale factors -- `scnb` on the
Lennard-Jones, `scee` on the Coulomb -- to the pairs three covalent bonds apart.

  BioSpring excludes only SPRUNG pairs from the non-bonded terms and scales nothing, so without
this its 1-4 pairs are live at full strength. AMBER's are not, and its bonded parameters were
fitted with those factors in place. Measured on 075.RNAHairpin, the 983 1-4 pairs carry +891
kJ.mol-1 of Lennard-Jones unscaled against AMBER's own +449 on the same pairs -- the error is
exactly the missing factor of two -- and on 073's ubiquitin +3069 against +1536. With the factor
on, the model's 1-4 Lennard-Jones lands within 0.7 % of AMBER's on the RNA and 0.1 % on the
protein, in absolute kJ.mol-1.

  Off by default, because enabling it changes the energy of every network already built. No `.nc`
has to be rebuilt: the 1-4 pairs are derived at setup from the spring mesh itself, by perceiving
which of its springs are covalent bonds. That perception needs the mesh to BE a bond-and-angle
graph -- which a `--rigidbody` mesh is and a `--cutoff` mesh is not -- and it says so: a perceived
bond or angle pair that carries no spring is counted and reported as a warning. Measured against
the force fields the generators read, the derivation is exact: 419 of 419 bonds, 1165 of 1165
1-2/1-3 pairs and 983 of 983 1-4 pairs on 075 (amber14/RNA.OL3), 1236 of 1236 bonds and 3261 of
3261 1-4 pairs on 073 (amber99sb), with nothing invented and nothing missed.

  It costs a second scan per candidate pair, over its own CSR, next to the spring exclusion's.
Measured over 3000 steps:

  | | 072.GK | 073.ubiquitin | 034.VirusCA |
  |---|---:|---:|---:|
  | CPU, steps/s off then on | 223.2 / 205.2 | 477.9 / 442.2 | 453.8 / 451.1 |
  | GPU, steps/s off then on | 404.5 / 329.3 | 488.0 / 395.7 | 1748 / 1762 |

  8 % on a core and 19 % on the device for an all-atom network, and nothing at all on a
coarse-grain one, which has no covalent bond to perceive: the list is empty, the kernels get a
null pointer and the test is one comparison.

* **scale14.steric = 0.5** *(dimensionless factor, float)* AMBER's `1/scnb`. Used only when
`scale14.enable` is set.

* **scale14.coulomb = 0.8333** *(dimensionless factor, float)* AMBER's `1/scee`, i.e. 1/1.2. Used
only when `scale14.enable` is set.
---
* **electrostaticgrid.enable = 0** *(boolean)* Enable the APBS electrostatic potential grid.
Renamed from `potentialgrid`, which said only that it was a grid: `densitygrid` below is one
too. An .msp still using the old name keeps working and prints one deprecation warning per
obsolete group.
* **electrostaticgrid.path = ""** *(string)* Name of the APBS potential grid file in OpenDX format.
* **electrostaticgrid.scale = 1** *(dimensionless factor, float)* Multiplier applied to electrostatic
forces derived from the grid.

`electrostaticgrid` and `coulomb` are independent and may be combined: a receptor built with
`pdb2spn --charge 0.0` contributes only through the map, while explicitly charged particles
(ions, say) feel the map *and* each other pairwise, with no double counting. Pairwise cost
scales with the number of *charged* particles, not the total. Before this was fixed, enabling
the grid alone silently produced no electrostatics at all: its gate was Coulomb's own flag,
so the DX file was never even read.
---
* **densitygrid.enable = 0** *(boolean)* Enable density grid.
* **densitygrid.path = ""** *(string)* Name of the density grid file in OpenDX format.
* **densitygrid.scale = 1** *(dimensionless factor, float)* Multiplier applied to forces derived
from the density grid (e.g. a SAXS/cryoEM-derived envelope). Independent from steric.gridscale
and electrostaticgrid.scale.

Dynamic Hydrogen Bonds
----------------------

A hydrogen bond here is a Morse well in the donor-acceptor **heavy-atom** distance, weighted at
both ends by how well the bond lines up with the direction each side offers -- the donor's
hydrogen standing in as "away from the antecedent", the acceptor's lone pairs standing in the
same way (`forcefield/shared/hbond_shared.h`). There is no explicit hydrogen in the model, and
the pairs are **not** fixed: they are re-assigned from the neighbour grid as the run goes, so a
bond can break and another form. For springs that cannot break, see `pdb2spn --static-hbond`
instead.

The well's own three parameters are compiled in rather than `.msp` settings: De = **16.7
kJ.mol-1**, re = **2.90 A**, a = **1.34 A-1** (`forcefield/ForceField.h`). `hbond.scale` is the
only handle on the depth.

* **hbond.enable = 0** *(boolean)* Enable the dynamic hydrogen-bond term.
* **hbond.scale = 1.0** *(dimensionless factor, float)* Multiplier on the well depth and on the
force. `1.0` means the 16.7 kJ.mol-1 above.
* **hbond.cutoff = 7.0** *(Angstroms, float)* Sizes the cells of this term's own neighbour grid.
It is the practical range of the Morse well -- under 1 % of the depth is left beyond ~7 A at the
shipped De/re/a -- and not a long-range, electrostatics-style cutoff.
* **hbond.path = ""** *(string)* The donor/acceptor table, a `.hbond` file: see
`data/reducerules/ProteinDonorAcceptor.hbond`, `DNADonorAcceptor.hbond` and
`RNADonorAcceptor.hbond`. It classifies each atom by its **group name**, the first column of the
`.grp` the network was reduced with, so there is nothing to translate in the `.msp`. Empty leaves
every atom neither donor nor acceptor, and the term with nothing to do.
* **hbond.log = ""** *(string)* Where to write the list of bonds actually held, one line each, at
every sample step. Empty writes nothing. A total energy and a bond count cannot say WHICH bonds
are held, and that is the only thing that settles a disagreement about pairing: a fully extended
chain reports a perfectly healthy count made entirely of contacts between neighbouring residues.

### Isostericity (nucleic acids)

A **reduction**, in [0, 1], on a pair of nucleotides whose base pair does not sit where a
Watson-Crick pair sits *with respect to its own backbone* -- the C1'-C1' distance and the two
glycosidic angles. It exists because the pairing term provably cannot see the difference: the
distance and both angles are preserved exactly when the partner is turned 180 deg about the
donor-acceptor axis, measured at **-0.00 +/- 0.01 kJ.mol-1** over 49 bonds, so without this a
mirrored pair scores like the real one.

* **hbond.isostericity = 0** *(boolean)* Enable the reduction. Off by default, because it is a
factor on an already-calibrated depth: every existing `.msp` keeps the trajectory it had.
* **hbond.isostericitytwosided = 0** *(boolean)* Whether the **distance** acceptance is a window
(`1`, both sides) or a floor (`0`, the default, long distances unpenalised). A window is wrong for
anything dynamic and measurably so: two bases on their way to pairing are FURTHER apart than the
centre, so the Gaussian switches the well off during the approach -- the factor is 0.157 at a
C1'-C1' of 12.5 A, 7e-4 at 14.5 A, and zero at the 18.6 A a free tRNA starts from. The parasites
this term exists to reject are SHORTER than a real pair, not longer:

  | form | native/parasite | approach kept |
  |---|---:|---:|
  | window (`= 1`) | 7.3 x | 5.3 % |
  | **floor (`= 0`, default)** | **7.0 x** | **100.0 %** |

  The window is kept only to reproduce the measurements that were made with it.
* **hbond.isostericitydistance = 10.46** *(Angstroms, float)* The Watson-Crick C1'-C1' centre,
measured on 8PKL over 1121 strict pairs (10.46 +/- 0.59 A).
* **hbond.isostericitylambda = 54.9** *(degrees, float)* The glycosidic-angle centre, same
measurement (54.9 +/- 5.6 deg).
* **hbond.isostericitydistancewidth = 1.5** *(Angstroms, float)* Gaussian width on the distance.
* **hbond.isostericitylambdawidth = 18.0** *(degrees, float)* Gaussian width on the angle.

  The two widths are **not** the standard deviations above. They were scanned against the native,
mirror and non-canonical populations at once, and came out wider, so that a real pair keeps 84 %
of its well while its mirror keeps 4.7 %.

Peptide Bond Formation
----------------------

The **second** stage of making a covalent bond during a run: once two residues have been held in
the attack conformation long enough, the spring topology is REWRITTEN so the new bond exists and
cannot come undone. The first stage -- bringing the nucleophile onto the electrophile -- is the
hydrogen-bond term above with an out-of-plane lobe site, and nothing here duplicates it. What
this adds is the irreversible part, which no potential well can express: a bond that has formed
is a different molecule, not a deeper minimum of the same one.

* **peptidebond.enable = 0** *(boolean)* Enable the topology rewrite. Off, and with no rules
path: forming a covalent bond is never something an `.msp` should start doing because it was
upgraded.
* **peptidebond.path = ""** *(string)* The `.rbody` the network was built from. The new springs
are created by the SAME rule that would have created them had the bond been there from the start,
so the product is the model's own peptide and not a special case. Required: with `enable` on and
this empty, setup dies.
* **peptidebond.group = PSI** *(string)* Which rule describes the new bond -- the suffix after the
residue code, so `PSI` picks up `A_PSI`, `R_PSI`, ...: {CA, C, O, +N, +H, +CA}, the peptide plane
ahead of a residue.
* **peptidebond.bond = C:N** *(string)* Which PAIR inside that group is the bond being made, as
`<electrophile>:<nucleophile>`. The rule itself cannot say: it is six atoms sprung pairwise, with
no distinguished pair among the fifteen.
* **peptidebond.stiffness = 650** *(kJ.mol-1.A-2, float)* What the new springs get.
* **peptidebond.distance = 3.2** *(Angstroms, float)* Below this the nucleophile counts as
attacking.
* **peptidebond.weight = 0.5** *(dimensionless, 0..1, float)* The angular factor both ends must
reach together.
* **peptidebond.dwell = 200** *(steps, integer)* How many CONSECUTIVE steps both conditions must
hold. The gate is a duration, not an instant: one step too far apart, or one step past the angle,
erases the count. The frame block on standard output reports the two extrema it watches --
closest approach against `distance`, best weight against `weight` -- and they are not necessarily
the same pair.
* **peptidebond.ramp = 2000** *(steps, integer)* Over how many steps the new springs are drawn
from the length they were BORN at -- the attack distance, so no energy appears out of nothing --
to the length an ideal trans peptide plane gives them. `0` leaves them where they were born, which
is a 2.9 A "peptide bond".
* **peptidebond.leaving = ""** *(string)* The names an atom may have to be the LEAVING GROUP,
space-separated (e.g. `OXT O3'`), empty by default. It is looked for among the atoms BONDED to the
electrophile, wherever they live, and what breaks is every spring between it and the
electrophile's residue.

  That one rule covers both chemistries. A free acid's OXT is inside that residue, so all its
springs go and it departs whole -- the water that leaves. An aminoacyl-tRNA's O3' belongs to the
tRNA, so only the ester C-O goes and the oxygen stays on its ribose: the peptide transfers to the
other tRNA and this one is left deacylated, which is peptidyl transfer. Breaking is both ways --
a spring is also what keeps a pair out of the steric and Coulomb terms, so what pushes the group
away afterwards is ordinary sterics.
* **peptidebond.sites = ""** *(string)* The `.psite` file saying WHERE a carbonyl carbon may be
attacked from: its antecedent, and the plane that splits its two faces. Its own file and its own
reader on purpose -- a reactive site is not a hydrogen bond, and filing it in the `.hbond` table
meant either file could silently overwrite the other's direction on an atom both named. Empty
falls back to whatever the `.hbond` table tagged, which is what every configuration written
before this setting existed relies on.
* **peptidebond.naming = ""** *(string)* The `.grp` the network was REDUCED with, when it was.
`amber.grp` renames every atom per residue -- C becomes AC for alanine, RC for arginine -- so on a
reduced network the names in `bond` and in the `.rbody` rules find nothing without it. Same table
and same purpose as `pdb2spn --rigidbody`'s own naming translation.
* **peptidebond.log = ""** *(string)* Parsed, printed in the configuration dump, and **read by
nothing** -- like `spring.cutoff` above. The per-step view of the gate goes to the standard-output
frame block instead.

Implicit Membrane (IMPALA)
-----------------

Particles with their transfer energies preconfigured via pdb2spn (using the forcefield and reducerules options) can interact with an implicit membrane according to the IMPALA model.

* **impala.enable = 0** *(boolean)* Enable IMPALA membrane interaction.
* **impala.scale = 1.0** *(dimensionless factor, float)* Multiplier applied to IMPALA forces
(particle transfer energies, in kJ.mol-1, come from the .ff file, see pdb2spn's forcefield/reducerule
options).
---
* **insertionvector.enable = 0** *(boolean)* Enable Insertion Vector.
* **insertionvector.vector = 0 0** *(int int)* Set the IDs of the two particles defining the insertion vector.

RigidBody
---------

Switch from flexible spring network dynamics mode to rigid body dynamics mode.

* **rigidbody.enable = 0** *(boolean)* Enable Rigid Body mode.
* **rigidbody.enablesampling = 0** *(boolean)* Enable Automatic Sampling of insertion into the implicit membrane.
* **rigidbody.enablemontecarlo = 0** *(boolean)* Enable Monte Carlo Rigid Body to run exploration of random steps in conformational space.
* **rigidbody.montecarlo_translation_norm = 0.1** *(float)* Magnitude of translation in angstroms (Å) for the Monte Carlo rigid body.
* **rigidbody.montecarlo_rotation_norm = 0.1** *(float)* Angle of rotation in degrees (°) for the Monte Carlo rigid body.
* **rigidbody.montecarlo_temperature = 298.1** *(float)* Temperature in Kelvin (K) for the Monte Carlo simulations rigid body.

Hydrophobicity (experimental)
--------------

Add a pseudo-hydrophobicity interaction for multimeric assembly into a rigid body in an implicit membrane.

This is a **two-body** term between particles, `E = -h_i h_j L exp(-r/L)`, and it
is not the IMPALA implicit-membrane term. The two read different columns of the
`.ff` file and must not be confused:

| `.ff` column | `.nc` variable | term | unit |
|---|---|---|---|
| 6, `transferIMP` | `hydrophobicityscale` | IMPALA, one-body against a membrane | kJ.mol-1.A-2 |
| 7, `Hydrophobicity` | `hydrophobicity` | this one, two-body between particles | see `scale` below |

* **hydrophobicity.enable = 0** *(boolean)* Enable the pairwise hydrophobic interaction.
* **hydrophobicity.scale = 1.0** *(dimensionless factor, float)* Multiplier applied to
hydrophobicity forces. The per-particle `h` comes from the seventh column of the `.ff`
file; the product `h_i h_j` is in kJ.mol-1.A-1. The product form means `h` cannot be a
signed hydrophobicity scale -- two hydrophilic particles would both be negative, their
product positive, and they would attract each other as if they were oil -- so a rectified
scale is required and a particle that is not hydrophobic gets exactly 0.
* **hydrophobicity.cutoff = 15.0** *(Angstroms, float)* Cutoff distance for the pair search.
* **hydrophobicity.decaylength = 10.0** *(Angstroms, float)* The distance over which the
attraction decays. A property of the solvent that you pick, in the same way
`coulomb.dielectric` is, rather than a universal constant. The default is the decay length
Israelachvili & Pashley measured between two **macroscopic** hydrophobic surfaces
(Nature 300:341, 1982, published as `22 exp(-D/10)` mJ.m-2 with D in Angstrom), which is
where this law comes from; later work puts the range anywhere between 3 and 10 A depending
on the system, and two coarse-grained beads are not two macroscopic plates. See example
058.Hydrophobicity.


Probe
-----

The probe is a charged entity used to explore and characterize binding sites through electrostatic interactions.

* **probe.enable = 0** *(boolean)* Enable probe.
* **probe.enableelectrostatic = 0** *(boolean)* Enable probe electrostatic interaction.
* **probe.enablesteric = 0** *(boolean)* Enable probe steric interaction.
* **probe.x = 1.0** *(Å, float)* Initial x position of the probe.
* **probe.y = 1.0** *(Å, float)* Initial y position of the probe.
* **probe.z = 1.0** *(Å, float)* Initial z position of the probe.
* **probe.mass = 1.0** *(Da, float)* Mass of the probe.
* **probe.epsilon = 1.0** *(kJ.mol-1, float)* Set the probe's interaction energy.
* **probe.radius = 1.0** *(Å, float)* Sets the probe's radius
* **probe.charge = 0.0** *(e, float)* Sets the probe's charge


Automatic Constraints (NOT IMPLEMENTED -- do not put these in a file)
--------------------------------------------------------------------

There is no `constraint` setting group in the code. `Configuration` registers none, so
`SafeConfigurationReader` kills the run on the first `constraint.*` line it reads
(`invalid parameter "constraint.enable"`). The four parameters this section used to list --
`enable`, `src`, `dest` and `scale`, the last in Da.A.fs-2 -- were a design note, never a feature,
and an `.msp` carrying them does not start.

To pull one group of particles toward another, steer it interactively over MDDriver/IMD (see
`doc/User_Manual.md`), optionally freezing the partner at build time with `pdb2spn --static`.


    
## References
[1]: Jurrus E, Engel D, Star K, et al. Improvements to the APBS biomolecular solvation software suite. Protein Sci. 2018;27(1):112-128. doi:10.1002/pro.3280  
[2]: Ducarme P, Rahman M, Brasseur R. IMPALA: a simple restraint field to simulate the biological membrane in molecular structure studies. Proteins. 1998;30(4):357-371. 