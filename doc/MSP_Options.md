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

* **simulation.timestep = 1** *(fs, float)* Timestep in fs used in this run
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
* **xtctrajectory.enable = 0** *(boolean)* Enables energies logging in xtc format.
* **xtctrajectory.frequency = 100** *(integer)* Frequence of energy logging.
* **xtctrajectory.path = ""** *(string)* Name of the xtc energies log.

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
* **dihedralnucleicbackbone.enable = 1** *(boolean)* Same, for a nucleotide's alpha..zeta.
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

A pre-computed electrostatic potential grid from APBS[1] can be use 
in some case, as well as an implicit membrane potential [2]. Details 
for these features are given in the appropriate examples. **(TODO: define 
which examples)**

* **steric.enable = 0** *(boolean)* Enable steric interaction.
* **steric.mode = linear** *(string)* Type of steric interaction Can be *linear, lennard-jones-8-6Lewitt, lennard-jones-8-6Zacharias, lennard-jones-12-6Amber*.
  `lennard-jones-8-6Amber` is the former spelling of `lennard-jones-12-6Amber`, still
  accepted: it is translated to the current name with a warning, so existing `.msp`
  files keep running the force field they ask for. Prefer the current name in new files.
* **steric.gridscale = 1** *(dimensionless factor, float)* Multiplier applied to steric forces
(steric stiffness in the *linear* mode is a fixed kJ.mol-1.A-2 constant in the code; the
Lennard-Jones modes use each particle's `epsilon`, in kJ.mol-1, from the .ff file).
* **steric.cutoff = 1** *(Angstroms, float)* Cutoff distance for steric calculation.
---
* **coulomb.enable = 0** *(boolean)* Enables Coulomb interaction.
* **coulomb.scale = 1.0** *(dimensionless factor, float)* Multiplier applied to electrostatic
forces (charges, in elementary charge units *e*, come from the .ff file).
* **coulomb.cutoff = 16.0** *(Angstroms, float)* Cutoff distance for Coulomb.
calculation
* **coulomb.dielectric = 1.0** *(dimensionless, float)* Relative dielectric constant used in
Coulomb's equation.
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


Automatic Constraints Parameters (not available yet)
--------------------------------

You can define automatic constraints to push a group of particles toward another
without manual interaction.

* **constraint.enable = 0** *(boolean)* Enables constraint for this run
* **constraint.src = ""** *(string)* The name of the first selection of atom
* **constraint.dest = ""** *(string)* The name of the second selection of atom
* **constraint.scale = 1.0** *(Da.A.fs-2, float)* Force module used for the constraint


    
## References
[1]: Jurrus E, Engel D, Star K, et al. Improvements to the APBS biomolecular solvation software suite. Protein Sci. 2018;27(1):112-128. doi:10.1002/pro.3280  
[2]: Ducarme P, Rahman M, Brasseur R. IMPALA: a simple restraint field to simulate the biological membrane in molecular structure studies. Proteins. 1998;30(4):357-371. 