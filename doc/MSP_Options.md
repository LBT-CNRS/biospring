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
* **simulation.neighborskin = 0** *(distance unit, float)* Extra margin added to the steric,
electrostatic and hydrophobic cutoffs when building their neighbor grids. When greater than
zero, the grid is only rebuilt once a particle has moved more than this margin since the last
rebuild, instead of every step, which reduces the cost of neighbor search. `0` (the default)
rebuilds the grid every step, which is always correct but can be slower for large systems.
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
* **dihedralphi.enable = 1** *(boolean)* Runtime debug on/off for phi (backbone) dihedral ghost
springs, independently of psi/omega/chi. Only meaningful if the topology was actually built with
`-dihedralbackbone` (or `-dihedral`) in the first place -- a family not built has no springs to
enable/disable either way. Defaults to enabled so an `.msp` written before this setting existed
keeps the same behaviour.
* **dihedralpsi.enable = 1** *(boolean)* Same as `dihedralphi.enable`, for the psi axis.
* **dihedralomega.enable = 1** *(boolean)* Same as `dihedralphi.enable`, for the omega
(peptide-bond) axis.
* **dihedralchi.enable = 1** *(boolean)* Same as `dihedralphi.enable`, for every side-chain
chi1-4 dihedral (the SIDECHAIN family in the `.bi.ff`).
* **dihedralplanarity.enable = 1** *(boolean)* Same as `dihedralphi.enable`, for the PLANARITY
impropers.

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

**These settings isolate a family's contribution; they do not undo the model.** Turning every
one of them off does *not* reproduce a topology built without the corresponding `pdb2spn`
flags: a family not requested at build time never gets a spring, a ghost particle or a NetCDF
entry, so nothing at runtime can bring it back -- and what you get instead is the rigid body
still carrying every ghost particle the disabled families created.

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

Note that `pdb2spn`'s own `-stiffness` default is still 1.0, which is not a usable value for a
rigid-body mesh: pass it explicitly.

Finally, `--rigidbody` takes every mesh spring's rest length from the **input structure's own
distances** (`RigidBodyBuilder` passes -1.0, meaning "use the current one"). Whatever geometry
the input has is thereby declared to be equilibrium. A torsion about a single bond does not
care, but a ring coordinate does: on an unminimised B-DNA the furanose pucker was frozen at
the input's value, putting nu1 31.6 deg off AMBER while raising `--stiffness` made it *worse*
(42.5 deg at k = 8000). Relax the structure under AMBER's bonded terms before building the
`.nc`, or the numbers above do not apply.
* **dihedral.tangentialonly = 0** *(boolean)* Project each ghost ring's reaction onto the
tangential direction about its own axis before it reaches the real atoms, so a torsion pushes a
substituent only *around* that axis -- which is exactly what AMBER's dihedral force does
(`F` is along `r_ij x r_jk`, hence perpendicular to both the axis and the i-j-k plane). What the
projection drops carries no torque about the axis at all, so the torque is preserved exactly;
the reaction on the two axis atoms is then balanced against zero rather than against the ghost
totals, which is what keeps the discarded part from simply reappearing there.

Measured on ubiquitin, this is not a small correction: 96.4 % of the force the rings apply
(99.8 % median) is radial or axial and merely deforms the rigid-body mesh. With the projection,
spring energy after 20000 steps falls from 442.90 to 0.80 kJ/mol, and per-atom agreement with
AMBER's own dihedral forces goes from a correlation of 0.000 to 0.650. Because the leak does not
scale with `--stiffness` while the mesh's resistance does, it is what forces a stiff mesh: with
the projection, stiffness can drop from 8000 to 500 and the timestep rise from 1.0 to 4.0 fs,
with *better* geometry at equal simulated time (CA-RMSD 1.60 -> 1.48 A over 20 ps).

Off by default: it changes the forces of an existing model, so it is opted into rather than
imposed. **Known limit**: the torque is right in total and in direction, but it is delivered
concentrated -- on 936 of 937 axes the ring pushes fewer substituents than AMBER (1.35 against
4.92 on average), so the force on the atoms it does push is about 2.4x too large.

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
* **potentialgrid.enable = 0** *(boolean)* Enable APBS potential grid.
* **potentialgrid.path = ""** *(string)* Name of the APBS potential grid file in OpenDX format.
* **potentialgrid.scale = 1** *(dimensionless factor, float)* Multiplier applied to electrostatic
forces derived from the potential grid.
---
* **densitygrid.enable = 0** *(boolean)* Enable density grid.
* **densitygrid.path = ""** *(string)* Name of the density grid file in OpenDX format.
* **densitygrid.scale = 1** *(dimensionless factor, float)* Multiplier applied to forces derived
from the density grid (e.g. a SAXS/cryoEM-derived envelope). Independent from steric.gridscale
and potentialgrid.scale.

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

* **hydrophobicity.enable = 0** *(boolean)* Enable Hydrophobicity interaction.
* **hydrophobicity.scale = 1.0** *(dimensionless factor, float)* Multiplier applied to
Hydrophobicity forces (per-particle hydrophobicity/transfer scale, in kJ.mol-1, comes from the
.ff file).
* **hydrophobicity.cutoff = 15.0** *(Angstroms, float)* Cutoff distance for Hydrophobicity.


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