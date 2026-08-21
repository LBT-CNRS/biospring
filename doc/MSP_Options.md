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
impropers that keep aromatic rings and the His hub flat.

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