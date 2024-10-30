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
* **spring.scale = 1.0** *(float)* Factor applied on the spring forces.
* **spring.cutoff = 1.0** *(float)* Factor applied on the spring forces.
---
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
* **steric.mode = linear** *(string)* Type of steric interaction Can be *linear, lennard-jones-8-6Lewitt, lennard-jones-8-6Zacharias, lennard-jones-8-6Amber*
* **steric.gridscale = 1** *(float)* Factor applied on steric forces.
* **steric.cutoff = 1** *(Angstroms, float)* Cutoff distance for steric calculation.
---
* **coulomb.enable = 0** *(boolean)* Enables Coulomb interaction.
* **coulomb.scale = 1.0** *(float)* Factor applied on electrostatic forces.
* **coulomb.cutoff = 16.0** *(Angstroms, float)* Cutoff distance for Coulomb.
calculation
* **coulomb.dielectric = 1.0** *(float)* Dielectric constant used in Coulomb equation.
---
* **potentialgrid.enable = 0** *(boolean)* Enable APBS potential grid.
* **potentialgrid.path = ""** *(string)* Name of the APBS potential grid file in OpenDX format.
* **potentialgrid.scale = 1** *(float)* Factor applied on electrostatic forces.
---
* **densitygrid.enable = 0** *(boolean)* Enable density grid.
* **densitygrid.path = ""** *(string)* Name of the density grid file in OpenDX format.

Implicit Membrane (IMPALA)
-----------------

Particles with their transfer energies preconfigured via pdb2spn (using the forcefield and reducerules options) can interact with an implicit membrane according to the IMPALA model.

* **impala.enable = 0** *(boolean)* Enable IMPALA membrane interaction.
* **impala.scale = 1.0** *(float)* Factor applied on IMPALA forces.
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
* **hydrophobicity.scale = 1.0** *(float)* Factor applied on Hydrophobicity forces.
* **hydrophobicity.cutoff = 15.0** *(Angstroms, float)* Cutoff distance for Hydrophobicity.


Probe
-----

The probe is a charged entity used to explore and characterize binding sites through electrostatic interactions.

* **probe.enable = 0** *(boolean)* Enable probe.
* **probe.enableelectrostatic = 0** *(boolean)* Enable probe electrostatic interaction.
* **probe.enablesteric = 0** *(boolean)* Enable probe steric interaction.
* **probe.x = 1.0** *(float)* Initial x position of the probe.
* **probe.y = 1.0** *(float)* Initial y position of the probe.
* **probe.z = 1.0** *(float)* Initial z position of the probe.
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