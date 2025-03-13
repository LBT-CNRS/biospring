# Introduction

BioSpring is a molecular simulation software based on a spring network, expecially designed for interactive simulations. The main goal of BioSpring is to allow user to interactively explore the biomechanical properties biomolecule using spring network model at different scales (all atom, coarsed grain, CA, ...).

As the input of BioSpring is a .nc file (NetCDF format), you need to use specific tools to convert .pdb file or .pqr file  to .nc file. These tools allow also user to provide rules to convert all atom representation (AA) to coarse-grain (CG) representation or CA represention (CA), before creating a spring network according a distance cutoff.


## BioSpring Theory


### Biomolecular structure and simulation

Experimentations (RMN or XRay) allow researchers to get the particle coordinates of a biomolecular structure, which are stored in a text file. A commonly used file format is the .pdb file, provided by the Protein Data Bank [http://www.rcsb.org/pdb/home/home.do](http://www.rcsb.org/pdb/home/home.do).

An example of this file can be found in the example directory. The particle coordinates of a biomolecular system are used as the initial state of molecular simulation used to study the dynamic behavior of biomolecules.


### Spring network model

We use the particle coordinates of the biomolecule to build a spring network between these particles:
When two particles `i` and `j` are closer than an arbitrary distance (usually between 7 and 15 Angstroms), a spring is created between these two particles, and the initial distance between these two particles is the equilibrium distance `Eij` for this spring.

Moreover, to control the rigidity and flexibility of the network, a stiffness parameter `k` for all spring is used (usually 0,6 kcal.mol<sup>-1</sup>.A<sup>-2</sup>).

Many papers show that spring network models reproduce the biomechanical behavior of biomolecules.

### Classical newtownian simulation

During a classical newtownian simulation to observe the behavior of this spring network after perturbation, we compute all the forces applied to the particles which are linked by a spring.

If the distance of particle `i` and `j` `Dij` is greater than initial equilibrium distance `Eij`, an attractive force is applied to the particle `i` and `j` and if the distance of particle `i` and `j` `Dij` is lower than initial equilibrium distance `Eij`, `i` and `j` are closer than initial equilibrium distance `Eij`, an repulsive force is applied to the particle `i` and `j`.

	Fij = k (Dij - Eij)

At each timestep t of the simulation, we first compute for each particle the sum of forces `Fi` applied on particle `i`, according to the previous formula.

We then compute for each particle its acceleration according to forces `Fi` and its mass `Mi`. The acceleration allow us to compute the new velocity `Vi` and the new position `Pi` for each particle `i`.

	Ai = Fi/Mi
	Vit = Vit-1 + Ai*dt
	Pit = Pit-1 + Vi*dt

### Steric and electrostatic interactions

As a particle represents an atom, we can consider it as a sphere. To avoid interpenetration between sphere, we use linear or Lennard-Jones (or van der Waals) potential, which necessites to compute the distances between all the particles.

Moreover particles are charged and are subject to electrostatic interactions, simulated using coulomb interaction (attractive force between opposite charge, repulsive in the other case), which requires also to compute the distances between all the particles. To avoid quadratic complexity of distance computing, we use a regular 3D grid to keep in real time the particle neighborhood at each step of the simulation.



## BioSpring input files description

Biospring is a simulation program which takes two files as input. The first one describes the structure and the biophysical properties of the system you want to study using the NetCDF binary format to describe this data (`.nc`), and the second one a text file which describes the parameters and settings for the spring network simulation (`.msp`).


### NetCDL file format, .cdl and .nc

We chose the NetCDF format (see [http://www.unidata.ucar.edu/software/netcdf/][netcdf]) which is a file format to describe the structure and the biophysical properties of the system because it is well adapted to describe array such as positions, charges, radius, ... , of an atoms set, and allows to quickly load a huge quantity of data.

There is two NetCDF formats, the human-readable one with `.cdl` extension, and the binaray one with `.nc`, used by BioSpring.

To convert human-readable file `model.cdl` to the binary file `model.nc`, the ncgen tool can be use as in this command line :

	ncgen -b -o model.nc model.cdl


### The BioSpring simulation parameter file .msp

Example of `.msp` file for BioSpring:

	#Molecular Simulation Parameter file

	#in fs
	simulation.nbsteps = -1
	simulation.timestep = 0.1
	simulation.samplerate = 10000

	spring.enable = 1
	spring.scale = 25

	#in Da.fs-1
	viscosity.enable = 1
	viscosity.value = 1



## Launching BioSpring


### without interaction

After generating model.nc file and setting simulation parameters param.msp and you can launch BioSpring using this command:

    biospring -s model.nc -c param.msp

Be carefull that you have set constraints and output parameters to get results. Indeed, without constraints, there is no reason that the spring
network comes out the initial equilibriumn, and without output parameters, you won't be able to get any results.

Or run the Docker image:

	docker run --init -v ./:/data ghcr.io/lbt-cnrs/biospring \
		biospring -s /data/model.nc -c /data/param.msp

### with interaction via VMD using MDDriver

BioSpring was especially designed to study the biomechanical properties of a molecule by interactively manipulating the molecule and see how it reacts to the user constraints.

BioSpring is also used to quickly build complex biostructural system, to provide for example pertinent starting points of Gromacs or NAMD simulation, or to interactively evaluate docking complex candidates.

We advice you to try BioSpring with MDDriver and VMD to understand how this interactive analysis could be used in your case.

#### Launching BioSpring with MDDriver support

To run BioSpring with VMD and MDDriver, start by launching the BioSpring tools with these MDDriver parameters.

	biospring -c model.nc -s param.msp --wait --port 8888 --log model.log

Or run the Docker image:

	docker run --rm -p 8888:8888 --init -v ./:/data ghcr.io/lbt-cnrs/biospring \
		biospring -s model.nc -c param.msp --wait --port 8888

Or use an alias for convenience:  
You can simplify the command by creating an alias and run the image with a shorter command:

	alias biospring-run='docker run --rm -p 8888:8888 --init -v ./:/data ghcr.io/lbt-cnrs/biospring biospring'

	biospring-run -s model.nc -c param.msp --wait --port 8888

Or reuse the Makefile in the example directory:  
Copy the file `example/base.mk` and some of other Makefile in example subdir as template for your own system (see README in example directory for usage)

#### Lauching VMD
Start the VMD program (double click on its icon or launch from the command line).

  - Open the description of the molecular system in the pdb file model.pdb that you have used for the simulation :
  - Connect VMD to the running simulation, specifying host name and port number on which the simulation was launched, then click connect.


#### Interact with the scene, using VMD tools

In a low Virtual Reality (VR) Context (no haptic devices, only the Mouse available) go to the menu Mouse>Force and choose a selection mode (atom, residue, fragment). Then you can pick atoms (residues or fragments) in the rendering window and influence the running simulation by applying additional user-defined forces.

For a higher VR context (with a haptic device for example), first lauch the VRPN server managing your VR device, then define a `.vmdsensor` file to configure your device (you might have to restart VMD when you change this file), and finally you can use the Graphics>Tools menu to define the interaction mode of your device with the simulation.



