# This example needs --rigidbody/--bondedinteraction and the nucleic data
# files (RNAAtomRigidGroups.rbody + NucleicAtomBonded.bi.ff), none of which
# are in a released BioSpring build yet: no Docker image has them, so there
# is no Makefile.docker here (unlike the older example/ subdirectories) --
# a local build is required.

PDB2SPN = pdb2spn
BIOSPRING = biospring

# BS_DATA_DIR points at the forcefield/reduce-rule data bundled with a local
# BioSpring install: source <biospring_install_path>/bin/BioSpring-config.sh
# first to set BioSpring_DATA_PATH (see doc/BUILDING.md).
BS_DATA_DIR = $(BioSpring_DATA_PATH)

RIGIDBODY = $(BS_DATA_DIR)/reducerules/RNAAtomRigidGroups.rbody
BONDED    = $(BS_DATA_DIR)/reducerules/NucleicAtomBonded.bi.ff

PDB = rna_hairpin.pdb

GRP       = $(BS_DATA_DIR)/reducerules/amber.rna.grp
FF        = $(BS_DATA_DIR)/forcefield/amber.rna.nbi.ff

# STIFFNESS is what holds the furanose rings, which carry no bonded rule of
# their own (see the README): the higher it is, the closer the pucker stays
# to the input. A stiffer mesh needs a shorter timestep -- the two are tied,
# dt scaling as 1/sqrt(stiffness) -- hence model.msp's 0.3 fs.
STIFFNESS = 3000.0

COMMON = --grp $(GRP) --nonbondedinteraction $(FF) \
         --rigidbody $(RIGIDBODY) --stiffness $(STIFFNESS)

prep : rna.nc

# The full bonded model: real AMBER bonds, valence angles and torsions on
# top of the articulated rigid-body mesh.
rna.nc : $(PDB)
	$(PDB2SPN) -s $(PDB) -o $@ $(COMMON) --bondedinteraction $(BONDED) \
	           --dihedral --stretching --bending

run : rna.nc
	$(BIOSPRING) -nc rna.nc -msp model.msp

# Point energy, timestep 0, one step: nothing moves, so these numbers depend
# only on the model and the input conformation. This is the table the README
# quotes.
point : rna.nc
	$(BIOSPRING) -nc rna.nc -msp point.msp 2>&1 \
	  | sed 's/\x1b\[[0-9;]*m//g' \
	  | grep -E "Stretch energy|Bend energy|Dihedral energy"

clean :
	rm -rf rna.nc rna.nc.pdb traj.pdb energies.csv mddriver.log

.PHONY : prep run point clean
