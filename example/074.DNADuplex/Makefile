# This example needs --rigidbody/--bondedinteraction and the nucleic data
# files (DNAAtomRigidGroups.rbody + NucleicAtomBonded.bi.ff), none of which
# are in a released BioSpring build yet: no Docker image has them, so there
# is no Makefile.docker here (unlike the older example/ subdirectories) --
# a local build is required.

PDB2SPN = pdb2spn
BIOSPRING = biospring

# BS_DATA_DIR points at the forcefield/reduce-rule data bundled with a local
# BioSpring install: source <biospring_install_path>/bin/BioSpring-config.sh
# first to set BioSpring_DATA_PATH (see doc/BUILDING.md).
BS_DATA_DIR = $(BioSpring_DATA_PATH)

RIGIDBODY = $(BS_DATA_DIR)/reducerules/DNAAtomRigidGroups.rbody
BONDED    = $(BS_DATA_DIR)/reducerules/NucleicAtomBonded.bi.ff

PDB = dna_duplex.pdb

# NO --grp and NO --nonbondedinteraction, and that is not an oversight:
# amber.grp maps the 20 amino acids only, and amber.nbi.ff has no nucleic
# entries, so there is nothing to map a nucleotide onto yet. This example is
# therefore the BONDED model alone -- see the README's "What this cannot
# show yet".
#
# STIFFNESS is what holds the furanose rings, which carry no bonded rule of
# their own (see the README): 3000 keeps them within a few degrees of the
# input pucker where the protein examples' 250 lets them relax. A stiffer
# mesh needs a shorter timestep, hence model.msp's 0.3 fs.
STIFFNESS = 3000.0

COMMON = --rigidbody $(RIGIDBODY) --stiffness $(STIFFNESS)

prep : dna.nc

# The full bonded model: real AMBER bonds, valence angles and torsions on
# top of the articulated rigid-body mesh.
dna.nc : $(PDB)
	$(PDB2SPN) -s $(PDB) -o $@ $(COMMON) --bondedinteraction $(BONDED) \
	           --dihedral --stretching --bending

run : dna.nc
	$(BIOSPRING) -nc dna.nc -msp model.msp

# Point energy, timestep 0, one step: nothing moves, so these numbers depend
# only on the model and the input conformation. This is the table the README
# quotes.
point : dna.nc
	$(BIOSPRING) -nc dna.nc -msp point.msp 2>&1 \
	  | sed 's/\x1b\[[0-9;]*m//g' \
	  | grep -E "Stretch energy|Bend energy|Dihedral energy"

clean :
	rm -rf dna.nc dna.nc.pdb traj.pdb energies.csv mddriver.log

.PHONY : prep run point clean
