# This example needs the --rigidbody/--bondedinteraction options and
# data/reducerules/ProteinAtomRigidGroups.rbody + ProteinAtomBonded.bi.ff,
# none of which are part of a released BioSpring build yet: no Docker image
# has them, so there is no Makefile.docker here (unlike the other example/
# subdirectories) -- a local build is required.

PDB2SPN = pdb2spn
BIOSPRING = biospring

# BS_DATA_DIR points at the forcefield/reduce-rule data bundled with a local
# BioSpring install: source <biospring_install_path>/bin/BioSpring-config.sh
# first to set BioSpring_DATA_PATH (see doc/BUILDING.md).
BS_DATA_DIR = $(BioSpring_DATA_PATH)

RIGIDBODY = $(BS_DATA_DIR)/reducerules/ProteinAtomRigidGroups.rbody
BONDED    = $(BS_DATA_DIR)/reducerules/ProteinAtomBonded.bi.ff
GRP       = $(BS_DATA_DIR)/reducerules/amber.grp
FF        = $(BS_DATA_DIR)/forcefield/amber.nbi.ff

STIFFNESS = 250.0
PDB = ubiquitin.pdb

# One .nc PER STAGE -- and that is the whole point of this example. The
# stages CANNOT be reproduced by toggling one .nc at runtime, because
# --stretching and --bending retune --rigidbody's springs IN PLACE: the
# uniform-stiffness model is overwritten at build time and no .msp key can
# bring it back. See the README, and doc/MSP_Options.md.
COMMON = --grp $(GRP) --nonbondedinteraction $(FF) --rigidbody $(RIGIDBODY) --stiffness $(STIFFNESS)

IMDLOG = mddriver.log
IMDOPT = -wait -port 3000 -debug 0 -log $(IMDLOG)

prep : stage1.nc stage2.nc stage3.nc

# Stage 1 -- articulated rigid body only. Every spring is the uniform
# --stiffness value; the .rbody cliques decide what is rigid and what stays
# free to hinge (phi/psi/chi). No force-field energy of any kind.
stage1.nc : $(PDB)
	$(PDB2SPN) -s $(PDB) -o $@ $(COMMON)

# Stage 2 -- + real dihedral wells. Bonds and angles STAY at the uniform
# value: only the torsions become real AMBER, so the difference against
# stage 1 is exactly what the dihedral model contributes and nothing else.
stage2.nc : $(PDB)
	$(PDB2SPN) -s $(PDB) -o $@ $(COMMON) --bondedinteraction $(BONDED) --dihedral

# Stage 3 -- + real bonds and angles. --bending REQUIRES --stretching (its
# 1-3 conversion needs STRETCH's real bond lengths), so these two arrive
# together; pdb2spn refuses --bending alone rather than silently ignoring it.
stage3.nc : $(PDB)
	$(PDB2SPN) -s $(PDB) -o $@ $(COMMON) --bondedinteraction $(BONDED) --dihedral --stretching --bending

run-stage1 : stage1.nc
	$(BIOSPRING) -nc stage1.nc -msp model.msp $(IMDOPT)

run-stage2 : stage2.nc
	$(BIOSPRING) -nc stage2.nc -msp model.msp $(IMDOPT)

run-stage3 : stage3.nc
	$(BIOSPRING) -nc stage3.nc -msp model.msp $(IMDOPT)

run : run-stage3

# Point energies of the three stages side by side, on the SAME input
# conformation (timestep 0, one step: nothing moves, so the three numbers
# are directly comparable). This is the table the README quotes.
compare : prep
	@printf "%-10s %12s %12s %12s %12s\n" stage spring stretch bend dihedral
	@for s in 1 2 3; do \
	    $(BIOSPRING) -nc stage$$s.nc -msp point.msp 2>&1 \
	      | sed 's/\x1b\[[0-9;]*m//g' \
	      | awk -v s=$$s '/Spring energy/{sp=$$4} /Stretch energy/{st=$$4} \
	                      /Bend energy/{be=$$4} /Dihedral energy/{di=$$4} \
	                      END{printf "%-10s %12s %12s %12s %12s\n", s, sp, st, be, di}'; \
	done

# PDB of stage 3's network with CONECT records for every spring, for
# visualization. No cutoff/rigidbody is given here, so no new spring is
# created at this step -- it only re-emits what stage3.nc already holds.
stage3.nc.pdb : stage3.nc
	$(PDB2SPN) -s stage3.nc -o $@ --pdbconect --stiffness $(STIFFNESS)

vis : stage3.nc.pdb
	vmd -e model.vmd

clean :
	rm -rf stage1.nc stage2.nc stage3.nc stage3.nc.pdb $(IMDLOG)

.PHONY : prep run run-stage1 run-stage2 run-stage3 compare vis clean
