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

# 8000 kJ/mol/A^2, not the 250 this example used to carry. The mesh is now
# the ONLY thing holding bonds and angles, so it has to be stiff enough to
# hold them against the dihedral rings: measured on this same ubiquitin,
# quenched 6 ps, AMBER's own judgement of the resulting geometry --
#
#     --stiffness     bonds    angles      dt max
#     250          29211.52  26949.30      2.0 fs
#     2500           968.37    729.06      0.5 fs
#     8000           447.81    344.91      1.0 fs
#     25000          424.96    287.52      0.5 fs
#
# in kJ/mol, against 443.69 and 278.83 for the input structure itself.
# 8000 is where the geometry is back to its input quality while a 1 fs step
# is still stable (dt scales as 1/sqrt(k), verified exactly).
STIFFNESS = 8000.0
PDB = ubiquitin.pdb

# One .nc PER STAGE -- and that is the whole point of this example. The two
# stages CANNOT be reproduced by toggling one .nc at runtime: a dihedral
# family not requested at build time never gets a spring, a particle or a
# NetCDF entry, so no .msp key can bring it back. See the README, and
# doc/MSP_Options.md.
COMMON = --grp $(GRP) --nonbondedinteraction $(FF) --rigidbody $(RIGIDBODY) --stiffness $(STIFFNESS)

IMDLOG = mddriver.log
IMDOPT = -wait -port 3000 -debug 0 -log $(IMDLOG)

prep : stage1.nc stage2.nc

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

run-stage1 : stage1.nc
	$(BIOSPRING) -nc stage1.nc -msp model.msp $(IMDOPT)

run-stage2 : stage2.nc
	$(BIOSPRING) -nc stage2.nc -msp model.msp $(IMDOPT)

run : run-stage2

# Point energies of the two stages side by side, on the SAME input
# conformation (timestep 0, one step: nothing moves, so the numbers are
# directly comparable). This is the table the README quotes.
compare : prep
	@printf "%-10s %12s %12s\n" stage spring dihedral
	@for s in 1 2; do \
	    $(BIOSPRING) -nc stage$$s.nc -msp point.msp 2>&1 \
	      | sed 's/\x1b\[[0-9;]*m//g' \
	      | awk -v s=$$s '/Spring energy/{sp=$$4} /Dihedral energy/{di=$$4} \
	                      END{printf "%-10s %12s %12s\n", s, sp, di}'; \
	done

# PDB of stage 2's network with CONECT records for every spring, for
# visualization. No cutoff/rigidbody is given here, so no new spring is
# created at this step -- it only re-emits what stage2.nc already holds.
stage2.nc.pdb : stage2.nc
	$(PDB2SPN) -s stage2.nc -o $@ --pdbconect --stiffness $(STIFFNESS)

vis : stage2.nc.pdb
	vmd -e model.vmd

clean :
	rm -rf stage1.nc stage2.nc stage2.nc.pdb $(IMDLOG)

.PHONY : prep run run-stage1 run-stage2 compare vis clean
