# This example needs the --rigidbody option, data/reducerules/
# ProteinAtomRigidGroups.rbody, the hbond mechanism and its .hbond table,
# none of which are yet part of any released BioSpring build: no Docker
# image has them yet, so there is no Makefile.docker here (unlike the other
# example/ subdirectories) -- a local build is required.

PDB2SPN = pdb2spn
BIOSPRING = biospring

# BS_DATA_DIR points at the forcefield/reduce-rule data bundled with a local
# BioSpring install: source <biospring_install_path>/bin/BioSpring-config.sh
# first to set BioSpring_DATA_PATH (see doc/BUILDING.md). Requires a build
# that includes data/reducerules/ProteinAtomRigidGroups.rbody.
BS_DATA_DIR = $(BioSpring_DATA_PATH)

# ProteinAtomRigidGroups.rbody: springs-only rotamer/hinge mechanism (all
# pairwise springs within each named rigid-body group -- CA tetrahedron,
# peptide plane, chi1..4 pivots, rings, hydroxyl/thiol rotations -- instead
# of a distance cutoff). See that file's header for the full convention.
RIGIDBODY = $(BS_DATA_DIR)/reducerules/ProteinAtomRigidGroups.rbody

# amber.grp/amber.ff: a genuine one-atom-per-rule identity mapping (not a
# real coarse-grain reduction) used here only to assign real per-atom AMBER
# masses/charges. Since it renames particles (e.g. N -> EN for Glu), pdb2spn
# reuses this same --grp file as --rigidbody's naming-translation table so
# it can still find atoms under their original name.
GRP = $(BS_DATA_DIR)/reducerules/amber.grp
FF = $(BS_DATA_DIR)/forcefield/amber.ff

# amber.hbond: donor/acceptor classification for the dynamic hydrogen-bond
# mechanism, keyed by amber.grp's renamed types (matches GRP/FF above). A
# local copy, not referenced through BS_DATA_DIR: model.msp's hbond.path
# points at it directly (hbond.enable/scale/cutoff/path all live in
# model.msp, see doc/MSP_Options.md), so no separate command-line option is
# needed at run/run_now/gdb below.

# Rigid-body springs represent real bonds/valence angles, not a soft elastic
# network: a much higher stiffness than a cutoff-based model, validated
# stable at the timestep set in model.msp (see that file's header).
STIFFNESS = 8000.0

# A synthetic 29-residue construct -- 10-residue helix + 10-residue loop +
# 9-residue beta hairpin -- built from scratch (ideal secondary-structure
# phi/psi via NeRF, side chains transplanted from real ubiquitin residues of
# the matching type, see doc/BUILDING_helix_loop_sheet.md-equivalent notes in
# the README) and supplied here already FULLY UNFOLDED (a single extended
# conformation, no helix/turn/sheet register at all). The point of this
# example is to pull the unfolded chain into shape interactively and watch
# whether --hbond finds and holds the right backbone hydrogen bonds as it
# refolds -- see the README for what to expect (and its real limits).
PDB = helix_loop_sheet.pdb

NC = model.nc
NC_PDB = model.nc.pdb
MSP = model.msp

# MDDriver option
# -wait, if the simulation server has to wait for client connection before starting (no value: presence = wait)
# -port port
# -debug debug level : 0 no log - 2 verbose log
# -log logfilename
IMDLOG = mddriver.log
IMDOPT = -wait -port 3000 -debug 0 -log $(IMDLOG)

prep : $(NC) $(NC_PDB)

$(NC) : $(PDB)
	$(PDB2SPN) -s $(PDB) -o $(NC) --grp $(GRP) --ff $(FF) --rigidbody $(RIGIDBODY) --stiffness $(STIFFNESS)

# PDB of the network with CONECT records for every rigid-body spring, for
# visualization (model.vmd loads this). Uses pdb2spn, not editspn: editspn
# has no --pdbconect equivalent. Passing --stiffness again here just
# restates the existing springs' own value (read from $(NC)); no cutoff/
# rigidbody is given, so no new spring is created at this step.
$(NC_PDB) : $(NC)
	$(PDB2SPN) -s $(NC) -o $(NC_PDB) --pdbconect --stiffness $(STIFFNESS)

run : $(NC) $(MSP)
	$(BIOSPRING) -nc $(NC) -msp $(MSP) $(IMDOPT)

run_now : $(NC) $(MSP)
	$(BIOSPRING) -nc $(NC) -msp $(MSP) -port 3000 -debug 0 -log $(IMDLOG)

vis : $(NC_PDB)
	vmd -e model.vmd

gdb : $(NC) $(MSP)
	gdb --args $(BIOSPRING) -nc $(NC) -msp $(MSP) $(IMDOPT)

clean :
	rm -rf $(NC) $(NC_PDB) $(IMDLOG)
