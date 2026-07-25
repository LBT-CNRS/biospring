# ------------------------------------------------------------------------------
# Example-specific configuration
# ------------------------------------------------------------------------------
RAW_PDB = ../model_in_membrane.pdb

# PDB2SPN INPUTS
FF = $(BS_DATA_DIR)/forcefield/amber.ff
GRP =  $(BS_DATA_DIR)/reducerules/amber.grp
CUTOFF = 5.0
STIFFNESS = 30.0

# Prepare the input files for the simulation
prep : pdb2spn

# ------------------------------------------------------------------------------
# Shared variables
# ------------------------------------------------------------------------------
PDB2SPN = pdb2spn
BIOSPRING = biospring
MARTINIZE = martinize2

# BS_DATA_DIR points at the forcefield/reduce-rule data bundled with a local
# BioSpring install: source <biospring_install_path>/bin/BioSpring-config.sh
# first to set BioSpring_DATA_PATH (see doc/BUILDING.md).
BS_DATA_DIR = $(BioSpring_DATA_PATH)

PDB2SPN_INPUT = pdb2spn_input.pdb

CDL = model.cdl
PDB = model.pdb
NC = model.nc # USED AS BIOSPRING INPUT
MSP = param.msp

# ------------------------------------------------------------------------------
# Classic (local) preparation and execution -- requires biospring/pdb2spn/
# martinize2 installed locally (and BioSpring_DATA_PATH set, see above).
# ------------------------------------------------------------------------------

# MARTINIZE TOOL
# INPUT: model_in_membrane.pdb (1BXW.pdb positioned in the membrane)
# OUTPUT: pdb2spn_input_cg.pdb (PDB2SPN INPUT)
martinize:
	@echo "Copy PDB input for Martinize"
	cp $(RAW_PDB) raw.pdb
	@echo "Run Martinize"
	$(MARTINIZE) -f raw.pdb \
                 -x $(PDB2SPN_INPUT) \
                 -ff martini30b32 \
                 -resid input \
                 -ignh
	rm raw.pdb

# PDB2SPN TOOL
# INPUT: pdb2spn_input_cg.pdb (Martinize output, CG PDB file)
# OUTPUT:
#  - model.cdl: Readable NetCDF file of the structure, used to verify the
#               structure informations
#
#  - model.nc: Binary  NetCDF file of the structure, used as BIOSPRING input
#              and as input structure for visualization tools like UnityMol
#              (in "BioSpring" version of UnityMol - March 2025)
#
#  - model.pdb: PDB file of the structure, used as input structure for
#               visualization tools (like UnityMol or VMD) and to determine the
#               insertion vector particle pair (see doc/MSP_Options.md)
pdb2spn:
	@if [ ! -f $(PDB2SPN_INPUT) ]; then \
		if [ ! -f $(RAW_PDB) ]; then \
			echo "Error: RAW_PDB file $(RAW_PDB) not found!"; \
			exit 1; \
		fi; \
		cp $(RAW_PDB) $(PDB2SPN_INPUT); \
	fi

	@echo "Run pdb2spn"
	$(PDB2SPN) -s $(PDB2SPN_INPUT) \
			   -o ${CDL} ${NC} $(PDB) \
			   --ff $(FF) \
			   --grp $(GRP) \
			   --cutoff $(CUTOFF) \
			   --stiffness $(STIFFNESS)

clean:
	rm -f \#*#
	rm -f $(PDB2SPN_INPUT)
	rm -f model.cdl model.nc model.pdb

pause:
	@read -p "Press enter to continue..." key

IMDLOG = mddriver.log
IMDOPT = -wait -port 8888 -debug 0 -log $(IMDLOG)

run : $(NC) $(MSP)
	$(BIOSPRING) -nc $(NC) -msp $(MSP) $(IMDOPT)

run_now : $(NC) $(MSP)
	$(BIOSPRING) -nc $(NC) -msp $(MSP) -port 8888 -debug 0 -log $(IMDLOG)
