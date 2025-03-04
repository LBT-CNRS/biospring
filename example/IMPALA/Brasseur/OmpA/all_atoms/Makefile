IMAGE = ghcr.io/lbt-cnrs/biospring
PDB2SPN = pdb2spn
BIOSPRING = biospring

BS_DATA_DIR = /biospring/share/biospring/data
FORCEFIELD_DIR = $(BS_DATA_DIR)/forcefield
REDUCERULES_DIR = $(BS_DATA_DIR)/reducerules


# PDB2SPN INPUTS
FF = $(FORCEFIELD_DIR)/amber.ff
GRP =  $(REDUCERULES_DIR)/amber.grp
RAW_PDB = ../model_in_membrane.pdb

# PDB2SPN OUTPUTS
CDL = model.cdl
PDB = model.pdb
NC = model.nc # USED AS BIOSPRING INPUT

# BIOSPRING INPUT PARAMETERS
MSP = param.msp


# PDB2SPN TOOL
# INPUT: model_in_membrane.pdb (1BXW.pdb positioned in the membrane)
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
	@echo "Copy PDB input for Martinize"
	cp $(RAW_PDB) .
	@echo "Run pdb2spn"
	docker run --rm --init -v $(PWD):/data $(IMAGE) \
	$(PDB2SPN) -s model_in_membrane.pdb \
			   -o ${CDL} ${NC} $(PDB) \
			   --grp $(GRP) \
			   --ff $(FF) \
			   --cutoff 4.0 \
			   --stiffness 30.0
	rm $(notdir $(RAW_PDB))

	
clean:
	rm -f \#*#
	rm -f model.cdl model.nc model.pdb

pause:
	@read -p "Press enter to continue..." key

# Run the docker image in interactive mode if you need to run the tools manually
interactive:
	docker run --rm -it -v $(PWD):/data $(IMAGE)

# Prepare the input files for the simulation
prep : pdb2spn

# Run the BioSpring simulation
# NOTE: The option --sasa-classifier biospring is important to use IMPALA
#       to at least compute once the particles surface areas. Otherwire you'll
#       have a warning message like "!! WARNING: freesasa_classifier_radius: 
#       radius not found fo res: ARG name: RBB"
# 
# NOTE2: By launching biospring, all the parameters used are displayed in the 
#        terminal. You can modify the parameters in the MSP file as you wish 
#        and restart biospring or tune certain parameters from UnityMol during 
#        an interactive simulation (for example IMPALA scale or Input Force, or 
#        even viscosity).
run : 
	@echo "Run BioSpring"
	docker run --network host --rm --init -v $(PWD):/data $(IMAGE) \
	$(BIOSPRING) -s $(NC) \
				 -c $(MSP) \
				 --wait --port 3000 \
				 --sasa-classifier biospring
