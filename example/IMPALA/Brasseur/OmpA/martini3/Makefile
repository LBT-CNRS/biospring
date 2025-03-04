RED    = \033[0;31m
GREEN  = \033[0;32m
YELLOW = \033[0;33m
BLUE   = \033[0;34m
RESET  = \033[0m
# ---------------------------------------------------------------------

IMAGE = ghcr.io/lbt-cnrs/biospring
PDB2SPN = pdb2spn
BIOSPRING = biospring
MARTINIZE = martinize2

BS_DATA_DIR = /biospring/share/biospring/data
FORCEFIELD_DIR = $(BS_DATA_DIR)/martini/forcefield
REDUCERULES_DIR = $(BS_DATA_DIR)/martini/reducerules

# PDB2SPN INPUTS
FF = $(FORCEFIELD_DIR)/martini_octanol.ff
GRP =  $(REDUCERULES_DIR)/martini.grp
RAW_PDB = ../model_in_membrane.pdb

# PDB2SPN OUTPUTS
CDL = model.cdl
PDB = model.pdb
NC = model.nc # USED AS BIOSPRING INPUT

# BIOSPRING INPUT PARAMETERS
MSP = param.msp

# MARTINIZE TOOL
# INPUT: model_in_membrane.pdb (1BXW.pdb positioned in the membrane)
# OUTPUT: pdb2spn_input_cg.pdb (PDB2SPN INPUT)
martinize:
	@echo "Copy PDB input for Martinize"
	cp $(RAW_PDB) .
	@echo "Run Martinize"
	docker run --rm -v $(PWD):/data $(IMAGE) \
    $(MARTINIZE) -f model_in_membrane.pdb \
                 -x pdb2spn_input_cg.pdb \
                 -ff martini30b32 \
                 -resid input \
                 -ignh
	rm $(notdir $(RAW_PDB))

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
	@echo "Run pdb2spn"
	docker run --rm --init -v $(PWD):/data $(IMAGE) \
	$(PDB2SPN) -s pdb2spn_input_cg.pdb \
			   -o ${CDL} ${NC} $(PDB) \
			   --ff $(FF) \
			   --grp $(GRP) \
			   --cutoff 5.0 \
			   --stiffness 30.0

	
clean:
	rm -f \#*#
	rm -f pdb2spn_input_cg.pdb
	rm -f model.cdl model.nc model.pdb

pause:
	@read -p "Press enter to continue..." key

# Run the docker image in interactive mode if you need to run the tools manually
interactive:
	docker run --rm -it -v $(PWD):/data $(IMAGE)

# Prepare the input files for the simulation
prep : martinize pause pdb2spn

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

expose_data:
	@echo "Host machine: $(GREEN)$(shell hostname)$(RESET)"
	@echo "On the client side, you can retrive the data with the following command:"
	@echo "\n    $(GREEN)make get_data HOSTNAME=$(shell hostname)$(RESET)\n"
	@docker run --network host --rm --init -v $(PWD):/data $(IMAGE) python -m http.server 4000 --directory /data

get_data :
	@if [ -z "$(HOSTNAME)" ]; then \
		echo "$(RED)Usage: make get_data HOSTNAME=<hostname>$(RESET)"; \
		exit 1; \
	fi

	@echo "$(GREEN)Retrieving data from http://$(HOSTNAME):4000/model.* ...$(RESET)"
	@wget -q -r -np -nH  --accept-regex='model.*' http://$(HOSTNAME):4000/
	@rm index.html

