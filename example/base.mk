
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

# PDB2SPN INPUTS
PDB2SPN_INPUT = pdb2spn_input.pdb

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
	cp $(RAW_PDB) raw.pdb
	@echo "Run Martinize"
	docker run --rm -v $(PWD):/data $(IMAGE) \
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
	docker run --rm --init -v $(PWD):/data $(IMAGE) \
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
run:
	@echo "Starting BioSpring in detached mode..."
	@CONTAINER_ID=$$(docker run -d --rm -p 8888:8888 --name biospring-container --init -v $(PWD):/data $(IMAGE) $(BIOSPRING) -s $(NC) -c $(MSP) --wait --port 8888 --sasa-classifier biospring); \
	if [ -z "$$CONTAINER_ID" ]; then \
		echo "$(RED)❌ Failed to start the container.$(RESET)"; \
		exit 1; \
	fi; \
	HOSTNAME=$$(hostname); \
	echo "$(GREEN)✅ BioSpring container started successfully!$(RESET)"; \
	echo ""; \
	echo "$(YELLOW)⚠️  Before continuing, please note that the BioSpring program is running now in idle mode within the container, waiting for a connection. By proceeding, we will reattach to the container and display BioSpring's logs in real-time. Afterward, you can use the following information to connect:$(RESET)"; \
	echo ""; \
	echo "  - On the local machine: $(GREEN)localhost:8888$(RESET) or $(GREEN)127.0.0.1:8888$(RESET)"; \
	echo "  - On the network: $(GREEN)$$HOSTNAME:8888$(RESET)"; \
	trap 'docker ps -q --filter "id=$$CONTAINER_ID" | grep -q . && docker stop $$CONTAINER_ID || true; \
		echo "$(YELLOW)⚠️ Program interrupted manually (SIGINT/SIGTERM).$(RESET)"; exit 130' INT TERM; \
	while true; do \
		read -p "🟢 Continue (c) or Kill the container (k)? " choice; \
		case $$choice in \
			[cC]) \
				echo "$(GREEN)Showing logs and attaching to the container...$(RESET)"; \
				docker logs -f $$CONTAINER_ID & \
				DOCKER_LOGS_PID=$$!; \
				docker attach $$CONTAINER_ID; \
				break;; \
			[kK]) \
				echo "$(RED)Killing the container...$(RESET)"; \
				if docker ps -q --filter "id=$$CONTAINER_ID" | grep -q .; then \
					docker stop $$CONTAINER_ID; \
				fi; \
				exit 0;; \
			*) \
				echo "$(YELLOW)Invalid input. Enter 'c' to continue or 'k' to kill.$(RESET)";; \
		esac; \
	done

# Run the BioSpring simulation without waiting for the connection
run_now:
	@echo "Run BioSpring NOOOW !!!!"
	docker run -p 8888:8888 --rm --init -v $(PWD):/data $(IMAGE) \
	$(BIOSPRING) -s $(NC) \
				 -c $(MSP) \
				 --sasa-classifier biospring

# If the client is on the same network but not on the same machine,
# you can expose the data to the client side before running the simulation.
# Run the command "make expose_data" on the host machine and on the client side
# run the command "make get_data HOSTNAME=<hostname>"
expose_data:
	@echo "Host machine: $(GREEN)$(shell hostname)$(RESET)"
	@echo "On the client side, you can retrive the data with the following command:"
	@echo "\n    $(GREEN)make get_data HOSTNAME=$(shell hostname)$(RESET)\n"
	@docker run -p 4000:4000 --rm --init -v $(PWD):/data $(IMAGE) python -m http.server 4000 --directory /data

# Retrieve the data from the host machine
get_data :
	@if [ -z "$(HOSTNAME)" ]; then \
		echo "$(RED)Usage: make get_data HOSTNAME=<hostname>$(RESET)"; \
		exit 1; \
	fi

	@echo "$(GREEN)Retrieving data from http://$(HOSTNAME):4000/model.* ...$(RESET)"
	@wget -q -r -np -nH  --accept-regex='model.*' http://$(HOSTNAME):4000/
	@rm index.html

pause:
	@read -p "Press enter to continue..." key

# Run the docker image in interactive mode if you need to run the tools manually
interactive:
	docker run --rm -it -v $(PWD):/data $(IMAGE)


usage:
	@echo "------------------------------------------------------------------------------"
	@echo "Usage:"
	@echo ""
	@echo "Step 1: Server side"
	@echo "  make prep"
	@echo "  make expose_data"
	@echo ""
	@echo "Step 2: Client side"
	@echo "  make get_data"
	@echo ""
	@echo "Step 3: Server side"
	@echo "  CTRL + C to stop the server exposing the data"
	@echo "  make run / make run_now"
	@echo "    - Check the line in the shell and note the HOSTNAME and PORT:"
	@echo "       MDDriver >      Interactive MD bind to /[HOSTNAME]/[PORT]"
	@echo ""
	@echo "Step 4: Client side"
	@echo "  open model.nc or model.pdb in VMD or UnityMol"
	@echo "  Connect to the server using the HOSTNAME and PORT"
	@echo ""
	@echo "Step 5: Server side"
	@echo "  CTRL + C to stop the server running the simulation"
	@echo "  make clean"
	@echo "------------------------------------------------------------------------------"