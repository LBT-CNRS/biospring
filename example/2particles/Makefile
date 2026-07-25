PDB2SPN = pdb2spn
BIOSPRING = biospring
NCGEN = ncgen

PDB = model.pdb
CDL = model.cdl
NC = model.nc
MSP = param.msp

# model.vmd expects a PDB with CONECT records under this exact name.
PDB_CONNECT = model.connect.pdb

# No reduction here (no --grp/--ff): each atom of model.pdb becomes one
# particle directly, matching the "2 particles" system.
CUTOFF = 5.0

# MDDriver option
# -wait, if the simulation server has to wait for client connection before starting (no value: presence = wait)
# -port port
# -debug debug level : 0 no log - 2 verbose log
# -log logfilename
IMDLOG = mddriver.log
IMDOPT = -wait -port 3000 -debug 0 -log $(IMDLOG)

prep : $(NC) $(PDB_CONNECT)

$(CDL) : $(PDB)
	$(PDB2SPN) -s $(PDB) -o $(CDL) --cutoff $(CUTOFF)

$(NC) : $(CDL)
	$(NCGEN) -b -o $(NC) $(CDL)

$(PDB_CONNECT) : $(PDB)
	$(PDB2SPN) -s $(PDB) -o $(PDB_CONNECT) --cutoff $(CUTOFF) --pdbconect

run : $(NC) $(MSP)
	$(BIOSPRING) -nc $(NC) -msp $(MSP) $(IMDOPT)

run_now : $(NC) $(MSP)
	$(BIOSPRING) -nc $(NC) -msp $(MSP) -port 3000 -debug 0 -log $(IMDLOG)

gdb : $(NC) $(MSP)
	gdb --args $(BIOSPRING) -nc $(NC) -msp $(MSP) $(IMDOPT)

clean :
	rm -rf $(CDL)
	rm -rf $(NC)
	rm -rf $(PDB_CONNECT)
	rm -rf $(IMDLOG)

# ------------------------------------------------------------------------------
# Docker-based execution -- use when biospring is not installed locally.
# $(NC)/$(MSP) must already exist (see the classic prep targets above, or run
# them once inside the container with `make interactive`).
# ------------------------------------------------------------------------------
RED    = \033[0;31m
GREEN  = \033[0;32m
YELLOW = \033[0;33m
RESET  = \033[0m

IMAGE = ghcr.io/lbt-cnrs/biospring

run_docker: $(NC) $(MSP)
	@echo "Starting BioSpring in detached mode..."
	@CONTAINER_ID=$$(docker run -d --rm -p 8888:8888 --name biospring-container --init -v $(PWD):/data $(IMAGE) $(BIOSPRING) -s $(NC) -c $(MSP) --wait --port 8888); \
	if [ -z "$$CONTAINER_ID" ]; then \
		echo "$(RED)Failed to start the container.$(RESET)"; \
		exit 1; \
	fi; \
	HOSTNAME=$$(hostname); \
	echo "$(GREEN)BioSpring container started successfully!$(RESET)"; \
	echo ""; \
	echo "  - On the local machine: $(GREEN)localhost:8888$(RESET) or $(GREEN)127.0.0.1:8888$(RESET)"; \
	echo "  - On the network: $(GREEN)$$HOSTNAME:8888$(RESET)"; \
	trap 'docker ps -q --filter "id=$$CONTAINER_ID" | grep -q . && docker stop $$CONTAINER_ID || true; exit 130' INT TERM; \
	while true; do \
		read -p "Continue (c) or Kill the container (k)? " choice; \
		case $$choice in \
			[cC]) \
				docker logs -f $$CONTAINER_ID & \
				docker attach $$CONTAINER_ID; \
				break;; \
			[kK]) \
				if docker ps -q --filter "id=$$CONTAINER_ID" | grep -q .; then \
					docker stop $$CONTAINER_ID; \
				fi; \
				exit 0;; \
			*) \
				echo "$(YELLOW)Invalid input. Enter 'c' to continue or 'k' to kill.$(RESET)";; \
		esac; \
	done

run_now_docker: $(NC) $(MSP)
	docker run -p 8888:8888 --rm --init -v $(PWD):/data $(IMAGE) \
	$(BIOSPRING) -s $(NC) -c $(MSP)

# Run the docker image in interactive mode if you need to run the tools manually
interactive:
	docker run --rm -it -v $(PWD):/data $(IMAGE)
