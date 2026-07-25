BASE_MK := $(shell d="$(CURDIR)"; while [ "$$d" != "/" ] && [ ! -f "$$d/base.mk" ]; do d=$$(dirname "$$d"); done; [ -f "$$d/base.mk" ] && echo "$$d/base.mk")

# Include base.mk if found
ifneq ($(wildcard $(BASE_MK)),)
    include $(BASE_MK)
    # $(info Included: $(BASE_MK))
endif

# RAW PDB INPUT
RAW_PDB = ../1AF6_in_membrane.pdb

# PDB2SPN INPUTS
FF = $(BS_DATA_DIR)/forcefield/amber.ff
GRP =  $(BS_DATA_DIR)/reducerules/amber.grp
CUTOFF = 5.0
STIFFNESS = 30.0

# Prepare the input files for the simulation
prep : pdb2spn

# ------------------------------------------------------------------------------
# Usage: run the command make usage
# ------------------------------------------------------------------------------