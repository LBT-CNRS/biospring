BASE_MK := $(shell echo $(PWD) | sed -E 's|(.*example).*|\1/base.mk|')

# Inclure base.mk si trouvé
ifneq ($(wildcard $(BASE_MK)),)
    include $(BASE_MK)
    # $(info Included: $(BASE_MK))
endif

# RAW PDB INPUT
RAW_PDB = ../1AF6_in_membrane.pdb

# PDB2SPN INPUTS
FF = $(BS_DATA_DIR)/martini/forcefield/martini_octanol.ff
GRP =  $(BS_DATA_DIR)/martini/reducerules/martini.grp
CUTOFF = 5.0
STIFFNESS = 30.0

# Prepare the input files for the simulation
prep : martinize pause pdb2spn

# ------------------------------------------------------------------------------
# Usage: run the command make usage
# ------------------------------------------------------------------------------