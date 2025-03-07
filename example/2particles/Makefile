BASE_MK := $(shell echo $(PWD) | sed -E 's|(.*example).*|\1/base.mk|')

# Inclure base.mk si trouvé
ifneq ($(wildcard $(BASE_MK)),)
    include $(BASE_MK)
    # $(info Included: $(BASE_MK))
endif


# ------------------------------------------------------------------------------
# Usage: make run / make run_now
# ------------------------------------------------------------------------------