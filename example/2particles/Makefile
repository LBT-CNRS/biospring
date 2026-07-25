BASE_MK := $(shell d="$(CURDIR)"; while [ "$$d" != "/" ] && [ ! -f "$$d/base.mk" ]; do d=$$(dirname "$$d"); done; [ -f "$$d/base.mk" ] && echo "$$d/base.mk")

# Include base.mk if found
ifneq ($(wildcard $(BASE_MK)),)
    include $(BASE_MK)
    # $(info Included: $(BASE_MK))
endif


# ------------------------------------------------------------------------------
# Usage: make run / make run_now
# ------------------------------------------------------------------------------