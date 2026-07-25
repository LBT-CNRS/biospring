NC = model.nc # USED AS BIOSPRING INPUT
MSP = param.msp

BIOSPRING = biospring

# ------------------------------------------------------------------------------
# Classic (local) execution -- requires biospring built/installed locally.
# model.nc/param.msp must already exist: this example's model.nc depends on
# forcefield/reduce-rule data bundled in the docker image, so there is no
# local prep here -- see Makefile.docker for that step.
# ------------------------------------------------------------------------------
IMDLOG = mddriver.log
IMDOPT = -wait -port 8888 -debug 0 -log $(IMDLOG)

run : $(NC) $(MSP)
	$(BIOSPRING) -nc $(NC) -msp $(MSP) $(IMDOPT)

run_now : $(NC) $(MSP)
	$(BIOSPRING) -nc $(NC) -msp $(MSP) -port 8888 -debug 0 -log $(IMDLOG)
