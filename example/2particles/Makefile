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
# Docker-based execution, if biospring is not installed locally: see
# Makefile.docker (make -f Makefile.docker run / make -f Makefile.docker usage)
# ------------------------------------------------------------------------------
