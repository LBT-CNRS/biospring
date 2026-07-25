
PDB2PDBCONECT = pdb2pdbconect
PDB2CDL = pdb2cdl
PDB2PQR = pdb2pqr
BIOSPRING = biospring

NCGEN=ncgen

GRP = model.grp
FF = model.ff
PDB = model.pdb
PQR = model.pqr
MSP = model.msp
MSPSPRING = model_spring.msp
MSPSPRINGSTERIC = model_spring_steric.msp
MSPSPRINGSTERICELECTROSTATIC = model_spring_steric_electrostatic.msp

PDB_CONECT = model.conect.pdb


CDL = model.cdl
NC = model.nc


# MDDriver option
# -wait 0 or 1, if the simulation server has to wait for client connection before starting
# -port port port
# -debub debug level : 0 no log - 2 verbose log
# -log logfilename
IMDLOG=mddriver.log
IMDOPT=-wait 1 -port 3000 -debug 0 -log $(IMDLOG)


prep :  $(PDB_CONECT) $(PDB) $(PQR) $(CDL) $(NC)


$(PQR) :
	$(PDB2PQR) -pdb $(PDB) -pqr $(PQR)  -grp $(GRP) -ff $(FF)


$(PDB_CONECT) :
	$(PDB2PDBCONECT) -pdb $(PDB) -pdbconect $(PDB_CONECT) -grp $(GRP) -ff $(FF)

$(CDL) : $(PDB) $(GRP) $(FF)
	$(PDB2CDL) -pdb $(PDB) -cdl $(CDL) -grp $(GRP) -ff $(FF) -springcutoff 9.0

$(NC) : $(CDL)
	$(NCGEN) -b -o $(NC) $(CDL)


run : $(NC) $(MSP)
	$(BIOSPRING)  -nc $(NC) -msp  $(MSP) $(IMDOPT)

run_spring : $(NC) $(MSP)
	$(BIOSPRING)  -nc $(NC) -msp  $(MSP) $(IMDOPT)

run_spring_steric : $(NC) $(MSPSPRINGSTERIC)
	$(BIOSPRING)  -nc $(NC) -msp  $(MSPSPRINGSTERIC) $(IMDOPT)

run_spring_steric_electrostatic : $(NC) $(MSPSPRINGSTERICELECTROSTATIC)
	$(BIOSPRING)  -nc $(NC) -msp  $(MSPSPRINGSTERICELECTROSTATIC) $(IMDOPT)


runopencl : $(NC) $(MSP)
	$(BIOSPRING)  -nc $(NC) -msp  $(MSP) $(IMDOPT) -opencl

gdb : $(NC) $(MSP)
	gdb --args $(BIOSPRING)  -nc $(NC) -msp  $(MSP) $(IMDOPT)


clean :
	rm -rf 	$(PDB_CONECT)
	rm -rf 	$(PQR)
	rm -rf $(NC)
	rm -rf $(IMDLOG)
