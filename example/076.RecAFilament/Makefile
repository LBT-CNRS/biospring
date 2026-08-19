# This example needs --rigidbody/--bondedinteraction and the nucleic data
# files, none of which are in a released BioSpring build yet: a local build
# is required, and there is no Makefile.docker.

PDB2SPN = pdb2spn
BIOSPRING = biospring

# BS_DATA_DIR points at the forcefield/reduce-rule data bundled with a local
# BioSpring install: source <biospring_install_path>/bin/BioSpring-config.sh
# first to set BioSpring_DATA_PATH (see doc/BUILDING.md).
BS_DATA_DIR = $(BioSpring_DATA_PATH)

PDB = reca_dna.pdb

# THE POINT OF THIS EXAMPLE: one system, two chemistries, so the protein and
# the nucleic parameter sets have to be used TOGETHER. pdb2spn takes one
# --grp and one --nonbondedinteraction, so they are concatenated here.
#
# That works because the type names cannot collide: amber.grp prefixes every
# protein type with an UPPERCASE residue letter (ACA for ALA's CA) and
# amber.dna.grp uses lowercase and digits (dC5, 1C1). The two files were
# built to be concatenated -- this example is what checks that they are.
GRP_PARTS = $(BS_DATA_DIR)/reducerules/amber.grp \
            $(BS_DATA_DIR)/reducerules/amber.dna.grp
FF_PARTS  = $(BS_DATA_DIR)/forcefield/amber.nbi.ff \
            $(BS_DATA_DIR)/forcefield/amber.dna.nbi.ff

GRP = complex.grp
FF  = complex.nbi.ff

# Likewise for the rigid-body and bonded rules: the protein and nucleic
# files describe disjoint residues, so their union is simply both.
RBODY_PARTS  = $(BS_DATA_DIR)/reducerules/ProteinAtomRigidGroups.rbody \
               $(BS_DATA_DIR)/reducerules/DNAAtomRigidGroups.rbody
BONDED_PARTS = $(BS_DATA_DIR)/reducerules/ProteinAtomBonded.bi.ff \
               $(BS_DATA_DIR)/reducerules/NucleicAtomBonded.bi.ff

RIGIDBODY = complex.rbody
BONDED    = complex.bi.ff

# 250 kJ/mol/A^2 is the protein examples' validated value. The furanose
# rings are held by this mesh alone and relax somewhat at 250 (see
# 074.DNADuplex's README) -- the compromise is deliberate here: this
# example is about the two chemistries coexisting, and raising the
# stiffness for 13 nucleotides would stiffen 1188 protein residues with
# them.
STIFFNESS = 250.0

IMDLOG = mddriver.log
IMDOPT = -wait -port 3000 -debug 0 -log $(IMDLOG)

prep : reca.nc

$(GRP) : $(GRP_PARTS)
	cat $(GRP_PARTS) > $@

$(FF) : $(FF_PARTS)
	cat $(FF_PARTS) > $@

$(RIGIDBODY) : $(RBODY_PARTS)
	cat $(RBODY_PARTS) > $@

$(BONDED) : $(BONDED_PARTS)
	cat $(BONDED_PARTS) > $@

reca.nc : $(PDB) $(GRP) $(FF) $(RIGIDBODY) $(BONDED)
	$(PDB2SPN) -s $(PDB) -o $@ --grp $(GRP) --nonbondedinteraction $(FF) \
	           --rigidbody $(RIGIDBODY) --stiffness $(STIFFNESS) \
	           --bondedinteraction $(BONDED) --dihedral --stretching --bending

run : reca.nc
	$(BIOSPRING) -nc reca.nc -msp model.msp $(IMDOPT)

point : reca.nc
	$(BIOSPRING) -nc reca.nc -msp point.msp 2>&1 \
	  | sed 's/\x1b\[[0-9;]*m//g' \
	  | grep -E "Stretch energy|Bend energy|Dihedral energy"

clean :
	rm -rf reca.nc reca.nc.pdb $(GRP) $(FF) $(RIGIDBODY) $(BONDED) \
	       traj.pdb energies.csv $(IMDLOG)

.PHONY : prep run point clean
