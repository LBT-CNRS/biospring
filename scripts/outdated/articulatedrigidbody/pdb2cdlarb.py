import sys
import argparse
import os

from atomicstructure import AtomicStructure
from rigidbodystructure import RigidBodyStructure
from springnetwork import SpringNetwork

def usage() :
	print 'pdb2cdlarb allows to create a spring network from a PDB file to simulate articulated rigid body based on rotameres degrees of freedom.'
	print '\n'

	print 'usage: python pdb2cdlarb [option] -pdb pdbname -cdl cdlname'
	print '\n'

	print 'positional arguments:'
	print '\t-pdb pdbname                   input pdb filename of dynamic structure'
	print '\t-cdl cdlname                   output cdl filename'

	print '\n'

	print 'optional arguments:'
#	print '\t-ff amber94                    forcefield directory name to use in the ff directory'
	print '\t-pdbout pdb                    output pdb filename for visualization purpose'
	print '\t-rotamere rotamerename         list of couples of atomname in pdb that are rotamere'
	print '\t--help                         display this help and exit'

def main(argv):
	description='pdb2cdlarb allows to create a spring network from a PDB file to simulate articulated rigid body based on rotameres degrees of freedom.'
	#parser = argparse.ArgumentParser(usage=usage(), description=description)
	parser = argparse.ArgumentParser(description=description)

	parser.add_argument('-pdb', action="store", help='foo help', dest='pdb')
	parser.add_argument('-cdl', action="store", help='foo help', dest='cdl')

	parser.add_argument('-ff', action="store", help='foo help', dest='ff')
	parser.add_argument('-pdbout', action="store", help='foo help', dest='pdbout')
	parser.add_argument('-rotamere', action="store", help='foo help', dest='rotamere')

	results = parser.parse_args()
	parser.print_usage()

	ff='amber94'

	if results.pdb==None or results.cdl==None :
		usage()
		sys.exit(1)
	else :
		s = AtomicStructure(ff)
		if results.rotamere!=None :
			s.parseRotamere(results.rotamere)
		s.buildFromPDB(results.pdb)
		#print s

		rb=RigidBodyStructure()
		rb.buildFromStructure(s)
		print rb

		sn=SpringNetwork()
		sn.buildFromRigidBodyStructure(rb)
		#print sn

		sn.toCDL(results.cdl)

		if results.pdbout!=None :
			sn.toPDB(results.pdbout)


if __name__ == "__main__":
	main(sys.argv[1:])







