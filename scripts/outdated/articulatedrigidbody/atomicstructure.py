import gromacs
import gromacs.fileformats as fileformats
import os

from forcefield import ForceField
from bond import Bond
from particle import Particle


class AtomicStructure () :
	def __init__(self,ff):
		self.atoms={}
		self.bonds={}


		self.aaabr2letter={}
		self.aaabr2letter["ALA"]='A'
		self.aaabr2letter["CYS"]='C'
		self.aaabr2letter["ASP"]='D'
		self.aaabr2letter["GLU"]='E'
		self.aaabr2letter["PHE"]='F'
		self.aaabr2letter["GLY"]='G'
		self.aaabr2letter["HIS"]='H'
		self.aaabr2letter["ILE"]='I'
		self.aaabr2letter["LYS"]='K'
		self.aaabr2letter["LEU"]='L'
		self.aaabr2letter["MET"]='M'
		self.aaabr2letter["ASN"]='N'
		self.aaabr2letter["PYL"]='O'
		self.aaabr2letter["PRO"]='P'
		self.aaabr2letter["GLN"]='Q'
		self.aaabr2letter["ARG"]='R'
		self.aaabr2letter["SER"]='S'
		self.aaabr2letter["THR"]='T'
		self.aaabr2letter["SEC"]='U'
		self.aaabr2letter["VAL"]='V'
		self.aaabr2letter["TRP"]='W'
		self.aaabr2letter["TYR"]='Y'

		self.rotameres=[]
		self.ff=ff
		self.defaultRotamere()


	def parseGRO(self, grofilename) :
		grofile = open(grofilename)
		lines = grofile.readlines()
		#title=lines[0]
		#nbatoms=int(lines[1])

		for line in lines[2:len(lines)-1] :
			resnumstr=str(line[0:5])
			resnamestr=str(line[5:10])
			atomnamestr=str(line[10:15])
			atomidstr=str(line[15:20])
			xstr=str(line[20:28])
			ystr=str(line[28:36])
			zstr=str(line[36:44])
			print resnumstr+','+resnamestr+','+atomnamestr+','+atomidstr+','+xstr+','+ystr+','+zstr
			atomid	= int(atomidstr)
			self.atoms[atomid].pos=[float(xstr)*10.0,float(ystr)*10.0,float(zstr)*10.0]

	def parseNameToType(self, n2tfilename) :
		n2tfile = open(n2tfilename)
		lines = n2tfile.readlines()
		type2name={}
		for line in lines :
			li=line.strip()
			if not li.startswith("#"):
				values=li.split('\t')
				if len (values)==2 :
					ntype=str(values[1])
					nname=str(values[0])
					type2name[ntype]=nname
		n2tfile.close()
		return	 type2name

	def defaultRotamere(self) :
		self.rotameres.append(('CA','C'))
		self.rotameres.append(('C','CA'))
		self.rotameres.append(('CB','CA'))
		self.rotameres.append(('CA','CB'))
		self.rotameres.append(('CA','N'))
		self.rotameres.append(('N','CA'))

	def parseRotamere(self, rotamerefilename) :
		rotamerefile = open(rotamerefilename)
		lines = rotamerefile.readlines()
		self.rotameres=[]
		for line in lines :
			li=line.strip()
			if not li.startswith("#"):
				values=li.split('\t')
				if len (values)==2 :
					a1=str(values[1])
					a2=str(values[0])
					self.rotameres.append((a1,a2))
					self.rotameres.append((a2,a1))

	def applyForceField(self, fffilename) :
		ff=ForceField()
		ff.parseForceField(fffilename)
		type2name=self.parseNameToType('ff'+os.path.sep+self.ff+os.path.sep+'nametotype.txt')
		for atom in self.atoms.values():
			pp=ff.particleproperties[type2name[atom.ptype]]
			atom.charge=pp.charge
			atom.radius=pp.radius
			atom.epsilon=pp.epsilon
			atom.mass=pp.mass
			atom.imp=pp.imp
			atom.hydro=pp.hydro

	def buildFromPDB(self, pdbfilename) :
		pdbfilelocal=os.path.basename(pdbfilename)
		pdb2gmxcmd = gromacs.tools.Pdb2gmx('ignh', f=[pdbfilename], p= pdbfilelocal +'.top', o=pdbfilelocal +'.gro',water='None', ff='amber94')
		pdb2gmxcmd.run()
		self.buildFromTopology(pdbfilelocal +'.top')
		self.parseGRO(pdbfilelocal +'.gro')

		self.applyForceField('ff'+os.path.sep+self.ff+os.path.sep+'ff.txt')



	def buildFromTopology(self, topofilename) :
		itp = fileformats.ITP(topofilename)
		databatoms=itp.sections['header'].sections['moleculetype'].sections['atoms'].data

		for particle in databatoms :
			(pid, ptype, resid, resname, patom, pidseq, pcharge, pmass, dummy) = particle
			p=Particle(pid, ptype, patom, resid)
			self.atoms[pid]=p

		databonds=itp.sections['header'].sections['moleculetype'].sections['bonds'].data
		for bond in databonds :
			(id1, id2, dummy1,dummy2) = bond
			if  id1<id2 :
				b=Bond(id1,id2)
				b.rotamere=self.RotamereRules(b)
				self.bonds[(id1,id2)]=b
				self.atoms[id1].addBond(b)
				self.atoms[id2].addBond(b)

			else:
				b=Bond(id2,id1)
				b.rotamere=self.RotamereRules(b)
				self.bonds[(id2,id1)]=b
				self.atoms[id1].addBond(b)
				self.atoms[id2].addBond(b)



	def RotamereRules(self, bond):
		p1=self.atoms[bond.id1]
		p2=self.atoms[bond.id2]

		if (p1.patom, p2.patom) in self.rotameres and p1.pparentid==p2.pparentid :
			return True

		return False;

	def findNext(self, atomfilter):
		for atom in self.atoms.values() :
			if atom.visited == False :
				if atomfilter==None:
					return atom
				elif atom.patom in atomfilter :
					return atom
		return None

	def AAABRtoAALetter(self, aaabr) :
		return self.aaabr2letter[aaabr]

	def __str__(self) :
		retstr=""
		for atom in self.atoms.values() :
			retstr=retstr+ str(atom)+'\n'

		for bond in self.bonds.values() :
			retstr=retstr+ str(bond)+'\n'
		return retstr