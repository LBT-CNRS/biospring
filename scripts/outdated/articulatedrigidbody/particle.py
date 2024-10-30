import math
from particleproperties import ParticleProperties

class Particle(ParticleProperties) :
	def __init__(self, pid, ptype, patom, pparentid):
		ParticleProperties.__init__(self, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		self.pid = pid
		self.ptype =ptype
		self.patom =patom
		self.pparentid=pparentid
		self.bonds={}
		self.pos=[0.0,0.0,0.0]

		self.bodyids=[]
		self.visited=False
		self.ghost=False

	def addBond(self,bond) :
		if bond.id1==self.pid:
			self.bonds[bond.id2]=bond
		else :
			self.bonds[bond.id1]=bond
	def addBodyId(self, bodyid) :
		self.bodyids.append(bodyid)

	def setVisited(self):
		self.visited=True

	def __str__(self) :
		strrep= '('+ 'pid :' + str(self.pid) + ','+' ghost :'+str(self.ghost)+ ','+'ptype :' + str(self.ptype) + ','+ 'patom :' + str(self.patom) + ','+' pparentid :' + str(self.pparentid) + ','+ ' pos : '+str(self.pos)+ ','+ ' charge : '+str(self.charge)+ ','+ ' radius : '+str(self.radius) + ','+ ' epsilon : '+str(self.epsilon) + ' mass : '+str(self.mass) + ','+ ' imp : '+str(self.imp) + ','+ ' hydro : '+str(self.hydro)+')'
		strrep= strrep +"=>"
		for bond in self.bonds.values() :
			if bond.id1==self.pid :
				neighbour=bond.id2
			else :
				neighbour=bond.id1
			strrep=strrep+' '+str(neighbour)
			if bond.rotamere :
				strrep=strrep+'*'
		return strrep

	@classmethod
	def distance(self, p1,p2) :
		dist = math.sqrt( (p2.pos[0] - p1.pos[0])**2 + (p2.pos[1] - p1.pos[1])**2+ (p2.pos[2] - p1.pos[2])**2)
		return dist