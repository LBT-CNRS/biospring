import copy

class RigidBody() :
	def __init__(self, bodyid) :
		self.bodyid=bodyid
		self.particles={}
	def addParticle(self, atom, ghost) :
		newatom=copy.deepcopy(atom)
		newatom.ghost=ghost
		newatom.pparentid=self.bodyid
		self.particles[newatom.pid]=newatom
	def __str__(self) :
		strrep= str(self.bodyid) + ':'
		for particle in self.particles.values() :
			strrep=strrep+ '\t'
			strrep=strrep+ str(particle)+'\n'
		return strrep