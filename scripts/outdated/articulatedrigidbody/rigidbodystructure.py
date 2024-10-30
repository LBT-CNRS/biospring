from rigidbody import RigidBody
from joint import Joint

class RigidBodyStructure() :
	def __init__(self) :
		self.bodies={}
		self.joints=[]


	def buildFromStructure(self,structure) :
		atomfilter=['CA', 'C', 'N', 'O','H']
		nextatom=structure.findNext(atomfilter)
		bodyid=0
		while  nextatom!=None :
			self.bodies[bodyid]=RigidBody(bodyid)
			self.RecursiveRigidBodyId(nextatom,structure,bodyid)
			bodyid=bodyid+1
			nextatom=structure.findNext(atomfilter)

		atomfilter=None
		nextatom=structure.findNext(atomfilter)
		while  nextatom!=None :
			self.bodies[bodyid]=RigidBody(bodyid)
			self.RecursiveRigidBodyId(nextatom,structure, bodyid)
			bodyid=bodyid+1
			nextatom=structure.findNext(atomfilter)

		for bond in structure.bonds.values() :
			if bond.rotamere :
				atom1 = structure.atoms[bond.id1]
				atom2 = structure.atoms[bond.id2]
				bodies1= atom1.bodyids
				bodies2= atom2.bodyids
				interbodies=list(set(bodies1) & set(bodies2))
				body1=self.bodies[interbodies[0]]
				body2=self.bodies[interbodies[1]]
				interparticles=list(set(body1.particles.keys()) & set(body2.particles.keys()))

				if body1.particles[interparticles[0]].ghost==False :
					p1=body1.particles[interparticles[0]]
				else :
					p1=body1.particles[interparticles[1]]

				if body2.particles[interparticles[0]].ghost==False :
					p2=body2.particles[interparticles[0]]
				else:
					p2=body2.particles[interparticles[1]]

				self.joints.append(Joint(self.bodies[interbodies[0]],p1, self.bodies[interbodies[1]], p2))

	def RecursiveRigidBodyId(self, atom,structure,bodyid):
		if not atom.visited :
			atom.setVisited()
			atom.addBodyId(bodyid)
			self.bodies[bodyid].addParticle(atom, False)
			atombonds=atom.bonds.values()
			atomids=atom.bonds.keys()
			for i in range(len(atomids)) :
				if not atombonds[i].rotamere :
					self.RecursiveRigidBodyId(structure.atoms[atomids[i]],structure, bodyid)
				else :
					structure.atoms[atomids[i]].addBodyId(bodyid)
					self.bodies[bodyid].addParticle(structure.atoms[atomids[i]], True)


	def __str__(self) :
		retstr=""
		for body in self.bodies.values() :
			retstr=retstr+str(body)+'\n'

		for joint in self.joints :
			retstr=retstr+str(joint)+'\n'
		return retstr