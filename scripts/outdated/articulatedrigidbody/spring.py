from particle import Particle

class Spring():
	def __init__(self,p1,p2, stiffness) :
		self.p1=p1
		self.p2=p2
		self.equilibrium=Particle.distance(p1, p2)
		self.stiffness=stiffness
	def __str__(self) :
		return 'Spring('+str(self.p1.pid)+','+str(self.p2.pid)+','+str(self.equilibrium)+','+str(self.stiffness)+','+')'