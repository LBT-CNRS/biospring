class Joint():
	def __init__(self, b1, p1, b2, p2) :
		self.b1=b1
		self.p1=p1

		self.b2=b2
		self.p2=p2

	def __str__(self) :
		strrep= "(Body id1 " + str(self.b1.bodyid) + " Particle id1 " + str(self.p1.pid) + ","+ "Body id2 " + str(self.b2.bodyid) + " Particle id2 " + str(self.p2.pid)+")"
		return strrep