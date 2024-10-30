from particleproperties import ParticleProperties

class ForceField() :
	def __init__(self):
		self.particleproperties={}
	def parseForceField(self, filename) :
		fffile = open(filename)
		lines = fffile.readlines()
		for line in lines :
			li=line.strip()
			if not li.startswith("#"):
				values=li.split('\t')
				name=str(values[0])
				charge=float(values[1])
				radius=float(values[2])
				epsilon =float(values[3])
				mass = float(values[4])
				imp = float(values[5])
				hydro = 0.0
				pp=ParticleProperties(charge,radius, epsilon, mass,imp, hydro)
				self.particleproperties[name]=pp
