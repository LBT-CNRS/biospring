import os

from spring import Spring
from string import Template



class SpringNetwork():
	def __init__(self) :
		self.particles=[]
		self.particleindexes={}
		self.springs=[]

	def buildFromRigidBodyStructure(self, rbs) :
		index=1
		for rb in rbs.bodies.values() :
			for particle in rb.particles.values() :
				self.particleindexes[(rb.bodyid, particle.pid)]=index
				self.particles.append(particle)
				index=index+1

		for rb in rbs.bodies.values() :
			self.computeSpringFromRigidBody(rb)

		for joint in rbs.joints :
			self.computeSpringFromJoint(joint)

	def computeSpringFromRigidBody(self, rigidbody):
		particles=rigidbody.particles.values()
		for i in range(0,len(particles)) :
			for j in range(i+1, len(particles)) :
				stiffness=100.0
				self.springs.append(Spring(particles[i],particles[j], stiffness))

	def computeSpringFromJoint(self, joint):
		body1=joint.b1
		body2=joint.b2
		p1ghost=None
		p2ghost=None
		for particle in body1.particles.values() :
			if particle.ghost and joint.p2.pid==particle.pid:
				p1ghost=particle


		for particle in body2.particles.values() :
			if particle.ghost and joint.p1.pid==particle.pid:
				p2ghost=particle

		print joint.p1
		print p2ghost

		print joint.p2
		print p1ghost


		stiffness=1000.0
		if not p1ghost==None and not p2ghost==None :
			self.springs.append(Spring(joint.p1,p2ghost, stiffness))
			self.springs.append(Spring(joint.p2,p1ghost, stiffness))
		else:
			print "Error connecting rigid body with spring!"

	def __str__(self) :
		retstr=""
		for spring in self.springs :
			retstr=retstr+str(spring)+'\n'
		return retstr

	def toPDB(self, pdbfilename) :
		pdbfile = open(pdbfilename, "w")
		for p in self.particles :
			index=self.particleindexes[(p.pparentid, p.pid)]
			line="%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n" % ("ATOM", index, p.patom, " ", p.pparentid, "A", p.pparentid," ", p.pos[0], p.pos[1], p.pos[2], 0.0, 0.0, p.patom[0],"  " )
			pdbfile.writelines(line)
		for s in self.springs :
			id1=self.particleindexes[(s.p1.pparentid, s.p1.pid)]
			id2=self.particleindexes[(s.p2.pparentid, s.p2.pid)]
			line="%-6s%5d%5d \n" % ("CONECT",  id1, id2)
			pdbfile.writelines(line)
		pdbfile.close()

	def toCDL(self, cdlfilename) :
		cdlfile = open(cdlfilename, "w")
		cdltemplatefile=open("template"+os.path.sep+"biospring.template.cdl", "r")
		lines=cdltemplatefile.readlines()
		strtemplate=""
		for line in lines :
			strtemplate=strtemplate+line
		template = Template(strtemplate)
		templatedict=self.buildTemplateDictionnary()

		resultstring=template.substitute(templatedict)
		cdlfile.write(resultstring)
		#for p in self.particles :
		#for s in self.springs :
		cdlfile.close()
		cdltemplatefile.close()

	def buildTemplateDictionnary(self) :
		templatedict={}
		templatedict["particlesnumber"]=""
		templatedict["springsnumber"]=""
		templatedict["particleids"]=""
		templatedict["coordinates"]=""
		templatedict["charges"]=""
		templatedict["radii"]=""
		templatedict["epsilon"]=""
		templatedict["mass"]=""
		templatedict["surfaceaccessibility"]=""
		templatedict["hydrophobicityscale"]=""
		templatedict["particlenames"]=""
		templatedict["resnames"]=""
		templatedict["resids"]=""
		templatedict["chainnames"]=""
		templatedict["dynamicstate"]=""
		templatedict["springs"]=""
		templatedict["springsstiffness"]=""
		templatedict["springsequilibrium"]=""

		templatedict["particlesnumber"]=str(len(self.particles))
		templatedict["springsnumber"]=str(len(self.springs))

		index=0
		for particle in self.particles :
			termcharacter=',\n'
			if index== (len(self.particles)-1) :
				termcharacter=''

			templatedict["particleids"]=templatedict["particleids"]+str(index)+termcharacter
			templatedict["coordinates"]=templatedict["coordinates"]+str(particle.pos[0])+','+str(particle.pos[1])+','+str(particle.pos[2])+termcharacter
			templatedict["charges"]=templatedict["charges"]+str(particle.charge)+termcharacter
			templatedict["radii"]=templatedict["radii"]+str(particle.radius)+termcharacter
			templatedict["epsilon"]=templatedict["epsilon"]+str(particle.epsilon)+termcharacter
			templatedict["mass"]=templatedict["mass"]+str(particle.mass)+termcharacter
			templatedict["surfaceaccessibility"]=templatedict["surfaceaccessibility"]+str(0.0)+termcharacter
			templatedict["hydrophobicityscale"]=templatedict["hydrophobicityscale"]+str(particle.hydro)+termcharacter
			templatedict["particlenames"]=templatedict["particlenames"]+"\""+str(particle.patom)+"\""+termcharacter
			templatedict["resnames"]=templatedict["resnames"]+"\""+str(particle.pparentid)+"\""+termcharacter
			templatedict["resids"]=templatedict["resids"]+str(particle.pparentid)+termcharacter
			templatedict["chainnames"]=templatedict["chainnames"]+"\""+str('A')+"\""+termcharacter
			templatedict["dynamicstate"]=templatedict["dynamicstate"]+str(1)+termcharacter
			index=index+1

		index=0
		for spring in self.springs :
			termcharacter=',\n'
			if index == (len(self.springs)-1) :
				termcharacter=''
			pindex1=self.particleindexes[(spring.p1.pparentid, spring.p1.pid)]-1
			pindex2=self.particleindexes[(spring.p2.pparentid, spring.p2.pid)]-1
			templatedict["springs"]=templatedict["springs"]+str(pindex1)+','+str(pindex2)+termcharacter
			templatedict["springsstiffness"]=templatedict["springsstiffness"]+str(spring.stiffness)+termcharacter
			templatedict["springsequilibrium"]=templatedict["springsequilibrium"]+str(spring.equilibrium)+termcharacter
			index=index+1

		return templatedict