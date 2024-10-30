
class Bond () :
	def __init__(self, id1,id2):
		self.id1=id1
		self.id2=id2
		self.rotamere=False

	def __str__(self) :
		strrep= '('+ 'id1 :' + str(self.id1) + ','+ 'id2 :' + str(self.id2) + ','+ 'rotamere :' + str(self.rotamere)+')'
		return strrep
