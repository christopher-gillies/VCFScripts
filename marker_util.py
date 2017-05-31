import re
class Marker:

	def __init__(self,location):
		match = re.match("(chr)?([^:]+):(\d+):([ACTG]+):([ACTG]+)",location)
		if match is None:
			raise ValueError("{0} not formatted correctly. Should be chr:pos:ref:alt".format(location))
		self.chrom = match.group(2)
		self.pos = int(match.group(3))
		self.ref = match.group(4)
		self.alt = match.group(5)
	
	def key(self):
		return self.__str__()
		
	@staticmethod
	def format(chrom,pos,ref,alt):
		return "{0}:{1}:{2}:{3}".format(chrom,pos,ref,alt)	
	
	def __str__(self):
		return Marker.format(self.chrom,self.pos,self.ref,self.alt)
	
	def __repr__(self):
		return self.__str__()