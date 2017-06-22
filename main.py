#import matlab.engine
from sympy import *

class CRN:
	def __init__(self, s, r, ip):
		self.species = s
		self.reactions = r
		self.initialPopulations = ip

	def __repr__(self):
		return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"

	def __str__(self):
		return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"


class Reaction:
	def __init__(self, r, p, ra):
		self.reactants = r
		self.products = p
		self.reactionrate = ra

	def __repr__(self):
		return "" + ' + '.join(["".join(x) for x in self.reactants]) + " ->{" + str(self.reactionrate) + "} " + ' + '.join(["".join(y) for y in self.products])

	def __str__(self):
		return "" + ' + '.join(["".join(x) for x in self.reactants]) + " ->{" + str(self.reactionrate) + "} " + ' + '.join(["".join(y) for y in self.products])

def exampleCRN():
	reaction1 = Reaction(['A'],['2B','A'], 1.0)
	reaction2 = Reaction(['2A','2C'],['C'], 1.0)
	reaction3 = Reaction(['C'],[], 1.0)
	return CRN(['A','B','C'], [reaction1,reaction2,reaction3], [5,2,1])


#def propensities(reations):
#	for reaction in reactions:


if __name__ == "__main__":
	crn = exampleCRN()
	test = symbols(['a','b','c']) * symbols('2')
	print test
	#reactants = crn.reactions.reactants
	#reactants = [ x.reactants for x in crn.reactions]


	#eng = matlab.engine.start_matlab()
	#ret = eng.symbolicLNA(crn.species,[ x.reactants for x in crn.reactions], [ x.products for x in crn.reactions], [ x.reactionrate for x in crn.reactions])
	

	#inspecies, inreactants, inproducts, inrates