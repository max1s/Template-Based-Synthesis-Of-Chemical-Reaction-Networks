#import matlab.engine
from sympy import *
from sympy import Matrix

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
	reaction1 = Reaction(['A','x_1'],['x_1','x_1'], 'k_1')
	reaction2 = Reaction(['x_1','x_2'],['x_2','x_2'], 'k_2')
	reaction3 = Reaction(['x_2'],['B'], 'k_3')
	return CRN(['A','B','x_1', 'x_2'], [reaction1,reaction2,reaction3], [5,2,1])


#def propensities(reations):
#	for reaction in reactions:

def propensity(reactions):
	propensities = []
	for reaction in reactions:
		propensity = symbols(str(reaction.reactionrate))
		for reactant in reaction.reactants:
			propensity *= symbols(reactant)
		propensities.append(propensity)
	return propensities

def netReactionChange(species, reactions):
	reactionChange = [] #ReferenceFrame('N')
	for reaction in reactions:
		netChange = []
		for specie in species:
			speciesChange = 0
			for reactant in reaction.reactants:
				if specie == reactant:
					speciesChange -= 1 
			for product in reaction.products:
				if specie == product:
					speciesChange += 1 
			netChange.append(speciesChange)
		reactionChange.append(Matrix(1, len(netChange), netChange))
	return reactionChange

#def flowFunction(propensities, reactionChange):
#	flow = []
#	for n in range(0, len(propensities)):
		#if (sum(reactionChange[n]) != 0):
#		flow.append(reactionChange[n]*propensities[n]) 
#	return Matrix(1, len(flow), flow)

def flowFunction(propensities, reactionChange):
	summationFunction = reactionChange[0]*propensities[0]
	for n in range(1, len(propensities)):
		summationFunction += reactionChange[n]*propensities[n]
	return summationFunction


def calculateCovariances(flow, species):
	print flow.shape
	print species.shape
	print flow.jacobian(species)

if __name__ == "__main__":
	crn = exampleCRN()
	print propensity(crn.reactions)
	print netReactionChange(crn.species, crn.reactions)
	#print propensity(crn.reactions)
	#print netReactionChange(crn.species, crn.reactions)
	calculateCovariances(flowFunction(propensity(crn.reactions), netReactionChange(crn.species, crn.reactions)), Matrix(1, len(crn.species), crn.species))
	#reactants = crn.reactions.reactants
	#reactants = [ x.reactants for x in crn.reactions]


	#eng = matlab.engine.start_matlab()
	#ret = eng.symbolicLNA(crn.species,[ x.reactants for x in crn.reactions], [ x.products for x in crn.reactions], [ x.reactionrate for x in crn.reactions])
	

	#inspecies, inreactants, inproducts, inrates