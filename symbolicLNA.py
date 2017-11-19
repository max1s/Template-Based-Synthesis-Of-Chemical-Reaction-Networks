#import matlab.engine
from sympy import *
from sympy import Matrix
import ipdb
from sympy import init_printing
import iSATParser


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

class stoichPair:
	def __init__(self, x, y):
		self.variable = symbols(x)
		self.coefficient = symbols(y)



class ReactionSketch:
	def __init__(self, r, opr, p, opp, coeff, ra, isop):
		self.reactants = r
		self.products = p
		self.optionalReactants = opr
		self.optionalProducts = opp
		self.coefficients = coeff
		self.reactionrate = ra
		self.isOptional = isop


	def __repr__(self):
		return "" + ' + '.join(["".join(x) for x in self.reactants]) + " ->{" + str(self.reactionrate) + "} " + ' + '.join(["".join(y) for y in self.products])

	def __str__(self):
		return "" + ' + '.join(["".join(x) for x in self.reactants]) + " ->{" + str(self.reactionrate) + "} " + ' + '.join(["".join(y) for y in self.products])

def AMParametricExample():
	#reaction1 = ReactionSketch([X], [B], [X, X], )
	pass
class OptionalReaction:
	def __init__(self, r, p, ra):
		self.reactants = r
		self.products = p
		self.reactionrate = ra

	def __repr__(self):
		return "" + ' + '.join(["".join(x) for x in self.reactants]) + " ->{" + str(self.reactionrate) + "} " + ' + '.join(["".join(y) for y in self.products])

	def __str__(self):
		return "" + ' + '.join(["".join(x) for x in self.reactants]) + " ->{" + str(self.reactionrate) + "} " + ' + '.join(["".join(y) for y in self.products])

class CRNSketch:
	def __init__(self, cs, os, r, opr, params):
		self.species = s
		self.optionalSpecies = os
		self.reactions = r
		self.optionalReactions = os
		self.parameters = params


	def __repr__(self):
		return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"

	def __str__(self):
		return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"


def exampleCRN():
	reaction1 = Reaction(['A','x_1'],['x_1','x_1'], 'k_1')
	reaction2 = Reaction(['x_1','x_2'],['x_2','x_2'], 'k_2')
	reaction3 = Reaction(['x_2'],['B'], 'k_3')
	return CRN(['A','B','x_1', 'x_2'], [reaction1,reaction2,reaction3], [5,2,1])

def AMExample():
	reaction1 = Reaction(['X','Y'], ['X','B'], 'k_1')
	reaction2 = Reaction(['Y','X'], ['Y','B'], 'k_2')
	reaction3 = Reaction(['X','B'], ['X','X'], 'k_3')
	reaction4 = Reaction(['Y', 'B'], ['Y','Y'], 'k_4')
	return CRN(['X', 'Y', 'B'], [reaction1,reaction2,reaction3,reaction4], [5,3,1])


def exampleParametricCRN():
	reaction1 = ReactionSketch(['A'],['B'], ['K'], [], [2,1], 'k_1', 0)
	reaction2 = ReactionSketch(['B'],['B'], ['K'], ['K'], [2,1], 'k_1', 0)
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

#alpha1=k1*SF*(lambda1A*A + lambda1B*B)*(c10+c11*K);
#alpha2=k2*SF*(c50 + c51*(lambda2A*A + lambda2B*B) )*(c31*K + c32*(K^2));
#alpha3=k3*SF;

def parametricPropensity(paramCRN):
	propensities = []
	for (reaction, i) in zip(paramCRN.reactions, range(len(paramCRN.species))):
		reactionRate = symbols(str(reaction.reactionrate))
		for reactant in reaction.reactants:
			propensity *= symbols(reactant)
		for optionalReactant in reaction.optionalReactant:
			propensity *= symbols('(' + optionalReactant +'A' + ' + ' + optionalReactant + 'B' + ')')

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


def jacobian(flow, species):
	return flow.jacobian(species)


def g(propensities, reactionChange):
	G = zeros(len(reactionChange.row(0)), len(reactionChange.row(0)))
	for i in range(len(propensities)):
		G += transpose(reactionChange.row(i)) * reactionChange.row(i) * propensities[i]
	return G

def generateCovarianceMatrix(speciesVector):
	mat = eye(len(speciesVector))
	for (m, i) in zip(speciesVector, range(len(speciesVector))):
		for (n, j) in zip(speciesVector, range(len(speciesVector))):
			if (m == n):
				mat[i,j] = 'cov' + m 
			else:
				mat[i, j] = 'cov' + n + m

	for x in range(len(speciesVector)):
		for y in range(len(speciesVector)):
			mat[x,y] = mat[y,x]
	#pprint(mat)
	return mat


	# G = reactionChange[1].transpose 	
	#G=v1'*v1*alpha1
	#G=G+(v2'*v2*alpha2)
	#G=G+(v3'*v3*alpha3)

   #C=[covX covXY ;
   #		covXY covY ]

   #dCovdt=J*C+C*(J')+G

if __name__ == "__main__":
	exampleParametricCRN()
	init_printing() 
	crn = exampleCRN()
	props = propensity(crn.reactions)
	nrc = netReactionChange(crn.species, crn.reactions)
	#props = propensity(crn.reactions)
	flow = flowFunction(props,nrc)
	#A = Matrix([crn.species])
	J = jacobian(flow, Matrix([crn.species]))
	#print covar
	G = g(Matrix(props), Matrix(nrc))
	C = generateCovarianceMatrix(crn.species)

	dCovdt = J*C + C*transpose(J)  + G

	#iSATParser.constructiSATFile()
	print crn
	print "propensities:"
	print props

	pprint(dCovdt)
	#reactants = crn.reactions.reactants
	#reactants = [ x.reactants for x in crn.reactions]
	

	#eng = matlab.engine.start_matlab()
	#ret = eng.symbolicLNA(crn.species,[ x.reactants for x in crn.reactions], [ x.products for x in crn.reactions], [ x.reactionrate for x in crn.reactions])
	

	#inspecies, inreactants, inproducts, inrates