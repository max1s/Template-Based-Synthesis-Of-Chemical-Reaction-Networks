# import matlab.engine
from sympy import *
from sympy import Matrix
import ipdb
from sympy import init_printing
import iSATParser
import itertools


class CRN:
    def __init__(self, s, r, ip):
        self.species = s
        self.reactions = r
        self.initialPopulations = ip

    def __repr__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"

    def __str__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"


class RateConstant:
    def __init__(self, name, min, max):
       self.name = name
       self.min = min
       self.max = max

    def __repr__(self):
         return self.name

    def __str__(self):
         return "Rate constant %s <= %s <= %s"  % (self.min, self.name, self.max)


class Reaction:
    def __init__(self, r, p, ra):
        self.reactants = r
        self.products = p
        self.reactionrate = ra


class LambdaChoice:
    def __init__(self, species, choiceNumber):
        self.choiceNo = choiceNumber
        self.species = species
        self.lambdas = [symbols("lam" + str(sp) + str(choiceNumber)) for sp in species]

    def constructChoice(self):
        return "(" + '+'.join([str(sp) + '*' + str(l) for sp, l in zip(self.species, self.lambdas)]) + ")"


    def contains(self, var):
        x = ''
        for lam in self.lambdas:
            if str(var) in str(lam):
                x += str(lam)
        return x

class Choice:
    def __init__(self, choiceName, minNumber, maxNumber):
        self.choiceName = choiceName
        self.choice = ['c' + str(choiceName) + str(x) for x in range(minNumber, maxNumber)]
        self.minNumber = minNumber
        self.maxNumber = maxNumber

    def constructChoice(self):
        chain = ""
        if (self.minNumber == 0):
            chain += self.choice[0]
        if (self.maxNumber > 0):
            chain += " + " + self.choice[1]
        if (self.maxNumber > 1):
            chain += (" + ").join([str(choice) + "*" + str(self.choiceName) + "^" + str(x) for choice, x in
                                   zip(self.choice, range(2, self.maxNumber))])
        return chain


class Species:
    def __init__(self, species, coefficient):
        self.species = species
        self.coefficient = coefficient

    def specRep(self):
        if isinstance(self.species, LambdaChoice):
            return str(self.coefficient) + "*" + self.species.constructChoice()
        elif isinstance(self.species, Choice):
            return str(self.coefficient) + "*" + " ( " + str(self.species.constructChoice()) + " ) "
        else:
            return str(self.coefficient) + "*" + str(self.species)

    def constructPropensity(self):
        if (isinstance(self.coefficient, int)):
            return self.specRep()
        else:
            raise NotImplementedError


class ReactionSketch:
    def __init__(self, r, opr, p, opp, ra, isop):
        self.reactants = r
        self.products = p
        self.lambdaReactants = opr
        self.lambdaProducts = opp
        self.reactionrate = ra
        self.isOptional = isop

    def __repr__(self):
        return "" + ' + '.join(["".join(x) for x in self.reactants]) + " ->{" + str(
            self.reactionrate) + "} " + ' + '.join(["".join(y) for y in self.products])

    def __str__(self):
        return "" + ' + '.join(["".join(x) for x in self.reactants]) + " ->{" + str(
            self.reactionrate) + "} " + ' + '.join(["".join(y) for y in self.products])



class OptionalReaction:
    def __init__(self, r, p, ra):
        self.reactants = r
        self.products = p
        self.reactionrate = ra

    def __repr__(self):
        return "" + ' + '.join(["".join(x) for x in self.reactants]) + " ->{" + str(
            self.reactionrate) + "} " + ' + '.join(["".join(y) for y in self.products])

    def __str__(self):
        return "" + ' + '.join(["".join(x) for x in self.reactants]) + " ->{" + str(
            self.reactionrate) + "} " + ' + '.join(["".join(y) for y in self.products])


class CRNSketch:
    def __init__(self, cs, r, opr):
        self.species = cs
        self.reactions = r
        self.optionalReactions = opr

    def __repr__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"

    def __str__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"

    def getSpecies(self):
        x = []
        for y in self.reactions:
            for react in y.reactants:
                if react not in x:
                    x.append(react.specRep())
            for react in y.products:
                if react not in x:
                    x.append(react.specRep())
        return x



def parametricPropensity(paramCRN):
    propensities = []
    for reaction in paramCRN.reactions:
        propensity = symbols(str(reaction.reactionrate.name))
        for reactant in reaction.reactants:
            propensity *= sympify(reactant.constructPropensity())
        propensities.append(propensity)
    return propensities


def parametricNetReactionChange(crn):
    reactionChange = []
    for reaction in crn.reactions:
        netChange = [''] * len(crn.species)
        for reactant in reaction.reactants:
            add(crn.species, netChange, reactant, '-')
        for product in reaction.products:
            if len(netChange) is 0:
                #netChange.append(product.specRep())
                add(crn.species, netChange, reactant, '+')
            else:
                #netChange.append("+" + product.specRep())
                add(crn.species, netChange, reactant, '+')
        netChange = [0 if n is '' else n for n in netChange]
        reactionChange.append(sympify(netChange))
    return reactionChange


def add(species, list,  fragment, part):
    for spec, i in zip(species, range(len(species))):
        if str(spec) in fragment.specRep():
            if "lam" in fragment.specRep():
                list[i] += part + fragment.species.contains(spec) + "*" + str(spec)
            else:
                list[i] += part + fragment.specRep()
    return list


def parametricFlow(propensities, reactionChange):
    container = [reactionChange[0] * propensities[0]]
    for n in range(1, len(propensities)):
        container.append(reactionChange[n] * propensities[n])

    return container



def parametricG(propensities, reactionChange):
    G = zeros(max(len(reactionChange.col(0)), len(reactionChange.row(0))),
              max(len(reactionChange.col(0)), len(reactionChange.row(0))))
    for i in range(len(propensities)):
        G += transpose(reactionChange.row(i)) * reactionChange.row(i) * propensities[i]
    return G


def generateCovarianceMatrix(speciesVector):
    mat = eye(len(speciesVector))
    for (m, i) in zip(speciesVector, range(len(speciesVector))):
        for (n, j) in zip(speciesVector, range(len(speciesVector))):
            if (m == n):
                mat[i, j] = 'cov' + m
            else:
                mat[i, j] = 'cov' + n + m

    for x in range(len(speciesVector)):
        for y in range(len(speciesVector)):
            mat[x, y] = mat[y, x]
    # pprint(mat)
    return mat

def generateAllTokens(crn):
    t = [x.free_symbols for x in sympify(crn.getSpecies())]
    a = reduce(lambda x, y : x |y, t)

    print a
def exampleParametricCRN():
    X = symbols('X')
    Y = symbols('Y')
    B = symbols('B')

    reaction1 = Reaction([Species(LambdaChoice([X, Y], 1), 1), Species(Y, 1)], [Species(X, 1), Species(B, 1)], RateConstant('k_1', 1, 2))
    reaction2 = Reaction([Species(LambdaChoice([X, Y], 2), 1), Species(Choice(X, 0, 2), 2)], [Species(Y, 1), Species(B, 1)], RateConstant('k_1', 1, 2))
    reaction3 = Reaction([Species(X, 1), Species(B, 1)], [Species(X, 1), Species(X, 1)], RateConstant('k_1', 1, 2))
    reaction4 = Reaction([Species(X, 1), Species(B, 1)], [Species(X, 1), Species(X, 1)], RateConstant('k_1', 1, 2))

    crn = CRNSketch([X, Y, B], [reaction1, reaction2, reaction3], [reaction4])

    prp = (parametricPropensity(crn))
    nrc = (parametricNetReactionChange(crn))
    dSpeciesdt = parametricFlow(prp, Matrix(nrc))
    J = Matrix(dSpeciesdt).jacobian([X, Y, B])
    G = parametricG(Matrix(prp), Matrix(nrc))

    C = generateCovarianceMatrix(['X', 'Y', 'B'])

    dCovdt = J * C + C * transpose(J)
    #pprint(dCovdt)
    generateAllTokens(crn)
    # generateAllTokens()




if __name__ == "__main__":
    exampleParametricCRN()


    # def exampleCRN():
    #     reaction1 = Reaction(['A', 'x_1'], ['x_1', 'x_1'], 'k_1')
    #     reaction2 = Reaction(['x_1', 'x_2'], ['x_2', 'x_2'], 'k_2')
    #     reaction3 = Reaction(['x_2'], ['B'], 'k_3')
    #     return CRN(['A', 'B', 'x_1', 'x_2'], [reaction1, reaction2, reaction3], [5, 2, 1])
    #
    #
    # def AMExample():
    #     reaction1 = Reaction(['X', 'Y'], ['X', 'B'], 'k_1')
    #     reaction2 = Reaction(['Y', 'X'], ['Y', 'B'], 'k_2')
    #     reaction3 = Reaction(['X', 'B'], ['X', 'X'], 'k_3')
    #     reaction4 = Reaction(['Y', 'B'], ['Y', 'Y'], 'k_4')
    #     return CRN(['X', 'Y', 'B'], [reaction1, reaction2, reaction3, reaction4], [5, 3, 1])

        # alpha1=k1*SF*(lambda1A*A + lambda1B*B)*(c10+c11*K);
    # alpha2=k2*SF*(c50 + c51*(lambda2A*A + lambda2B*B) )*(c31*K + c32*(K^2));
    # alpha3=k3*SF;

    # def flowFunction(propensities, reactionChange):
    #	flow = []
    #	for n in range(0, len(propensities)):
    # if (sum(reactionChange[n]) != 0):
    #		flow.append(reactionChange[n]*propensities[n])
    #	return Matrix(1, len(flow), flow)

    # def g(propensities, reactionChange):
    #     G = zeros(max(len(reactionChange.col(0)), len(reactionChange.row(0))),
    #               max(len(reactionChange.col(0)), len(reactionChange.row(0))))
    #     for i in range(len(propensities)):
    #         G += transpose(reactionChange.row(i)) * reactionChange.row(i) * propensities[i]
    #     return G
    # G = reactionChange[1].transpose
    # G=v1'*v1*alpha1
    # G=G+(v2'*v2*alpha2)
    # G=G+(v3'*v3*alpha3)

    # C=[covX covXY ;
    #		covXY covY ]

    # dCovdt=J*C+C*(J')+G


    # exampleParametricCRN()
    # init_printing()
    # crn = exampleCRN()
    # props = propensity(crn.reactions)
    # nrc = netReactionChange(crn.species, crn.reactions)
    # props = propensity(crn.reactions)

    # flow = flowFunction(props,nrc)
    # print(flow)
    # print(type(flow))
    # #A = Matrix([crn.species])
    # J = jacobian(flow, Matrix([crn.species]))
    # #print covar
    # G = g(Matrix(props), Matrix(nrc))
    # C = generateCovarianceMatrix(crn.species)
    #
    # dCovdt = J*C + C*transpose(J)  + G
    #
    # #iSATParser.constructiSATFile()
    # print crn
    # print "propensities:"
    # print props
    #
    # pprint(dCovdt)
    # reactants = crn.reactions.reactants
    # reactants = [ x.reactants for x in crn.reactions]


    # eng = matlab.engine.start_matlab()
    # ret = eng.symbolicLNA(crn.species,[ x.reactants for x in crn.reactions], [ x.products for x in crn.reactions], [ x.reactionrate for x in crn.reactions])


    # inspecies, inreactants, inproducts, inrates

    # def propensity(reactions):
    #     propensities = []
    #     for reaction in reactions:
    #         propensity = symbols(str(reaction.reactionrate))
    #         for reactant in reaction.reactants:
    #             propensity *= symbols(reactant)
    #         propensities.append(propensity)
    #     return propensities
    #
    #
    # def netReactionChange(species, reactions):
    #     reactionChange = []  # ReferenceFrame('N')
    #     for reaction in reactions:
    #         netChange = []
    #         for specie in species:
    #             speciesChange = 0
    #             for reactant in reaction.reactants:
    #                 if specie == reactant:
    #                     speciesChange -= 1
    #             for product in reaction.products:
    #                 if specie == product:
    #                     speciesChange += 1
    #             netChange.append(speciesChange)
    #         reactionChange.append(Matrix(1, len(netChange), netChange))
    #     return reactionChange

    # elif(len(self.coefficient) is 2):
    # 	if(self.coefficient[-1] is 2):
    # 		return "(" + self.coefficient[0]*self.specRep() + self.coefficient[1]*symbols(self.specRep() + "^2") + ")"
    # 	else:
    # 		return "(" + self.coefficient[0]  + self.coefficient[1]*self.specRep() + ")"
    # elif(len(self.coefficient) is 3):
    # 	if(self.coefficient[-1] is 3):
    # 		return  "(" + self.coefficient[0]*self.specRep() + self.coefficient[1]*symbols(self.specRep() + "^2") + self.coefficient[2]*symbols(self.specRep() + "^3") +  ")"
    # 	elif(self.coefficient[-1] is 2):
    # 		return  "(" + self.coefficient[0] + self.coefficient[1]*self.specRep() + self.coefficient[2]*symbols(self.specRep() + "^2") + ")"
    #

    # def flowFunction(propensities, reactionChange):
    #     summationFunction = reactionChange[0] * propensities[0]
    #     for n in range(1, len(propensities)):
    #         summationFunction += reactionChange[n] * propensities[n]
    #     return summationFunction