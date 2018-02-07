# import matlab.engine
from sympy import *
from sympy import Matrix
import sys
from sympy import init_printing
import itertools
from six import string_types
import iSATParser
from functools import reduce


class InputSpecies:
    def __init__(self, name, concentration, ode):
        self.name = name
        self.concentration = concentration
        self.ode = ode

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
        if self.minNumber == 0:
            chain += self.choice[0]
        if self.maxNumber > 0:
            chain += " + " + self.choice[1]
        if self.maxNumber > 1:
            chain += " + ".join([str(choice) + "*" + str(self.choiceName) + "^" + str(x) for choice, x in
                                 zip(self.choice, list(range(2, self.maxNumber)))])
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
        if isinstance(self.coefficient, int):
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
        self.t = symbols('t')

    def __repr__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"

    def __str__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"

    def getSpecies(self):
        x = []

        all_reactions = self.reactions[:]
        all_reactions.extend(self.optionalReactions)

        for y in all_reactions:
            for react in y.reactants:
                if react not in x:
                    x.append(react.specRep())
            for react in y.products:
                if react not in x:
                    x.append(react.specRep())
        return x

    def getRawSpecies(self):
        x = []
        for y in self.reactions:
            for react in y.reactants:
                if react not in x:
                    x.append(react)
            for react in y.products:
                if react not in x:
                    x.append(react)
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
            add(crn.species, netChange, product, '+')

        netChange = [0 if n is '' else n for n in netChange]
        reactionChange.append(sympify(netChange))
    return reactionChange


def add(species, stoichList,  fragment, part):
    for spec, i in zip(species, list(range(len(species)))):
        if str(spec) in fragment.specRep():
            if "lam" in fragment.specRep():
                stoichList[i] += part + fragment.species.contains(spec) + "*" + str(spec)
            else:
                stoichList[i] += part + fragment.specRep()
    return stoichList


def parametricFlow(propensities, reactionChange):
    container = [reactionChange[0] * propensities[0]]
    for n in range(1, len(propensities)):
        container.append(reactionChange[n] * propensities[n])

    return container

def michmenton(S, Vmax, v, Km):
    return v * ( Vmax / (Km + S))

def hill(L, n, Ka, k):
    return k * ( L^n / (Ka^n + L^n) )

def michaelisMentonFlow(species, Vmax, v, Km):
    m = Matrix(1, len(species))
    for spec, i in zip(species, len(species)):
        m[1, i] = michmenton(spec, Vmax, v[i], Km)
    return m

def hillKineticsFlow(species, Ka, k, n):
    m = Matrix(1, len(species))
    for spec, i in zip(species, len(species)):
        m[1, i] = hill(spec, n, Ka, k[i])
    return m


def parametricG(propensities, reactionChange):
    G = zeros(max(len(reactionChange.col(0)), len(reactionChange.row(0))),
              max(len(reactionChange.col(0)), len(reactionChange.row(0))))
    for i in range(len(propensities)):
        G += transpose(reactionChange.row(i)) * reactionChange.row(i) * propensities[i]
    return G


def generateCovarianceMatrix(speciesVector):
    mat = eye(len(speciesVector))
    for (m, i) in zip(speciesVector, list(range(len(speciesVector)))):
        for (n, j) in zip(speciesVector, list(range(len(speciesVector)))):
            if m == n:
                mat[i, j] = 'cov' + str(m)
            else:
                mat[i, j] = 'cov' + str(n) + str(m)

    for x in range(len(speciesVector)):
        for y in range(len(speciesVector)):
            mat[x, y] = mat[y, x]
    # pprint(mat)
    return mat

def generateAllTokens(crn, derivatives, C = set()):
    sym = [x.free_symbols for x in sympify(crn.getSpecies())]
    a = reduce(lambda x, y : x |y, sym)
    b = set()
    if len(C) is not 0:
        b = C.free_symbols
    if derivatives is None:
        return a | b
    else:
        return a | b | sympify(derivatives)



def derivative(species, flowdict, crn):

    # define function for each state variable
    functions = {}
    function_reverse = {}
    for variable in flowdict:
        new_function = Function(variable.name)(crn.t)
        functions[variable] = new_function
        function_reverse[new_function] = variable

    function_flowdict = {}
    constants = []
    for variable in flowdict:
        if flowdict[variable] is None:
            constants.append(variable.subs(functions))
        else:
            function_flowdict[variable.subs(functions)] = flowdict[variable].subs(functions)

    x1 = function_flowdict[functions[species]] # first derivative of species
    x2 = Derivative(x1, crn.t).doit() # second derivative

    # substitute in to replace first derivatives
    for function in function_flowdict:
        derivative_string = "Derivative(" + str(function) + ", t)"
        x2 = x2.subs(derivative_string, function_flowdict[function])

    for constant in constants:
        derivative_string = "Derivative(" + str(constant) + ", t)"
        x2 = x2.subs(derivative_string, 0)

    # replace functions with the corresponding symbols
    x2 = x2.subs(function_reverse)

    return simplify(x2)


#rate,ratemax, constant
def flowDictionary(crn, species, isLNA, derivatives, kinetics='massaction', firstConstant='2', secondConstant='2'):
    if not derivatives:
        derivatives = set()

    if isLNA:
        a = dict.fromkeys(generateAllTokens(crn, set(derivatives), generateCovarianceMatrix(species)))
    else:
        a = dict.fromkeys(generateAllTokens(crn, set(derivatives)))

    if kinetics == 'massaction':
        prp = (parametricPropensity(crn))
        nrc = (parametricNetReactionChange(crn))
        dSpeciesdt = parametricFlow(prp, Matrix(nrc))
    elif kinetics == 'hill':
        dSpeciesdt = hillKineticsFlow(species,firstConstant, [y.reactionrate for y in x for x in crn.reactions], secondConstant )
    elif kinetics == 'michaelis-menton':
        dSpeciesdt = michaelisMentonFlow(species, firstConstant, [y.reactionrate for y in x for x in crn.reactions], secondConstant )
    for sp, i in zip(species, list(range(len(species)))):
        if isinstance(sp, str):
            a[symbols(sp)] = dSpeciesdt[i]
        else:
            a[sp] = dSpeciesdt[i]

    for der in derivatives:
        a[symbols(der)] = derivative(symbols('X'), a, crn)
    jmat = [x for x in species]
    J = Matrix(dSpeciesdt).jacobian(jmat)
    G = parametricG(Matrix(prp), Matrix(nrc))
    C = generateCovarianceMatrix(species)
    dCovdt = J * C + C * transpose(J)
    for i in range(C.cols*C.rows):
        a[C[i]] = dCovdt[i]
    for key in a:
        if a[key] is None and not isinstance(a[key], str):
            a[key] = 0
    return a



def hillFlowDictionary(crn, species, isLNA, derivatives):
    a = dict.fromkeys(generateAllTokens(crn, set(derivatives), generateCovarianceMatrix(species))) if isLNA else dict.fromkeys(generateAllTokens(crn, set.add(derivatives)))
    prp = (parametricPropensity(crn))
    nrc = (parametricNetReactionChange(crn))
    dSpeciesdt = parametricFlow(prp, Matrix(nrc))
    for sp, i in zip(species, list(range(len(species)))):
        if isinstance(sp, str):
            a[symbols(sp)] = dSpeciesdt[i]
        else:
            a[sp] = dSpeciesdt[i]

    jmat = [x for x in species]
    J = Matrix(dSpeciesdt).jacobian(jmat)
    G = parametricG(Matrix(prp), Matrix(nrc))
    C = generateCovarianceMatrix(species)
    dCovdt = J * C + C * transpose(J)
    for i in range(C.cols*C.rows):
        a[C[i]] = dCovdt[i]
    for key in a:
        if a[key] is None and not isinstance(a[key], str):
            a[key] = 0
    return a

def intDictionary(crn, species, covariance, flowdict):
    #getInt(crn)
    a = dict.fromkeys(list(flowdict.keys()))
    t = [x for x in crn.getRawSpecies()]
    for reaction in crn.reactions:
        for reactant in reaction.reactants:
            if isinstance(reactant.species, LambdaChoice):
                for x in reactant.species.lambdas:
                    a[x] = len(reactant.species.species)
            elif isinstance(reactant.species, Choice):
                for x in reactant.species.choice:
                    a[symbols(x)] = reactant.species.maxNumber
            else:
                a[reactant.species] = reactant.coefficient

        for product in reaction.products:
            if isinstance(product.species, LambdaChoice):
                for x in product.species.lambdas:
                    a[x] = len(product.species.species)
            elif isinstance(product.species, Choice):
                for x in product.species.choice:
                    a[symbols(x)] = product.species.maxNumber
            else:
                a[product.species] = product.coefficient

    i = 0
    for spec in species:
        i = max(a[spec], i)
    i = i**2 + 1
    for co in covariance:
        a[co] = i
    return a






def exampleParametricCRN():
    X = symbols('X')
    Y = symbols('Y')
    B = symbols('B')


    reaction1 = Reaction([Species(LambdaChoice([X, Y], 1), 1), Species(Y, 1)], [Species(X, 1), Species(B, 1)], RateConstant('k_1', 1, 2))
    reaction2 = Reaction([Species(LambdaChoice([X, Y], 2), 1), Species(Choice(X, 0, 2), 2)], [Species(Y, 1), Species(B, 1)], RateConstant('k_1', 1, 2))
    reaction3 = Reaction([Species(X, 1), Species(B, 1)], [Species(X, 1), Species(X, 1)], RateConstant('k_1', 1, 2))
    reaction4 = Reaction([Species(X, 1), Species(B, 1)], [Species(X, 1), Species(X, 1)], RateConstant('k_1', 1, 2))

    crn = CRNSketch([X, Y, B], [reaction1, reaction2, reaction3], [reaction4])


    #pprint(dCovdt)
    isLNA = False
    derivatives = set(['dXdt']) # set(['dXdt'])
    flow = flowDictionary(crn, [X, Y, B], isLNA, derivatives)
    ints = intDictionary(crn, [X, Y, B], generateCovarianceMatrix([X, Y, B]), flow)
    specification = [(0, 'X = 0'), (0.5, 'X = 0.5'), (1, 'X = 0') ]
    #file = iSATParser(flow, ints, specification)

    rate_constants = {}
    for reaction in crn.reactions:
        rate = reaction.reactionrate
        if str(rate) not in list(rate_constants.keys()):
            rate_constants[str(rate)] = rate
    rate_constants = list(rate_constants.values())

    d,i,m,p,t = iSATParser.constructSpecification(specification, flow, ints, '', rate_constants=rate_constants)
    spec = iSATParser.constructISAT(d,i,m,p,t)
    print(spec)





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