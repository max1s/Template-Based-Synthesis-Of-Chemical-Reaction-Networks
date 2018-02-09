from sympy import *
import itertools
from six import string_types
import iSATParser
from functools import reduce


class InputSpecies:
    def __init__(self, name, ode, initial_value=None):
        self.name = name  # a string, to be printed
        self.symbol = Symbol(name)  # a sympy symbol, to be used when constructing formulae
        self.initial_value = initial_value
        self.ode = ode

    def __str__(self):
        return self.name


class Species:
    def __init__(self, name, initial_value=None):
        self.name = name
        self.symbol = Symbol(name)
        self.initial_value = initial_value

    def __str__(self):
        return self.name


class Term:
    # Represents conjunction of a species (or InputSpecies) with a stoichiometric coefficient
    def __init__(self, species, coefficient):
        self.species = species
        self.coefficient = coefficient

    def specRep(self):
        if isinstance(self.species, LambdaChoice):
            return str(self.coefficient) + "*" + self.species.constructChoice()
        elif isinstance(self.species, Choice):
            return str(self.coefficient) + "*" + " ( " + str(self.species.constructChoice()) + " ) "
        else:
            return str(self.coefficient) + "*" + str(self.species.name)

    def constructPropensity(self):
        if isinstance(self.coefficient, int):
            return self.specRep()
        else:
            raise NotImplementedError


class RateConstant:
    def __init__(self, name, minimum, maximum):
        self.name = name
        self.min = minimum
        self.max = maximum

    def __repr__(self):
        return self.name

    def __str__(self):
        return "Rate constant %s <= %s <= %s" % (self.min, self.name, self.max)


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
    def __init__(self, cs, r, opr, input_species=None):
        self.species = cs
        self.reactions = r
        self.optionalReactions = opr
        self.input_species = input_species

        self.t = symbols('t')

    def __repr__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"

    def __str__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"

    def getSpecies(self):
        # Construct list of only those species (or InputSpecies) that participate in a reaction
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

    def getRateConstants(self):
        rate_constants = {}
        for reaction in self.reactions:
            rate = reaction.reactionrate
            if str(rate) not in list(rate_constants.keys()):
                rate_constants[str(rate)] = rate
        rate_constants = list(rate_constants.values())
        return rate_constants

def parametricPropensity(paramCRN):
    # Returns a list: each element is a sympy expression corresponding to the propensity of the n'th reaction
    propensities = []
    for reaction in paramCRN.reactions:
        propensity = symbols(str(reaction.reactionrate.name))
        for reactant in reaction.reactants:
            propensity *= sympify(reactant.constructPropensity())
        propensities.append(propensity)
    return propensities


def parametricNetReactionChange(crn):
    # Returns a 2D list: change[reaction_index][species_index] is a string representing the stoichiometry change

    change = []
    for reaction in crn.reactions:
        netChange = ['0'] * len(crn.species)
        for reactant in reaction.reactants:
            add_stoichiometry_change(crn.species, netChange, reactant, '-')
        for product in reaction.products:
            add_stoichiometry_change(crn.species, netChange, product, '+')

        change.append(sympify(netChange))
    return change


def add_stoichiometry_change(species, stoichiometry_change, fragment, sign):
    for i, sp in enumerate(species):
        if str(sp) in fragment.specRep():
            if "lam" in fragment.specRep():
                new_term = " %s%s * %s" % (sign, fragment.coefficient, fragment.species.contains(sp))
            else:
                new_term = "%s%s" % (sign, fragment.coefficient)
            stoichiometry_change[i] = " + ".join([stoichiometry_change[i], new_term])
    return stoichiometry_change


def parametricFlow(propensities, reactionChange):
    return Matrix(reactionChange).transpose() * Matrix(propensities)


def michmenton(S, Vmax, v, Km):
    return v * (Vmax / (Km + S))


def hill(L, n, Ka, k):
    return k * (L ^ n / (Ka ^ n + L ^ n))


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

def generateCovarianceMatrix(speciesVector):
    # TODO: make sure we ignore InputSpecies
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
    return mat


def generateAllTokens(crn, C=set()):
    sym = [x.free_symbols for x in sympify(crn.getSpecies())]
    a = reduce(lambda x, y: x | y, sym)
    b = set()
    if len(C) is not 0:
        b = C.free_symbols
    return a | b


def derivative(derivatives, flowdict, crn):
    # define function for each state variable
    funcs = {}
    function_reverse = {}
    for variable in flowdict:
        new_function = Function(variable.name)(crn.t)
        funcs[variable] = new_function
        function_reverse[new_function] = variable

    function_flowdict = {}
    constants = []
    for variable in flowdict:
        if flowdict[variable] is None:
            constants.append(variable.subs(funcs))
        else:
            function_flowdict[variable.subs(funcs)] = flowdict[variable].subs(funcs)

    results = {}
    for d in derivatives:

        if d["is_variance"]:
            species = symbols("cov" + d["variable"])
        else:
            species = symbols(d["variable"])

        order = d["order"]
        name = d["name"]

        x1 = function_flowdict[funcs[species]]  # first derivative of species

        xn = x1
        for i in range(order):
            xn = Derivative(xn, crn.t).doit()

            # substitute in to replace first derivatives
            for func in function_flowdict:
                derivative_string = "Derivative(" + str(func) + ", t)"
                xn = xn.subs(derivative_string, function_flowdict[func])

            for constant in constants:
                derivative_string = "Derivative(" + str(constant) + ", t)"
                xn = xn.subs(derivative_string, 0)

        # replace functions with the corresponding symbols
        xn = xn.subs(function_reverse)
        results[symbols(name)] = simplify(xn)

    return results


def flowDictionary(crn, species, isLNA, derivatives, kinetics='massaction', firstConstant='2', secondConstant='2'):
    if not derivatives:
        derivatives = set()

    if isLNA:
        a = dict.fromkeys(generateAllTokens(crn, generateCovarianceMatrix(species)))
    else:
        a = dict.fromkeys(generateAllTokens(crn))

    if kinetics == 'massaction':
        prp = (parametricPropensity(crn))
        nrc = (parametricNetReactionChange(crn))
        dSpeciesdt = parametricFlow(prp, nrc)
    elif kinetics == 'hill':
        dSpeciesdt = hillKineticsFlow(species, firstConstant, [y.reactionrate for y in x for x in crn.reactions],
                                      secondConstant)
    elif kinetics == 'michaelis-menton':
        dSpeciesdt = michaelisMentonFlow(species, firstConstant, [y.reactionrate for y in x for x in crn.reactions],
                                         secondConstant)
    for sp, i in zip(species, list(range(len(species)))):
        if isinstance(sp, str):
            a[symbols(sp)] = dSpeciesdt[i]
        else:
            a[sp.symbol] = dSpeciesdt[i]

    if isLNA:
        jmat = [sp.symbol for sp in species]
        J = Matrix(dSpeciesdt).jacobian(jmat)
        C = generateCovarianceMatrix(species)
        dCovdt = J * C + C * transpose(J)
        for i in range(C.cols * C.rows):
            a[C[i]] = dCovdt[i]

    a.update(derivative(derivatives, a, crn))

    for sp in crn.input_species:
        a[sp.symbol] = sp.ode

    for key in a:
        if a[key] is None and not isinstance(a[key], str):
            a[key] = 0
    return a


def intDictionary(crn, species, isLNA, flowdict):
    a = dict.fromkeys(list(flowdict.keys()))
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
    i = i ** 2 + 1

    if isLNA:
        for co in generateCovarianceMatrix(species):
            a[co] = i
    return a


def exampleParametricCRN():
    X = Species('X')
    Y = Species('Y')
    B = Species('B')

    reaction1 = Reaction([Term(LambdaChoice([X, Y], 1), 1), Term(Y, 1)], [Term(X, 1), Term(B, 1)],
                         RateConstant('k_1', 1, 2))
    reaction2 = Reaction([Term(LambdaChoice([X, Y], 2), 1), Term(Choice(X, 0, 2), 2)],
                         [Term(Y, 1), Term(B, 1)], RateConstant('k_2', 1, 2))
    reaction3 = Reaction([Term(X, 1), Term(B, 1)], [Term(X, 1), Term(X, 1)], RateConstant('k_3', 1, 2))
    reaction4 = Reaction([Term(X, 1), Term(B, 1)], [Term(X, 1), Term(X, 1)], RateConstant('k_4', 1, 2))

    input1 = InputSpecies("Input1", sympify("0.1*t + 54.2735055776743*exp(-(0.04*t - 2.81375654916915)**2) + 35.5555607722356/(1.04836341039216e+15*(1/t)**10.0 + 1)"), 15)
    reaction5 = Reaction([Term(input1, 1)], [Term(B, 1)], RateConstant('k_input', 1, 2))

    isLNA = False
    derivatives = [{"variable": 'X', "order": 1, "is_variance": False, "name": "X_dot"},
                   {"variable": 'X', "order": 2, "is_variance": False, "name": "X_dot_dot"},
                   {"variable": 'X', "order": 2, "is_variance": True, "name": "covX_dot_dot"}]
    derivatives = []
    specification = [(0, 'X = 0'), (0.5, 'X = 0.5'), (1, 'X = 0')]
    # file = iSATParser(flow, ints, specification)

    crn = CRNSketch([X, Y, B], [reaction1, reaction2, reaction3, reaction5], [reaction4], [input1])

    flow = flowDictionary(crn, [X, Y, B], isLNA, derivatives)
    ints = intDictionary(crn, [X, Y, B], isLNA, flow)

    d, i, m, p, t = iSATParser.constructSpecification(specification, flow, ints, '', rate_constants=crn.getRateConstants())
    spec = iSATParser.constructISAT(d, i, m, p, t)
    print(spec)


if __name__ == "__main__":
    exampleParametricCRN()
