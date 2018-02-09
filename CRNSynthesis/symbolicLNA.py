from sympy import *
import itertools
from six import string_types
from functools import reduce


class InputSpecies:
    def __init__(self, name, ode, initial_value=None):
        self.name = name  # a string, to be printed
        self.symbol = Symbol(name)  # a sympy symbol, to be used when constructing formulae
        self.initial_value = initial_value
        self.ode = ode

    def __str__(self):
        return self.name

    def iSATDefinition(self):
        return "\tfloat[%s, %s] %s;\n" % (0, 100, self.name) # TODO: set min/max better

    def iSATInitialization(self):
        return "\t%s = %s;\n" % (self.name, self.initial_value)

    def get_species(self):
        return [self]

    def get_real_species(self):
        return []

class Species:
    def __init__(self, name, initial_value=None, initial_min=None, initial_max=None):
        self.name = name
        self.symbol = Symbol(name)
        self.initial_value = initial_value

        self.initial_min = initial_min
        self.initial_max = initial_max

    def __str__(self):
        return self.name

    def iSATDefinition(self):
        return "\tfloat[0, %s] %s;\n" % (10, self.name)

    def iSATInitialization(self):
        if self.initial_value:
            return "\t%s = %s;\n" % (self.name, self.initial_value)

        terms = []
        if self.initial_min:
            terms.append("(%s >= %s)" % (self.name, self.initial_min))

        if self.initial_max:
            terms.append("(%s <= %s)" % (self.name, self.initial_max))

        if len(terms) > 0:
            return "\t" + " and ".join(terms) + ";\n"
        return ""

    def get_species(self):
        return [self]

    def get_real_species(self):
        return [self]

class Term:
    # Represents conjunction of a species (or InputSpecies) with a stoichiometric coefficient
    def __init__(self, species, coefficient):
        self.species = species
        self.coefficient = coefficient

    def specRep(self):
        coefficient = self.coefficient
        if isinstance(coefficient, Choice):
            coefficient = coefficient.symbol

        if isinstance(self.species, LambdaChoice):
            return str(coefficient) + "*" + self.species.constructChoice()
        else:
            return str(coefficient) + "*" + str(self.species.name)

    def constructPropensity(self):
        if not isinstance(self.coefficient, int) and not isinstance(self.coefficient, Choice):
            raise NotImplementedError

        coefficient = self.coefficient
        if isinstance(coefficient, Choice):
            coefficient = coefficient.symbol

        if isinstance(self.species, LambdaChoice):
            return self.species.constructChoice() + "**" + str(coefficient)
        else:
            return str(self.species.name) + "**" + str(coefficient)

    def get_species(self):
        return self.species.get_species()

    def get_real_species(self):
        return self.species.get_real_species()


class RateConstant:
    def __init__(self, name, minimum=0.1, maximum=10):
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

    def get_species(self):
        return sum([x.get_species() for x in self.species], [])  # flatten list-of-lists into list

    def get_real_species(self):
        return sum([x.get_real_species() for x in self.species], [])  # flatten list-of-lists into list

    def contains(self, variable):
        x = ''
        for lam in self.lambdas:
            if str(variable) in str(lam):
                x += str(lam)
        return x

    def format_constraint(self):
        clauses = []
        for active_value in self.lambdas:
            # generate term in which only element i is on
            subclauses = []
            for lam in self.lambdas:
                if lam == active_value:
                    subclauses.append("(%s = 1)" % lam)
                else:
                    subclauses.append("(%s = 0)" % lam)

            clause = "(" + " and ".join(subclauses) + ")"
            clauses.append(clause)

        return "\t" + " or ".join(clauses) + ";\n"

    def iSATDefinition(self):
        declarations = ["\tfloat[0, 1] %s;" % str(lam) for lam in self.lambdas]
        return "\n".join(declarations)


class Choice:
    def __init__(self, choiceNumber, minValue, maxValue):
        self.choiceNumber = choiceNumber
        self.name = str('c' + str(self.choiceNumber))
        self.symbol = Symbol(self.name)
        self.minValue = minValue
        self.maxValue = maxValue

    def format_constraint(self):
        clauses = ["(%s = %s)" % (self.name, x) for x in range(self.minValue, self.maxValue+1)]
        return "\t" + " or ".join(clauses) + ";\n"

    def iSATDefinition(self):
        return "\tfloat[%s, %s] %s;\n" % (self.minValue, self.maxValue, self.name)

    def __str__(self):
        return 'c' + str(self.choiceNumber)


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
    def __init__(self, r, opr, input_species=None):
        self.reactions = r
        self.optionalReactions = opr
        self.input_species = input_species

        self.species = self.getSpecies()
        self.real_species = self.getSpecies(include_inputs=False)
        self.input_species = list(set(self.species)-set(self.real_species))

        self.choice_variables = set()
        self.lambda_variables = set()
        self.record_decision_variables()

        self.t = symbols('t')

    def __repr__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"

    def __str__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.reactions]) + "\n]"

    def record_decision_variables(self):
        for reaction in self.reactions:
            terms = reaction.reactants[:]
            terms.extend(reaction.products)

            for term in terms:
                if isinstance(term.species, LambdaChoice):
                    self.lambda_variables.add(term.species)

                if isinstance(term.coefficient, Choice):
                    self.choice_variables.add(term.coefficient)

    def generateAllTokens(self, C=set()):

        species_strings = [] # entries will be strings representing caluses that appear reactants/products

        all_reactions = self.reactions[:]
        all_reactions.extend(self.optionalReactions)

        for y in all_reactions:
            reactants_or_products = y.reactants[:]
            reactants_or_products.extend(y.products)
            for react in reactants_or_products:
                if react not in species_strings:
                    species_strings.append(react.specRep())


        sym = [x.free_symbols for x in sympify(species_strings)]
        a = reduce(lambda x, y: x | y, sym)
        b = set()
        if len(C) is not 0:
            b = C.free_symbols
        return a | b



    def getSpecies(self, include_inputs=True):
        # Construct list of only those species (or InputSpecies) that participate in a reaction
        x = set()

        all_reactions = self.reactions[:]
        all_reactions.extend(self.optionalReactions)

        for y in all_reactions:
            reactants_or_products = y.reactants[:]
            reactants_or_products.extend(y.products)
            for sp in reactants_or_products:
                if sp not in x:
                    if include_inputs:
                        x = x.union(sp.get_species())
                    else:
                        x = x.union(sp.get_real_species())

        return list(x)

    def getRateConstants(self):
        rate_constants = {}
        for reaction in self.reactions:
            rate = reaction.reactionrate
            if str(rate) not in list(rate_constants.keys()):
                rate_constants[str(rate)] = rate
        rate_constants = list(rate_constants.values())
        return rate_constants

    def flow(self, isLNA, derivatives, kinetics='massaction', firstConstant='2', secondConstant='2'):
        if not derivatives:
            derivatives = set()

        if isLNA:
            a = dict.fromkeys(self.generateAllTokens(generateCovarianceMatrix(self.real_species)))
        else:
            a = dict.fromkeys(self.generateAllTokens())

        if kinetics == 'massaction':
            propensities = self.parametricPropensity()
            stoichiometry_change = self.parametricNetReactionChange()
            dSpeciesdt = Matrix(stoichiometry_change).transpose() * Matrix(propensities)
        elif kinetics == 'hill':
            dSpeciesdt = hillKineticsFlow(self.species, firstConstant, [y.reactionrate for y in x for x in self.reactions],
                                          secondConstant)
        elif kinetics == 'michaelis-menton':
            dSpeciesdt = michaelisMentonFlow(self.species, firstConstant, [y.reactionrate for y in x for x in self.reactions],
                                             secondConstant)
        for i, sp in enumerate(self.real_species):
            if isinstance(sp, str):
                a[symbols(sp)] = dSpeciesdt[i]
            else:
                a[sp.symbol] = dSpeciesdt[i]

        if isLNA:
            jmat = [sp.symbol for sp in self.species]
            J = Matrix(dSpeciesdt).jacobian(jmat)
            C = generateCovarianceMatrix(self.species)
            dCovdt = J * C + C * transpose(J)
            for i in range(C.cols * C.rows):
                a[C[i]] = dCovdt[i]

        a.update(derivative(derivatives, a))

        for sp in self.input_species:
            a[sp.symbol] = sp.ode

        constants_to_remove = []
        for key in a:
            if a[key] is None and not isinstance(a[key], str):
                a[key] = 0
            if str(key) not in [str(sp) for sp in self.species]:
                constants_to_remove.append(key)

        # remove constant keys from flowDict, as they are handled separately when output generated
        for key in constants_to_remove:
            a.pop(key, None)

        return a

    def parametricPropensity(self):
        # Returns a list: each element is a sympy expression corresponding to the propensity of the n'th reaction
        propensities = []
        for reaction in self.reactions:
            propensity = symbols(str(reaction.reactionrate.name))
            for reactant in reaction.reactants:
                propensity *= sympify(reactant.constructPropensity())
            propensities.append(propensity)
        return propensities

    def parametricNetReactionChange(self):
        # Returns a 2D list: change[reaction_index][species_index] is a string representing the stoichiometry change

        change = []
        for reaction in self.reactions:
            netChange = ['0'] * len(self.real_species)
            for reactant in reaction.reactants:
                add_stoichiometry_change(self.real_species, netChange, reactant, '-')
            for product in reaction.products:
                add_stoichiometry_change(self.real_species, netChange, product, '+')

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


# TODO: these have same bug as parametric flow
# should just be modifying propensity calculation
def michaelisMentonFlow(species, Vmax, v, Km):
    m = Matrix(1, len(species))
    for spec, i in zip(species, len(species)):
        m[1, i] = v[i] * (Vmax / (Km + spec))
    return m


def hillKineticsFlow(species, Ka, k, n):
    m = Matrix(1, len(species))
    for spec, i in zip(species, len(species)):
        m[1, i] = k[i] * (spec ^ n / (Ka ** n + spec ** n))
    return m


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
    return mat


def derivative(derivatives, flowdict):
    # define function for each state variable
    funcs = {}
    function_reverse = {}
    for variable in flowdict:
        new_function = Function(variable.name)(Symbol('t'))
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
            xn = Derivative(xn, Symbol('t')).doit()

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
