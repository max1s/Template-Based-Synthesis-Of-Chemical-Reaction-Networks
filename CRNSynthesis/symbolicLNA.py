from sympy import *
import itertools
from six import string_types
from functools import reduce
from operator import mul


class InputSpecies:
    """
    This class represents an Input Species - that is, a species whose concentration follows a fixed trajectory,
    rather than being determined by the kinetics of reactions in which it is a reactant or product.
    """

    def __init__(self, name, ode, initial_value=None):
        """
        :param name: name of the species (string)
        :param ode: time-derivative of the cocnentration, expressed in terms of 't' (SymPy expression)
        :param initial_value: initial value of the concentration (float)

        """

        self.name = name  # a string, to be printed
        self.symbol = Symbol(name)  # a sympy symbol, to be used when constructing formulae
        self.initial_value = initial_value
        self.ode = ode

    def __str__(self):
        return self.name

    def iSATDefinition(self):
        """
        Returns string that defines variable in iSAT (.hys) format.
        """

        return "\tfloat[%s, %s] %s;\n" % (0, 100, self.name) # TODO: set min/max better

    def iSATInitialization(self):
        """
        Returns string that assigns initial value to variable in iSAT (.hys) format.
        """
        return "\t%s = %s;\n" % (self.name, self.initial_value)

    def dRealDefinition(self):
        """
        Returns string that defines variable in dReach (.drh) format.
        """
        return "\t[%s, %s] %s;\n" % (0, 10, self.name) # TODO: set min/max better

    def dRealInitialization(self):
        """
        Returns string that assigns initial value to variable in dReach (.drh) format.
        """
        return "\t(%s = %s)\n" % (self.name, self.initial_value)

    def get_species(self):
        """
        Returns an list containing self, as this InputSpecies is a species.
        """
        return [self]

    def get_real_species(self):
        """
        Returns an empty list, as this InputSpecies is a not a real (non-input) species.
        """

        return []

class Species:
    """
    This class represents an Species - that is, a species whose concentration is determined by the kinetics of reactions
     in which it is a reactant or product.
    """

    def __init__(self, name, initial_value=None, initial_min=None, initial_max=None):
        """
        This constructor requires the speification of either the initial value of the concentration,
         or upper and lower limits on the initial cocentration.

        :param name:  name of the species (string)
        :param initial_value:  initial value of concentration for this species (float)
        :param initial_min:  lower limit of possible initial values of concentration for this species (float)
        :param initial_max:  upper limit of possible initial value of concentration for this species (float)
        """

        self.name = name
        self.symbol = Symbol(name)
        self.initial_value = initial_value

        self.initial_min = initial_min
        self.initial_max = initial_max

    def __str__(self):
        return self.name

    def iSATDefinition(self):
        """
        Returns string that defines variable in iSAT (.hys) format.
        """
        return "\tfloat[0, %s] %s;\n" % (1, self.name)

    def iSATInitialization(self):
        """
        Returns string that assigns initial value to variable in iSAT (.hys) format.
        """
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

    def dRealDefinition(self):
        """
        Returns string that defines variable in dReach (.drh) format.
        """
        return "\t[0, %s] %s;\n" % (10, self.name)

    def dRealInitialization(self):
        """
        Returns string that assigns initial value to variable in dReach (.drh) format.
        """
        if self.initial_value:
            return "\t(%s = %s)\n" % (self.name, self.initial_value)

        terms = []
        if self.initial_min:
            terms.append(" (%s >= %s)" % (self.name, self.initial_min))

        if self.initial_max:
            terms.append(" (%s <= %s)" % (self.name, self.initial_max))

        if len(terms) > 0:
            return "\t ( and" + " ".join(terms) + ")\n"
        return ""


    def get_species(self):
        """
        Returns an list containing self, as this Species is a species.
        """
        return [self]

    def get_real_species(self):
        """
        Returns an list containing self, as this InputSpecies is a real (non-input) species.
        """
        return [self]

class Term:
    """
    This class represents the pairing of a species (or InputSpecies) with a stoichiometric coefficient.
    The species may be either a specific species, or a LambdaChoice object representing a choice of species.
    The stoichiometric coefficient may be either an integer or a Choice object representing a choice of values.
    """

    def __init__(self, species, coefficient):
        """

        :param species: a Species, InputSpecies, or LambdaChoice object
        :param coefficient: an Integer or Choice object representing the stoichiometric coefficient
        """

        self.species = species
        self.coefficient = coefficient

    def __str__(self):
        return "%s%s" % (self.coefficient, self.species)

    def specRep(self):
        coefficient = self.coefficient
        if isinstance(coefficient, Choice):
            coefficient_expression = 0
            for value in list(range(coefficient.minValue, coefficient.maxValue+1)):
                coefficient_expression += value * symbols(coefficient.name + "_" + str(value))
            coefficient = coefficient_expression

        if isinstance(self.species, LambdaChoice):
            return str(coefficient * sympify(self.species.constructChoice()))
        else:
            return str(coefficient * sympify(self.species.name))

    def constructPropensity(self):
        if not isinstance(self.coefficient, int) and not isinstance(self.coefficient, Choice):
            raise NotImplementedError

        if isinstance(self.species, LambdaChoice):
            term_to_exponentiate = self.species.constructChoice()
        else:
            term_to_exponentiate = str(self.species.name)

        coefficient = self.coefficient

        if isinstance(coefficient, Choice):
            terms =[]
            for value in list(range(coefficient.minValue, coefficient.maxValue+1)):
                terms.append("(%s ** %s)*%s_%s" % (term_to_exponentiate, value, coefficient.name, value))
            return " + ".join(terms)

        else:
            return term_to_exponentiate + "**" + str(coefficient)


    def get_species(self):
        return self.species.get_species()

    def get_real_species(self):
        return self.species.get_real_species()


class RateConstant:
    def __init__(self, name, minimum=0.1, maximum=10):
        self.name = name
        self.min = minimum
        self.max = maximum
        self.symbol = symbols(name)

    def __repr__(self):
        return self.name

    def __str__(self):
        return "Rate constant %s <= %s <= %s" % (self.min, self.name, self.max)


class Reaction(object):
    def __init__(self, r, p, ra):
        self.reactants = r

        for i, r in enumerate(self.reactants):
            if isinstance(r, tuple):
                self.reactants[i] = Term(r[0], r[1])

        self.products = p
        for i, p in enumerate(self.products):
            if isinstance(p, tuple):
                self.products[i] = Term(p[0], p[1])

        self.reactionrate = ra
        self.is_optional = False
        self.variable_name = "" # Boolean variable indicating whether optional reaction occurs

    def __repr__(self):
        return "" + ' + '.join([repr(x) for x in self.reactants]) + " ->{" + repr(self.reactionrate) + "} " + ' + '.join([repr(y) for y in self.products])

    def __str__(self):
        return "" + ' + '.join([str(x) for x in self.reactants]) + " ->{" + str(self.reactionrate) + "} " + ' + '.join([str(y) for y in self.products])

    def get_rate_constants(self):
        return [self.reactionrate]

    def get_propensity(self):
        propensity = symbols(str(self.reactionrate.name))
        for reactant in self.reactants:
            propensity *= sympify(reactant.constructPropensity())
        if self.is_optional:
            propensity *= sympify(self.variable_name)
        return propensity

class ArbitraryRateReaction(Reaction):
        def __init__(self, r, p, ra, constants):
            super(ArbitraryRateReaction, self).__init__(r, p, ra)
            self.reactionrate = sympify(ra) # rate is a sympy formula, rather than constant
            self.constants = constants

        def get_rate_constants(self):
            return self.constants

        def get_propensity(self):
            return self.reactionrate


class HillActivationReaction(Reaction):
    def __init__(self, r, p, Kmax, Ka, n):
        super(HillActivationReaction, self).__init__(r, p, Kmax)
        self.Ka = Ka
        self.n = n

    def get_rate_constants(self):
        return [self.reactionrate, self.Ka, self.n]

    def get_propensity(self):
        species = getSpeciesFromTerm(self.reactants[0])
        return self.reactionrate.symbol * (species ** self.n.symbol / (self.Ka.symbol ** self.n.symbol + species ** self.n.symbol))

class HillRepressionReaction(Reaction):
    def __init__(self, r, p, Kmax, Ka, n):
        super(HillRepressionReaction, self).__init__(r, p, Kmax)
        self.Ka = Ka
        self.n = n

    def get_rate_constants(self):
        return [self.reactionrate, self.Ka, self.n]

    def get_propensity(self):
        species = getSpeciesFromTerm(self.reactants[0])
        return self.reactionrate.symbol / (self.Ka.symbol ** self.n.symbol + species ** self.n.symbol)


class MichaelisMentenReaction(Reaction):
    def __init__(self, r, p, Kmax, Km):
        super(MichaelisMentenReaction, self).__init__(r, p, Kmax)
        self.Km = Km

    def get_rate_constants(self):
        return [self.reactionrate, self.Km]

    def get_propensity(self):
        species = getSpeciesFromTerm(self.reactants[0])
        return species * (self.reactionrate.symbol / (self.Km.symbol + species))


def getSpeciesFromTerm(initial_term):
    # This returns something like specRep without the stoichiometric coefficients

    terms = []

    if isinstance(initial_term, TermChoice):
        for i, term in enumerate(initial_term.possible_terms):
            if isinstance(term.species, LambdaChoice):
                terms.append(sympify("%s%s * (%s)" % (initial_term.base_variable_name, i, term.species.constructChoice())))
            else:
                terms.append(sympify("%s%s * (%s)" % (initial_term.base_variable_name, i, term.species.name)))

        return sum(terms)

    elif isinstance(initial_term.species, LambdaChoice):
        return sympify(initial_term.species.constructChoice())
    else:
        return sympify(initial_term.species.name)


class LambdaChoice:
    def __init__(self, species, choiceNumber):
        self.choiceNo = choiceNumber
        self.species = species
        self.lambdas = [symbols("lam" + str(sp) + str(choiceNumber)) for sp in species]

    def __str__(self):
        return "(" + " or ".join([str(x) for x in self.species]) + ")"

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

    def iSATformat_constraint(self):
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

    def dRealformat_constraint(self):
        clauses = []
        for active_value in self.lambdas:
            # generate term in which only element i is on
            subclauses = []
            for lam in self.lambdas:
                if lam == active_value:
                    subclauses.append("(%s = 1)" % lam)
                else:
                    subclauses.append("(%s = 0)" % lam)

            clause = "( and " + "".join(subclauses) + ")"
            clauses.append(clause)

        return "\t" + " ( or " + "".join(clauses) + ")\n"

    def iSATDefinition(self):
        """
        Returns string that defines variable in iSAT (.hys) format.
        """

        declarations = ["\tfloat [0, 1] %s;" % str(lam) for lam in self.lambdas]
        return "\n".join(declarations)

    def dRealDefinition(self):
        """
        Returns string that defines variable in dReach (.drh) format.
        """
        declarations = ["\t[0, 1] %s;" % str(lam) for lam in self.lambdas]
        return "\n".join(declarations)


class Choice:
    def __init__(self, choiceNumber, minValue, maxValue):
        self.choiceNumber = choiceNumber
        self.name = str('c' + str(self.choiceNumber))
        self.symbol = Symbol(self.name)
        self.minValue = minValue
        self.maxValue = maxValue

    def iSATformat_constraint(self):
        clauses = []

        for i in list(range(self.minValue, self.maxValue + 1)):
            terms = []
            for j in list(range(self.minValue, self.maxValue + 1)):
                terms.append("(%s_%s = %s)" % (self.name, j, int(i==j)))
            clauses.append(" and ".join(terms))

        return "\t" + " or ".join([" (%s) " % clause for clause in clauses]) + ";\n"

    def dRealformat_constraint(self):
        clauses = []

        for i in list(range(self.minValue, self.maxValue + 1)):
            terms = []
            for j in list(range(self.minValue, self.maxValue + 1)):
                terms.append("(%s_%s = %s)" % (self.name, j, int(i==j)))
            clauses.append(" ( and " + "".join(terms) + ")")

        alternatives = " ".join([" (%s) " % clause for clause in clauses])
        return "\t ( or %s) " % alternatives

    def iSATDefinition(self):
        """
        Returns string that defines variable in iSAT (.hys) format.
        """

        s = ""
        for value in list(range(self.minValue, self.maxValue + 1)):
            s += "\tfloat [0, 1] %s_%s;\n" % (self.name, value)
        return s

    def dRealDefinition(self):
        """
        Returns string that defines variable in dReach (.drh) format.
        """
        s = ""
        for value in list(range(self.minValue, self.maxValue + 1)):
            s += "\t[0, 1] %s_%s;\n" % (self.name, value)
        return s

    def __str__(self):
        return self.name


class TermChoice:
    def __init__(self, termChoiceNumber, choices):
        self.termChoiceNumber = termChoiceNumber
        self.possible_terms = choices # each choice is a term
        self.base_variable_name = "tc_%s_" % self.termChoiceNumber

        for i, term in enumerate(self.possible_terms):
            if isinstance(term, tuple):
                self.possible_terms[i] = Term(term[0], term[1])

    def __str__(self):
        return "(" + " or ".join([str(x) for x in self.possible_terms]) + ")"

    def constructPropensity(self):
        terms = []
        for i, t in enumerate(self.possible_terms):
            terms.append("%s%s * (%s)" % (self.base_variable_name, i, t.constructPropensity()))
        return " + ".join(terms)

    def specRep(self):
        substrings = []
        for i, term in enumerate(self.possible_terms):
            substrings.append("%s%s * (%s)" % (self.base_variable_name, i, term.specRep()))

        return " + ".join(substrings)

    def list_decision_variables(self):
        return [self.base_variable_name + str(i) for i in list(range(0, len(self.possible_terms)))]

    def get_species(self):
        x = set()
        for choice in self.possible_terms:
            x = x.union(choice.get_species())
        return list(x)

    def get_real_species(self):
        x = set()
        for choice in self.possible_terms:
            x = x.union(choice.get_real_species())
        return list(x)

    def iSATDefinition(self):
        """
        Returns string that defines variable in iSAT (.hys) format.
        """

        declarations = []
        for i, term in enumerate(self.possible_terms):
            declarations.append("\tfloat[0, 1] %s%s;" % (self.base_variable_name, i))
        return "\n".join(declarations)

    def dRealDefinition(self):
        """
        Returns string that defines variable in dReach (.drh) format.
        """
        declarations = []
        for i, term in enumerate(self.possible_terms):
            declarations.append("\t [0, 1] %s%s;" % (self.base_variable_name, i))
        return "\n".join(declarations)

    def iSATformat_constraint(self):
        clauses = []
        for active_value in list(range(len(self.possible_terms))):

                # generate term in which only element i is on
                subclauses = []
                for i in list(range(len(self.possible_terms))):
                    if i == active_value:
                        subclauses.append("(%s%s = 1)" % (self.base_variable_name, i))
                    else:
                        subclauses.append("(%s%s = 0)" % (self.base_variable_name, i))

                clause = "(" + " and ".join(subclauses) + ")"
                clauses.append(clause)

        return "\t" + " or ".join(clauses) + ";\n"

    def dRealformat_constraint(self):
        clauses = []
        for active_value in list(range(len(self.possible_terms))):

                # generate term in which only element i is on
                subclauses = []
                for i in list(range(len(self.possible_terms))):
                    if i == active_value:
                        subclauses.append("(%s%s = 1)" % (self.base_variable_name, i))
                    else:
                        subclauses.append("(%s%s = 0)" % (self.base_variable_name, i))

                clause = "( and " + "".join(subclauses) + ")"
                clauses.append(clause)

        return "\t" + "( or " + "".join(clauses) + ") \n"



class CRNSketch:
    """
    Represents a sketch of a Chemical Reaction Network.
    """

    def __init__(self, reactions, optional_reactions, input_species=None):
        self.reactions = reactions
        self.optionalReactions = optional_reactions

        self.all_reactions = self.reactions[:]
        self.all_reactions.extend(self.optionalReactions)

        self.input_species = input_species

        self.species = self.getSpecies()
        self.real_species = self.getSpecies(include_inputs=False)
        self.input_species = list(set(self.species)-set(self.real_species))
        self.derivatives = set()

        self.choice_variables = set()
        self.lambda_variables = set()
        self.joint_choice_variables = set()
        self.record_decision_variables()

        self.t = symbols('t')

        for i, reaction in enumerate(self.optionalReactions):
            reaction.is_optional = True
            reaction.variable_name = "o" + str(i)

    def __repr__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.all_reactions]) + "\n]"

    def __str__(self):
        return "[" + '\n' + '\n'.join([str(x) for x in self.all_reactions]) + "\n]"

    def record_decision_variables(self):
        """
        Iterates over all reactions in self.all_reactions, populating the lists of decision variables
        ``self.joint_choice_variables``, ``self.lambda_variables`` and ``self.choice_variables``
        """

        for reaction in self.all_reactions:
            terms = reaction.reactants[:]
            terms.extend(reaction.products)

            choice_terms = []
            for term in terms:
                if isinstance(term, TermChoice):
                    self.joint_choice_variables.add(term)
                    choice_terms.extend(term.possible_terms)

            terms.extend(choice_terms)
            for term in terms:
                if isinstance(term, TermChoice):
                    continue
                if isinstance(term.species, LambdaChoice):
                    self.lambda_variables.add(term.species)

                if isinstance(term.coefficient, Choice):
                    self.choice_variables.add(term.coefficient)

    def generateAllTokens(self, covariance_matrix=None):
        """
        Returns a dictionary of all free_symbols that will appear in the system's rate equations.
        :param covariance_matrix: covariance matrix (needed only if LNA approximation is used)
        """

        sym = set()
        species_strings = [] # entries will be strings representing caluses that appear reactants/products
        for y in self.all_reactions:
            reactants_or_products = y.reactants[:]
            reactants_or_products.extend(y.products)
            for react in reactants_or_products:
                new_symbols = sympify(react.specRep()).free_symbols
                sym = sym.union(new_symbols)

        if covariance_matrix is not None:
            sym = sym.union(covariance_matrix.free_symbols)

        return sym

    def getSpecies(self, include_inputs=True):
        """
        Returns list of only those species (or InputSpecies) that participate in one or more reactions as either a
        reactant or product

        :param include_inputs: if True, include InputSpecies in list.
        """

        x = set()
        for y in self.all_reactions:
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
        """
        Generate list of all rate constants (needed when generating iSAT/dReach files)
        :return:
        """
        rate_constants = {}
        for reaction in self.all_reactions:
            for rate in reaction.get_rate_constants():
                if str(rate) not in list(rate_constants.keys()):
                    rate_constants[str(rate)] = rate
        rate_constants = list(rate_constants.values())
        return rate_constants

    def flow(self, isLNA, derivatives):
        """
        Returns a dictionary in which keys are species-names, and values are SympPy expressions for their time derivatives
        :param isLNA: if True, use the Linear Noise Approximation (LNA); if False use the Rate Reaction Equations (RRE)
        :param derivatives: list specifying which derivatives appear in the specification (each element is a dict)
        """

        if not derivatives:
            derivatives = set()
        self.derivatives = derivatives

        if isLNA:
            a = dict.fromkeys(self.generateAllTokens(generateCovarianceMatrix(self.real_species)))
        else:
            a = dict.fromkeys(self.generateAllTokens())

        propensities = self.parametricPropensity()
        stoichiometry_change = self.parametricNetReactionChange()
        dSpeciesdt = Matrix(stoichiometry_change).transpose() * Matrix(propensities)

        for i, sp in enumerate(self.real_species):
            if isinstance(sp, str):
                a[symbols(sp)] = dSpeciesdt[i]
            else:
                a[sp.symbol] = dSpeciesdt[i]

        if isLNA:
            jmat = [sp.symbol for sp in self.real_species]
            J = Matrix(dSpeciesdt).jacobian(jmat)
            C = generateCovarianceMatrix(self.real_species)
            dCovdt = J * C + C * transpose(J)
            for i in range(C.cols * C.rows):
                a[C[i]] = dCovdt[i]

        derivative_expressions = compute_derivatives(derivatives, a)
        a.update(derivative_expressions)

        for sp in self.input_species:
            a[sp.symbol] = sp.ode

        constants_to_remove = []
        for key in a:
            if a[key] is None and not isinstance(a[key], str):
                a[key] = 0
            if str(key) not in [str(sp) for sp in self.species] and not str(key).startswith("cov"):
                constants_to_remove.append(key)

        # remove constant keys from flowDict, as they are handled separately when output generated
        for key in constants_to_remove:
            if str(key) not in [x["name"] for x in derivatives]:
                a.pop(key, None)

        return a

    def get_cost(self):
        """
        Construct a SymPy expression for the topology cost in terms of the decision variables.

        The cost is defined as:
         3*(total number of species that participate in a reaction) +
          6*(sum of stoichiometric coefficients of all reactants) +
          5*(sum of stoichiometric coefficients of all products) +

        :return:
        """

        # Weighted sum of stoichiometric coefficients for each reaction
        change = []
        for reaction in self.all_reactions:
            netChange = ['0'] * len(self.real_species)
            for reactant in reaction.reactants:
                add_stoichiometry_change(self.real_species, netChange, reactant, '6*')

            for product in reaction.products:
                add_stoichiometry_change(self.real_species, netChange, product, '5*')

            change.append(sympify(netChange))

        # if a species is optional (i.e. only occurs in lamdaChoice), find all the lamda variables dning to its selections

        cost = sum([sum(x) for x in change])

        cost += 3 * len(self.real_species) # initially include all species, even optional ones that are not used

        # get list of all decision-variable constants
        constants = set()
        for cv in self.lambda_variables:
            for lam in cv.lambdas:
                constants.add(lam)
        for jcv in self.joint_choice_variables:
            for i in list(range(len(jcv.possible_terms))):
                constants.add("%s%s" % (jcv.base_variable_name, i))

        for species in self.real_species:
            species_index = self.real_species.index(species)
            stoich_sum = sum(Matrix(change)[species_index, :])

            # replace all constants in sum by zero: if result is non-zero the species is not optional
            for constant in constants:
                stoich_sum = stoich_sum.subs(sympify(constant), 0)
            if stoich_sum != 0:
                break

            # replace all non-constants by one and subtract constant to get condition under which species is absent
            stoich_prod = reduce(mul, Matrix(change)[species_index,:], 1)
            symbols = stoich_prod.free_symbols - constants
            for s in symbols:
                stoich_prod = stoich_prod.subs(s, 1)

            for integer in stoich_prod.atoms(S(2),0):
                stoich_prod = stoich_prod.subs(integer, 1)  # replace coefficients with 1
            stoich_prod = stoich_prod.expand()
            for integer in stoich_prod.atoms(S(2),0):
                stoich_prod = stoich_prod.subs(integer, 0)  # remove constant term

            cost -= 3*(1-stoich_prod)

        return cost.simplify()

    def parametricPropensity(self):
        """
        Returns a list of reaction propensities. The n'th element is a SymPy expression for the propensity of the n'th reaction
        """

        propensities = []
        for reaction in self.all_reactions:
            propensities.append(reaction.get_propensity())
        return propensities

    def parametricNetReactionChange(self):
        """
        Construct a matrix recording the change in stoichiometry of each species in each reaction.
        Returns a 2D list: change[reaction_index][species_index] is a string representing the stoichiometry change

        :return:
        """

        change = []
        for reaction in self.all_reactions:
            netChange = ['0'] * len(self.real_species)
            for reactant in reaction.reactants:
                add_stoichiometry_change(self.real_species, netChange, reactant, '-')
            for product in reaction.products:
                add_stoichiometry_change(self.real_species, netChange, product, '+')

            change.append(sympify(netChange))
        return change

    def simplify_expression(self, full_expression):
        f = full_expression.expand()

        for cv in self.choice_variables:
            for i in list(range(cv.minValue, cv.maxValue + 1)):
                for j in list(range(cv.minValue, cv.maxValue + 1)):
                    if i != j:
                        expression = "%s_%s * %s_%s" % (cv.name, i, cv.name, j)
                        f = f.subs(sympify(expression), 0)

        for lc in self.lambda_variables:
            for active_value in lc.lambdas:
                for lam in lc.lambdas:
                    if lam != active_value:
                        expression = "%s * %s" % (active_value, lam)
                        f = f.subs(sympify(expression), 0)

        for jc in self.joint_choice_variables:
            for active_value in list(range(len(jc.possible_terms))):
                for i in list(range(len(jc.possible_terms))):
                    if i != active_value:
                        expression = "%s%s * %s%s" % (jc.base_variable_name, i, jc.base_variable_name, active_value)
                        f = f.subs(sympify(expression), 0)

        return f


def add_stoichiometry_change(species, stoichiometry_change, fragment, sign):
    """
    Update the matrix of stoichiometry changes for each reaction. This allows the matrix to be constructed incrementally,
    by iterating over the reactants and products of each reaction.

    :param species: list of species
    :param stoichiometry_change: matrix recording the change in stoichiometry of each species in each reaction
    :param fragment: a reactant or product term (which may be a choice)
    :param sign: "+" or "-" (string)
    :return: udated stoichiometry_change matrix
    """
    for i, sp in enumerate(species):
        if str(sp) in fragment.specRep():
            if "lam" in fragment.specRep():
                new_term = "%s(%s)" % (sign, str(sympify(fragment.specRep()).expand().coeff(sp)))
            elif "tc_" in fragment.specRep():
                new_term = "%s(%s)" % (sign, str(sympify(fragment.specRep()).expand().coeff(sp)))
            else:
                new_term = "%s(%s)" % (sign, str(sympify(fragment.specRep()).expand().coeff(sp)))
            stoichiometry_change[i] = " + ".join([stoichiometry_change[i], new_term])
    return stoichiometry_change


def generateCovarianceMatrix(species_names):
    """
    Returns a SympPy Matrix that is a covariance matrix. Each entry is a SymPY Symbol whose name is determined by
     species_names.

    For example::

        >> generateCovarianceMatrix(['A', 'B', 'C'])
        >>   Matrix([
            [ covA, covAB, covAC],
            [covAB,  covB, covBC],
            [covAC, covBC,  covC]])

    :param species_names: list of names of species
    """
    mat = eye(len(species_names))

    for i, m in enumerate(species_names):
        for j, n in enumerate(species_names):
            if m == n:
                mat[i, j] = 'cov' + str(m)
            else:
                mat[i, j] = 'cov' + str(n) + str(m)

    for x in range(len(species_names)):
        for y in range(len(species_names)):
            mat[x, y] = mat[y, x]
    return mat


def compute_derivatives(derivatives, flow):
    """
    Returns a dictionary containing SymPy expressions for the requested derivatives of state variables.

    :param derivatives: list specifying which derivatives appear in the specification (each element is a dict)
    :param flow: dictionary in which keys are species-names, and values are SympPy expressions for their time derivatives
    """
    # define function for each state variable
    funcs = {}
    function_reverse = {}
    for variable in flow:
        new_function = Function(variable.name)(Symbol('t'))
        funcs[variable] = new_function
        function_reverse[new_function] = variable

    function_flowdict = {}
    constants = []
    for variable in flow:
        if flow[variable] is None:
            constants.append(variable.subs(funcs))
        else:
            function_flowdict[variable.subs(funcs)] = flow[variable].subs(funcs)

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
        for i in range(order - 1):
            xn = Derivative(xn, Symbol('t')).doit()

            # substitute in to replace first derivatives
            for func in function_flowdict:
                derivative_string = "Derivative(" + str(func) + ", t)"
                xn = xn.subs(sympify(derivative_string), function_flowdict[func])

            for constant in constants:
                derivative_string = "Derivative(" + str(constant) + ", t)"
                xn = xn.subs(sympify(derivative_string), 0)

        # replace functions with the corresponding symbols
        xn = xn.subs(function_reverse)
        results[symbols(name)] = simplify(xn)

    return results
