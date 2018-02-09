class Declaration:
    def __init__(self, crn, numModes):
        self.crn = crn
        self.numModes = numModes

    def constructiSAT(self):
        s = "\nDECL \n"

        s += "define MAX_TIME = 1;\n"  # TODO: set this sensibly

        s += "\t-- declare time variables\n"
        s += "\tfloat [0, MAX_TIME] time;\n"
        s += "\tfloat [0, MAX_TIME] delta_time;\n\n"

        if len(self.crn.real_species) > 0:
            s += "\n\t-- Define State Variables\n"
        for d in self.crn.real_species:
            s += d.iSATDefinition()

        if len(self.crn.input_species) > 0:
            s += "\n\t-- Define Input Variables\n"
        for d in self.crn.input_species:
            s += d.iSATDefinition()

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t-- Lambda Variables\n"
        for lam in self.crn.lambda_variables:
            s += lam.iSATDefinition() + "\n"

        if len(self.crn.choice_variables) > 0:
            s += "\n\t-- Choice Variables\n"
        for c in self.crn.choice_variables:
            s += c.iSATDefinition() + "\n"

        if len(self.crn.optionalReactions) > 0:
            s += "\n\t-- Optional Reaction Variables\n"
        for optional_reaction in self.crn.optionalReactions:
            s += "\tfloat [0, 1] %s;\n" % optional_reaction.variable_name

        if len(self.crn.getRateConstants()) > 0:
            s += "\n\t-- Rate constants\n"
        for rate in self.crn.getRateConstants():
            s += "\tfloat [%s, %s] %s;\n" % (rate.min, rate.max, rate.name)

        s += "\n\t--Define modes\n"
        for i in range(1, self.numModes + 1):
            s += "\tboole mode_%s;\n" % i

        return s

    def constructdReal(self):
        raise NotImplementedError


class Transition:
    def __init__(self, crn, flows, numModes):
        self.crn = crn

        self.flow = flows
        self.numModes = numModes

    def constructiSAT(self):
        s = "\n\nTRANS \n\n"

        s += "\t-- time constraint\n"
        s += "\ttime' = time + delta_time;\n\n"


        s += "\t-- No state change without time consumption.\n"

        terms = ["(%s' = %s)" % (x.name, x.name) for x in self.crn.real_species]
        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = ["(%s' = %s)" % (x.name, x.name) for x in self.crn.choice_variables]
        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = ["(%s' = %s)" % (x.name, x.name) for x in self.crn.getRateConstants()]
        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = ["(%s' = %s)" % (x.variable_name, x.variable_name) for x in self.crn.optionalReactions]
        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = []
        for lambda_choice in self.crn.lambda_variables:
            terms.extend(["(%s' = %s)" % (x, x) for x in lambda_choice.lambdas])
        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)
        s += "\n"

        if len(self.crn.getRateConstants()) > 0:
            s += "\n\t-- Rate constants are fixed\n"
        for rate in self.crn.getRateConstants():
            s += "\t(d.%s/d.time = 0);\n" % rate.name

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t-- Lambda variables are fixed\n"
        for lam in self.crn.lambda_variables:
            s += "".join(["\t(d.%s/d.time = 0);\n" % x.name for x in lam.lambdas])

        if len(self.crn.choice_variables) > 0:
            s += "\n\t-- Choice variables are fixed\n"
        for c in self.crn.choice_variables:
            s += "\t(d.%s/d.time = 0);\n" % c.name

        if len(self.crn.optionalReactions) > 0:
            s += "\n\t-- Optional reaction variables are fixed\n"
        for c in self.crn.optionalReactions:
            s += "\t(d.%s/d.time = 0);\n" % c.variable_name

        mode_list = ["mode_" + str(x) for x in range(1, self.numModes + 1)]
        modes_string = " or ".join(mode_list)

        s += "\n\n\t-- Flows\n"
        s += ''.join(['\t(%s) -> %s;\n' % (modes_string, x.constructiSAT()) for x in self.flow])

        return s


class Initial:
    def __init__(self, crn, numModes):
        self.crn = crn
        self.numModes = numModes

    def constructiSAT(self):
        s = "\nINIT\n"

        s += "\ttime = 0;\n\n"

        # mode exclusion
        mode_list = ["mode_%s" % x for x in range(1, self.numModes + 1)]
        s += "\t-- cannot be in two modes at the same time. We start in mode_1.\n"
        s += "\tmode_1 = 1;\n"
        s += "\t" + " + ".join(mode_list) + " = 1;\n\n"

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t-- Integer encoding of lambda variables\n"
        for lam in self.crn.lambda_variables:
            s += lam.format_constraint()

        if len(self.crn.choice_variables) > 0:
            s += "\n\t-- Integer encoding of choice variables\n"
        for c in self.crn.choice_variables:
            s += c.format_constraint()

        if len(self.crn.choice_variables) > 0:
            s += "\n\t-- Integer encoding of optional reaction variables\n"
        for optional_reaction in self.crn.optionalReactions:
            s += "\t(%s = 0) or (%s = 1);\n" % (optional_reaction.variable_name, optional_reaction.variable_name)



        if len(self.crn.real_species) > 0:
            s += "\n\t-- Limits on initial conditions\n"
        for sp in self.crn.real_species:
            s += sp.iSATInitialization()

        if len(self.crn.input_species) > 0:
            s += "\n\t-- Initial conditions of inputs\n"
        for sp in self.crn.input_species:
            s += sp.iSATInitialization()


        return s

    def constructdReal(self):
        raise NotImplementedError


class Flow:
    def __init__(self, var, t, fl):
        self.variable = var
        self.time = t
        self.flow = fl

    def constructiSAT(self):
        # Python represents powers as a**b, whereas iSAT uses a^b
        flow = str(self.flow).replace('**', '^')

        s = ""
        s += ("\t(d." + str(self.variable) + "/d." + str(self.time) + " = " + flow + ")")
        return s

    def constructdReal(self):
        raise NotImplementedError


class Mode:
    def __init__(self, m, inv):
        self.modeName = m
        self.invariants = inv

    def constructiSAT(self):

        s = "\n\n"

        for invariant in self.invariants:
            s += "\t-- transition into mode %s\n" % self.modeName
            if invariant[0] is not None:
                # constraint on mode start time
                s += "\t mode_%s -> (time >= %s);\n" % (self.modeName, invariant[0])
            # constraint imposed by mode on state
            s += "\t mode_%s -> (%s);\n" % (self.modeName, invariant[1])

        return s

    def constructdReal(self):
        raise NotImplementedError


class Post:
    def __init__(self, t, specs, mod):
        self.time = t
        self.specification = specs
        self.mode = mod

    def constructiSAT(self):
        s = ""
        s += "\nTARGET \n"
        s += "\tmode_" + str(self.mode) + ' ' + "and (time <" + str(self.time) + ") "
        s.join('and (' + str(x) + ')' for x in self.specification)
        s += ";"
        return s

    def constructdReal(self):
        raise NotImplementedError


class Specification:
    def __init__(self, d, i, m, p, add):
        self.declaration = d
        self.initial = i
        self.modes = m
        self.post = p
        self.additionalFlowConstraints = add


class SpecificationPart:
    def __init__(self, t, l):
        self.time = t
        self.logic = l


def MTLConverter(specification, flow, maxtime=1):
    post = Post(maxtime, [], 1)
    modes = []
    timeList = [maxtime]
    noOfModes = 1
    for specificationPart in specification:
        if specificationPart[0] is maxtime:
            post.specification += specificationPart
        else:
            if specificationPart[0] not in timeList:
                modes.append(Mode(noOfModes, [specificationPart]))
                noOfModes = noOfModes + 1
                timeList.append(specificationPart[0])
            else:
                modes[noOfModes].invariants += specificationPart
    return modes, post


def constructISAT(crn, specification, flow, costFunction=''):
    m_flow = [Flow(x, 'time', y) for x, y in flow.items()]
    m_specification = [SpecificationPart(x, y) for x, y in specification]
    numModes = min(1, len(specification))

    d = Declaration(crn, numModes).constructiSAT()
    i = Initial(crn, numModes).constructiSAT()
    t = Transition(crn, m_flow, numModes).constructiSAT()

    flow, post = MTLConverter(specification, m_flow, 1)
    f = [x.constructiSAT() for x in flow]
    p = post.constructiSAT()

    return d + i + t + ''.join(f) + p


def constructdReal(crn, specification, flow, costFunction=''):
    pass
