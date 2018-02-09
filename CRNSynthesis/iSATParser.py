class Constant:
    def __init__(self, name, value):
        self.constantName = name
        self.constantValue = value


class DeclType:
    def __init__(self, name, minparam, maxparam, declaration_type):
        self.name = name
        self.minimumParameter = minparam

        self.maximumParameter = maxparam
        if not maxparam:
            self.maximumParameter = 10

        self.type = declaration_type

    def constructiSAT(self):
        return "%s [%s, %s] %s;" % (self.type, self.minimumParameter, self.maximumParameter, self.name)

    def constructdReal(self):
        pass


class Declaration:
    def __init__(self, crn, decltype, decConstants, numModes, reactionRates):
        self.crn = crn
        self.declarationOfParameter = decltype
        self.declarationOfConstants = decConstants
        self.numModes = numModes
        self.reactionRates = reactionRates

    def constructiSAT(self):
        s = "\nDECL \n"

        s += "define MAX_TIME = 1;\n"  # TODO: set this sensibly

        s += "\t-- declare time variables\n"
        s += "\tfloat [0, MAX_TIME] time;\n"
        s += "\tfloat [0, MAX_TIME] delta_time;\n\n"

        if self.declarationOfConstants is not 0:
            s += "\t-- Define constants\n"
            for constant in self.declarationOfConstants:
                s += "\tdefine " + constant.constantName + ' = ' + constant.constantValue + ';\n'

        if len(self.declarationOfParameter) > 0:
            s += "\n\t-- Define State Variables\n"
        for d in self.declarationOfParameter:
            s += "\t" + d.constructiSAT() + '\n'

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t-- Lambda Variables\n"
        for lam in self.crn.lambda_variables:
            s += lam.iSATDefinition() + "\n"

        if len(self.reactionRates) > 0:
            s += "\n\t-- Rate constants\n"
        for rate in self.reactionRates:
            s += "\tfloat [%s, %s] %s;\n" % (rate.min, rate.max, rate.name)

        s += "\n\t--Define modes\n"
        for i in range(1, self.numModes + 1):
            s += "\tboole mode_%s;\n" % i

        return s

    def constructdReal(self):
        raise NotImplementedError


class Transition:
    def __init__(self, crn, decltype, decConstants, reactionRates, flows, numModes):
        self.crn = crn
        self.declarationOfParameter = decltype
        self.declarationOfConstants = decConstants
        self.reactionRates = reactionRates

        self.flow = flows
        self.numModes = numModes

    def constructiSAT(self):
        s = "\n\nTRANS \n\n"

        s += "\t-- time constraint\n"
        s += "\ttime' = time + delta_time;\n\n"

        terms = ["(%s' = %s)" % (x.name, x.name) for x in self.declarationOfParameter]
        s += "\t-- No state change without time consumption.\n"
        s += "\t(delta_time = 0) -> (%s);" % " and ".join(terms)
        s += "\n\n"

        if self.declarationOfConstants is not 0:
            s += "\t-- Constants are fixed\n"
            for constant in self.declarationOfConstants:
                s += "\t(d.%s/d.time = 0);\n" % constant.constantName

        if len(self.reactionRates) > 0:
            s += "\n\t-- Rate constants are fixed\n"
        for rate in self.reactionRates:
            s += "\t(d.%s/d.time = 0);\n" % rate.name

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t-- Lambda variables are fixed\n"
        for lam in self.crn.lambda_variables:
            s += "".join(["\t(d.%s/d.time = 0);\n" % x.name for x in lam.lambdas])

        mode_list = ["mode_" + str(x) for x in range(1, self.numModes + 1)]
        modes_string = " or ".join(mode_list)

        s += "\n\n\t-- Flows\n"
        s += ''.join(['\t(%s) -> %s;\n' % (modes_string, x.constructiSAT()) for x in self.flow])

        return s


class InitVals:
    def __init__(self, sp, init, sgn):
        self.species = sp
        self.initial = init
        self.sign = sgn


class IntegerConstraint:
    def __init__(self, var, minInt, maxInt):
        self.variableName = var
        self.minimumInteger = minInt
        self.maximumInteger = maxInt

    def constructiSAT(self):
        s = "\t("
        for i in range(self.minimumInteger, self.maximumInteger - 1):
            s += "(" + self.variableName + " = " + str(i) + ") or "
        s += "(" + self.variableName + " = " + str(self.maximumInteger) + "); \n"
        return s

    def constructdReal(self):
        raise NotImplementedError


class Initial:
    def __init__(self, crn, spinitvalpair, numModes, integerC=None, costC=None):
        self.crn = crn
        self.speciesInitialValuePair = spinitvalpair
        self.integerConstraints = integerC
        self.costConstraints = costC
        self.numModes = numModes

    def constructiSAT(self):
        s = "\nINIT\n"

        s += "\ttime = 0;\n\n"

        if self.speciesInitialValuePair is not 0:
            for pair in self.speciesInitialValuePair:
                s += (pair.species + " " + pair.sign + " " + pair.initial + ';\n')
        s.join(x + ';\n' for x in self.integerConstraints) if self.integerConstraints is not None else 0
        s.join(self.costConstraints)

        # mode exclusion
        mode_list = ["mode_%s" % x for x in range(1, self.numModes + 1)]
        s += "\t-- cannot be in two modes at the same time. We start in mode_1.\n"
        s += "\tmode_1 = 1;\n"
        s += "\t" + " + ".join(mode_list) + " = 1;\n\n"

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t-- Integer encoding of lambda variables\n"
        for lam in self.crn.lambda_variables:
            s += lam.format_constraint()


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


def constructSpecification(crn, specification, flow, declaration, costFunction, integerConstraints=None, constants=None,
                           initialValues=None, rate_constants=None):
    m_flow = [Flow(x, 'time', y) for x, y in flow.items()]
    m_specification = [SpecificationPart(x, y) for x, y in specification]
    m_integerConstraints = [IntegerConstraint(x, y.min, y.max) for x, y in
                            integerConstraints] if integerConstraints is not None else 0

    #m_decltypes = [DeclType(x, 0, y, 'float') for x, y in declaration.items()]

    m_contants = [Constant(x, y) for x, y in constants] if constants is not None else 0

    numModes = len(specification)
    d = Declaration(crn, [], m_contants, numModes, rate_constants)
    i = Initial(crn, m_integerConstraints, numModes, None, costFunction)
    t = Transition(crn, [], m_contants, rate_constants, m_flow, numModes)
    m, p = MTLConverter(specification, m_flow, 1)

    return d, i, m, p, t


def constructISAT(decl, initial, flow, post, transitionStart):
    d = decl.constructiSAT()
    i = initial.constructiSAT()

    t = transitionStart.constructiSAT()
    f = [x.constructiSAT() for x in flow]
    p = post.constructiSAT()
    return d + i + t + ''.join(f) + p


def constructdReal(decl, initial, flow, post):
    pass
