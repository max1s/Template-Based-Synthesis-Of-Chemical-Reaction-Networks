
class Constant:
    def __init__(self, name, value):
        self.constantName = name
        self.constantValue = value

class DeclType:
    def __init__(self, name, minparam, maxparam, type):
        self.name = name
        self.minimumParameter = minparam
        self.maximumParameter = maxparam
        self.type = type

    def constructiSAT(self):
        s = ''
        s += self.type + ' [' + str(self.minimumParameter) + ", " + str(self.maximumParameter) + '] ' + str(self.name) + ';'
        return s

    def constructdReal(self):
        pass

class Declaration:
    def __init__(self, decltype, decConstants, numModes, reactionRates):
        self.declarationOfParameter = decltype
        self.declarationOfConstants = decConstants
        self.numModes = numModes
        self.reactionRates = reactionRates

    def constructiSAT(self):
        s = "\nDECL \n"

        s += "define MAX_TIME = 1;\n" # TODO: set this sensibly

        s += "\t-- declare time variables\n"
        s += "\tfloat [0, MAX_TIME] time;\n"
        s += "\tfloat [0, MAX_TIME] delta_time;\n\n"

        if self.declarationOfConstants is not 0:
            for constant in self.declarationOfConstants:
                s += "\tdefine " + constant.constantName + ' = ' + constant.constantValue + ';\n'
            s += '\n'

        for d in self.declarationOfParameter:
            s += "\t" + d.constructiSAT() + '\n'

        s += "\n\t-- Rate constants\n"
        for rate in self.reactionRates:
            s += "\tdefine float %s [%s, %s];\n" % (rate.name, rate.min, rate.max)

        s += "\n"
        for i in range(1, self.numModes+1):
            s += "\tboole mode_%s;\n" % i

        return s

    def constructdReal(self):
        raise NotImplementedError


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
            for i in range(minInt, maxInt-1):
                s += "(" + self.variableName + " = " + str(i) + ") or "
            s += "(" + self.variableName + " = " + str(maxInt) + "); \n"
            return s

        def constructdReal(self):
            raise NotImplementedError



class Initial:
    def __init__(self, spinitvalpair, numModes, integerC=None, costC=None):
        self.speciesInitialValuePair = spinitvalpair
        self.integerConstraints = integerC
        self.costConstraints = costC
        self.numModes = numModes

    def constructiSAT(self):
        s = "\nINIT\n"

        s += "\ttime = 0;\n\n"

        if self.speciesInitialValuePair is not 0:
            for pair in self.speciesInitialValuePair:
                s += (pair.species + " " + pair.sign + " " +  pair.initial + ';\n')
        s.join(x + ';\n' for x in self.integerConstraints) if self.integerConstraints is not None else 0
        s.join(self.costConstraints)

        # mode exclusion
        mode_list = ["mode_%s" % x for x in range(1, self.numModes+1)]
        s += "\t-- cannot be in two modes at the same time. We start in mode_1.\n"
        s += "\tmode_1 = 1;\n"
        s += "\t" + " + ".join(mode_list) + " = 1;\n\n"
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
    def __init__(self, m, inv, fl):
        self.modeName = m
        self.invariants = inv
        self.flow = fl

    def constructiSAT(self):
        s = "\nTRANS \n "

        s += "\t-- time constraint\n"
        s += "\ttime' = time + delta_time;\n\n"

        for invariant in self.invariants:
             if invariant[0] is not None:
                    # constraint on mode start time
                   s += "\t mode_%s -> (time >= %s);\n" % (self.modeName, invariant[0])
             # constraint imposed by mode on state
             s += "\t mode_%s -> (%s);\n" % (self.modeName, invariant[1])

        s += ''.join(['\tmode_' + str(self.modeName) + ' -> ' + x.constructiSAT() + ';\n' for x in self.flow])
        return s

    def constructdReal(self):
        raise NotImplementedError

class Post:
    def __init__(self, t, specs, mod):
        self.time = t
        self.specification =  specs
        self.mode = mod

    def constructiSAT(self):
        s = ""
        s +="\nTARGET \n"
        s += "\t" + str(self.mode) + ' ' + "and (time <" + str(self.time)+ ") "
        s.join('and (' + str(x) + ')' for x in self.specification)
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
            if(specificationPart[0] not in timeList):
                modes.append(Mode(noOfModes, [specificationPart], flow))
                noOfModes = noOfModes + 1
                timeList.append(specificationPart[0])
            else:
                modes[noOfModes].invariants += specificationPart
    return modes, post


def constructSpecification(specification, flow, declaration, costFunction, integerConstraints=None, constants=None, initialValues=None, rate_constants=None):
    m_flow = [Flow(x, 'time', y) for x,y in flow.iteritems()]
    m_specification = [SpecificationPart(x, y) for x,y in specification]
    m_integerConstraints = [IntegerConstraint(x, y.min, y.max) for x,y in integerConstraints] if integerConstraints is not None else 0
    m_decltypes = [DeclType(x, 0, y, 'float') for x,y in declaration.iteritems()]


    m_contants = [Constant(x, y) for x,y in constants] if constants is not None else 0

    numModes = len(specification)
    d = Declaration(m_decltypes, m_contants, numModes, rate_constants)
    i = Initial(m_integerConstraints, numModes, None, costFunction)
    m,p = MTLConverter(specification, m_flow, 1)

    return d,i,m,p

def constructISAT(decl, initial, flow, post):
    d = decl.constructiSAT()
    i = initial.constructiSAT()
    f =[x.constructiSAT() for x in flow]
    p = post.constructiSAT()
    return d + i + ''.join(f) + p


def constructdReal(decl, initial, flow, post):
    pass




