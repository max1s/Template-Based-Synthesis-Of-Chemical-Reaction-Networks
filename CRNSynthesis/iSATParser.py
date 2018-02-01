
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
    def __init__(self, decltype, decConstants, numModes):
        self.declarationOfParameter = decltype
        self.declarationOfConstants = decConstants
        self.numModes = numModes

    def constructiSAT(self):
        s = "DECL \n"
        if self.declarationOfConstants is not 0:
            for constant in self.declarationOfConstants:
                "define " + constant.constantName + ' = ' + constant.constantValue + ';\n'
            s += '\n'

        for d in self.declarationOfParameter:
            s += d.constructiSAT() + '\n'

        for i in range(1, self.numModes+1):
            s += "boole mode_%s;\n" % i

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
            s = "("
            for i in range(minInt, maxInt-1):
                s += "(" + self.variableName + " = " + str(i) + ") or "
            s += "(" + self.variableName + " = " + str(maxInt) + "); \n"
            return s

        def constructdReal(self):
            raise NotImplementedError



class Initial:
    def __init__(self, spinitvalpair, integerC=None, costC=None):
        self.speciesInitialValuePair = spinitvalpair
        self.integerConstraints = integerC
        self.costConstraints = costC

    def constructiSAT(self):
        s = "INIT\n"
        if self.speciesInitialValuePair is not 0:
            for pair in self.speciesInitialValuePair:
                s += (pair.species + " " + pair.sign + " " +  pair.initial + ';\n')
        s.join(x + ';\n' for x in self.integerConstraints) if self.integerConstraints is not None else 0
        s.join(self.costConstraints)
        return s

    def constructdReal(self):
        raise NotImplementedError

class Flow:
    def __init__(self, var, t, fl):
        self.variable = var
        self.time = t
        self.flow = fl

    def constructiSAT(self):
        s = ""
        s += ("(" + str(self.variable) + "/" + str(self.time) + " = " + str(self.flow) + ")")
        return s

    def constructdReal(self):
        raise NotImplementedError

class Mode:
    def __init__(self, m, inv, fl):
        self.modeName = m
        self.invariants = inv
        self.flow = fl

    def constructiSAT(self):
        s = "TRANS \n "
        s += ''.join(['mode_' + str(self.modeName) + ' -> ' + str(x) + ';\n' for x in self.invariants])
        s += ''.join(['mode_' + str(self.modeName) + ' -> ' + x.constructiSAT() + ';\n' for x in self.flow])
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
        s +="TARGET \n"
        s +=  str(self.mode) + ' ' + "and (time <" + str(self.time)+ ") "
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


def constructSpecification(specification, flow, declaration, costFunction, integerConstraints=None, constants=None, initialValues=None):
    m_flow = [Flow(x, 't', y) for x,y in flow.iteritems()]
    m_specification = [SpecificationPart(x, y) for x,y in specification]
    m_integerConstraints = [IntegerConstraint(x, y.min, y.max) for x,y in integerConstraints] if integerConstraints is not None else 0
    m_decltypes = [DeclType(x, 0, y, 'float') for x,y in declaration.iteritems()]
    m_contants = [Constant(x, y) for x,y in constants] if constants is not None else 0

    d = Declaration(m_decltypes, m_contants, len(specification))
    i = Initial(m_integerConstraints, None, costFunction)
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




