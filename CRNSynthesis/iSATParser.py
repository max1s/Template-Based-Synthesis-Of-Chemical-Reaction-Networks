
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
        s += self.type + ' [' + self.minimumParameter + ", " + self.maximumParameter + ']' + self.name + ';'
        return s

    def constructdReal(self):
        pass

class Declaration:
    def __init__(self, decltype, decConstants):
        self.declarationOfParameter = decltype
        self.declarationOfConstants = decConstants

    def constructiSAT(self):
        s = "DECL \n"
        for constant in self.declarationOfConstants:
            "define " + constant.constantName + ' = ' + constant.constantValue + ';\n'
        s += '\n'
        for d in self.declarationOfParameter:
            s += d.constructiSAT() + '\n'
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
        for pair in self.speciesInitialValuePair:
            s += (pair.species + " " + pair.sign + " " +  pair.initial + ';\n')
        s.join(x + ';\n' for x in self.integerConstraints)
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
        s.join("(" + self.variable + "/" + self.time + " = " + self.flow + ")")
        return s

    def constructdReal(self):
        raise NotImplementedError

class Mode:
    def __init__(self, m, inv, fl):
        self.modeName = m
        self.invariants = inv
        self.flow = fl

    def constructiSAT(self):
        s = "TRANS\n"
        s.join((self.modeName + ' -> ' + x + ';\n' for x in self.invariants))
        s.join((self.modeName + ' -> ' + x.constructiSAT() + ';\n' for x in self.flow))
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
        s +=  self.mode + ' ' + "and (time <" + self.time + ") "
        s.join('and (' + x + ')' for x in self.specification)
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
        if specificationPart.time is maxtime:
            post.specification += specificationPart
        else:
            if(specificationPart.time not in timeList):
                noOfModes = noOfModes + 1
                modes += Mode(noOfModes, [specificationPart], flow)
                timeList += specificationPart.time
            else:
                modes[noOfModes].invariants += specificationPart
    return modes, post


def constructSpecification(specification, flow, declaration, integerConstraints=None, constants=None, initialValues=None, costFunction):
    m_flow = [Flow(x, 't', y) for x,y in flow]
    m_specification = [SpecificationPart(x, y) for x,y in specification]
    m_integerConstraints = [IntegerConstraint(x, y.min, y.max) for x,y in integerConstraints]
    m_decltypes = [DeclType(x, y[0], y[1], z) for x,y,z in declaration]
    m_contants = [Constant(x, y) for x,y in constants]

    d = Declaration(m_decltypes, m_contants)
    i = Initial(m_integerConstraints, None, costFunction)
    m,p = MTLConverter(specification, m_flow, 1)

    return d,i,m,p

def constructISAT(decl, initial, flow, post):
    d = decl.constructiSAT()
    i = initial.constructiSAT()
    f =[x.constructiSAT() for x in flow]
    p = post.constructiSAT()

    return d + i + flatten(f) + p




def constructdReal(decl, initial, flow, post):
    pass




