def addConstant(var, val):
    return "define " + var + ' = ' + str(val)

def addSpecies(species, parameters):
    return "float  [" + str(parameters[0]) + "," + str(parameters[1]) + "] " + species


def addStochConstant(stochConstant, parameters):
    decl = "float  [" + str(parameters[0]) + "," + str(parameters[1]) + "] " + stochConstant
    #constraints = '((' + parameters[0]  + ['((' + x + '= 1) or ( = 2));' for x in range(int(parameters[0]) + 1, int(parameters[1]) + 1)]
def addEquations(species, eq, file):
    pass

def constructiSATFile(reactionSketch, conditions, otherConstraints):
    pass