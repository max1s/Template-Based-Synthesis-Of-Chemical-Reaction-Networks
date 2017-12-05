from CRNSynthesis.symbolicLNA import *
from sympy import init_printing, Matrix, transpose, pprint


def exampleCRN():
    reaction1 = Reaction(['A', 'x_1'], ['x_1', 'x_1'], 'k_1')
    reaction2 = Reaction(['x_1', 'x_2'], ['x_2', 'x_2'], 'k_2')
    reaction3 = Reaction(['x_2'], ['B'], 'k_3')
    return CRN(['A', 'B', 'x_1', 'x_2'], [reaction1, reaction2, reaction3], [5, 2, 1])


def AMExample():
    reaction1 = Reaction(['X', 'Y'], ['X', 'B'], 'k_1')
    reaction2 = Reaction(['Y', 'X'], ['Y', 'B'], 'k_2')
    reaction3 = Reaction(['X', 'B'], ['X', 'X'], 'k_3')
    reaction4 = Reaction(['Y', 'B'], ['Y', 'Y'], 'k_4')
    return CRN(['X', 'Y', 'B'], [reaction1, reaction2, reaction3, reaction4], [5, 3, 1])

def exampleParametricCRN():

    X = symbols('X')
    Y = symbols('Y')
    B = symbols('B')


    reaction1 = Reaction([Species(LambdaChoice([X, Y], 1), 1), Species(Y, 1)], [Species(X, 1), Species(B, 1)], 'k_1')
    reaction2 = Reaction([Species(LambdaChoice([X,Y], 2), 1), Species(Choice(X, 0, 2), 2)], [Species(Y, 1), Species(B, 1)], 'k_2')
    reaction3 = Reaction([Species(X, 1), Species(B, 1)], [Species(X, 1), Species(X, 1)], 'k_3')
    reaction4 = Reaction([Species(X, 1), Species(B, 1)], [Species(X, 1), Species(X, 1)], 'k_3')

    crn = CRNSketch([X, Y, B], [reaction1,reaction2,reaction3], [reaction4])

    prp = (parametricPropensity(crn))
    nrc = (parametricNetReactionChange(crn))
    dSpeciesdt = parametricFlow(prp, nrc)

    J = Matrix(dSpeciesdt).jacobian([X, Y, B])
    pprint(J)
    pprint(Matrix(prp))
    pprint(Matrix(nrc))
    G = parametricG(Matrix(prp), Matrix(nrc))

    C = generateCovarianceMatrix(['X','Y','B'])

    print (J * C  + C * transpose(J)).shape
    print G.shape

    dCovdt = J * C  + C * transpose(J)
    pprint (dCovdt)
    # nrc = netReactionChange(crn.species, crn.reactions)
    # #props = propensity(crn.reactions)
    # flow = flowFunction(props,nrc)
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
    # self.reactants = r
    # self.products = p
    # self.lambdaReactants = opr
    # self.lambdaProducts = opp
    # self.reactionrate = ra
    # self.isOptional = isop




def printCRNDetails(crn):
    #crn = exampleCRN()

    #crn = AMExample()
    print "CRN:"
    print crn

    props = propensity(crn.reactions)
    print "\nPropensities:"
    print props

    nrc = netReactionChange(crn.species, crn.reactions)
    flow = flowFunction(props, nrc)
    J = flow.jacobian(Matrix([crn.species]))
    G = g(Matrix(props), Matrix(nrc))
    C = generateCovarianceMatrix(crn.species)

    dCovdt = J * C + C * transpose(J) + G
    print "\ndCovdt:"
    init_printing()
    pprint(dCovdt)

    print "\n\n\n"
    # iSATParser.constructiSATFile()




if __name__ == "__main__":
    print "Example CRN:"
    printCRNDetails(exampleCRN())

    #print "AM Example:"
    #printCRNDetails(AMExample())

    print "New Parametric CRN Example:"
    exampleParametricCRN()


