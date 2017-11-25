from CRNSynthesis.symbolicLNA import *
from sympy import init_printing, Matrix, transpose, pprint


def exampleParametricCRN():
    reaction1 = ReactionSketch(['A'], ['B'], ['K'], [], [2, 1], 'k_1', 0)
    reaction2 = ReactionSketch(['B'], ['B'], ['K'], ['K'], [2, 1], 'k_1', 0)

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

if __name__ == "__main__":
    exampleParametricCRN()
    init_printing()
    crn = exampleCRN()
    props = propensity(crn.reactions)
    nrc = netReactionChange(crn.species, crn.reactions)
    # props = propensity(crn.reactions)
    flow = flowFunction(props, nrc)
    # A = Matrix([crn.species])
    J = jacobian(flow, Matrix([crn.species]))
    # print covar
    G = g(Matrix(props), Matrix(nrc))
    C = generateCovarianceMatrix(crn.species)

    dCovdt = J * C + C * transpose(J) + G

    # iSATParser.constructiSATFile()
    print crn
    print "propensities:"
    print props

    pprint(dCovdt)
    # reactants = crn.reactions.reactants
    # reactants = [ x.reactants for x in crn.reactions]


    # eng = matlab.engine.start_matlab()
    # ret = eng.symbolicLNA(crn.species,[ x.reactants for x in crn.reactions], [ x.products for x in crn.reactions], [ x.reactionrate for x in crn.reactions])


    # inspecies, inreactants, inproducts, inrates