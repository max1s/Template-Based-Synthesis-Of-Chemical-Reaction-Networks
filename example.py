from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from sympy import init_printing, Matrix, transpose, pprint

def exampleCRN():
    A = Species('A')
    B = Species('B')
    x1 = Species('x_1')
    x2 = Species('x_2')

    reaction1 = Reaction([Term(A,1), Term(x1,1)], [Term(x1,1), Term(x1,1)], RateConstant('k_1', 1, 2))
    reaction2 = Reaction([Term(x1,1), Term(x2,1)], [Term(x2,2)], RateConstant('k_2', 1, 2))
    reaction3 = Reaction([Term(x2,1)], [Term(B,1)], RateConstant('k_3', 1, 2))

    crn = CRNSketch([reaction1, reaction2, reaction3], [], [])
    isLNA = False
    derivatives = []
    specification = []

    flow = crn.flow(isLNA, derivatives)
    return iSATParser.constructISAT(crn, specification, flow, costFunction='')


def AMExample():
    X = Species('X')
    Y = Species('Y')
    B = Species('B')


    reaction1 = Reaction([Term(X,1), Term(Y,1)], [Term(X,1), Term(B,1)], RateConstant('k_1', 0.1, 10))
    reaction2 = Reaction([Term(Y,1), Term(X,1)], [Term(Y,1), Term(B,1)], RateConstant('k_2', 0.1, 10))
    reaction3 = Reaction([Term(X,1), Term(B,1)], [Term(X,1), Term(X,1)], RateConstant('k_3', 0.1, 10))
    reaction4 = Reaction([Term(Y,1), Term(B,1)], [Term(Y,1), Term(Y,1)], RateConstant('k_4', 0.1, 10))

    crn = CRNSketch([reaction1, reaction2, reaction3, reaction4], [], [])
    isLNA = False
    derivatives = []
    specification = []

    flow = crn.flow(isLNA, derivatives)
    return iSATParser.constructISAT(crn, specification, flow, costFunction='')


def exampleParametricCRN():

    X = Species('X')
    Y = Species('Y')
    B = Species('B')

    k1 = RateConstant('k_1', 0.1, 10)
    k2 = RateConstant('k_2', 0.1, 10)
    k3 = RateConstant('k_3', 0.1, 10)


    reaction1 = Reaction([Term(LambdaChoice([X, Y], 1), 1), Term(Y, 1)], [Term(X, 1), Term(B, 1)], k1)
    reaction2 = Reaction([Term(LambdaChoice([X,Y], 2), 1), Term(X, Choice(0, 0, 2))], [Term(Y, 1), Term(B, 1)], k2)
    reaction3 = Reaction([Term(X, 1), Term(B, 1)], [Term(X, 1), Term(X, 1)], k3)
    reaction4 = Reaction([Term(X, 1), Term(B, 1)], [Term(X, 1), Term(X, 1)], k3)

    crn = CRNSketch([reaction1,reaction2,reaction3], [reaction4], [])

    isLNA = False
    derivatives = []
    specification = []

    flow = crn.flow(isLNA, derivatives)
    return iSATParser.constructISAT(crn, specification, flow, costFunction='')


def exampleParametricCRN_complete():
    X = Species('X', initial_max=5)
    Y = Species('Y', initial_value=12)
    B = Species('B')

    reaction1 = Reaction([(LambdaChoice([X, Y], 1), 1), (Y, 1)], [(X, 1), (B, 1)],
                         RateConstant('k_1', 1, 2))
    reaction2 = Reaction([(LambdaChoice([X, Y], 2), 1), (X, Choice(0, 0, 3))],
                         [(Y, 1), (B, 1)], RateConstant('k_2', 1, 2))
    reaction3 = Reaction([(X, 1), (B, 1)], [(X, 1), (X, 1)], RateConstant('k_3', 1, 2))
    reaction4 = Reaction([(X, 1), (B, 1)], [(X, 1), (X, 1)], RateConstant('k_4', 1, 2))

    input1 = InputSpecies("Input1", sympify("0.1*t + 54.2735055776743*exp(-(0.04*t - 2.81375654916915)**2) + 35.5555607722356/(1.04836341039216e+15*(1/t)**10.0 + 1)"), 15)
    reaction5 = Reaction([Term(input1, 1)], [Term(B, 1)], RateConstant('k_input', 1, 2))

    isLNA = False
    derivatives = [{"variable": 'X', "order": 1, "is_variance": False, "name": "X_dot"},
                   {"variable": 'X', "order": 2, "is_variance": False, "name": "X_dot_dot"},
                   {"variable": 'X', "order": 2, "is_variance": True, "name": "covX_dot_dot"}]
    derivatives = []
    specification = [(0, 'X = 0'), (0.5, 'X = 0.5'), (1, 'X = 0')]

    crn = CRNSketch([reaction1, reaction2, reaction3, reaction5], [reaction4], [input1])
    flow = crn.flow(isLNA, derivatives)
    return iSATParser.constructISAT(crn, specification, flow, costFunction='')



if __name__ == "__main__":
    # print exampleCRN()
    # print AMExample()
    # print exampleParametricCRN()
    print exampleParametricCRN_complete()

