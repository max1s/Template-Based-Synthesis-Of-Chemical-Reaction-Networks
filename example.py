from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCaller
from CRNSynthesis.solverCaller import SolverCallerDReal
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


def exampleLotkaOscillator(solver='iSAT'):
    A = Species('A')
    B = Species('B')

    k1 = RateConstant('k_1', 0.01, 2)
    k2 = RateConstant('k_2', 0.01, 2)
    k3 = RateConstant('k_3', 0.01, 2)

    reaction1 = Reaction([Term(A, 1)],[Term(A, 2)], k1)
    reaction2 = Reaction([Term(A,1),Term(B,1)], [Term(B, 2)], k2)
    reaction3 = Reaction([Term(B,1)],[], k3)

    crn = CRNSketch([reaction1, reaction2, reaction3], [], [])

    isLNA = False
    derivatives = [{"variable": 'A', "order": 1, "is_variance": False, "name": "A_dot"},
                   {"variable": 'B', "order": 1, "is_variance": False, "name": "B_dot"}]

    specification = [(0.1, 'A_dot = 0'), (0.2, 'A_dot < 0 or A_dot > 0'), (0.3, 'A_dot = 0')]

    flow = crn.flow(isLNA, derivatives)
    if solver is 'iSAT':
        return iSATParser.constructISAT(crn, specification, flow)
    elif solver is 'dReal':
        return iSATParser.constructdReal(crn, specification, flow)


def exampleToggleSwitch():
    lacl = Species('lacl')
    clts = Species('clts')

    empty_1 = Species('empty_1')
    empty_2 = Species('empty_2')
    empty_3 = Species('empty_3')
    empty_4 = Species('empty_4')

    a_2 = RateConstant('a_2', 0.1, 2)
    a_4 = RateConstant('a_4', 0.1, 2)

    alpha_1 = RateConstant('alpha_1', 0.1, 5)
    alpha_2 = RateConstant('alpha_2', 0.1, 5)


    beta = RateConstant('beta', 0.1, 5)
    gamma = RateConstant('gamma', 0.1, 5)

   #a_1 = RateConstant(str(alpha_1.symbol / 1 + (clts.symbol ** beta.symbol)), 0, 5)
   # a_3 = RateConstant(str(alpha_2.symbol / 1 + (lacl.symbol ** gamma.symbol)), 0, 5)

    a_1 = RateConstant('a_1', 0, 5)
    a_3 = RateConstant('a_3', 0 ,5)
    input1 = InputSpecies("Input1", sympify(
        "0.1*t + 54.2735055776743*exp(-(0.04*t - 2.81375654916915)**2) + 35.5555607722356/(1.04836341039216e+15*(1/t)**10.0 + 1)"),
                          15)


    reaction1 = Reaction([Term(empty_1, 1)], [Term(lacl, 1)], a_1)
    reaction2 = Reaction([Term(lacl, 1)], [Term(empty_2,1)], a_2)
    reaction3 = Reaction([Term(empty_3, 1)],[Term(clts, 1)], a_3)
    reaction4 = Reaction([Term(clts, 1)], [Term(empty_4, 1)], a_4)

    crn = CRNSketch([reaction1, reaction2, reaction3, reaction4], [], [input1])

    isLNA = False
    derivatives = []
    specification = []

    flow = crn.flow(isLNA, derivatives)

    return iSATParser.constructdReal(crn, specification, flow)


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
    return iSATParser.constructdReal(crn, specification, flow)


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
    #derivatives = []
    specification = [(0, 'X = 0'), (0.5, 'X = 0.5'), (1, 'X = 0')]

    crn = CRNSketch([reaction1, reaction2, reaction3, reaction5], [reaction4], [input1])
    flow = crn.flow(isLNA, derivatives)
    crn.get_cost()

    return iSATParser.constructdReal(crn, specification, flow, costFunction='')


def exampleJointAlternative():

    X = Species('X', initial_max=5)
    Y = Species('Y', initial_value=12)
    B = Species('B')

    reaction1 = Reaction([TermChoice(0, [(X, 2), (Y, 3)])], [(B, 1)], RateConstant('k_1', 1, 2))

    isLNA = False
    derivatives = []
    specification = [(0, 'X = 0'), (0.5, 'X = 0.5'), (1, 'X = 0')]

    crn = CRNSketch([reaction1], [], [])
    flow = crn.flow(isLNA, derivatives)
    return flow, iSATParser.constructISAT(crn, specification, flow, costFunction='')


def bellshape_example():
    sc = SolverCaller("./bellshape.hys")
    result_files = sc.single_synthesis(cost=60)

    print("\n\n\n\n")
    for file_name in result_files:
        print("\n\n" +  sc.getCRNValues(file_name))


def hill_function_example():
    X = Species('X')
    Y = Species('Y')

    vmax = RateConstant('V_max', 0.1, 10)
    km = RateConstant('Km', 0.1, 10)
    n = RateConstant('n', 0.1, 10)

    #reaction = HillReaction([(X, 2)], [(Y,3)], vmax, km, n)
    #reaction = HillReaction([ TermChoice(0, [(X, 2), (Y, 2)])], [(Y,4)], vmax, km, n)

    reaction = HillReaction([ TermChoice(0, [(LambdaChoice([X, Y], 1), 2), (Y, 4)])], [(Y,5)], vmax, km, n)
    reaction = HillReaction([ (LambdaChoice([X, Y], 1), 2), (Y, 4)], [(Y,5)], vmax, km, n)

    derivatives = []

    crn = CRNSketch([reaction], [], [])
    flow = crn.flow(False, derivatives)
    crn.get_cost()

    return iSATParser.constructISAT(crn, [], flow, costFunction='')

def michaelis_menten_example():
    X = Species('X')
    Y = Species('Y')

    vmax = RateConstant('V_max', 0.1, 10)
    km = RateConstant('K_m', 0.1, 10)

    reaction = MichaelisMentenReaction([(X, 2)], [(Y,3)], vmax, km)

    derivatives = []

    crn = CRNSketch([reaction], [], [])
    flow = crn.flow(False, derivatives)
    crn.get_cost()

    return iSATParser.constructISAT(crn, [], flow, costFunction='')


def mixedMMExample():
    A = Species('A')
    g = Species('g')
    Ag = Species('Ag')
    Out = Species('Out')

    vAmax = RateConstant('VA_max', 0.5, 5)
    vBmax = RateConstant('VB_max', 0.5, 5)
    vConstant = RateConstant('vConstant', 0.5, 0.5)


    reaction1 = MichaelisMentenReaction([(A, 1), (g, 1)], [(Ag, 1), (A, 1)], vAmax, vConstant)
    reaction2 = MichaelisMentenReaction([(Ag, 1)], [(A, 1), (g,1)], vBmax, vConstant)
    reaction3 = Reaction([Term(Ag,1)], [Term(Out,1)] , RateConstant('k1', 0.1, 5))
    reaction4 = Reaction([Term(Out,1)], [], RateConstant('k2', 0.1, 5))

    isLNA = True
    derivatives = [{"variable": 'Out', "order": 1, "is_variance": False, "name": "dOut"},
                   {"variable": 'Out', "order": 2, "is_variance": False, "name": "ddOut"}]

    #Here specification should be when second derivative is 0 we want variance of O < certain threshold
    #specification = [(0, 'O = 0'), (0.5, 'O < 0.5'), (1, 'O > 1')]
    specification = [('', '', 'ddOut = 0'), ('', 'covAg > 10'), ('', 'Out > 1')]
    crn = CRNSketch([reaction1,reaction2,reaction3,reaction4],[],[])
    flow = crn.flow(isLNA, derivatives)
    return flow, iSATParser.constructdReal(crn, specification, flow)

def complete_process():

    flow, problem_string = exampleJointAlternative()
    with open("./test.hys", "w") as f:
        f.write(problem_string)

    sc = SolverCaller("./test.hys")
    result_file = sc.single_synthesis(cost=60, precision=0.1)
    param_values = sc.getCRNValues(result_file)

    print("\n\nSpecific CRN identified is:\n")
    print(sc.get_parametrised_flow(flow, param_values))

def dReal_process(filename):
    problem_string = exampleParametricCRN()
    print problem_string
    with open("./" + filename, "w") as f:
        f.write(problem_string)

    sc = SolverCallerDReal(model_path="./" + filename)
    result_file = sc.single_synthesis(precision=0.1)
    param_values = sc.getCRNValues(result_file[0])
    print param_values
    quit()
    #print("\n\nSpecific CRN identified is:\n")
    #print(sc.get_parametrised_flow(flow, param_values))


def complete_MM():

    flow, problem_string = mixedMMExample()
    print problem_string
    with open("./testMM.hys", "w") as f:
        f.write(problem_string)

    #sc = SolverCaller("./testMM.hys")
    #result_file = sc.single_synthesis(cost=100, precision=0.1)
    #param_values = sc.getCRNValues(result_file)

    #print("\n\nSpecific CRN identified is:\n")
    #print(sc.get_parametrised_flow(flow, param_values))


if __name__ == "__main__":
    # print(exampleCRN())
    # print(AMExample())
    # print(exampleParametricCRN())
     #print(exampleParametricCRN_complete())
     #print(exampleJointAlternative())
    #dReal_process('test.drh', )
    print exampleToggleSwitch()
    # print(complete_process())
    #complete_MM()
    # print(hill_function_example())
    #f,i = mixedMMExample()
    # print i
    # SolverCaller()