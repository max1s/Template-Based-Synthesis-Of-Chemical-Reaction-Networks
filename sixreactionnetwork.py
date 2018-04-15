from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCaller
from CRNSynthesis.solverCaller import SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint



def sixReactionNetwork():

    input1 = InputSpecies("Input1", sympify("0.1*t + 54.2735055776743*exp(-(0.04*t - 2.81375654916915)**2) + 35.5555607722356/(1.04836341039216e+15*(1/t)**10.0 + 1)"), 15)


    POne = Species('POne', initial_max=5)
    PTwo = Species('PTwo', initial_max=5)
    PThree = Species('PThree', initial_max=5)

    k1 = RateConstant('k_1', 0, 5)
    k2 = RateConstant('k_2', 0, 5)
    k3 = RateConstant('k_3', 0, 5)
    k4 = RateConstant('k_4', 0, 5)
    k5 = RateConstant('k_5', 0, 5)
    k6 = RateConstant('k_6', 0, 5)

    reaction1 = Reaction([Term(POne, 1)], [Term(PTwo, 1)], k1)
    reaction2 = Reaction([Term(PTwo, 1)], [Term(PThree, 1)], k2)
    reaction3 = Reaction([Term(POne, 1)], [Term(PThree, 1)], k3)
    reaction4 = Reaction([Term(PThree, 1)], [Term(PTwo, 1)], k4)
    reaction5 = Reaction([Term(PThree, 1)], [Term(POne, 1)], k5)
    reaction6 = Reaction([Term(PTwo, 1)], [Term(POne, 1)], k6)



    #bellshape
    specification = [(0, 'PThree_dot = 0'), (0.5, 'PThree_dot = 0.5'), (1, 'PThree_dot = 0')]

    isLNA = False
    derivatives = [{"variable": 'PThree', "order": 1, "is_variance": False, "name": "PThree_dot"},
                   {"variable": 'PThree', "order": 2, "is_variance": False, "name": "PThree_dot_dot"}]

    crn = CRNSketch([reaction1, reaction2, reaction3, reaction4, reaction5, reaction6], [], [input1])
    flow = crn.flow(isLNA, derivatives)
    crn.get_cost()

    return flow, iSATParser.constructdReal(crn, specification, flow)


def complete_process():

    flow, problem_string = sixReactionNetwork()
    with open("./sixreactionnetwork.hys", "w") as f:
        f.write(problem_string)

    sc = SolverCaller("./sixreactionnetwork.hys")
    result_file = sc.single_synthesis(cost=60, precision=0.1)
    print result_file
    #param_values = sc.getCRNValues(result_file)

    #print("\n\nSpecific CRN identified is:\n")
    #print(sc.get_parametrised_flow(flow, param_values))

def dReal_process(filename):
    problem_string = sixReactionNetwork()
    print problem_string
    with open("./" + filename, "w") as f:
        f.write(problem_string)

    sc = SolverCallerDReal(model_path="./" + filename)
    result_file = sc.single_synthesis(precision=0.1)
    param_values = sc.getCRNValues(result_file[0])
    print param_values
    quit()

if __name__ == "__main__":
    a,b = sixReactionNetwork()
    print b
   #print(complete_process())
    #complete_MM()
    # SolverCaller()