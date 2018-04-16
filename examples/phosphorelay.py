from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCaller
from CRNSynthesis.solverCaller import SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint



def exampleParametricCRN_complete():
    L1 = Species('L1')
    L2 = Species('L2')
    B = Species('B')
    L1p = Species('L1p')
    L2p = Species('L2p')
    L3p = Species('L3p')

    k1 = RateConstant('k_1', 0.1, 5)
    k2 = RateConstant('k_2', 0.1, 5)
    k3 = RateConstant('k_3', 0.1, 5)
    k4 = RateConstant('k_4', 0.1, 5)

    reaction1 = Reaction([Term(L1, 1), Term(B, 1)], [Term(B,1), Term(L1p, 1)], k1)
    reaction2 = Reaction([Term(L2, 1), Term(L1p, 1)], [Term(L1, 1), Term(L2p, 1)], k2)
    reaction3 = Reaction([Term(L2p, 1), Term(L3, 1)], [Term(L2, 1), Term(L3p, 1)], k3)
    reaction4 = Reaction([Term(L3p, 1)], [Term(L3, 1), k4])

    isLNA = True
    derivatives = [{"variable": 'L3p', "order": 1, "is_variance": False, "name": "L3p_dot"},
                   {"variable": 'L3p', "order": 2, "is_variance": False, "name": "L3p_dot_dot"}
                   ]

    specification = [(0, 'L3p_dot_dot > 0'), (0.5, 'L3p_dot_dot = 0'), (1, 'L3p_dot_dot < 0')]

    crn = CRNSketch([reaction1, reaction2, reaction3, reaction4], [], [])
    flow = crn.flow(isLNA, derivatives)
    crn.get_cost()

    return iSATParser.constructdReal(crn, specification, flow, costFunction='')


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