from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCaller
from CRNSynthesis.solverCaller import SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint



def bellshape():
    A = Species('A', initial_max=0.1)
    B = Species('B', initial_value=1)
    K = Species('K', initial_value=1)
    o = Species('o', initial_value=1)

    k1 = RateConstant('k_1', 0.1, 5)
    k2 = RateConstant('k_2', 0.1, 5)
    k3 = RateConstant('k_3', 0.1, 5)

    reaction1 = Reaction([(LambdaChoice([A, B], 1), 1), Term(K, Choice(0, 0, 2))], [Term(K, Choice(1, 0, 2))], k1)
    reaction2 = Reaction([(LambdaChoice([A, B], 2), 1), Term(K, Choice(2, 0, 2)) ], [(LambdaChoice([A, B], 2), 1), Term(K, Choice(3, 0, 2))], k2)
    reaction3 = Reaction([Term(o, 1)], [Term(K, Choice(4,1,2))], k3)

    isLNA = False
    derivatives = [{"variable": 'K', "order": 1, "is_variance": False, "name": "K_dot"}]
    #derivatives = []
    #bellshape
    specification = [(0, 'K_dot = 0'), (0.5, 'K_dot = 0.5'), (1, 'K_dot = 0')]


    crn = CRNSketch([reaction1, reaction2], [reaction3], [])
    flow = crn.flow(isLNA, derivatives)
    crn.get_cost()

    return iSATParser.constructdReal(crn, specification, flow, costFunction='')


def complete_process():

    flow, problem_string = bellshape()
    with open("./bellshape.hys", "w") as f:
        f.write(problem_string)

    sc = SolverCaller("./test.hys")
    result_file = sc.single_synthesis(cost=60, precision=0.1)
    param_values = sc.getCRNValues(result_file)

    print("\n\nSpecific CRN identified is:\n")
    print(sc.get_parametrised_flow(flow, param_values))

def dReal_process(filename):
    problem_string = bellshape()
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
    print complete_process()
