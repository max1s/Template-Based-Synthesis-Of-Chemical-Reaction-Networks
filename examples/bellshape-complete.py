from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCaller
from sympy import init_printing, Matrix, transpose, pprint


def formCRN():
    K = Species('K')
    A = Species('A')
    B = Species('B')

    lam1 = LambdaChoice([A, B], 1)
    lam2 = LambdaChoice([A, B], 2)
    c1 = Choice(1, 0, 2)
    c2 = Choice(2, 0, 2)
    c3 = Choice(3, 0, 2)
    c4 = Choice(4, 0, 2)
    c5 = Choice(5, 0, 1)
    c6 = Choice(6, 0, 1)
    c7 = Choice(7, 1, 2)

    reaction1 = Reaction([(lam1, 1), (K, c1)], [(K, c2)], RateConstant('k_1', 0, 1))
    reaction2 = Reaction([(lam2, c5), (K, c3)], [(lam2, c6), (K, c4)], RateConstant('k_2', 0, 1))
    reaction3 = Reaction([], [TermChoice(1, [(lam2, 1), (K, c7)])], RateConstant('k_3', 1, 1))

    return CRNSketch([reaction1, reaction2], [reaction3], [])

derivatives = [{"variable": 'K', "order": 1, "is_variance": False, "name": "dK_dt"}]
specification = [('', 'dK_dt >= 0', '((K > 0.3) and (dK_dt < 0.001))'), ('', 'dK_dt <= 0', '((K >= 0) and (K < 0.1))')]

crn = formCRN()
flow = crn.flow(False, derivatives)
hys = iSATParser.constructISAT(crn, specification, flow)

with open('bellshape.hys', 'w') as file:
    file.write(hys)

sc = SolverCaller("./bellshape.hys", isat_path="../isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt")

result_files = sc.single_synthesis(cost=0)

for file_name in result_files:
    print("\n\n")
    # print(sc.getCRNValues(file_name))

    vals, all_vals = sc.getCRNValues(file_name)
    initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

    print("Initial Conditions", initial_conditions)
    print("Flow:", parametrised_flow)

