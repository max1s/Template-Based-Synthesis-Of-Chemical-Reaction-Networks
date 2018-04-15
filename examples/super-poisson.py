from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt

def formCRN():
    A = Species('A')
    B = Species('B')

    c1 = Choice(1, 0, 2)
    c2 = Choice(2, 0, 2)
    c3 = Choice(3, 0, 2)

    lam1 = LambdaChoice([A, B], 1)
    lam2 = LambdaChoice([A, B], 1)

    r1 = Reaction([], [(A, c1), (lam1, c2)], RateConstant('k_1', 0, 1))
    r2 = Reaction([(A, 1)], [(lam2, c3)], RateConstant('k_2', 0, 1))

    return CRNSketch([r1, r2], [], [])

derivatives = []

specification = [('', '(covA > A)', '')]

crn = formCRN()
flow = crn.flow(True, derivatives)
hys = iSATParser.constructISAT(crn, specification, flow)
with open('super-poisson.hys', 'w') as file:
    file.write(hys)

specification_dreal = [('', '(covA > A)', '')]
drh = iSATParser.constructdReal(crn, specification_dreal, flow)
with open('super-poisson.drh', 'w') as file:
    file.write(drh)


# Try to solve using iSAT
sc = SolverCallerISAT("./super-poisson.hys", isat_path="../isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt")

result_files = sc.single_synthesis(cost=0)

for file_name in result_files:
    print("\n\n")
    # print(sc.getCRNValues(file_name))

    vals, all_vals = sc.getCRNValues(file_name)
    initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

    print("Initial Conditions", initial_conditions)
    print("Flow:", parametrised_flow)

    t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow)
    print("\n\n")
    print(variable_names)
    print(sol)
    savetxt(file_name + "-simulation.csv", sol, delimiter=",")


# Try to solve using dReal
sc = SolverCallerDReal("./super-poisson.drh", dreal_path="../dReal-3.16.09.01-linux/bin/dReach")
result_files = sc.single_synthesis(cost=0)
