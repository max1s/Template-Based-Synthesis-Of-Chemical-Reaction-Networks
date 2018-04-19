from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt

def formCRN():

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


    return CRNSketch([reaction1, reaction2, reaction3, reaction4, reaction5, reaction6], [], [input1])

derivatives = [{"variable": 'PThree', "order": 1, "is_variance": False, "name": "PThree_dot"},
               {"variable": 'PThree', "order": 2, "is_variance": False, "name": "PThree_dot_dot"}]
#specification = [('', 'PThree_dot >= 0', '((PThree > 0.1) and (PThree_dot < 0.001))'), ('', 'PThree_dot <= 0', '')]
specification = [('', '(PThree_dot >= 0)', '(PThree_dot = 0)' ), ('', '(PThree_dot < 0)', '(PThree_dot < 0)')]
crn = formCRN()
flow = crn.flow(False, derivatives)
hys = iSATParser.constructISAT(crn, specification, flow)
with open('sixreaction.hys', 'w') as file:
    file.write(hys)

specification_dreal = [('', 'PThree_dot >= 0', '(and (PThree > 0.1) (PThree_dot < 0.001))'), ('', 'PThree_dot <= 0', '(and (PThree >= 0)(PThree < 0.1))')]
drh = iSATParser.constructdReal(crn, specification_dreal, flow)
with open('sixreaction.drh', 'w') as file:
    file.write(drh)


# Try to solve using iSAT
sc = SolverCallerISAT("./sixreaction.hys", isat_path="../isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt")

result_files = sc.single_synthesis(cost=0)

for file_name in result_files:
    print("\n\n")
    # print(sc.getCRNValues(file_name))

    vals, all_vals = sc.getCRNValues(file_name)
    initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

    print("Initial Conditions", initial_conditions)
    print("Flow:", parametrised_flow)

    t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                   plot_name=file_name + "-simulation.png")
    print("\n\n")
    print(variable_names)
    print(sol)
    savetxt(file_name + "-simulation.csv", sol, delimiter=",")


# Try to solve using dReal
sc = SolverCallerDReal("./sixreaction.drh", dreal_path="../dReal-3.16.09.01-linux/bin/dReach")
result_files = sc.single_synthesis(cost=0)
