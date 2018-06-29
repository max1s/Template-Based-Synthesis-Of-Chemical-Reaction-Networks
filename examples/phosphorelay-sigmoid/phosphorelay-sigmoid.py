from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt

def form_crn():
    B = Species('B')
    L1 = Species('L1')
    L1p = Species('L1p')
    L2 = Species('L2')
    L2p = Species('L2p')
    L3 = Species('L3')
    L3p = Species('L3p')

    r1 = Reaction([(B, 1), (L1, 1)], [(B, 1), (L1p, 1)], RateConstant('k_1', 0, 1))
    r2 = Reaction([(L2, 1), (L1p, 1)], [(L2p, 1), (L1, 1)], RateConstant('k_2', 0, 1))
    r3 = Reaction([(L3, 1), (L2p, 1)], [(L3p, 1), (L2, 1)], RateConstant('k_3', 0, 1))
    r4 = Reaction([(L3p, 1)], [(L3, 1)], RateConstant('k_4', 0, 1))
    r5 = Reaction([], [(B, 1)], RateConstant('k_5', 1, 1))

    return CRNSketch([r1, r2, r3, r4, r5], [], [])

def synthesize_with_isat(crn):
    derivatives = [{"variable": 'L3p', "order": 1, "is_variance": False, "name": "L3p_dot"},
                   {"variable": 'L3p', "order": 2, "is_variance": False, "name": "L3p_dot_dot"}]
    flow = crn.flow(False, derivatives)

    specification = [('', '(L3p_dot >= 0) and (L3p_dot_dot >= 0)', '(L3p_dot_dot < 0.001)'),
                     ('', '(L3p_dot >= 0) and (L3p_dot_dot <= 0)', '(L3p > 100)')]

    hys = iSATParser.constructISAT(crn, specification, flow)
    with open('sigmoid.hys', 'w') as file:
        file.write(hys)

    sc = SolverCallerISAT("./sigmoid.hys", isat_path="../isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt")

    result_files = sc.single_synthesis(cost=0)


    for file_name in result_files:
        vals, all_vals = sc.getCRNValues(file_name)
        initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

        print("\n\nInitial Conditions", initial_conditions)
        print("Flow:", parametrised_flow)

        t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow)
        print("\n\n")
        print(variable_names)
        print(sol)
        savetxt(file_name + "-simulation.csv", sol, delimiter=",")


def synthesize_with_dreal(crn):
    derivatives = [{"variable": 'L3p', "order": 1, "is_variance": False, "name": "L3p_dot"},
                   {"variable": 'L3p', "order": 2, "is_variance": False, "name": "L3p_dot_dot"}]
    flow = crn.flow(False, derivatives)

    specification_dreal = [('', '(and (L3p_dot >= 0)(L3p_dot_dot >= 0))', '(L3p_dot_dot < 0.001)'),
                           ('', '(and (L3p_dot >= 0)(L3p_dot_dot <= 0))', '(L3p > 100)')]
    drh = iSATParser.constructdReal(crn, specification_dreal, flow)
    with open('sigmoid.drh', 'w') as file:
        file.write(drh)

    sc = SolverCallerDReal("./sigmoid.drh", dreal_path="../dReal-3.16.09.01-linux/bin/dReach")
    result_files = sc.single_synthesis(cost=0)

    for file_name in result_files:
        vals, all_vals = sc.getCRNValues(file_name)
        initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

        print("\n\nInitial Conditions", initial_conditions)
        print("Flow:", parametrised_flow)

        t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                       plot_name=file_name + "-simulationdreal.png")
        print("\n\n")
        print(variable_names)
        print(sol)
        savetxt(file_name + "-simulationdreal.csv", sol, delimiter=",")

if __name__ == "__main__":
    crn = form_crn()

    synthesize_with_isat(crn)
    synthesize_with_dreal(crn)
