from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt
from numpy import linspace

def form_crn():
    B = Species('B', initial_value=0)
    L1 = Species('L1', initial_value=0.33)
    L1p = Species('L1p', initial_value=0)
    L2 = Species('L2', initial_value=0.33)
    L2p = Species('L2p', initial_value=0)
    L3 = Species('L3', initial_value=0.33)
    L3p = Species('L3p', initial_value=0)

    r1 = Reaction([(B, 1), (L1, 1)], [(B, 1), (L1p, 1)], RateConstant('k_1', 0.42, 0.42))
    r2 = Reaction([(L2, 1), (L1p, 1)], [(L2p, 1), (L1, 1)], RateConstant('k_2', 0.47, 0.47))
    r3 = Reaction([(L3, 1), (L2p, 1)], [(L3p, 1), (L2, 1)], RateConstant('k_3', 0.52, 0.52))
    r4 = Reaction([(L3p, 1)], [(L3, 1)], RateConstant('k_4', 0.162, 0.162))
    r5 = Reaction([], [(B, 1)], RateConstant('k_5', 1, 1))

    return CRNSketch([r1, r2, r3, r4, r5], [], [])

def synthesize_with_isat(crn):
    derivatives = [{"variable": 'L3p', "order": 1, "is_variance": False, "name": "L3p_dot"},
                   {"variable": 'L3p', "order": 2, "is_variance": False, "name": "L3p_dot_dot"}]
    flow = crn.flow(False, derivatives)

    #specification = [('', '(L3p_dot >= 0) and (L3p_dot_dot >= 0)', '(L3p_dot_dot < 0.001)'),
    #                 ('', '(L3p_dot >= 0) and (L3p_dot_dot <= 0)', '(L3p > 100)')]

    specification = [('','','')]

    flow = crn.flow(False, derivatives)
    hys = iSATParser.constructISAT(crn, specification, flow, max_time=100)
    with open('phosphorelay.hys', 'w') as file:
        file.write(hys)


    sc = SolverCallerISAT("./phosphorelay.hys",
                          isat_path="../isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt")

    result_files = sc.single_synthesis(cost=0, precision=0.01, msw=0.05)
    # result_files = ["./results/sixreactionstar_0_0.01-isat.txt"]

    for file_name in result_files:
        print("\n\n")
        # print(sc.getCRNValues(file_name))

        vals, all_vals = sc.getCRNValues(file_name)
        initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

        print("Initial Conditions", initial_conditions)
        print("Flow:", parametrised_flow)

        t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                       plot_name=file_name + "-simulation.png",
                                                       t=linspace(0, 100, 1000))
        print("\n\n")
        print(variable_names)
        print(sol)
        savetxt(file_name + "-simulation.csv", sol, delimiter=",")


def synthesize_with_dreal(crn):
    derivatives = [{"variable": 'L3p', "order": 1, "is_variance": False, "name": "L3p_dot"},
                   {"variable": 'L3p', "order": 2, "is_variance": False, "name": "L3p_dot_dot"}]
    flow = crn.flow(False, derivatives)

    #specification_dreal = [('', '(and (L3p_dot >= 0)(L3p_dot_dot >= 0))', '(L3p_dot_dot < 0.001)'),
    #                       ('', '(and (L3p_dot >= 0)(L3p_dot_dot <= 0))', '(L3p > 100)')]
    specification_dreal = [('','','')]
    drh = iSATParser.constructdReal(crn, specification_dreal, flow, max_time=100)
    with open('phosphorelay.drh', 'w') as file:
        file.write(drh)

    sc = SolverCallerDReal("./phosphorelay.drh", dreal_path="../dReal-3.16.09.01-linux/bin/dReach")
    result_files = sc.single_synthesis(cost=0, precision=0.1)


    for file_name in result_files:
            vals, all_vals = sc.getCRNValues('./phosphorelay_0_0.smt2.proof')
            initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

            print("Initial Conditions", initial_conditions)
            print("Flow:", parametrised_flow)
            t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                           plot_name=file_name + "-simulationdreal.png", t = linspace(0, 10, 100))
            print("\n\n")
            print(variable_names)
            print(sol)
            savetxt(file_name + "-simulationdreal.csv", sol, delimiter=",")
if __name__ == "__main__":
    crn = form_crn()

    #synthesize_with_isat(crn)
    synthesize_with_dreal(crn)
