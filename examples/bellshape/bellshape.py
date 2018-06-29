from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt
from numpy import linspace

def form_crn():
    K = Species('K', initial_value=0.001)
    A = Species('A', initial_value=0.5)
    B = Species('B', initial_value=0.001)

    lam1 = LambdaChoice([A, B], 1)
    lam2 = LambdaChoice([A, B], 2)
    c1 = Choice(1, 0, 2)
    c2 = Choice(2, 0, 2)
    c3 = Choice(3, 0, 2)
    c4 = Choice(4, 0, 2)
    c5 = Choice(5, 0, 1)
    c6 = Choice(6, 0, 2)
    c7 = Choice(7, 1, 2)

    reaction1 = Reaction([(A, 1), (K, 1)], [(K, 2)], RateConstant('k_1', 0.55, 0.57))
    reaction2 = Reaction([(B, 1), (K, 1)], [(B, 2)], RateConstant('k_2', 0.42, 0.44))
    reaction3 = Reaction([], [TermChoice(1, [(lam2, 1), (K, c7)])], RateConstant('k_3', 0, 0))

    return CRNSketch([reaction1, reaction2], [], [])

def synthesize_with_isat(crn):

    derivatives = [{"variable": 'K', "order": 1, "is_variance": False, "name": "dK_dt"}]
    specification = [('','','')]
    specification = [('(K < 0.1)', 'dK_dt >= 0', '(dK_dt = 0) and (K > 0.4)'), ('', 'dK_dt <= 0', '')]

    flow = crn.flow(False, derivatives)
    hys = iSATParser.constructISAT(crn, specification, flow, max_time=100)
    with open('bellshape.hys', 'w') as file:
        file.write(hys)


    sc = SolverCallerISAT("./bellshape.hys",
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
    derivatives = [{"variable": 'K', "order": 1, "is_variance": False, "name": "dK_dt"}]
    flow = crn.flow(False, derivatives)


    specification_dreal = [('K < 0.1', 'dK_dt >= 0', '(and (K > 0.4) (dK_dt = 0) )'), ('','dK_dt <= 0', '')]
    #specification_dreal  = [('K < 0.1', '',  'K > 0.4'), ('','', '')]
    #('', 'dK_dt <= 0', '(and (K >= 0)(K < 0.1))')]

    #specification_dreal = []
    drh = iSATParser.constructdReal(crn, specification_dreal, flow, max_time=100)
    with open('bellshape.drh', 'w') as file:
        file.write(drh)

    sc = SolverCallerDReal("./bellshape.drh", dreal_path="../dReal-3.16.09.01-linux/bin/dReach")
    result_files = sc.single_synthesis(cost=0, precision=0.1)


    for file_name in result_files:
            vals, all_vals = sc.getCRNValues('./bellshape_1_0.smt2.proof')
            initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

            print("Initial Conditions", initial_conditions)
            print("Flow:", parametrised_flow)
            t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                           plot_name=file_name + "-simulationdreal.png", t = linspace(0, 100, 1000))
            print("\n\n")
            print(variable_names)
            print(sol)
            savetxt(file_name + "-simulationdreal.csv", sol, delimiter=",")


if __name__ == "__main__":
    crn = form_crn()

    synthesize_with_isat(crn)
    synthesize_with_dreal(crn)
