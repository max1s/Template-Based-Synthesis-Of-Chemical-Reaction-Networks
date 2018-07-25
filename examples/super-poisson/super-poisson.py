from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt
from numpy import linspace

def form_crn():
    A = Species('A', initial_value=0)
    B = Species('B', initial_value=0)
    lam1 = LambdaChoice([A, B], 1)


    r1 = Reaction([], [(A, 2)], RateConstant('k_1', 0.23, 0.23))
    r2 = Reaction([(A, 1)], [], RateConstant('k_2', 0.94, 0.94))

    return CRNSketch([r1, r2], [], [])

def synthesize_with_isat(crn):
    derivatives = []
    flow = crn.flow(True, derivatives)

    # specification_dreal = [('','cA > A','')]
    specification = [('','','')]

    flow = crn.flow(False, derivatives)
    hys = iSATParser.constructISAT(crn, specification, flow, max_time=100)
    with open('superpoisson.hys', 'w') as file:
        file.write(hys)


    sc = SolverCallerISAT("./superpoisson.hys",
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
    derivatives = []
    flow = crn.flow(True, derivatives)

    specification_dreal = [('','varA >= A','')]
    #specification_dreal = [('','','')]
    drh = iSATParser.constructdReal(crn, specification_dreal, flow, max_time=100)
    with open('superpoisson.drh', 'w') as file:
        file.write(drh)

    sc = SolverCallerDReal("./superpoisson.drh", dreal_path="../dReal-3.16.09.01-linux/bin/dReach")
    result_files = sc.single_synthesis(cost=0, precision=0.1)


    for file_name in result_files:
            vals, all_vals = sc.getCRNValues('./superpoisson_0_0.smt2.proof')
            initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

            print("Initial Conditions", initial_conditions)
            print("Flow:", parametrised_flow)
            t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                           plot_name=file_name + "-simulationdreal.png",
                                                           t = linspace(0, 10, 100), lna=True)
            print("\n\n")
            print(variable_names)
            print(sol)
            savetxt(file_name + "-simulationdreal.csv", sol, delimiter=",")
if __name__ == "__main__":
    crn = form_crn()

    #synthesize_with_isat(crn)
    synthesize_with_dreal(crn)
