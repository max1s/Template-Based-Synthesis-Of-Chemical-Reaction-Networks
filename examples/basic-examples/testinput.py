from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt
from numpy import linspace

def form_crn():
    rate = "-(-50 + t)/(50 * exp((-50 + t)**2/100))"
    input1 = InputSpecies("Input1", sympify(rate), initial_value=0.01)
    A = Species('A')
    reaction1 = Reaction([Term(input1, 1), Term(A, 1)],[Term(input1, 1)], RateConstant('inpt', 1.0, 1.0))

    return CRNSketch([reaction1], [], [input1])


def synthesize_with_isat(crn):
    isLNA = False
    derivatives = []
    #derivatives = [{"variable": 'cIts', "order": 1, "is_variance": False, "name": "dcIts"}]
    specification = []

    a = ''
    a.

    flow = crn.flow(isLNA, derivatives)

    hys = iSATParser.constructISAT(crn, specification, flow, max_time=100, other_constraints='',scale_factor=1)

    with open("./toggleswitch.hys", "w") as f:
        f.write(hys)

    sc = SolverCallerISAT("./toggleswitch.hys", isat_path="../isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt")

    result_files = sc.single_synthesis(cost=0, precision=0.01, msw=0.1)

    for file_name in result_files:
        vals, all_vals = sc.getCRNValues(file_name)
        initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

        print("\n\nInitial Conditions", initial_conditions)
        print("Flow:", parametrised_flow)

        t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                       plot_name=file_name + "-simulation.png", t = linspace(0, 100, 1000))
        print("\n\n")
        print(variable_names)
        print(sol)
        savetxt(file_name + "-simulation.csv", sol, delimiter=",")


def synthesize_with_dreal(crn):

    isLNA = False
    derivatives = []
    specification_dreal = [('','','')]

    flow = crn.flow(isLNA, derivatives)

    drh = iSATParser.constructdReal(crn, specification_dreal, flow)
    with open('testinput.drh', 'w') as file:
        file.write(drh)

    sc = SolverCallerDReal("./testinput.drh", dreal_path="../dReal-3.16.09.01-linux/bin/dReach")
    result_files = sc.single_synthesis(cost=0)

    for file_name in result_files:
            vals, all_vals = sc.getCRNValues('./testinput_1_0.smt2.proof')
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
    synthesize_with_dreal(crn)
