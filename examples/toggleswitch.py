from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt


def form_crn():
    lacl = Species('lacl')
    clts = Species('clts')

    empty_1 = Species('empty_1')
    empty_2 = Species('empty_2')
    empty_3 = Species('empty_3')
    empty_4 = Species('empty_4')

    a_2 = RateConstant('a_2', 0.1, 2)
    a_4 = RateConstant('a_4', 0.1, 2)

    alpha_1 = RateConstant('alpha_1', 0.1, 5)
    alpha_2 = RateConstant('alpha_2', 0.1, 5)


    beta = RateConstant('beta', 0.1, 5)
    gamma = RateConstant('gamma', 0.1, 5)

   #a_1 = RateConstant(str(alpha_1.symbol / 1 + (clts.symbol ** beta.symbol)), 0, 5)
   # a_3 = RateConstant(str(alpha_2.symbol / 1 + (lacl.symbol ** gamma.symbol)), 0, 5)

    a_1 = RateConstant('a_1', 0, 5)
    a_3 = RateConstant('a_3', 0 ,5)
    input1 = InputSpecies("Input1", sympify(
        "0.1*t + 54.2735055776743*exp(-(0.04*t - 2.81375654916915)**2) + 35.5555607722356/(1.04836341039216e+15*(1/t)**10.0 + 1)"),
                          15)


    reaction1 = Reaction([Term(empty_1, 1)], [Term(lacl, 1)], a_1)
    reaction2 = Reaction([Term(lacl, 1)], [Term(empty_2,1)], a_2)
    reaction3 = Reaction([Term(empty_3, 1)],[Term(clts, 1)], a_3)
    reaction4 = Reaction([Term(clts, 1)], [Term(empty_4, 1)], a_4)

    return CRNSketch([reaction1, reaction2, reaction3, reaction4], [], [input1])


def synthesize_with_isat(crn):
    isLNA = False
    derivatives = []
    specification = []

    flow = crn.flow(isLNA, derivatives)

    hys = iSATParser.constructdReal(crn, specification, flow)

    with open("./toggleswitch.hys", "w") as f:
        f.write(hys)

    sc = SolverCallerISAT("./toggleswitch.hys", isat_path="../isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt")

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

    isLNA = False
    derivatives = []
    specification_dreal = []

    flow = crn.flow(isLNA, derivatives)

    drh = iSATParser.constructdReal(crn, specification_dreal, flow)
    with open('toggleswitch.drh', 'w') as file:
        file.write(drh)

    sc = SolverCallerDReal("./toggleswitch.drh", dreal_path="../dReal-3.16.09.01-linux/bin/dReach")
    result_files = sc.single_synthesis(cost=0)

    for file_name in result_files:
            vals, all_vals = sc.getCRNValues('./sixreactionnetwork_1_0.smt2.proof')
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
