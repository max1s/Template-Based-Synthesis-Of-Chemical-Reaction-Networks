from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt


def form_crn():

    lacI = Species('lacI')
    cIts = Species('cIts')
    lacIStar = Species('lacIStar')
    IPTG = Species('IPTG')

    input1 = InputSpecies("Input1", sympify("(-355.556 * exp(-177.778 * (-0.5 + t) ** 2) * (-0.5 + t))"), initial_value=0.01)

    k_1 = RateConstant('k_1', 0.1, 2)
    k_2 = RateConstant('k_2', 0.1, 2)

    alpha_1 = RateConstant('alpha_1', 0.1, 2)
    alpha_2 = RateConstant('alpha_2', 0.1, 2)

    # Gamma and Beta have special meanings to sympy
    beta = RateConstant('beta_1', 2, 2)
    gamma = RateConstant('gamma_1', 2, 2)

    degredation_1 = RateConstant('degredation_1', 0.1, 2)
    degredation_2 = RateConstant('degredation_2', 0.1, 2)

    Ka = RateConstant('Ka', 1, 1)

    reaction1 = ArbitraryRateReaction([], [Term(lacI, 1)], 'alpha_2 / (1 + ((k_2 * lacI)/ (k_1 * IPTG) + k_2) ** gamma_1)', [alpha_2, k_1, k_2, gamma])
    reaction2 = HillRepressionReaction([Term(lacI, 1)], [Term(lacI, 1),Term(cIts, 1)], alpha_1, Ka, beta)

    reaction3 = Reaction([Term(lacI, 1)], [], degredation_1)
    reaction4 = Reaction([Term(cIts, 1)], [], degredation_2)

    return CRNSketch([reaction1, reaction2, reaction3, reaction4], [], [input1])


def synthesize_with_isat(crn):
    isLNA = False
    derivatives = []
    specification = []

    flow = crn.flow(isLNA, derivatives)

    hys = iSATParser.constructISAT(crn, specification, flow)

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
