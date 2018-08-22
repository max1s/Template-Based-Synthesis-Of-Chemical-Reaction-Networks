from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt
from numpy import linspace

def form_crn():
    #rate = sympify("SF*(-6*(1/(SF*t))^10.0/(SF*t*(49*(1/(SF*t))^10.0 + 1)^2) + 42*(1/(SF*t))^10.0/(SF*t*(39*(1/(SF*t))^10.0 + 1)^2))")
    rate =  sympify("SF*(-1.2175*10^20*(1/(SF*t))^10.0/(SF*t*(1.29*10^17*(1/(SF*t))^10.0 + 1)^2) + 74539337510769.1*(1/(SF*t))^10.0/(SF*t*(74539337510.7691*(1/(SF*t))^10.0 + 1)^2))")
    #rate = sympify("SF*(-7.48*10^21*(1 / (SF*t))^10.0/(SF*t*(7.1645*10^18*(1/(SF*t))^10.0 + 1)^2) + 6.851*10^17*(1 / (SF * t)) ^ 10.0 / (SF * t * (656575214694504 * (1 / (SF * t)) ^ 10.0 + 1) ^ 2))") #try three

    #rate = sympify("0")
    lacI = Species('lacI', initial_max=200, initial_min=0, initial_value=0)
    cIts = Species('cIts', initial_max=100, initial_min=0, initial_value=0)
    #lacIStar = Species('lacIStar', initial_max = 100, initial_min=0)
    # IPTG = Species('IPTG')

    input1 = InputSpecies("Input1",rate, initial_value=0)
    #k_1 = RateConstant('k_1', 0.2, 2)
    #k_2 = RateConstant('k_2', 0.2, 2)

    alpha_1 = RateConstant('alpha_1', 156.25, 156.25) #156.25
    alpha_2 = RateConstant('alpha_2', 15.6, 15.6) #15.6

    degredation_1 = RateConstant('degredation_1', 1, 1)
    degredation_2 = RateConstant('degredation_2', 1, 1)

    Ka = RateConstant('Ka', 1, 1)


    reaction1 = ArbitraryRateReaction([Term(input1, 1)], [Term(cIts, 1), Term(input1, 1)], ' (Input1)*alpha_2 + (1 - Input1) * alpha_2 / (1 + lacI ** 2)', [alpha_2])
    reaction2 = HillRepressionReaction([Term(cIts, 1)], [Term(lacI, 1), Term(cIts, 1)], alpha_1, Ka, 2)

    reaction3 = Reaction([Term(lacI, 1)], [], degredation_1)
    reaction4 = Reaction([Term(cIts, 1)], [], degredation_2)

    return CRNSketch([reaction1, reaction2, reaction3, reaction4], [], [input1])


def synthesize_with_isat(crn):
    isLNA = False
    derivatives = []
    #derivatives = [{"variable": 'cIts', "order": 1, "is_variance": False, "name": "dcIts"}]
    specification = [('','','')]

    flow = crn.flow(False, derivatives)
    hys = iSATParser.constructISAT(crn, specification, flow, max_time=100)
    with open('toggleswitch.hys', 'w') as file:
        file.write(hys)


    sc = SolverCallerISAT("./toggleswitch.hys",
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

    isLNA = False
    derivatives = []
    #specification_dreal = [('','','')]
    specification_dreal = [('','','(and (inputTime > 70)(cIts > 10))')]
    flow = crn.flow(False, derivatives)
    drh = iSATParser.constructdReal(crn, specification_dreal, flow, max_time=100)
    with open('toggleswitch.drh', 'w') as file:
        file.write(drh)

    sc = SolverCallerDReal("./toggleswitch.drh", dreal_path="../dReal-3.16.09.01-linux/bin/dReach")
    result_files = sc.single_synthesis(cost=0, precision=0.1)


    for file_name in result_files:
            vals, all_vals = sc.getCRNValues('./toggleswitch_0_0.smt2.proof')
            initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

            print("Initial Conditions", initial_conditions)
            print("Flow:", parametrised_flow)
            t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                           plot_name=file_name + "-simulationdreal.png",
                                                           t = linspace(0, 100, 200)) # mode_times=all_vals["time"])

            print("\n\n")
            print(variable_names)
            print(sol)
            savetxt(file_name + "-simulationdreal.csv", sol, delimiter=",")


if __name__ == "__main__":
    crn = form_crn()

    #synthesize_with_isat(crn)
    synthesize_with_dreal(crn)
