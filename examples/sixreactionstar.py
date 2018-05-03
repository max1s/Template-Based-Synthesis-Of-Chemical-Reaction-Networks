from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt
from numpy import linspace

def form_crn():
    input1 = InputSpecies("Input1", sympify("-(-150 + t)/(50 * exp((-150 + t)**2/100))"), initial_value=0.01)

    POne = Species('POne',initial_max=1)
    PTwo = Species('PTwo', initial_max=1)
    PThree = Species('PThree', initial_max=1)
    POneStar = Species('POneStar',initial_max=1)
    PTwoStar = Species('PTwoStar', initial_max=1)
    PThreeStar = Species('PThreeStar', initial_max=1)

    k1 = RateConstant('k_1', 0, 1)
    k2 = RateConstant('k_2', 0, 1)
    k3 = RateConstant('k_3', 0, 1)
    k4 = RateConstant('k_4', 0, 1)
    k5 = RateConstant('k_5', 0, 1)
    k6 = RateConstant('k_6', 0, 1)

    k7 = RateConstant('k_7', 0, 1)
    k8 = RateConstant('k_8', 0, 1)
    k9 = RateConstant('k_9', 0, 1)
    k10 = RateConstant('k_10', 0, 1)
    k11 = RateConstant('k_11', 0, 1)
    k12 = RateConstant('k_12', 0, 1)

    k13 = RateConstant('k_13', 0, 1)
    k14 = RateConstant('k_14', 0, 1)
    k15 = RateConstant('k_15', 0, 1)
    k16 = RateConstant('k_16', 0, 1)
    k17 = RateConstant('k_17', 0, 1)
    k18 = RateConstant('k_18', 0, 1)

    k19 = RateConstant('k_19', 0, 1)
    k20 = RateConstant('k_20', 0, 1)
    k21 = RateConstant('k_21', 0, 1)
    k22 = RateConstant('k_22', 0, 1)
    k23 = RateConstant('k_23', 0, 1)
    k24 = RateConstant('k_24', 0, 1)

    reactionI = Reaction([Term(input1, 1), Term(POne, 1)], [Term(POneStar, 1)], RateConstant('inpt', 1.0, 1.0))

    # POne              -> POneStar
    # POne + POneStar   -> POneStar + POneStar
    # POne + PTwoStar   -> POneStar + PTwoStar
    # POne + PThreeStar -> POneStar + PThreeStar

    reaction1 = Reaction([Term(POne, 1)],                      [Term(POneStar, 1)],                      k1)
    reaction2 = Reaction([Term(POne, 1), Term(POneStar, 1)],   [Term(POneStar, 1), Term(POneStar, 1)],   k2)
    reaction3 = Reaction([Term(POne, 1), Term(PTwoStar, 1)],   [Term(POneStar, 1), Term(PTwoStar, 1)],   k3)
    reaction4 = Reaction([Term(POne, 1), Term(PThreeStar, 1)], [Term(POneStar, 1), Term(PThreeStar, 1)], k4)


    # POneStar              -> POne
    # POneStar + POneStar   -> POne + POneStar
    # POneStar + PTwoStar   -> POne + PTwoStar
    # POneStar + PThreeStar -> POne + PThreeStar
    reaction5 = Reaction([Term(POneStar, 1)],                      [Term(POne, 1)],                      k5)
    reaction6 = Reaction([Term(POneStar, 1), Term(POneStar, 1)],   [Term(POne, 1), Term(POneStar, 1)],   k6)
    reaction7 = Reaction([Term(POneStar, 1), Term(PTwoStar, 1)],   [Term(POne, 1), Term(PTwoStar, 1)],   k7)
    reaction8 = Reaction([Term(POneStar, 1), Term(PThreeStar, 1)], [Term(POne, 1), Term(PThreeStar, 1)], k8)


    # PTwo              -> PTwoStar
    # PTwo + POneStar   -> PTwoStar + POneStar
    # PTwo + PTwoStar   -> PTwoStar + PTwoStar
    # PTwo + PThreeStar -> PTwoStar + PThreeStar
    reaction9 =  Reaction([Term(PTwo, 1)],                      [Term(PTwoStar, 1)],                       k9)
    reaction10 = Reaction([Term(PTwo, 1), Term(POneStar, 1)],   [Term(PTwoStar, 1), Term(POneStar, 1)],   k11)
    reaction11 = Reaction([Term(PTwo, 1), Term(PTwoStar, 1)],   [Term(PTwoStar, 1), Term(PTwoStar, 1)],   k10)
    reaction12 = Reaction([Term(PTwo, 1), Term(PThreeStar, 1)], [Term(PTwoStar, 1), Term(PThreeStar, 1)], k12)


    # PTwoStar              -> PTwo
    # PTwoStar + POneStar   -> PTwo + POneStar
    # PTwoStar + PTwoStar   -> PTwo + PTwoStar
    # PTwoStar + PThreeStar -> PTwo + PThreeStar
    reaction13 = Reaction([Term(PTwoStar, 1)],                      [Term(PTwo, 1)],                      k13)
    reaction14 = Reaction([Term(PTwoStar, 1), Term(POneStar, 1)],   [Term(PTwo, 1), Term(POneStar, 1)],   k15)
    reaction15 = Reaction([Term(PTwoStar, 1), Term(PTwoStar, 1)],   [Term(PTwo, 1), Term(PTwoStar, 1)],   k14)
    reaction16 = Reaction([Term(PTwoStar, 1), Term(PThreeStar, 1)], [Term(PTwo, 1), Term(PThreeStar, 1)], k16)

 
    # PThree              -> ThreeStar
    # PThree + POneStar   -> PThreeStar + POneStar
    # PThree + PTwoStar -> PThreeStar + PTwoStar
    # PThree + PThreeStar -> PThreeStar + PThreeStar
    reaction17 = Reaction([Term(PThree, 1)],                      [Term(PThreeStar, 1)],                      k17)
    reaction18 = Reaction([Term(PThree, 1), Term(POneStar, 1)],   [Term(PThreeStar, 1), Term(POneStar, 1)],   k19)
    reaction19 = Reaction([Term(PThree, 1), Term(PTwoStar, 1)],   [Term(PThreeStar, 1), Term(PTwoStar, 1)],   k20)
    reaction20 = Reaction([Term(PThree, 1), Term(PThreeStar, 1)], [Term(PThreeStar, 1), Term(PThreeStar, 1)], k18)


    # PThreeStar -> PThree
    # PThreeStar + POneStar -> PThree + POneStar
    # PThreeStar + PTwoStar -> PThree + PTwoStar
    # PThreeStar + PThreeStar -> PThree + PThreeStar
    reaction21 = Reaction([Term(PThreeStar, 1)],                      [Term(PThree, 1)],                      k21)
    reaction22 = Reaction([Term(PThreeStar, 1), Term(POneStar, 1)],   [Term(PThree, 1), Term(POneStar, 1)],   k23)
    reaction23 = Reaction([Term(PThreeStar, 1), Term(PTwoStar, 1)],   [Term(PThree, 1), Term(PTwoStar, 1)],   k24)
    reaction24 = Reaction([Term(PThreeStar, 1), Term(PThreeStar, 1)], [Term(PThree, 1), Term(PThreeStar, 1)], k22)


    return CRNSketch([reaction1, reaction2, reaction3, reaction4, reaction5, reaction6, reaction7, reaction8, reaction9, reaction10, reaction11,
                      reaction12, reaction13, reaction14, reaction15, reaction16, reaction17, reaction18, reaction19, reaction20, reaction21,
                      reaction22, reaction23, reaction24, reactionI], [], [input1])


def synthesize_with_isat(crn):
    #derivatives = []
    derivatives = [{"variable": 'PThreeStar', "order": 1, "is_variance": False, "name": "PThreeStar_dot"}]
    flow = crn.flow(False, derivatives)
    specification = [('','','')]

    specification =  [('','','(inputTime > 100) and (PThreeStar < 0.75)'), ('','','(PThreeStar > 0.9)')]  # try to capture an inverted-bellshape response

    hys = iSATParser.constructISAT(crn, specification, flow, max_time=350, other_constraints='',scale_factor=1)
    #with open('sixreactionstar.hys', 'w') as file:
    #    file.write(hys)

    sc = SolverCallerISAT("./sixreactionstar.hys", isat_path="../isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt")

    result_files = sc.single_synthesis(cost=0, precision=0.01, msw=0.05)
    #result_files = ["./results/sixreactionstar_0_0.01-isat.txt"]

    for file_name in result_files:
        print("\n\n")
        # print(sc.getCRNValues(file_name))

        vals, all_vals = sc.getCRNValues(file_name)
        initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

        print("Initial Conditions", initial_conditions)
        print("Flow:", parametrised_flow)

        t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                       plot_name=file_name + "-simulation.png", t = linspace(0, 350, 1000))
        print("\n\n")
        print(variable_names)
        print(sol)
        savetxt(file_name + "-simulation.csv", sol, delimiter=",")


def synthesize_with_dreal(crn):
    derivatives = [{"variable": 'PThreeStar', "order": 1, "is_variance": False, "name": "PThreeStar_dot"}]
    flow = crn.flow(False, derivatives)

    #specification_dreal = [('', '', 'PThree > 0.4 '), ('', '', 'PThree < 0.3')]
    #specification_dreal = [('','','(PThreeStar > 0.5) '),('','','(PThreeStar < 0.4)')]
    #specification_dreal = [('', 'PThree_dot >= 0', '(and (PThree > 0.3) (PThree_dot = 0))'), ('', 'PThree_dot <= 0', '(and (PThree >= 0)(PThree < 0.1))')]

    specification_dreal = [('','','inputTime > 100'), ('', '', 'PThreeStar < 0.75'), ('', '', 'POneStar > 0.9')]
    drh = iSATParser.constructdReal(crn, specification_dreal, flow, max_time=350, other_constraints='',scale_factor=1)
    with open('sixreactionstar.drh', 'w') as file:
        file.write(drh)

    sc = SolverCallerDReal("./sixreactionstar.drh", dreal_path="../dReal-3.16.09.01-linux/bin/dReach")
    result_files = sc.single_synthesis(cost=0, precision=0.1)

    for file_name in result_files:
            vals, all_vals = sc.getCRNValues('./sixreactionstar_1_0.smt2.proof')
            initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

            print("Initial Conditions", initial_conditions)
            print("Flow:", parametrised_flow)
            t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                           plot_name=file_name + "-simulationdreal.png", t = linspace(0, 10, 1000))
            print("\n\n")
            print(variable_names)
            print(sol)
            savetxt(file_name + "-simulationdreal.csv", sol, delimiter=",")


if __name__ == "__main__":
    crn = form_crn()
    synthesize_with_isat(crn)
    #synthesize_with_dreal(crn)
