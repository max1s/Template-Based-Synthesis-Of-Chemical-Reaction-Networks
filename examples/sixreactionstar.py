from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt
from numpy import linspace

def form_crn():
    input1 = InputSpecies("Input1", sympify("(-355.556 * exp(-177.778 * (-0.5 + t) ** 2) * (-0.5 + t))"), initial_value=0.01)

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

    reactionI = Reaction([Term(input1, 1)], [Term(POne, 1)], RateConstant('inpt', 1.0, 1.0))

    #POne ->POneStar
    #POne + POneStar -> POneStar + POneStar
    reaction1 = Reaction([Term(POne, 1)], [Term(POneStar, 1)], k1)
    reaction2 = Reaction([Term(POne, 1),Term(POneStar, 1)], [Term(POneStar, 1),Term(POneStar, 1)], k2)

    #POne + PTwoStar -> POneStar + PTwoStar
    #POne + PThreeStar -> POneStar + PThreeStar
    reaction3 = Reaction([Term(POne, 1),Term(PTwoStar, 1)], [Term(POneStar, 1),Term(PTwoStar, 1)], k3)
    reaction4 = Reaction([Term(POne, 1),Term(PThreeStar, 1)], [Term(POneStar, 1),Term(PThreeStar, 1)], k4)


    # POneStar -> POne
    # POneStar + POneStar -> POne + POneStar
    reaction5 = Reaction([Term(POneStar, 1)], [Term(POne, 1)], k5)
    reaction6 = Reaction([Term(POneStar, 1),Term(POneStar, 1)], [Term(POne, 1),Term(POneStar, 1)], k6)

    # POneStar + PTwoStar -> POne + PTwoStar
    # POneStar + PThreeStar -> POne + PThreeStar
    reaction7 = Reaction([Term(POneStar, 1),Term(PTwoStar, 1)], [Term(POne, 1),Term(PTwoStar, 1)], k7)
    reaction8 = Reaction([Term(POneStar, 1),Term(PThreeStar, 1)], [Term(POne, 1),Term(PThreeStar, 1)], k8)


    # PTwo -> PTwoStar
    # PTwo + POneStar -> PTwoStar + POneStar
    reaction9 = Reaction([Term(PTwo, 1)], [Term(PTwoStar, 1)], k9)
    reaction10 = Reaction([Term(PTwo, 1), Term(POneStar, 1)], [Term(PTwoStar, 1), Term(POneStar, 1)], k11)

    # PTwo + PTwoStar -> PTwoStar + PTwoStar
    # PTwo + PThreeStar -> PTwoStar + PThreeStar
    reaction11 = Reaction([Term(PTwo, 1),Term(PTwoStar, 1)], [Term(PTwoStar, 1),Term(PTwoStar, 1)], k10)
    reaction12 = Reaction([Term(PTwo, 1), Term(PThreeStar, 1)], [Term(PTwoStar, 1), Term(PThreeStar, 1)], k12)


    # PTwoStar -> PTwo
    # PTwoStar + POneStar -> PTwo + POneStar
    reaction13 = Reaction([Term(PTwoStar, 1)], [Term(PTwo, 1)], k13)
    reaction14 = Reaction([Term(PTwoStar, 1), Term(POneStar, 1)], [Term(PTwo, 1), Term(POneStar, 1)], k15)

    # PTwoStar + PTwoStar -> PTwo + PTwoStar
    # PTwoStar + PThreeStar -> PTwo + PThreeStar
    reaction15 = Reaction([Term(PTwoStar, 1),Term(PTwoStar, 1)], [Term(PTwo, 1),Term(PTwoStar, 1)], k14)
    reaction16 = Reaction([Term(PTwoStar, 1), Term(PThreeStar, 1)], [Term(PTwo, 1), Term(PThreeStar, 1)], k16)

 
    # PThree -> ThreeStar
    # PThree + POneStar -> PThreeStar + POneStar
    reaction17 = Reaction([Term(PThree, 1)], [Term(PThreeStar, 1)], k17)
    reaction18 = Reaction([Term(PThree, 1), Term(POneStar, 1)], [Term(PThreeStar, 1), Term(POneStar, 1)], k19)

    # PThree + PTwoStar -> PThreeStar + PTwoStar
    # PThree + PThreeStar -> PThreeStar + PThreeStar
    reaction19 = Reaction([Term(PThree, 1), Term(PTwoStar, 1)], [Term(PThreeStar, 1), Term(PTwoStar, 1)], k20)
    reaction20 = Reaction([Term(PThree, 1), Term(PThreeStar, 1)], [Term(PThreeStar, 1), Term(PThreeStar, 1)], k18)


    # PThreeStar -> PThree
    # PThreeStar + POneStar -> PThree + POneStar
    reaction21 = Reaction([Term(PThreeStar, 1)], [Term(PThree, 1)], k21)
    reaction22 = Reaction([Term(PThreeStar, 1), Term(POneStar, 1)], [Term(PThree, 1), Term(POneStar, 1)], k23)

    # PThreeStar + PTwoStar -> PThree + PTwoStar
    # PThreeStar + PThreeStar -> PThree + PThreeStar
    reaction23 = Reaction([Term(PThreeStar, 1), Term(PTwoStar, 1)], [Term(PThree, 1), Term(PTwoStar, 1)], k24)
    reaction24 = Reaction([Term(PThreeStar, 1), Term(PThreeStar, 1)], [Term(PThree, 1), Term(PThreeStar, 1)], k22)


    return CRNSketch([reaction1, reaction2, reaction3, reaction4, reaction5, reaction6, reaction7, reaction8, reaction9, reaction10, reaction11,
                      reaction12, reaction13, reaction14, reaction15, reaction16, reaction17, reaction18, reaction19, reaction20, reaction21,
                      reaction22, reaction23, reaction24, reactionI], [], [input1])


def synthesize_with_isat(crn):
    #derivatives = []
    derivatives = [{"variable": 'PThreeStar', "order": 1, "is_variance": False, "name": "PThreeStar_dot"}]
    flow = crn.flow(False, derivatives)
    specification = [('','','')]
    # specification = [('', 'PThree_dot >= 0', '((PThree > 0.1) and (PThree_dot < 0.001))'), ('', 'PThree_dot <= 0', '')]
    # specification = [('', '(K > 0.3) and (PThree_dot >= 0)', '(PThree_dot = 0)' ), ('', '(PThree_dot < 0)', '(K < 0.1) and (PThree_dot < 0)')]
    # specification = [('', '', 'PThree > 0.4 '), ('', '', 'PThree < 0.3')]
    #specification = [('','PThree_dot >= 0','PThree_dot = 0'),('','PThree_dot < 0','')]
    #specification = [('','PThreeStar_dot >= 0',' (PThreeStar_dot = 0) and (PThreeStar > 0.3) '),('','PThreeStar_dot <= 0','(PThreeStar < 0.2)')]
    #specification =  [('','','(PThree > 0.3) '),('','','(PThree < 0.2)')]
    hys = iSATParser.constructISAT(crn, specification, flow, max_time=10, other_constraints='PThree < 0.2; POne + POneStar = 1; PTwo + PTwoStar = 1; PThree + PThreeStar = 1;',scale_factor=1)
    with open('sixreactionstar.hys', 'w') as file:
        file.write(hys)

    sc = SolverCallerISAT("./sixreactionstar.hys", isat_path="../isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt")

    result_files = sc.single_synthesis(cost=0, precision=0.01, msw=0.05)

    for file_name in result_files:
        print("\n\n")
        # print(sc.getCRNValues(file_name))

        vals, all_vals = sc.getCRNValues(file_name)
        initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

        print("Initial Conditions", initial_conditions)
        print("Flow:", parametrised_flow)

        t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                       plot_name=file_name + "-simulation.png", t = linspace(0, 10, 1000))
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
    specification_dreal = [('','','')]
    drh = iSATParser.constructdReal(crn, specification_dreal, flow, max_time=10, other_constraints='(and (PThreeStar < 0.4) (POne + POneStar = 1) (PTwo + PTwoStar = 1) (PThree + PThreeStar = 1))',scale_factor=1)
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
    #synthesize_with_isat(crn)
    synthesize_with_dreal(crn)
