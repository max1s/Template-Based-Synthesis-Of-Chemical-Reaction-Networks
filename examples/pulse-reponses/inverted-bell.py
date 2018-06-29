from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt
from numpy import linspace


# This is based on topology 583

def parametrise_flow_365(flow):
    vals = {}

    sa = 0.01
    vals = {'SF': 1, 'k_1': sa, 'k_5': sa, 'k_9': sa, 'k_13': sa, 'k_17': sa, 'k_21': sa}
    vals['inpt'] = 1

    # rows are species doing effect, so go down cols
    vals['k_3'] = 0  # P2 activating P1
    vals['k_4'] = 0.87  # P3 activating P1
    vals['k_7'] = 0.94  # P2 inactivating P1
    vals['k_8'] = 0  # P3 inactivating P1

    vals['k_11'] = 0.03  # P1 activating P2
    vals['k_12'] = 0  # P3 activating P2
    vals['k_15'] = 0  # P1 inactivating P2
    vals['k_16'] = 0.43  # P3 inactivating P2

    vals['k_19'] = 0  # P1 activating P3
    vals['k_20'] = 0.43  # P2 activating P3
    vals['k_23'] = 0.85  # P1 inactivating P3
    vals['k_24'] = 0.99  # P2 inactivating P3

    # No self interactions
    for i in [2, 6, 9, 10, 13, 14, 15, 17, 18, 22]:
        vals['k_' + str(i)] = 0

    for x in flow:
        for v in vals:
            flow[x] = flow[x].subs(sympify(v), sympify(vals[v]))

    print(flow)
    return flow


def create_rate_constant(name, val):
    return RateConstant(name, val, val)


def form_crn():
    input1 = InputSpecies("Input1", sympify("-(-50 + t)/(50 * exp((-50 + t)**2/100))"), initial_value=0.1)

    POne = Species('POne', initial_max=1, initial_value=0.95)
    PTwo = Species('PTwo', initial_max=1, initial_value=0.05)
    PThree = Species('PThree', initial_max=1, initial_value=0.06)
    POneStar = Species('POneStar', initial_max=1, initial_value=0.05)
    PTwoStar = Species('PTwoStar', initial_max=1, initial_value=0.95)
    PThreeStar = Species('PThreeStar', initial_max =1, initial_value=0.94)

    sa = 0.01

    # vals = {'SF': 1, 'k_1': sa, 'k_5': sa, 'k_9': sa, 'k_13': sa, 'k_17': sa, 'k_21': sa}
    # create_rate_constant('k_1', sa)
    #[2, 6, 9, 10, 13, 14, 15, 17, 18, 22]:

    k1 = create_rate_constant('k_1', sa)
    k2 = create_rate_constant('k_2', 0)
    k3 = create_rate_constant('k_3', 0)
    k4 = create_rate_constant('k_4', 0)
    k5 = create_rate_constant('k_5', sa)
    k6 = create_rate_constant('k_6', 0)

    k7 = create_rate_constant('k_7', 1)
    k8 = create_rate_constant('k_8', 1)
    k9 = create_rate_constant('k_9', sa)
    k10 = create_rate_constant('k_10', 0)
    k11 = create_rate_constant('k_11', 0)
    k12 = create_rate_constant('k_12', 1)

    k13 = create_rate_constant('k_13', sa)
    k14 = create_rate_constant('k_14', 0)
    k15 = create_rate_constant('k_15', 1)
    k16 = create_rate_constant('k_16', 0)
    k17 = create_rate_constant('k_17', sa)
    k18 = create_rate_constant('k_18', 0)

    k19 = create_rate_constant('k_19', 0)
    k20 = create_rate_constant('k_20', 1)
    k21 = create_rate_constant('k_21', sa)
    k22 = create_rate_constant('k_22', 0)
    k23 = create_rate_constant('k_23', 1)
    k24 = create_rate_constant('k_24', 0)

    reactionI = Reaction([Term(input1, 1), Term(POne, 1)], [Term(POneStar, 1)], RateConstant('inpt', 1.0, 1.0))

    # POne              -> POneStar
    # POne + POneStar   -> POneStar + POneStar
    # POne + PTwoStar   -> POneStar + PTwoStar
    # POne + PThreeStar -> POneStar + PThreeStar

    reaction1 = Reaction([Term(POne, 1)], [Term(POneStar, 1)], k1)
    reaction2 = Reaction([Term(POne, 1), Term(POneStar, 1)], [Term(POneStar, 1), Term(POneStar, 1)], k2)
    reaction3 = Reaction([Term(POne, 1), Term(PTwoStar, 1)], [Term(POneStar, 1), Term(PTwoStar, 1)], k3)
    reaction4 = Reaction([Term(POne, 1), Term(PThreeStar, 1)], [Term(POneStar, 1), Term(PThreeStar, 1)], k4)

    # POneStar              -> POne
    # POneStar + POneStar   -> POne + POneStar
    # POneStar + PTwoStar   -> POne + PTwoStar
    # POneStar + PThreeStar -> POne + PThreeStar
    reaction5 = Reaction([Term(POneStar, 1)], [Term(POne, 1)], k5)
    reaction6 = Reaction([Term(POneStar, 1), Term(POneStar, 1)], [Term(POne, 1), Term(POneStar, 1)], k6)
    reaction7 = Reaction([Term(POneStar, 1), Term(PTwoStar, 1)], [Term(POne, 1), Term(PTwoStar, 1)], k7)
    reaction8 = Reaction([Term(POneStar, 1), Term(PThreeStar, 1)], [Term(POne, 1), Term(PThreeStar, 1)], k8)

    # PTwo              -> PTwoStar
    # PTwo + POneStar   -> PTwoStar + POneStar
    # PTwo + PTwoStar   -> PTwoStar + PTwoStar
    # PTwo + PThreeStar -> PTwoStar + PThreeStar
    reaction9 = Reaction([Term(PTwo, 1)], [Term(PTwoStar, 1)], k9)
    reaction10 = Reaction([Term(PTwo, 1), Term(POneStar, 1)], [Term(PTwoStar, 1), Term(POneStar, 1)], k11)
    reaction11 = Reaction([Term(PTwo, 1), Term(PTwoStar, 1)], [Term(PTwoStar, 1), Term(PTwoStar, 1)], k10)
    reaction12 = Reaction([Term(PTwo, 1), Term(PThreeStar, 1)], [Term(PTwoStar, 1), Term(PThreeStar, 1)], k12)

    # PTwoStar              -> PTwo
    # PTwoStar + POneStar   -> PTwo + POneStar
    # PTwoStar + PTwoStar   -> PTwo + PTwoStar
    # PTwoStar + PThreeStar -> PTwo + PThreeStar
    reaction13 = Reaction([Term(PTwoStar, 1)], [Term(PTwo, 1)], k13)
    reaction14 = Reaction([Term(PTwoStar, 1), Term(POneStar, 1)], [Term(PTwo, 1), Term(POneStar, 1)], k15)
    reaction15 = Reaction([Term(PTwoStar, 1), Term(PTwoStar, 1)], [Term(PTwo, 1), Term(PTwoStar, 1)], k14)
    reaction16 = Reaction([Term(PTwoStar, 1), Term(PThreeStar, 1)], [Term(PTwo, 1), Term(PThreeStar, 1)], k16)

    # PThree              -> ThreeStar
    # PThree + POneStar   -> PThreeStar + POneStar
    # PThree + PTwoStar -> PThreeStar + PTwoStar
    # PThree + PThreeStar -> PThreeStar + PThreeStar
    reaction17 = Reaction([Term(PThree, 1)], [Term(PThreeStar, 1)], k17)
    reaction18 = Reaction([Term(PThree, 1), Term(POneStar, 1)], [Term(PThreeStar, 1), Term(POneStar, 1)], k19)
    reaction19 = Reaction([Term(PThree, 1), Term(PTwoStar, 1)], [Term(PThreeStar, 1), Term(PTwoStar, 1)], k20)
    reaction20 = Reaction([Term(PThree, 1), Term(PThreeStar, 1)], [Term(PThreeStar, 1), Term(PThreeStar, 1)], k18)

    # PThreeStar -> PThree
    # PThreeStar + POneStar -> PThree + POneStar
    # PThreeStar + PTwoStar -> PThree + PTwoStar
    # PThreeStar + PThreeStar -> PThree + PThreeStar
    reaction21 = Reaction([Term(PThreeStar, 1)], [Term(PThree, 1)], k21)
    reaction22 = Reaction([Term(PThreeStar, 1), Term(POneStar, 1)], [Term(PThree, 1), Term(POneStar, 1)], k23)
    reaction23 = Reaction([Term(PThreeStar, 1), Term(PTwoStar, 1)], [Term(PThree, 1), Term(PTwoStar, 1)], k24)
    reaction24 = Reaction([Term(PThreeStar, 1), Term(PThreeStar, 1)], [Term(PThree, 1), Term(PThreeStar, 1)], k22)

    return CRNSketch(
        [reaction1, reaction2, reaction3, reaction4, reaction5, reaction6, reaction7, reaction8, reaction9, reaction10,
         reaction11,
         reaction12, reaction13, reaction14, reaction15, reaction16, reaction17, reaction18, reaction19, reaction20,
         reaction21,
         reaction22, reaction23, reaction24, reactionI], [], [input1])




def synthesize_with_dreal(crn):


    derivatives = [{"variable": 'PThreeStar', "order": 1, "is_variance": False, "name": "PThreeStar_dot"}]
    flow = crn.flow(False, derivatives)



    specification_dreal = [('', '', '(inputTime > 20)'),
                           ('', '(PThreeStar_dot < 0)', '(abs(PThreeStar_dot) < 0.1)'),
                           ('', '(PThreeStar_dot > 0)', '(abs(PThreeStar_dot) < 0.1)'),
                           ('', '(abs(PThreeStar_dot) < 0.1)', '(inputTime > 100)')]

    #flow = crn.flow(False, derivatives)
    #specification_dreal = [('', '', '')]

    drh = iSATParser.constructdReal(crn, specification_dreal, flow, max_time=350, other_constraints='', scale_factor=1)
    with open('inverted-bell.drh', 'w') as file:
        file.write(drh)

    sc = SolverCallerDReal("./inverted-bell.drh", dreal_path="../dReal-3.16.09.01-linux/bin/dReach")
    result_files = sc.single_synthesis(cost=0, precision=0.1, max_depth=len(specification_dreal))

    for file_name in result_files:
        vals, all_vals = sc.getCRNValues('./inverted-bell_3_0.smt2.proof')
        initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

        print("Initial Conditions", initial_conditions)
        print("Flow:", parametrised_flow)
        t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                       plot_name=file_name + "-simulationdreal.png",
                                                       t=linspace(0, 300, 1000))
        print("\n\n")
        print(variable_names)
        print(sol)
        savetxt(file_name + "-simulationdreal.csv", sol, delimiter=",")


if __name__ == "__main__":
    crn = form_crn()
    synthesize_with_dreal(crn)
