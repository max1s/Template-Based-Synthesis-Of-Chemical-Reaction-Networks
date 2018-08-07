from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCallerISAT, SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint
from numpy import savetxt
from numpy import linspace


# This is based on topology 362

def create_rate_constant(name, val):
    return RateConstant(name, 0, 0.1)

def zero(name, val):
    return RateConstant(name, 0, 0)

def form_crn():
    input1 = InputSpecies("Input1", sympify("-(-50 + t)/(50 * exp((-50 + t)**2/100))"), initial_value=0.1)

    POne = Species('POne', initial_max=1, initial_value=0.95)
    PTwo = Species('PTwo', initial_max=1, initial_value=0.05)
    PThree = Species('PThree', initial_max=1, initial_value=0.06)
    POneStar = Species('POneStar', initial_max=1, initial_value=0.05)
    PTwoStar = Species('PTwoStar', initial_max=1, initial_value=0.95)
    PThreeStar = Species('PThreeStar', initial_max =1, initial_value=0.94)

    sa = 0
    sd = 0

    # vals = {'SF': 1, 'k_1': sa, 'k_5': sa, 'k_9': sa, 'k_13': sa, 'k_17': sa, 'k_21': sa}
    # create_rate_constant('k_1', sa)
    #[2, 6, 9, 10, 13, 14, 15, 17, 18, 22]:

    k1 = create_rate_constant('k_1', sa)
    k2 = zero('k_2', 0)
    k3 = zero('k_3', 0)
    k4 = zero('k_4', 0)
    k5 = zero('k_5', sa)
    k6 = zero('k_6', 0)

    k7 = create_rate_constant('k_7', 1)
    k8 = zero('k_8', 1)
    k9 = zero('k_9', sa)
    k10 = zero('k_10', 0)
    k11 = create_rate_constant('k_11', 0)
    k12 = zero('k_12', 1)

    k13 = zero('k_13', sa)
    k14 = zero('k_14', 0)
    k15 = zero('k_15', 1)
    k16 = create_rate_constant('k_16', 0)
    k17 = zero('k_17', sa)
    k18 = zero('k_18', 0)

    k19 = zero('k_19', 0)
    k20 = create_rate_constant('k_20', 1)
    k21 = zero('k_21', sa)
    k22 = zero('k_22', 0)
    k23 = create_rate_constant('k_23', 1)
    k24 = zero('k_24', 0)

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
         reaction12, reaction13, reaction14, reaction15, reaction16, reaction17, reaction18, reaction19, reaction20,
         reaction21, reaction22, reaction23, reaction24, reactionI], 
        [], [input1])




def synthesize_with_dreal(crn, name, specification_dreal, result_file):


    derivatives = [{"variable": 'PThreeStar', "order": 1, "is_variance": False, "name": "PThreeStar_dot"}]
    flow = crn.flow(False, derivatives)


    #flow = crn.flow(False, derivatives)
    #specification_dreal = [('', '', '')]

    drh = iSATParser.constructdReal(crn, specification_dreal, flow, max_time=350, other_constraints='', scale_factor=1)
    with open(name, 'w') as file:
        file.write(drh)

    sc = SolverCallerDReal(name, dreal_path="/code/dReal-3.16.06.02-linux/bin/dReach")
    result_files = sc.single_synthesis(cost=0, precision=0.1, max_depth=len(specification_dreal))

    for file_name in [result_file]: # TODO: remove ths aegument

        print file_name

        # './inverted-bell_3_0.smt2.proof'
        vals, all_vals = sc.getCRNValues(file_name)
        initial_conditions, parametrised_flow = sc.get_full_solution(crn, flow, all_vals)

        print("Initial Conditions", initial_conditions)
        print("Flow:", parametrised_flow)
        t, sol, variable_names = sc.simulate_solutions(initial_conditions, parametrised_flow,
                                                       plot_name=file_name + "-simulationdreal.png",
                                                       t=linspace(0, 300, 1000),
                                                       mode_times=all_vals["time"])
        print("\n\n")
        print(variable_names)
        print(sol)
        savetxt(file_name + "-simulationdreal.csv", sol, delimiter=",")


if __name__ == "__main__":
    crn = form_crn()




    print "\n\n\n\nSynthesizing an inverted bellshape:"
    name = "./inverted-bellshape.drh"
    spec = [('', '', '(inputTime > 20)'),
                       ('', '(PThreeStar_dot < -0.1)', '(abs(PThreeStar_dot) < 0.1)'),
                       ('', '(PThreeStar_dot > 0.1)', '(abs(PThreeStar_dot) < 0.1)'),
                       ('', '(abs(PThreeStar_dot) < 0.1)', '(inputTime > 100)')]
                       # on previous iteration, constraints on gradient were /- 0 rather than +/- 0.1

    synthesize_with_dreal(crn, name, spec, 'inverted-bellshape_3_0.smt2.proof')




    # this will be unsat unless input starts at zero. But with rpecision 0.1, got spurious solution...
    print "\n\n\n\nSynthesizing a switch:"
    name = "./switch.drh"
    spec = [('', '', '(and (inputTime > 20)(inputTime < 50))'),
                               ('', '(PThreeStar_dot > 0)', '(and (PThreeStar > 0.5)(abs(PThreeStar_dot) <= 0.01))'),
                               ('', 'abs(PThreeStar_dot) < 0.01', 'inputTime > 100')]

    synthesize_with_dreal(crn, name, spec, 'switch_2_0.smt2.proof')



    


    # TODO: fix this - it is currently unsat
    print "\n\n\n\nSynthesizing a derivative response:"
    name = "./derivative.drh"
    spec = [('', '', 'inputTime > 20'), ('', 'PThreeStar_dot < 0', '(PThreeStar_dot=0)'), ('','PThreeStar_dot > 0','(PThreeStar_dot=0'), ('','PThreeStar_dot < 0','PThreeStar_dot = 0'), ('','PThreeStar_dot = 0','inputTime > 100')]

    synthesize_with_dreal(crn, name, spec, 'derivative_4_0.smt2.proof')





    print "\n\n\n\nSynthesizing a non-response:"
    name = "./constant.drh"
    spec = [('', '', '(and (inputTime > 5)(inputTime < 15)'),  # 'constant' mode starts before pulse, but after an equilibration time
                           ('', ' abs(PThreeStar_dot) < 0.1', 'inputTime > 100')]


    synthesize_with_dreal(crn, name, spec, 'constant_1_0.smt2.proof')







# TODO: constant [should give unsat]
# TODO: integrator and constant not implemented
