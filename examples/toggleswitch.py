from CRNSynthesis.symbolicLNA import *
from CRNSynthesis import iSATParser
from CRNSynthesis.solverCaller import SolverCaller
from CRNSynthesis.solverCaller import SolverCallerDReal
from sympy import init_printing, Matrix, transpose, pprint


def exampleToggleSwitch():
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

    crn = CRNSketch([reaction1, reaction2, reaction3, reaction4], [], [input1])

    isLNA = False
    derivatives = []
    specification = []

    flow = crn.flow(isLNA, derivatives)

    return iSATParser.constructdReal(crn, specification, flow)


def complete_process():

    flow, problem_string = toggleSwitch()
    with open("./toggleswitch.hys", "w") as f:
        f.write(problem_string)

    sc = SolverCaller("./test.hys")
    result_file = sc.single_synthesis(cost=60, precision=0.1)
    param_values = sc.getCRNValues(result_file)

    print("\n\nSpecific CRN identified is:\n")
    print(sc.get_parametrised_flow(flow, param_values))

def dReal_process(filename):
    problem_string = toggleSwitch()
    print problem_string
    with open("./" + filename, "w") as f:
        f.write(problem_string)

    sc = SolverCallerDReal(model_path="./" + filename)
    result_file = sc.single_synthesis(precision=0.1)
    param_values = sc.getCRNValues(result_file[0])
    print param_values
    quit()
    #print("\n\nSpecific CRN identified is:\n")
    #print(sc.get_parametrised_flow(flow, param_values))



if __name__ == "__main__":
    complete_process()
