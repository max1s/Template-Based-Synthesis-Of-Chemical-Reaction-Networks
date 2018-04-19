"""
This sub-module is responsible for calling iSAT or dReal, parsing their output, and substituting the extracted
parameter values into the ``CRNSketch`` to obtain a model that can be simulated.

It is also able to perform optimization by iteratively solving the SAT-ODE problem, decreasing the permitted cost
after each iteration.

"""

import subprocess as sub
import re
import os
from sympy import sympify
from scipy.integrate import odeint
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

class SolverCaller(object):
    """
    This abstract class is extended by ``SolverCallerISAT`` and ``SolverCallerDReal``, which are specialized to call the
    corresponding solvers.
    """
    def __init__(self, model_path="./bellshape.hys"):
        self.model_path = model_path

        self.results_folder = "/bellshaperesults"

        directory_name, file_name = os.path.split(model_path)
        self.model_name, _ = os.path.splitext(file_name)

        self.results_dir = os.path.join(directory_name, "results")

        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

    def single_synthesis(self, cost=20, precision=0.1, msw=0):
        """
        Call the solver once to synthesize a single system. Interpretation of precision and msw depends on which solver
        is used.

        :param cost: maxmimum permitted value for the cost (if 0, cost is ignored)
        :param precision:
        :param msw:
        :return:
        """
        return self.optimal_synthesis_decreasing_cost(max_cost=cost, min_cost=cost, precision=precision, msw=msw)

    def optimal_synthesis_decreasing_cost(self, max_cost=35, min_cost=10, precision=0.1):
        pass

    def simulate_solutions(self, initial_conditions, parametrised_flow, t=False, plot_name=""):
        """
        Numerically integrates the system using ``scipy.odeint`` to obtain a simulated time-course.
        Requires a specific initial_condition and flow dictionary in which parameters have been
        replaced by specific numerical values.

        :param initial_conditions: dictionary in which keys are species names (string) and values are corresponding concentration (float)
        :param parametrised_flow: dictionary in which keys are species names, and values are SympPy expressions for their time derivatives
        :param t: vector of times at which the system state will be calculated
        :param plot_name: name of file in which plot of results should be saved
        :return:
        """
        if not t:
            t = np.linspace(0, 1, 100)

        ic = []
        species_list = parametrised_flow.keys()
        for i, species in enumerate(species_list):
            ic.append(initial_conditions[str(species)])

        sol = odeint(self.gradient_function, ic, t, args=(parametrised_flow,species_list))
        variable_names = [str(x) for x in parametrised_flow]

        if plot_name:
            lines = plt.plot(sol)
            plt.legend(iter(lines), variable_names)
            plt.savefig(plot_name + "-simulation.png")
            plt.xlabel("Time")

        return t, sol, variable_names

    @staticmethod
    def gradient_function(X, t, flow, species_list):
        """
        Evaluates the time-derivative of the system so that it can be numerically integrated by
        ``scipy.odeint`` to obtain a simulated time-course.

        :param X: vector of current concentrations (in order given by species_list)
        :param t: current time (ignored, as not needed to evaluate derivatives)
        :param flow: dictionary in which keys are species-names, and values are SympPy expressions for their time derivatives
        :param species_list: list of species names (SymPy objects)
        :return:
        """
        vals = {"t": t}

        for i, species in enumerate(species_list):
            vals[species] = X[i]

        result = []

        for species in flow:
            result.append(flow[species].evalf(subs=vals))

        return result


class SolverCallerISAT(SolverCaller):
    """
    This class is responsible for calling iSAT to solve a SAT-ODE problem, parsing the result, and substituting the
    extracted parameter values into the ``CRNSketch`` to obtain a model that can be simulated
    """
    def __init__(self, model_path="./bellshape.hys", isat_path=""):
        """

        :param model_path: path to the .hys file containing the SAT-ODE problem to be solved
        :param isat_path: path to the iSAT binary
        """
        super(SolverCallerISAT, self).__init__(model_path)

        self.isat_path = isat_path
        if not isat_path:
            self.isat_path = "./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt"

    def optimal_synthesis_decreasing_cost(self, max_cost=35, min_cost=10, precision=0.1, msw=0):
        """
        Call iSAT repeatedly, decreasing the permitted cost by 1 between each iteration.

        :param max_cost: the maximum cost permitted on the first iteration
        :param min_cost: the maximum cost permitted on the final iteration
        :param precision: value of --prabs parameter to eb passed to iSAT
        :param msw: value of --msw parameter to eb passed to iSAT
        :return: list of file names, each containing the output from iSAT from one iteration
        """
        cost = max_cost
        result_file_names = []
        while cost >= min_cost:
            self.edit_cost(cost)
            result_file_name = self.call_solver(precision, cost, ' --ode-opts --continue-after-not-reaching-horizon', msw=msw)
            result_file_names.append(result_file_name)
            cost -= 1
        return result_file_names

    def edit_cost(self, cost):
        """
        Edit the model file to update the MAX_COST limit.
        :param cost: maximum permitted cost - 0 means no limit applied (float)
        """
        with open(self.model_path, 'r') as f:
            lines = f.read().split('\n')

            for line_number, line in enumerate(lines):
                if "define MAX_COST = " in line:
                    lines[line_number] = "define MAX_COST = %s; " % cost

                if "define NO_COST_LIMIT = " in line:
                    if cost == 0:
                        lines[line_number] = "define NO_COST_LIMIT = 1;"
                    else:
                        lines[line_number] = "define NO_COST_LIMIT = 0;"

        with open(self.model_path, 'w') as f:
            f.write('\n'.join(lines))

    def call_solver(self, precision, cost, otherPrams, max_depth=2, msw=0):
        """
        Call iSAT, and save its output to a file.

        :param precision: value of --prabs parameter to be passed to iSAT
        :param cost:  maximum value of cost [determines output file name]
        :param otherPrams: string containing other arguments to pass to iSAT
        :param max_depth: maximum unrolling depth for BMC
        :param msw: value of --msw (minimum splitting width) parameter to pass to iSAT
        :return: name of output file
        """

        if msw == 0:
            msw = precision * 5

        out_file = os.path.join(self.results_dir, "%s_%s_%s-isat.txt" % (self.model_name, cost, precision))
        command = "%s --i %s --prabs=%s --msw=%s --max-depth=%s %s " % \
                  (self.isat_path, self.model_path, precision, msw, max_depth, otherPrams)

        with open(out_file, "w") as f:
            print("Calling solver!\n " + command)
            sub.call(command.split(), stdout=f, stderr=sub.PIPE)

        return out_file


    def getCRNValues(self, file_path):
        """
        Parse the output of iSAT, and extract parameter values and initial conditions.

        Returns a ``constant_values`` dictionary, containing only values of that do not change over time,
        and a ``all_values`` dictionary that also contains the initial value of variables that change over time.

        :param file_path: path to the file containing iSAT output
        """

        p = re.compile(r"(.+?) \(.+?\):")
        p2 = re.compile(r".+?[\[\(](.+?),(.+?)[\]\)].+")

        constant_values = {}
        all_values = {} # includes state variables
        var_name = False
        with open(file_path, "r") as f:
            for line in f:

                if p.match(line):
                    var_name = p.match(line).groups()[0].strip()

                    if "solver" in var_name or "_trigger" in var_name:
                        var_name = False

                elif p2.match(line) and var_name:
                    # save this row's values
                    values = p2.match(line).groups()

                    if var_name not in constant_values.keys():
                        constant_values[var_name] = values
                        all_values[var_name] = values

                    elif constant_values[var_name] != values:
                        # if we've already recorded a different value, it's because value changes between modes
                        # it's not a constant parameter, so don't record it
                        constant_values.pop(var_name, None)

        return constant_values, all_values

    def get_full_solution(self, crn, flow, vals):
        """
        Use values extracted from iSAT output to construct initial conditions dictionary and replace parameters in flow
        dictionary with their numerical values.

        :param crn:
        :param flow:
        :param vals:
        :return:
        """
        initial_conditions = {}

        var_names = [str(var) for var in flow.keys()]

        for val in vals:
            if val in var_names:
                initial_conditions[val] = (float(vals[val][0]) + float(vals[val][1])) / 2
            else:
                for x in flow:
                    mean_val = (float(vals[val][0]) + float(vals[val][1])) / 2
                    flow[x] = flow[x].subs(val, mean_val)

        parametrised_flow = dict(flow)
        for x in crn.derivatives:
            derivative_symbol = sympify(x["name"])
            del parametrised_flow[derivative_symbol]
            # del initial_conditions[str(derivative_symbol)]

        return initial_conditions, parametrised_flow


class SolverCallerDReal(SolverCaller):

    def __init__(self, model_path="./bellshape.drh", dreal_path="/Users/maxtby/local/bin/dreach"):
        """

        :param model_path: path to the .drh file containing the SAT-ODE problem to be solved
        :param dreal_path:
        """
        super(SolverCallerDReal, self).__init__(model_path)
        self.dreal_path = dreal_path

    def optimal_synthesis_decreasing_cost(self, max_cost=35, min_cost=10, precision=0.1, msw=0):
        """
        Call dReal repeatedly, decreasing the permitted cost by 1 between each iteration.

        :param max_cost: the maximum cost permitted on the first iteration
        :param min_cost: the maximum cost permitted on the final iteration
        :param precision: value of --precision parameter to be passed to dReach
        :param msw: ignored
        :return: list of file names, each containing the output from dReach from one iteration
        """

        cost = max_cost
        result_file_names = []
        while cost >= min_cost:
            self.edit_cost(cost)
            result_file_name = self.call_solver(precision, cost, '')
            result_file_names.append(result_file_name)
            cost -= 1
        return result_file_names

    def edit_cost(self, cost):
        """
        Edit the model file to update the MAX_COST limit.
        :param cost: maximum permitted cost - 0 means no limit applied (float)
        """

        with open(self.model_path, 'r') as f:
            lines = f.read().split('\n')

            for line_number, line in enumerate(lines):
                if "define MAX_COST = " in line:
                    lines[line_number] = "#define MAX_COST %s" % cost

                if "define NO_COST_LIMIT = " in line:
                    if cost == 0:
                        lines[line_number] = "#define NO_COST_LIMIT 1"
                    else:
                        lines[line_number] = "#define NO_COST_LIMIT 0;"

        with open(self.model_path, 'w') as f:
            f.write('\n'.join(lines))

    def call_solver(self, precision, cost, otherPrams, max_depth=2):
        """
        Call dReach, and save its output to a file.

        :param precision: value of --precision parameter to be passed to dReach
        :param cost:  maximum value of cost [determines output file name]
        :param otherPrams: string containing other arguments to pass to dReach
        :param max_depth: maximum unrolling depth for BMC
        :return: name of output file
        """
        out_file = os.path.join(self.results_dir, "%s_%s_%s-dreal.txt" % (self.model_name, cost, precision))
        command = "%s -k %s %s --precision %s %s" % \
                  (self.dreal_path, max_depth, self.model_path, precision, otherPrams)

        with open(out_file, "w") as f:
            print("Calling solver!\n " + command)
            sub.call(command.split(), stdout=f, stderr=sub.PIPE)

        return out_file

    def getCRNValues(self, file_path):
        """
        Parse the output of dReach, and extract parameter values and initial conditions.

        Returns a ``constant_values`` dictionary, containing values of that do not change over time, and a ``all_values``
        dictionary that contains the initial value of variables that change over time.

        :param file_path: path to the file containing dReach output
        """

        results = ''
        with open(file_path) as f:
            results = f.read()

        return results

        # constant_values = {}
        # all_values = {}  # includes state variables
        #
        # for t in results["traces"][0]:
        #     var_name = t["key"].replace("_0_0", "")
        #
        #     interval = t["values"][0]["enclosure"]
        #
        #     single_value = True
        #
        #     for v in t["values"]:
        #         # if i[0] != interval[0] or i[1] != interval[1] :
        #         if v["enclosure"] != interval:
        #             single_value = False
        #             break
        #
        #     if single_value:
        #         constant_values[var_name] = interval
        #     all_values[var_name] = interval
        #
        # return constant_values, all_values
