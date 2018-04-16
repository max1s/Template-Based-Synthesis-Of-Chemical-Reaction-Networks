import subprocess as sub
import re
import io
import os
from os import walk
from sympy import sympify
from scipy.integrate import odeint
import numpy as np
import json

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

class SolverCaller(object):
    def __init__(self, model_path="./bellshape.hys"):
        self.model_path = model_path

        self.results_folder = "/bellshaperesults"

        directory_name, file_name = os.path.split(model_path)
        self.model_name, _ = os.path.splitext(file_name)

        self.results_dir = os.path.join(directory_name, "results")

        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

    def single_synthesis(self, cost=20, precision=0.1):
        return self.optimal_synthesis_decreasing_cost(max_cost=cost, min_cost=cost, precision=precision)

    def optimal_synthesis_decreasing_cost(self, max_cost=35, min_cost=10, precision=0.1):
        pass

    def simulate_solutions(self, initial_conditions, parametrised_flow, t=False, plot_name=""):
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
        vals = {"t": t}

        for i, species in enumerate(species_list):
            vals[species] = X[i]

        result = []

        for species in flow:
            result.append(flow[species].evalf(subs=vals))

        return result


class SolverCallerISAT(SolverCaller):
    def __init__(self, model_path="./bellshape.hys", isat_path=""):
        super(SolverCallerISAT, self).__init__(model_path)

        self.isat_path = isat_path
        if not isat_path:
            self.isat_path = "./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt"

    def optimal_synthesis_decreasing_cost(self, max_cost=35, min_cost=10, precision=0.1):
        cost = max_cost
        result_file_names = []
        while cost >= min_cost:
            self.edit_cost(cost)
            result_file_name = self.call_solver(precision, cost, ' --ode-opts --continue-after-not-reaching-horizon')
            result_file_names.append(result_file_name)
            cost -= 1
        return result_file_names

    def edit_cost(self, cost):
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

    def call_solver(self, precision, cost, otherPrams, max_depth=2):
        out_file = os.path.join(self.results_dir, "%s_%s_%s-isat.txt" % (self.model_name, cost, precision))
        command = "%s --i %s --prabs=%s --msw=%s --max-depth=%s %s " % \
                  (self.isat_path, self.model_path, precision, precision * 5, max_depth, otherPrams)

        with open(out_file, "w") as f:
            print("Calling solver!\n " + command)
            sub.call(command.split(), stdout=f, stderr=sub.PIPE)

        return out_file


    def getCRNValues(self, file_path):
        # Return a dictionary containing the ranges for each parameter as returned from the solver
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
        initial_conditions = {}

        var_names = [str(var) for var in flow.keys()]

        for val in vals:
            if val in var_names:
                initial_conditions[val] = (float(vals[val][0]) + float(vals[val][1])) / 2
            else:
                for x in flow:
                    mean_val = (float(vals[val][0]) + float(vals[val][1])) / 2
                    flow[x] = flow[x].subs(val, mean_val)
                    print "Now", flow[x]

        parametrised_flow = dict(flow)
        for x in crn.derivatives:
            derivative_symbol = sympify(x["name"])
            del parametrised_flow[derivative_symbol]
            # del initial_conditions[str(derivative_symbol)]

        return initial_conditions, parametrised_flow


class SolverCallerDReal(SolverCaller):

    def __init__(self, model_path="./bellshape.hys", dreal_path="/Users/maxtby/local/bin/dreach"):
        super(SolverCallerDReal, self).__init__(model_path)
        self.dreal_path = dreal_path

    def optimal_synthesis_decreasing_cost(self, max_cost=35, min_cost=10, precision=0.1):
        cost = max_cost
        result_file_names = []
        while cost >= min_cost:
            self.edit_cost(cost)
            result_file_name = self.call_solver(precision, cost, '')
            result_file_names.append(result_file_name)
            cost -= 1
        return result_file_names


    def edit_cost(self, cost):
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
        out_file = os.path.join(self.results_dir, "%s_%s_%s-dreal.txt" % (self.model_name, cost, precision))
        command = "%s -k %s %s --precision %s %s" % \
                  (self.dreal_path, max_depth, self.model_path, precision, otherPrams)

        with open(out_file, "w") as f:
            print("Calling solver!\n " + command)
            sub.call(command.split(), stdout=f, stderr=sub.PIPE)

        return out_file

    def getCRNValues(self, file_path):
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
