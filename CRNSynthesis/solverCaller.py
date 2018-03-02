import subprocess as sub
import re
import io
import os
from os import walk
from sympy import sympify

class SolverCaller():
    
    def __init__(self, model_path="./bellshape.hys", isat_path=""):

        self.isat_path = isat_path
        if not isat_path:
            self.isat_path = "./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt"

        self.results_folder = "/bellshaperesults"
        self.model_path = model_path

        directory_name, file_name = os.path.split(model_path)
        self.model_name, _ = os.path.splitext(file_name)

        self.results_dir = os.path.join(directory_name, "results")

        if not os.path.exists(self.results_dir):
            os.makedirs(self.results_dir)

    def single_synthesis(self, cost=20, precision=0.1):
        return self.optimal_synthesis_decreasing_cost(max_cost=cost, min_cost=cost, precision=precision)

    def optimal_synthesis_decreasing_cost(self, max_cost=35, min_cost=10, precision=0.1):
        cost = max_cost
        result_file_names = []
        while cost >= min_cost:
            self.editCost(cost)
            result_file_name = self.costCallSolver(precision, cost, ' --ode-opts --continue-after-not-reaching-horizon')
            result_file_names.append(result_file_name)
            cost -= 1
        return result_file_names


    # This function edits the cost within the " + self.model_name + " file such that we can automate the synthesis experiment
    # If the cost is 0 we comment out references to cost in the .hys file.
    def editCost(self, cost):
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

    def costCallSolver(self, precision, cost, otherPrams, max_depth=2):
        out_file = os.path.join(self.results_dir, "%s_%s_%s.txt" % (self.model_name, cost, precision))
        command = "%s --i %s --prabs=%s --msw=%s --max-depth=%s %s " % \
                  (self.isat_path, self.model_path, precision, precision * 5, max_depth, otherPrams)

        with open(out_file, "w") as f:
            print("Calling solver!\n " + command)
            sub.call(command.split(), stdout=f, stderr=sub.PIPE)

        return out_file


    def getCRNValues(self, file_path):
        # Return a dictionary containing the ranges for each parameter as returned from the solver
        p = re.compile(r"(.+?) \(.+?\):")
        p2 = re.compile(r".+?\[(.+?),(.+?)].+")

        var_values = {}
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

                    if var_name not in var_values.keys():
                        var_values[var_name] = values
                        all_values[var_name] = values

                    elif var_values[var_name] != values:
                        # if we've already recorded a different value, it's because value changes between modes
                        # it's not a constant parameter, so don't record it
                        var_values.pop(var_name, None)

        return var_values, all_values

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

        parametrised_flow = dict(flow)
        for x in crn.derivatives:
            derivative_symbol = sympify(x["name"])
            del parametrised_flow[derivative_symbol]
            del initial_conditions[str(derivative_symbol)]

        return initial_conditions, parametrised_flow

