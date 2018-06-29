from sympy import integrate, Symbol, sympify

"""

This sub-module is responsible for actually generating a iSAT (.hys) or dReach (.drh) file
that encodes the CRN Synthesis problem as an SMT-ODE problem.

The public interfaces are the ``constructISAT`` or ``constructdReal`` functions.
These will construct various objects to represent part of the problem
(``Declaration``, ``Initial``, ``Transition``, ``Post``),
call the ``constructISAT`` or ``constructdReal`` methods on them, assemble the results, and return a string.

"""

class Declaration:
    """
    This class contains all of the variable declarations necessary to construct a iSAT (.hys) or dReach (.drh) file
    encoding the CRN Synthesis problem as an SMT-ODE problem.
    """
    def __init__(self, crn, numModes, flows):
        """

        :param crn: a CRNSketch object
        :param numModes: the number of modes in the specification (int)
        :param flows: dictionary in which keys are species-names, and values are SympPy expressions for their time derivatives
        """
        self.crn = crn
        self.numModes = numModes
        self.flows = flows

    def constructiSAT(self, max_time=1, scale_factor=1):
        """
        Returns a string containing variable definitions in iSat (.hys) format.

        :param max_time: value for MAX_TIME by which goal condition must be reached
        :param scale_factor: scale factor for time (multiplicative factor for all rates)
        """
        s = "\nDECL \n"

        s += "define MAX_TIME = %s;\n" % max_time
        s += "\tdefine SF = %s;\n" % scale_factor

        s += "\t-- declare time variables\n"
        s += "\tfloat [0, MAX_TIME] time;\n"
        s += "\tfloat [0, MAX_TIME] delta_time;\n\n"
        if len(self.crn.input_species) > 0:
            s += "\tfloat [0, MAX_TIME] inputTime;\n\n"


        s += "\t-- declare cost variables\n"
        s += "\tdefine MAX_COST = 1000;\n"
        s += "\tdefine NO_COST_LIMIT = 0;\n"

        if len(self.crn.real_species) > 0:
            s += "\n\t-- Define State Variables\n"
        for d in self.crn.real_species:
            s += d.iSATDefinition()

        if len(self.crn.input_species) > 0:
            s += "\n\t-- Define Input Variables\n"
        for d in self.crn.input_species:
            s += d.iSATDefinition()

        if len(self.crn.derivatives) > 0:
            s += "\n\t-- Define Derivative Variables\n"
        for d in self.crn.derivatives:
            s += "\tfloat [-4, 4] %s;\n" % d["name"]

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t-- Lambda Variables\n"
        for lam in self.crn.lambda_variables:
            s += lam.iSATDefinition() + "\n"

        if len(self.crn.choice_variables) > 0:
            s += "\n\t-- Choice Variables\n"
        for c in self.crn.choice_variables:
            s += c.iSATDefinition() + "\n"

        if len(self.crn.joint_choice_variables) > 0:
            s += "\n\t-- Joint choice Variables\n"
        for c in self.crn.joint_choice_variables:
            s += c.iSATDefinition() + "\n"

        if len(self.crn.optionalReactions) > 0:
            s += "\n\t-- Optional Reaction Variables\n"
        for optional_reaction in self.crn.optionalReactions:
            s += "\tfloat [0, 1] %s;\n" % optional_reaction.variable_name

        if len(self.crn.getRateConstants()) > 0:
            s += "\n\t-- Rate constants\n"
        for rate in self.crn.getRateConstants():
            s += "\tfloat [%s, %s] %s;\n" % (rate.min, rate.max, rate.name)

        s += "\n\t--Define modes\n"
        for i in range(1, self.numModes + 1):
            s += "\tboole mode_%s;\n" % i

        return s

    def constructdReal(self, max_time=1, scale_factor=1):
        """
        Returns a string containing variable definitions in dReach (.drh) format.

        :param max_time: value for MAX_TIME by which goal condition must be reached
        :param scale_factor: scale factor for time (multiplicative factor for all rates)
        """

        s = "\n #define MAX_TIME %s\n" % max_time
        s += "\t#define SF %s\n\n" % scale_factor

        s += "\t// declare cost variables\n"
        s += "\t#define MAX_COST 1000\n"
        s += "\t#define NO_COST_LIMIT 0\n\n"

        s += "\t// define derivatives\n"
        for x in self.flows:
            f = x.constructdRealDerivativeDefinitions()
            if f:
                s += "\t%s\n" % f
        s += "\n"

        s += "\t// declare time variables\n"
        s += "\t[0, MAX_TIME] time;\n"
        if len(self.crn.input_species) > 0:
            s += "\t[0, MAX_TIME] inputTime;\n\n"

        if len(self.crn.real_species) > 0:
            s += "\n\t//Define State Variables\n"
        for d in self.crn.real_species:
            s += d.dRealDefinition()

        if len(self.crn.input_species) > 0:
            s += "\n\t// Define Input Variables\n"
        for d in self.crn.input_species:
            s += d.dRealDefinition()

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t// Lambda Variables\n"
        for lam in self.crn.lambda_variables:
            s += lam.dRealDefinition() + "\n"

        if len(self.crn.choice_variables) > 0:
            s += "\n\t// Choice Variables\n"
        for c in self.crn.choice_variables:
            s += c.dRealDefinition() + "\n"

        if len(self.crn.joint_choice_variables) > 0:
            s += "\n\t// Joint choice Variables\n"
        for c in self.crn.joint_choice_variables:
            s += c.dRealDefinition() + "\n"

        if len(self.crn.optionalReactions) > 0:
            s += "\n\t// Optional Reaction Variables\n"
        for optional_reaction in self.crn.optionalReactions:
            s += "\t[0, 1] %s;\n" % optional_reaction.variable_name

        if len(self.crn.getRateConstants()) > 0:
            s += "\n\t// Rate constants\n"
        for rate in self.crn.getRateConstants():
            s += "\t [%s, %s] %s;\n" % (rate.min, rate.max, rate.name)

        return s


class Transition:
    """
    This class contains all of the conditions describing mode transitions necessary to construct a iSAT (.hys)
     or dReach (.drh) file encoding the CRN Synthesis problem as an SMT-ODE problem.
    """

    def __init__(self, crn, flows, modes):
        """

        :param crn: a CRNSketch object
        :param flows: dictionary in which keys are species-names, and values are SympPy expressions for their time derivatives
        :param modes: list of tuples describing each mode. The first element is each tuple is ignored; the second encodes a condition that must hold at all times during that mode; the third encodes a condition that must hold at the time that the system transitions out of that mode.
        """

        self.crn = crn
        self.modes = modes
        self.flow = flows

    def constructiSAT(self):
        """
        Returns a string describing mode transitions in iSat (.hys) format.
        """
        s = "\n\nTRANS \n\n"

        s += "\t-- time constraint\n"
        s += "\ttime' = time + delta_time;\n\n"

        # s += "\t-- must progress through modes in order\n"
        # s += "\t mode_1' -> mode_1;\n"
        # for i in list(range(2, len(self.modes)+1)):
        #     s += "\tmode_%s' -> (mode_%s or mode_%s);\n" % (i, i, i-1)

        s += "\n\n\t-- invariant conditions during modes\n"
        for mode_index, mode in enumerate(self.modes):
            if mode[1]:
                s += "\tmode_%s  -> (%s);\n" % (mode_index+1, mode[1])
                s += "\tmode_%s'  -> (%s);\n" % (mode_index+1, mode[1])

        s += "\n\t-- jump conditions between modes\n"
        for mode_index, mode in enumerate(self.modes[:-1]):
            if len(mode) > 2 and mode[2]: # post-condition on mode
                s += "\t(mode_%s and mode_%s') -> (%s);\n" % (mode_index+1, mode_index+2, mode[2])

        s += "\n\t-- No state change without time consumption.\n"

        terms = ["(%s' = %s)" % (x.name, x.name) for x in self.crn.species]
        if len(self.crn.input_species) > 0:
            terms.append("(inputTime' = inputTime)")

        for c in self.crn.choice_variables:
            terms = ["(%s_%s' = %s_%s)" % (c.name, i, c.name, i) for i in list(range(c.minValue, c.maxValue + 1))]
            s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        if (len(terms) > 0):
            s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = ["(%s' = %s)" % (x.name, x.name) for x in self.crn.getRateConstants()]
        if (len(terms) > 0):
            s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = ["(%s' = %s)" % (x.variable_name, x.variable_name) for x in self.crn.optionalReactions]
        if (len(terms) > 0):
            s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = []
        for lambda_choice in self.crn.lambda_variables:
            terms.extend(["(%s' = %s)" % (x, x) for x in lambda_choice.lambdas])
        if(len(terms) > 0):
            s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)
            s += "\n"

        if len(self.crn.getRateConstants()) > 0:
            s += "\n\t-- Rate constants are fixed\n"
        for rate in self.crn.getRateConstants():
            s += "\t(d.%s/d.time = 0);\n" % rate.name
            s += "\t(%s = %s');\n" % (rate.name, rate.name)

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t-- Lambda variables are fixed\n"
        for lam in self.crn.lambda_variables:
            s += "".join(["\t(d.%s/d.time = 0);\n" % x.name for x in lam.lambdas])
            s += "\t" + " and ".join(["(%s = %s')" % (x.name, x.name) for x in lam.lambdas]) + ";\n"

        if len(self.crn.choice_variables) > 0:
            s += "\n\t-- Choice variables are fixed\n"
        for c in self.crn.choice_variables:
            values = list(range(c.minValue, c.maxValue + 1))
            for i in values:
                s += "\t(d.%s_%s/d.time = 0);\n" % (c.name, i)
            s += "\t" + " and ".join(["(%s_%s = %s_%s')" % (c.name, i, c.name, i) for i in values]) + ";\n"

        if len(self.crn.joint_choice_variables) > 0:
            s += "\n\t-- Joint choice variables are fixed\n"
        for c in self.crn.joint_choice_variables:
            for v in c.list_decision_variables():
                s += "\t(d.%s/d.time = 0);\n" % v
            s += "\t" + " and ".join(["(%s = %s')" % (v, v) for v in c.list_decision_variables()]) + ";\n"

        if len(self.crn.optionalReactions) > 0:
            s += "\n\t-- Optional reaction variables are fixed\n"
        for c in self.crn.optionalReactions:
            s += "\t(d.%s/d.time = 0);\n" % c.variable_name
            s += "\t(%s = %s');\n" % (c.variable_name, c.variable_name)

        numModes = max(1, len(self.modes))
        mode_list = ["mode_" + str(x) for x in range(1, numModes + 1)]
        modes_string = " or ".join(mode_list)

        s += "\n\n\t-- Flows\n"
        s += ''.join(['\t(%s) -> %s;\n' % (modes_string, x.constructiSAT()) for x in self.flow])
        if len(self.crn.input_species) > 0:
            s += '\t(%s) -> (d.inputTime/d.time  = 1);\n' % modes_string

        return s

    def constructdReal(self):
        """
        Returns a string describing mode transitions in dReach (.drh) format.
        """
        s = ''
        if len(self.modes) is not 0:
            for mode_index in range(0, len(self.modes)):
                mode = self.modes[mode_index]
                s += "{ mode %s;\n\n" % (mode_index+1)

                #mode invariants
                s += "invt: "
                s += "\n\n\t // invariant conditions during modes\n"

                if mode[1]:
                    s += "\t " + ''.join(mode[1]) + ";"

                #flow variables
                s += '\n\nflow: \n'
                if len(self.crn.getRateConstants()) > 0:
                    s += "\n\t // Rate constants are fixed\n"
                for rate in self.crn.getRateConstants():
                    s += "\td/dt[%s] = 0;\n" % rate.name

                if len(self.crn.lambda_variables) > 0:
                    s += "\n\t // Lambda variables are fixed\n"
                for lam in self.crn.lambda_variables:
                    s += "".join(["\td/dt[%s] = 0;\n" % x.name for x in lam.lambdas])

                if len(self.crn.choice_variables) > 0:
                    s += "\n\t // Choice variables are fixed\n"
                for c in self.crn.choice_variables:
                    values = list(range(c.minValue, c.maxValue + 1))
                    for i in values:
                        s += "\td/dt[%s_%s] = 0;\n" % (c.name, i)

                if len(self.crn.joint_choice_variables) > 0:
                    s += "\n\t // Joint choice variables are fixed\n"
                for c in self.crn.joint_choice_variables:
                    for v in c.list_decision_variables():
                        s += "\td/dt[%s] = 0;\n" % v


                if len(self.crn.optionalReactions) > 0:
                    s += "\n\t // Optional reaction variables are fixed\n"
                for c in self.crn.optionalReactions:
                    s += "\td/dt[%s] = 0;\n" % c.variable_name


                s += "\n\n\t// Flows\n"
                for x in self.flow:
                    f = x.constructdReal()
                    if f:
                        s += '\t%s\n' % (f)

                if len(self.crn.input_species) > 0:
                    s += "\td/dt[inputTime] = 1;"

                #mode jump
                terms = ["(%s' = %s)" % (x.name, x.name) for x in self.crn.real_species]

                for c in self.crn.choice_variables:
                    terms.extend(["(%s_%s' = %s_%s)" % (c.name, i, c.name, i) for i in list(range(c.minValue, c.maxValue + 1))])

                terms.extend(["(%s' = %s)" % (x.name, x.name) for x in self.crn.getRateConstants()])
                terms.extend(["(%s' = %s)" % (x.variable_name, x.variable_name) for x in self.crn.optionalReactions])

                for lambda_choice in self.crn.lambda_variables:
                    terms.extend(["(%s' = %s)" % (x, x) for x in lambda_choice.lambdas])




                s += "\n\njump: "

                if (mode_index + 1) < len(self.modes):
                    s += "\n\n\t // jump conditions during modes\n"
                    invariants = "\t ( and " + "".join(terms) + ");"
                    s += "%s ==> @%s %s\n" % (mode[2], mode_index+2, invariants)


                s += "\n\n }"

        else:
            s += "{mode 1; \n"
            s += "invt: true; \n"
            s += "flow: "
            s += "\t"

            if len(self.crn.getRateConstants()) > 0:
                s += "\n\t // Rate constants are fixed\n"
            for rate in self.crn.getRateConstants():
                s += "\td/dt[%s] = 0;\n" % rate.name

            if len(self.crn.lambda_variables) > 0:
                s += "\n\t // Lambda variables are fixed\n"
            for lam in self.crn.lambda_variables:
                s += "".join(["\td/dt[%s] = 0;\n" % x.name for x in lam.lambdas])

            if len(self.crn.choice_variables) > 0:
                s += "\n\t // Choice variables are fixed\n"
            for c in self.crn.choice_variables:
                values = list(range(c.minValue, c.maxValue + 1))
                for i in values:
                    s += "\td/dt[%s_%s] = 0;\n" % (c.name, i)

            if len(self.crn.joint_choice_variables) > 0:
                s += "\n\t // Joint choice variables are fixed\n"
            for c in self.crn.joint_choice_variables:
                for v in c.list_decision_variables():
                    s += "\td/dt[%s] = 0;\n" % v

            if len(self.crn.optionalReactions) > 0:
                s += "\n\t // Optional reaction variables are fixed\n"
            for c in self.crn.optionalReactions:
                s += "\td/dt[%s] = 0;\n" % c.variable_name

            s += "\n\n\t// Flows\n"
            for x in self.flow:
                f = x.constructdReal()
                if f:
                    s += '\t%s\n' % (f)

            # mode jump
            terms = ["(%s' = %s)" % (x.name, x.name) for x in self.crn.real_species]

            for c in self.crn.choice_variables:
                terms.extend(
                    ["(%s_%s' = %s_%s)" % (c.name, i, c.name, i) for i in list(range(c.minValue, c.maxValue + 1))])

            terms.extend(["(%s' = %s)" % (x.name, x.name) for x in self.crn.getRateConstants()])
            terms.extend(["(%s' = %s)" % (x.variable_name, x.variable_name) for x in self.crn.optionalReactions])

            for lambda_choice in self.crn.lambda_variables:
                terms.extend(["(%s' = %s)" % (x, x) for x in lambda_choice.lambdas])

            s += "\n\njump: \n\n }"


        return s

class Initial:
    """
    This class represents the initial condition for the SMT-ODE problem.
    """

    def __init__(self, crn, numModes, other_constraints):
        """

        :param crn: a CRNSketch object
        :param numModes: the number of modes in the specification (int)
        :param other_constraints: other user-constraints (string)
        """
        self.crn = crn
        self.numModes = numModes
        self.other_constraints = other_constraints

    def constructiSAT(self):
        """
        Returns a string describing initial conditions in iSAT (.hys) format.
        """
        s = "\nINIT\n"

        s += "\ttime = 0;\n\n"

        if len(self.crn.input_species) > 0:
            s += "\tinputTime = 0;\n\n"


        s += "\t-- cost condition\n"
        s += "\t(((%s) <= MAX_COST) or (NO_COST_LIMIT = 1));\n\n" % self.crn.get_cost()

        # mode exclusion
        mode_list = ["mode_%s" % x for x in range(1, self.numModes + 1)]
        s += "\t-- cannot be in two modes at the same time. We start in mode_1.\n"
        s += "\tmode_1 = 1;\n"
        s += "\t" + " + ".join(mode_list) + " = 1;\n\n"

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t-- Integer encoding of lambda variables\n"
        for lam in self.crn.lambda_variables:
            s += lam.iSATformat_constraint()

        if len(self.crn.choice_variables) > 0:
            s += "\n\t-- Integer encoding of choice variables\n"
        for c in self.crn.choice_variables:
            s += c.iSATformat_constraint()

        if len(self.crn.choice_variables) > 0:
            s += "\n\t-- Integer encoding of optional reaction variables\n"
        for optional_reaction in self.crn.optionalReactions:
            s += "\t(%s = 0) or (%s = 1);\n" % (optional_reaction.variable_name, optional_reaction.variable_name)

        if len(self.crn.joint_choice_variables) > 0:
            s += "\n\t-- Integer encoding of joint choice variables\n"
        for joint_choice_variables in self.crn.joint_choice_variables:
            s += joint_choice_variables.iSATformat_constraint()


        if len(self.crn.real_species) > 0:
            s += "\n\t-- Limits on initial conditions\n"
        for sp in self.crn.real_species:
            s += sp.iSATInitialization()

        if len(self.crn.input_species) > 0:
            s += "\n\t-- Initial conditions of inputs\n"
        for sp in self.crn.input_species:
            s += sp.iSATInitialization()

        if self.other_constraints:
            s += "\n\t-- Manually specified constraints\n"
            s += "\t%s\n" % self.other_constraints


        return s

    def constructdReal(self):
        """
        Returns a string describing initial conditions in dReach (.drh) format.
        """
        s = "\ninit:\n\n"

        s += " @1 (and\n"
        s += "\t // cost condition\n"
        s += "\t( (or ( (%s) <= MAX_COST) (NO_COST_LIMIT = 1)))\n\n" % self.crn.get_cost()

        if len(self.crn.input_species) > 0:
            s += "\t(and (inputTime = 0))\n\n"

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t// Integer encoding of lambda variables\n"
        for lam in self.crn.lambda_variables:
            s += lam.dRealformat_constraint()

        if len(self.crn.choice_variables) > 0:
            s += "\n\t // Integer encoding of choice variables\n"
        for c in self.crn.choice_variables:
            s += c.dRealformat_constraint()
        if len(self.crn.choice_variables) > 0:
            s += "\n\t // Integer encoding of optional reaction variables\n"
        for optional_reaction in self.crn.optionalReactions:
            s += "\t(or (%s = 0) (%s = 1))\n" % (optional_reaction.variable_name, optional_reaction.variable_name)

        if len(self.crn.joint_choice_variables) > 0:
            s += "\n\t // Integer encoding of joint choice variables\n"
        for joint_choice_variables in self.crn.joint_choice_variables:
            s += joint_choice_variables.dRealformat_constraint()


        if len(self.crn.real_species) > 0:
            s += "\n\t // Limits on initial conditions\n"
        for sp in self.crn.real_species:
            s += sp.dRealInitialization()

        if len(self.crn.input_species) > 0:
            s += "\n\t // Initial conditions of inputs\n"
        for sp in self.crn.input_species:
            s += sp.dRealInitialization()

        if self.other_constraints:
            s += "\n\t //Manually specified constraints\n"
            s += "\t%s\n" % self.other_constraints

        s += ");"
        return s


class Flow:
    """
    This class represents the dynamics of a single species.
    """

    def __init__(self, var, t, flow, crn):
        """

        :param var: variable representing the species (SymPy object?)
        :param t: string representing time
        :param flow: SymPy expression represnting time-derivative of variable
        :param crn: a CRNSketch object
        """
        self.variable = var
        self.time = t
        self.flow = flow.subs(sympify('t'), sympify('inputTime'))
        self.crn = crn

    def constructiSAT(self):
        """
        Returns a string representing the dynamics of a single species in iSAT (.hys) format.
        """

        # Python represents powers as a**b, whereas iSAT uses a^b
        flow = str(self.flow).replace('**', '^')

        derivative_names = [x["name"] for x in self.crn.derivatives]

        if str(self.variable) in derivative_names:
            return "\t(%s = %s)" % (self.variable, flow)
        else:
            return "\t(d.%s/d.%s  = %s)" % (self.variable, self.time, flow)

    def constructdReal(self):
        """
        Returns a string representing the dynamics of a single species in dReach (.drh) format.
        """
        flow = str(self.flow).replace('**', '^')

        derivative_names = [x["name"] for x in self.crn.derivatives]

        if str(self.variable) not in derivative_names:
            return "\td/dt[%s]  = %s;" % (self.variable, flow)

    def constructdRealDerivativeDefinitions(self):
        
        #flow = integrate(self.flow, Symbol('t')).doit()
        flow = self.flow

        flow = str(flow).replace('**', '^')

        derivative_names = [x["name"] for x in self.crn.derivatives]

        if str(self.variable) in derivative_names:
            return "\t#define %s  (%s)" % (self.variable, flow)


class Post:
    """
    This class represents the final target/goal condition that applies at the end-time of the SMT-ODE problem.
    """

    def __init__(self, t, modes):
        """

        :param t: string representing time
        :param modes: list of tuples describing each mode. The first element is each tuple is ignored; the second encodes a condition that must hold at all times during that mode; the third encodes a condition that must hold at the time that the system transitions out of that mode.
        """
        self.time = t
        self.modes = modes

    def constructiSAT(self):
        """
        Returns a string representing the final target/goal condition in iSAT (.hys) format.
        """
        if len(self.modes) > 0 and len(self.modes[-1]) > 2 and self.modes[-1][2]:
            post_condition = "and (" + self.modes[-1][2] + ")"
        else:
            post_condition = ""

        s = "\nTARGET \n"
        s += "\tmode_%s and (time <= %s) %s;\n" % (len(self.modes), self.time, post_condition)
        return s

    def constructdReal(self):
        """
        Returns a string representing the final target/goal condition in dReach (.drh) format.
        """

        if len(self.modes) > 0 and len(self.modes[-1]) > 2 and self.modes[-1][2]:
            post_condition = "and (" + self.modes[-1][2] + ")"
        else:
            post_condition = "true"

        s = "\n\ngoal: \n"
        s += "\t@%s %s;\n" % (len(self.modes), post_condition) if (len(self.modes) > 0) else  "\t@%s %s;\n" % (len(self.modes) + 1, post_condition)

        return s

def constructISAT(crn, modes, flow, other_constraints=False, scale_factor=1, max_time=1):
    """
    Returns a string in iSAT (.hys) format encoding the CRN Synthesis problem as an SMT-ODE problem.

    :param crn: a CRNSketch object
    :param modes: list of tuples describing each mode. The first element is each tuple is ignored; the second encodes a condition that must hold at all times during that mode; the third encodes a condition that must hold at the time that the system transitions out of that mode.
    :param flow: dictionary in which keys are species-names, and values are SympPy expressions for their time derivatives
    :param other_constraints: other user-constraints (string)
    """
    m_flow = [Flow(x, 'time', y, crn) for x, y in flow.items()]
    numModes = max(1, len(modes))

    d = Declaration(crn, numModes, []).constructiSAT(scale_factor=scale_factor, max_time=max_time)
    i = Initial(crn, numModes, other_constraints).constructiSAT()
    t = Transition(crn, m_flow, modes).constructiSAT()
    p = Post(max_time, modes).constructiSAT()  # TODO: set maxtime

    return d + i + t + p


def constructdReal(crn, modes, flow, other_constraints=False, scale_factor=1, max_time=1):
    """

    Returns a string in dReach (.drh) format encoding the CRN Synthesis problem as an SMT-ODE problem.

    :param crn: a CRNSketch object
    :param modes: list of tuples describing each mode. The first element is each tuple is ignored; the second encodes a condition that must hold at all times during that mode; the third encodes a condition that must hold at the time that the system transitions out of that mode.
    :param flow: dictionary in which keys are species-names, and values are SympPy expressions for their time derivatives
    :param other_constraints: other user-constraints (string)
    """
    m_flow = [Flow(x, 'time', y, crn) for x, y in flow.items()]
    numModes = max(1, len(modes))
    d = Declaration(crn, numModes, m_flow).constructdReal(scale_factor=scale_factor, max_time=max_time)
    i = Initial(crn, numModes, other_constraints).constructdReal()
    t = Transition(crn, m_flow, modes).constructdReal()
    p = Post(max_time, modes).constructdReal()
    return d + t + i + p
