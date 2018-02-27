class Declaration:
    def __init__(self, crn, numModes):
        self.crn = crn
        self.numModes = numModes

    def constructiSAT(self):
        s = "\nDECL \n"

        s += "define MAX_TIME = 1;\n"  # TODO: set this sensibly

        s += "\t-- declare time variables\n"
        s += "\tfloat [0, MAX_TIME] time;\n"
        s += "\tfloat [0, MAX_TIME] delta_time;\n\n"

        s += "\t-- declare cost variables\n"
        s += "\tdefine MAX_COST = 100;\n"
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
            s += "\tfloat [0, 10] %s;\n" % d["name"]

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

    def constructdReal(self):

        s = "\n #define MAX_TIME 1\n"  # TODO: set this sensibly

        s += "\t# declare cost variables\n"
        s += "\tdefine MAX_COST 100\n"
        s += "\tdefine NO_COST_LIMIT 0\n\n"

        s += "\t# declare time variables\n"
        s += "\t[0, MAX_TIME] time;\n"


        if len(self.crn.real_species) > 0:
            s += "\n\t #Define State Variables\n"
        for d in self.crn.real_species:
            s += d.dRealDefinition()

        if len(self.crn.input_species) > 0:
            s += "\n\t  #Define Input Variables\n"
        for d in self.crn.input_species:
            s += d.dRealDefinition()

        if len(self.crn.derivatives) > 0:
            s += "\n\t #Define Derivative Variables\n"
        for d in self.crn.derivatives:
            s += "\t[-4, 4] %s;\n" % d["name"]

        if len(self.crn.lambda_variables) > 0:
            s += "\n\t #Lambda Variables\n"
        for lam in self.crn.lambda_variables:
            s += lam.dRealDefinition() + "\n"

        if len(self.crn.choice_variables) > 0:
            s += "\n\t #Choice Variables\n"
        for c in self.crn.choice_variables:
            s += c.dRealDefinition() + "\n"

        if len(self.crn.joint_choice_variables) > 0:
            s += "\n\t #Joint choice Variables\n"
        for c in self.crn.joint_choice_variables:
            s += c.dRealDefinition() + "\n"

        if len(self.crn.optionalReactions) > 0:
            s += "\n\t #Optional Reaction Variables\n"
        for optional_reaction in self.crn.optionalReactions:
            s += "\t[0, 1] %s;\n" % optional_reaction.variable_name

        if len(self.crn.getRateConstants()) > 0:
            s += "\n\t #Rate constants\n"
        for rate in self.crn.getRateConstants():
            s += "\t [%s, %s] %s;\n" % (rate.min, rate.max, rate.name)

        return s


class Transition:
    def __init__(self, crn, flows, modes):
        self.crn = crn

        self.flow = flows
        self.modes = modes

    def constructiSAT(self):
        s = "\n\nTRANS \n\n"

        s += "\t-- time constraint\n"
        s += "\ttime' = time + delta_time;\n\n"

        s += "\t-- must progress through modes in order\n"
        s += "\t mode_1' -> mode_1;\n"
        for i in list(range(2, len(self.modes)+1)):
            s += "\tmode_%s' -> (mode_%s or mode_%s);\n" % (i, i, i-1)

        s += "\n\n\t-- invariant conditions during modes\n"
        for mode_index, mode in enumerate(self.modes):
            s += "\tmode_%s  -> (%s);\n" % (mode_index+1, mode[1])
            s += "\tmode_%s'  -> (%s);\n" % (mode_index+1, mode[1])

        s += "\n\t-- jump conditions between modes\n"
        for mode_index, mode in enumerate(self.modes[:-1]):
            if len(mode) > 2 and mode[2]: # post-condition on mode
                s += "\t(mode_%s and mode_%s') -> (%s);\n" % (mode_index+1, mode_index+2, mode[2])

        s += "\n\t-- No state change without time consumption.\n"

        terms = ["(%s' = %s)" % (x.name, x.name) for x in self.crn.real_species]
        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        for c in self.crn.choice_variables:
            terms = ["(%s_%s' = %s_%s)" % (c.name, i, c.name, i) for i in list(range(c.minValue, c.maxValue + 1))]
            s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = ["(%s' = %s)" % (x.name, x.name) for x in self.crn.getRateConstants()]
        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = ["(%s' = %s)" % (x.variable_name, x.variable_name) for x in self.crn.optionalReactions]
        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = []
        for lambda_choice in self.crn.lambda_variables:
            terms.extend(["(%s' = %s)" % (x, x) for x in lambda_choice.lambdas])
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

        return s

    def constructdReal(self):
        s = "\n\nTRANS \n\n"

        s += "\t-- time constraint\n"
        s += "\ttime' = time + delta_time;\n\n"

        s += "\t-- must progress through modes in order\n"
        s += "\t mode_1' -> mode_1;\n"
        for i in list(range(2, len(self.modes)+1)):
            s += "\tmode_%s' -> (mode_%s or mode_%s);\n" % (i, i, i-1)

        s += "\n\n\t-- invariant conditions during modes\n"
        for mode_index, mode in enumerate(self.modes):
            s += "\tmode_%s  -> (%s);\n" % (mode_index+1, mode[1])
            s += "\tmode_%s'  -> (%s);\n" % (mode_index+1, mode[1])

        s += "\n\t-- jump conditions between modes\n"
        for mode_index, mode in enumerate(self.modes[:-1]):
            if len(mode) > 2 and mode[2]: # post-condition on mode
                s += "\t(mode_%s and mode_%s') -> (%s);\n" % (mode_index+1, mode_index+2, mode[2])

        s += "\n\t-- No state change without time consumption.\n"

        terms = ["(%s' = %s)" % (x.name, x.name) for x in self.crn.real_species]
        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        for c in self.crn.choice_variables:
            terms = ["(%s_%s' = %s_%s)" % (c.name, i, c.name, i) for i in list(range(c.minValue, c.maxValue + 1))]
            s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = ["(%s' = %s)" % (x.name, x.name) for x in self.crn.getRateConstants()]
        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = ["(%s' = %s)" % (x.variable_name, x.variable_name) for x in self.crn.optionalReactions]
        s += "\t(delta_time = 0) -> (%s);\n" % " and ".join(terms)

        terms = []
        for lambda_choice in self.crn.lambda_variables:
            terms.extend(["(%s' = %s)" % (x, x) for x in lambda_choice.lambdas])
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

        return s

class Initial:
    def __init__(self, crn, numModes, other_constraints):
        self.crn = crn
        self.numModes = numModes
        self.other_constraints = other_constraints

    def constructiSAT(self):
        s = "\nINIT\n"

        s += "\ttime = 0;\n\n"

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
        s = "\ninit:\n\n"

        s += " @1 ( and (time = 0) \n"
        s += "\t #cost condition\n"
        s += "\t( (or ( (%s) <= MAX_COST) (NO_COST_LIMIT = 1)))\n\n" % self.crn.get_cost()


        if len(self.crn.lambda_variables) > 0:
            s += "\n\t#Integer encoding of lambda variables\n"
        for lam in self.crn.lambda_variables:
            s += lam.dRealformat_constraint()

        if len(self.crn.choice_variables) > 0:
            s += "\n\t #Integer encoding of choice variables\n"
        for c in self.crn.choice_variables:
            s += c.dRealformat_constraint()
        if len(self.crn.choice_variables) > 0:
            s += "\n\t #Integer encoding of optional reaction variables\n"
        for optional_reaction in self.crn.optionalReactions:
            s += "\t(or (%s = 0) (%s = 1))\n" % (optional_reaction.variable_name, optional_reaction.variable_name)

        if len(self.crn.joint_choice_variables) > 0:
            s += "\n\t #Integer encoding of joint choice variables\n"
        for joint_choice_variables in self.crn.joint_choice_variables:
            s += joint_choice_variables.dRealformat_constraint()


        if len(self.crn.real_species) > 0:
            s += "\n\t #Limits on initial conditions\n"
        for sp in self.crn.real_species:
            s += sp.dRealInitialization()

        if len(self.crn.input_species) > 0:
            s += "\n\t #Initial conditions of inputs\n"
        for sp in self.crn.input_species:
            s += sp.dRealInitialization()

        if self.other_constraints:
            s += "\n\t $Manually specified constraints\n"
            s += "\t%s\n" % self.other_constraints

        s += ");"
        return s


class Flow:
    def __init__(self, var, t, fl, crn):
        self.variable = var
        self.time = t
        self.flow = fl
        self.crn = crn

    def constructiSAT(self):
        # Python represents powers as a**b, whereas iSAT uses a^b
        flow = str(self.flow).replace('**', '^')

        derivative_names = [x["name"] for x in self.crn.derivatives]

        if str(self.variable) in derivative_names:
            return "\t(%s = %s)" % (self.variable, flow)
        else:
            return "\t(d.%s/d.%s  = %s)" % (self.variable, self.time, flow)

    def constructdReal(self):
        raise NotImplementedError


class Mode:
    def __init__(self, m, inv):
        self.modeName = m
        self.invariants = inv

    def constructiSAT(self):

        s = "\n\n"

        for invariant in self.invariants:
            s += "\t-- transition into mode %s\n" % self.modeName
            if invariant[0] is not None:
                # constraint on mode start time
                s += "\t mode_%s -> (time >= %s);\n" % (self.modeName, invariant[0])
            # constraint imposed by mode on state
            s += "\t mode_%s -> (%s);\n" % (self.modeName, invariant[1])

        return s

    def constructdReal(self):

        s = "\n\n"

        for invariant in self.invariants:
            s += "\t #transition into mode %s\n" % self.modeName
            if invariant[0] is not None:
                s += "\t mode_%s -> (time >= %s);\n" % (self.modeName, invariant[0])
            s += "\t mode_%s -> (%s);\n" % (self.modeName, invariant[1])

        return s


class Post:
    def __init__(self, t, modes):
        self.time = t
        self.modes = modes

    def constructiSAT(self):

        if len(self.modes) > 0 and len(self.modes[-1]) > 2 and self.modes[-1][2]:
            post_condition = "and (" + self.modes[-1][2] + ")"
        else:
            post_condition = ""

        s = "\nTARGET \n"
        s += "\tmode_%s and (time < %s) %s;\n" % (len(self.modes), self.time, post_condition)
        return s

    def constructdReal(self):
        raise NotImplementedError

def constructISAT(crn, modes, flow, other_constraints=False):
    m_flow = [Flow(x, 'time', y, crn) for x, y in flow.items()]
    numModes = max(1, len(modes))

    d = Declaration(crn, numModes).constructiSAT()
    i = Initial(crn, numModes, other_constraints).constructiSAT()
    t = Transition(crn, m_flow, modes).constructiSAT()
    p = Post(1, modes).constructiSAT()  # TODO: set maxtime

    return d + i + t + p


def constructdReal(crn, modes, flow, other_constraints=False):
    m_flow = [Flow(x, 'time', y, crn) for x, y in flow.items()]
    numModes = max(1, len(modes))

    d = Declaration(crn, numModes).constructdReal()
    i = Initial(crn, numModes, other_constraints).constructdReal()
    t = Transition(crn, m_flow, modes).constructdReal()
    p = Post(1, modes).constructdReal()
    print d + i + t + p
    quit()
    return d + i + t + p
