import subprocess as sub
import re
import io
import os
from os import walk


#The output from the solver are saved in text files in the folder /bellshaperesults. Using constructResultsSummary() these can be parsed and the results summarised.
# Using getRunTimeAndResult(file, name) we can also retrieve the candidate CRN solution from a specific file

#A procedure that automates calls to the solver for each bellshape model where each model increases the discrete statespace as described in the paper. 
def stateSpaceExperiment():
	files =[]
	cwd = os.getcwd()
	for (dirpath, dirnames, filenames) in walk(cwd + "/bellshapemodels"):
		files.extend(filenames)
#Each of these different models within the folder /bellshapemodels is tried at 3 precisions 10^-1 10^-3 and 10^-5		

	for f in files:
		callSolver(stateSpaceConstructForCommandLine(f, 0.1, "--continue-after-not-reaching-horizon")) #the option --continue-after-not-reaching-horizon is used as an optimization for the solver.
		callSolver(stateSpaceConstructForCommandLine(f, 0.001, "--continue-after-not-reaching-horizon"))
		callSolver(stateSpaceConstructForCommandLine(f, 0.00001, "--continue-after-not-reaching-horizon"))

#A procedure that calls the solver for costs starting at a maximum cost and decreasing to a minimum cost. In order to generate the results found within the paper
# We start at a cost of 35 down to 10.
def optimalsynthesisExperiment(maximumcost, minimumcost):
	cwd = os.getcwd()
	cost = maximumcost
	while(cost > 0 ):
		editCost("bellshape.hys", cost)
		callSolver(costConstructForCommandLine("bellshape.hys", 0.1, cost, "--continue-after-not-reaching-horizon"))
		cost -= 1

#Once the solver has generated all of the relevant text files and outputed them to the folder /bellshaperesultsexp1 or /bellshaperesultsexp2  we can 
# construct a summary of the runtimes and results.
def constructResultSummary(experiment):
	files =[]
	results = []
	cwd = os.getcwd()
	for (dirpath, dirnames, filenames) in walk(cwd + "/bellshaperesults" + experiment):
		files.extend(filenames)

	for f in files:
		file = open("/bellshapemodelsresults/" + f).read().split('\n')
		results.append(getRunTimeAndResult(file, f))
	return results


#We can also generate individual runtimes and results.
def getRunTimeAndResult(file, name):
	vals = getCRNValues(file)
	#We check if the bellshape has been synthesised correctly. Sometimes even if the result is UNSAT, due to the option "--continue-after-not-reaching-horizon" 
	#passed to the ODE solver, it will output a candidate solution.
	return  + "  " + file[-1] + "SAT" if vals['K'] > 0.3 else return file[-1] + "UNSAT"


#We call the solver on the command line and wait (using .communicate()) for its output.
def callSolver(executableString):
	p = sub.Popen(executableString.split(),stdout=sub.PIPE,stderr=sub.PIPE)
	output, errors = p.communicate()
	return output
	 

#We construct a string which will be executed on the command line for the statespaceexperiment. We also specify where the output will be piped to.
def stateSpaceConstructForCommandLine(model, precision, otherPrams):
	return "./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt --i + bellshapemodels/" + model + " --prabs=" + precision + " --msw=" + precision*5 + " " + otherPrams + " > /bellshaperesultsexp1/" + model + precision + ".txt"

#We construct a string which will be executed on the command line for the optimal synthesis experiment. We also specify where the output will be piped to.
def costConstructForCommandLine(model, precision, cost, otherPrams):
	return "./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt --i + bellshape.hys --prabs=" + precision + " --msw=" + precision*5 + " " + otherPrams + " > /bellshaperesultsexp2/" + model + cost + precision + ".txt"
	

def constructforCommandLine(model,precision,otherPrams):
	if(model != ""):
		if(precision != ""):
			return "./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt --i bellshape" + model + ".hys + --prabs=" + precision + " --msw=" + precision*5 + " " + otherPrams  
		else:
			return  "./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt --i bellshape" + model + ".hys + --prabs=0.01 --msw=0.1 " + otherPrams
	else:
		if(precision != ""):
			return "./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt --i bellshape.hys + --prabs=" + precision + " --msw=" + precision*5 + " " + otherPrams  
		else:
			return  "./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt --i bellshape.hys + --prabs=0.01 --msw=0.1 " + otherPrams


#This function edits the cost within the bellshape.hys file such that we can automate the synthesis experiment
#If the cost is 0 we comment out references to cost in the .hys file. 
def editCost(filename, i):
	f = open("bellshape.hys")
	file = f.read().split('\n')
		for linenumber,line in enumerate(s, 1):
			if "--cost" in line:
				if cost == 0:
					file[linenumber] = '--6*(c1c + c5c + c3c + 2*(1-c3c)) + 5*(c2 + c6 + c4 ) <= 0; --cost'
				else
					file[linenumber] =	'6*(c1c + c5c + c3c + 2*(1-c3c)) + 5*(c2 + c6 + c4 ) <= ' + cost +'; --cost'    
				if "--costfinal" in line:
				if cost == 0:
					file[linenumber] =  'mode_2 and (time <= MAX_TIME) and (K >= 0) and (K < 0.1); --and (6*(c1c + c5c + c3c + 2*(1-c3c)) + 5*(c2 + c6 + c4 ) <= 0); --costfinal'
				else
					file[linenumber] = 	'mode_2 and (time <= MAX_TIME) and (K >= 0) and (K < 0.1) and (6*(c1c + c5c + c3c + 2*(1-c3c)) + 5*(c2 + c6 + c4 ) <= ' + cost +'); --costfinal'
	f.write('\n'.join(file)).close()

#This function is used to construct CRN from the output of an individual file. 
def constructResults(file):
	x = constructCRN(getCRNValues(file))
	for y in x:
		print y
	print getRunTimeAndResult(file)



#This function builds a dictionary (vals) of synthesized parameters which are outputed by the solver. It scrapes the output file for the relevant parameters.
def getCRNValues(s):
	vals = {}

	storedLineNo = 0
	flag = 0

	for linenumber,line in enumerate(s, 1):

		if flag > 0:
			if flag == 1:
				vals['lambda1'] = 'A' if (re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[') == '0') else 'B'
			if flag == 2:
				vals['lambda2'] = 'A' if (re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[') == '0') else 'B'
			if flag == 3:
				vals['c1'] = '1' if (re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[') == '0') else '0'
			if flag == 4:
				vals['c2'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
			if flag == 5:
				vals['c3'] = '2' if (re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[') == '0') else '1'
			if flag == 6:
				vals['c4'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
			if flag == 7:
				vals['c5'] = '1' if (re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[') == '0') else '0'
			if flag == 8:
				vals['c6'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
			if flag == 9:
				vals['c7'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
			if flag == 10:
				vals['choice1'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[') 
			#We take the average of the parameter interval given for the rates.
			if flag == 11:
				rate1 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')) 
				rate2 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']'))
				vals['k1'] = str((rate1+rate2)/2)
			if flag == 12:
				rate1 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')) 
				rate2 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']'))
				vals['k2'] = str((rate1+rate2)/2)
			if flag == 13:
				rate1 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')) 
				rate2 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']')) 
				vals['k3'] = str((rate1+rate2)/2)
			flag = 0;

		if 'lam1 (float):' in line:
			storedLineNo = linenumber
			flag = 1
		if 'lam2 (float):' in line:
			storedLineNo = linenumber
			flag = 2
		if 'c1c (float):' in line:
			storedLineNo = linenumber
			flag = 3
		if 'c2 (float):' in line:
			storedLineNo = linenumber
			flag = 4
		if 'c3c (float):' in line:
			storedLineNo = linenumber
			flag = 5
		if 'c4 (float):' in line:
			storedLineNo = linenumber
			flag = 6
		if 'c5c (float):' in line:
			storedLineNo = linenumber
			flag = 7
		if 'c6 (float):' in line:
			storedLineNo = linenumber
			flag = 8
		if 'c7 (float):' in line:
			storedLineNo = linenumber
			flag = 9
		if 'choice1 (float):' in line:
			storedLineNo = linenumber
			flag = 10
		if 'k1 (float):' in line:
			storedLineNo = linenumber
			flag = 11
		if 'k2 (float):' in line:
			storedLineNo = linenumber
			flag = 12
		if 'k3 (float):' in line:
			storedLineNo = linenumber
			flag = 13
	return vals

#We can construct the CRN based on the parameter values derived from the previous function.
def constructCRN(x):
	reaction1 = x['lambda1'] 
	reaction1 += x['c1'] + "*K" if x['c1'] is not '0' else ""
	reaction1 += " ->{" + x['k1'] + "} "
	reaction1 += x['c2'] + "*K" if x['c2'] is not '0' else ""

	reaction2 = x['c5'] + "*" + x['lambda1'] if x['c5'] is not '0' else ""
	reaction2 += " + " + x['c3'] + "*K" if x['c3'] is not '0' else ""
	reaction2 += "->{" + x['k2'] + "} " 
	reaction2 += x['c6'] + "*" + x['lambda2'] if x['c6'] is not '0' else ""
	reaction2 += " + " + x['c4'] + "*K" if x['c4'] is not '0' else ""


	reaction3 = "->{" + x['k3'] + "} "
	reaction3 += x['c7'] + "*K" if x['choice1'] is '0' and x['c7'] is not '0'  else ""
	reaction3 += x['lambda2'] if x['choice1'] is '1' else ""
	reaction3 = "" if len(reaction3) < 20 else ""

	return [reaction1, reaction2, reaction3]




stateSpaceExperiment()





