
def precisionExperiment():
	files =[]
	cwd = os.getcwd()
	callSolver(autoConstructForCommandLine("superpoisson.hys", 0.001, '')) 
	callSolver(autoConstructForCommandLine("superpoisson10.hys", 0.001, ''))
	callSolver(autoConstructForCommandLine("superpoisson100.hys", 0.001, ''))

def constructResultSummary():
	files =[]
	results = []
	cwd = os.getcwd()
	for (dirpath, dirnames, filenames) in walk(cwd + "/phosphorelayresults"):
		files.extend(filenames)

	for f in files:
		file = open("/phosphorelayresults/" + f).read().split('\n')
		results.append(getRunTimeAndResult(file, f))
	return results



def getRunTimeAndResult(file, name):
	vals = getCRNValues(file)
	return  + "  " + file[-1]


def callSolver(executableString):
	p = sub.Popen(executableString.split(),stdout=sub.PIPE,stderr=sub.PIPE)
	output, errors = p.communicate()
	#print subprocess.check_output(['ls','-l'])
	return output
	 

def autoConstructForCommandLine(model, precision, otherPrams):
	return "./isat-ode-r2806-static-x86_64-generic-noSSE-stripped.txt --i" + model + " --prabs=" + precision + " --msw=" + precision*5 + " " + otherPrams + " > /phosphorelayresults/" + model + precision ".txt"


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


def constructforSpecificFile(file):
	return constructCRN(getCRNValues(file))

def getRunTimeAndResult(file, name):
	vals = getCRNValues(file)
	return  + "  " + file[-1] + "SAT"


def constructforSpecificFile(file):
	return constructCRN(getCRNValues(file))

#def intervalExperiment():
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
				vals['c1'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
			if flag == 4:
				vals['c2'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
			if flag == 5:
				vals['c3'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
			if flag == 6:
				rate1 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')) 
				rate2 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']'))
				vals['k1'] = str((rate1+rate2)/2)
			if flag == 7:
				rate1 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')) 
				rate2 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']'))
				vals['k2'] = str((rate1+rate2)/2)
			if flag == 8:
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
		if 'c1 (float):' in line:
			storedLineNo = linenumber
			flag = 3
		if 'c2 (float):' in line:
			storedLineNo = linenumber
			flag = 4
		if 'c3 (float):' in line:
			storedLineNo = linenumber
			flag = 5
		if 'k1 (float):' in line:
			storedLineNo = linenumber
			flag = 6
		if 'k2 (float):' in line:
			storedLineNo = linenumber
			flag = 7
	return vals

#$\tau_1: \to^{k_1} c_1 A + c_2 \lambda_1; \quad 
#\tau_2: A \to^{k_2} c_3 \lambda_2; $
def constructCRN(x):
	reaction1 = " ->{" + x['k1'] + "} " x['lambda1'] 
	reaction1 += x['c1'] + "*A" if x['c1'] is not '0' else ""
	reaction1 += " + " x['c2'] + "*" + x['lambda1'] if x['c2'] is not '0' else ""

	reaction2 = "A ->{" + x['k2'] + "} " 
	reaction2 +=  x['c3'] + "*" + x['lambda2'] if x['c3'] is not '0' else ""

	return [reaction1, reaction2,]
