#def precisionExperiment():


def precisionExperiment():
	files =[]
	cwd = os.getcwd()
	f = "phosphorelaywithcovariances.hys"
	f2 = "phosphorelaywithoutcovariances.hys"
		
	callSolver(autoConstructForCommandLine(f, 0.1, '')) 
	callSolver(autoConstructForCommandLine(f, 0.001, ''))
	callSolver(autoConstructForCommandLine(f, 0.00001, ''))

	callSolver(autoConstructForCommandLine(f2, 0.1, '')) 
	callSolver(autoConstructForCommandLine(f2, 0.001, ''))
	callSolver(autoConstructForCommandLine(f2, 0.00001, ''))



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



def getCRNValues(s):
	vals = {}

	storedLineNo = 0
	flag = 0

	for linenumber,line in enumerate(s, 1):

		if flag > 0:
			if flag == 1:
				rate1 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')) 
				rate2 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']'))
				vals['k1'] = str((rate1+rate2)/2)
			if flag == 2:
				rate1 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')) 
				rate2 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']'))
				vals['k2'] = str((rate1+rate2)/2)
			if flag == 3:
				rate1 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')) 
				rate2 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']')) 
				vals['k3'] = str((rate1+rate2)/2)
			if flag == 4:
				rate1 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')) 
				rate2 =  float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']')) 
				vals['k3'] = str((rate1+rate2)/2)
			flag = 0;

		if 'k1 (float):' in line:
			storedLineNo = linenumber
			flag = 1
		if 'k2 (float):' in line:
			storedLineNo = linenumber
			flag = 2
		if 'k3 (float):' in line:
			storedLineNo = linenumber
			flag = 3
		if 'k4 (float):' in line:
			storedLineNo = linenumber
			flag = 4
	return vals

def constructCRN(x):

	reaction1 = "L1 + B ->{" + x['k1'] + "}  B + L1p"
	reaction2 = "L2 + L1p ->{" + x['k2'] + "} L1 + L2p"
	reaction3 = "L2p + L3 ->{" + x['k3'] + "}  L2 + L3p"
	reaction4 = "L3p ->{" + x['k4'] + "} L3"

	return [reaction1, reaction2, reaction3, reaction4]
