import subprocess as sub
import re
import io
import os
from os import walk


# This function is used to construct CRN from the output of an individual file.
def constructResults(file):
    x = constructCRN(getCRNValues(file))
    for y in x:
        print y
    print getRunTimeAndResult(file)


# This function builds a dictionary (vals) of synthesized parameters which are outputed by the solver. It scrapes the output file for the relevant parameters.
def getCRNValues(s):
    vals = {}

    storedLineNo = 0
    flag = 0

    for linenumber, line in enumerate(s, 1):

        if flag > 0:
            if flag == 1:
                vals['lambda1'] = 'A' if (
                re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[') == '0') else 'B'
            if flag == 2:
                vals['lambda2'] = 'A' if (
                re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[') == '0') else 'B'
            if flag == 3:
                vals['c1'] = '1' if (
                re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[') == '0') else '0'
            if flag == 4:
                vals['c2'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
            if flag == 5:
                vals['c3'] = '2' if (
                re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[') == '0') else '1'
            if flag == 6:
                vals['c4'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
            if flag == 7:
                vals['c5'] = '1' if (
                re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[') == '0') else '0'
            if flag == 8:
                vals['c6'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
            if flag == 9:
                vals['c7'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
            if flag == 10:
                vals['choice1'] = re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('[')
            # We take the average of the parameter interval given for the rates.
            if flag == 11:
                rate1 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('['))
                rate2 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']'))
                vals['k1'] = str((rate1 + rate2) / 2)
            if flag == 12:
                rate1 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('['))
                rate2 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']'))
                vals['k2'] = str((rate1 + rate2) / 2)
            if flag == 13:
                rate1 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[0].strip('['))
                rate2 = float(re.search('\[.*\]', s[storedLineNo + 1]).group(0).split(',')[1].strip(']'))
                vals['k3'] = str((rate1 + rate2) / 2)
            flag = 0

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


# We can construct the CRN based on the parameter values derived from the previous function.
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
