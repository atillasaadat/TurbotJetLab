import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt

from IPython import embed

fullFile = "group4data/group4_full.txt"
midFile = "group4data/group4_mid.txt"

#['Time', 'Comp Inlet Press', 'Comp Exit Press', 'Turb Inlet Press', 'Turb Exit Press', 
#'Noss Exit Press', 'Fuel Flow', 'RPM', 'Thrust', 'Comp Inlet', 'Comp Exit Temp', 'Turb Inlet Temp', 'Turb Exit Temp', 'EGT']

#[nan, u'PSIG', u'PSIG', u'PSIG', u'PSIG', u'PSIG', u'Gal/Hour',
#u'RPM', u'Lbs.', u'\xb0C', u'\xb0C', u'\xb0C', u'\xb0C', u'\xb0C'],

R = 287.0 #J/(kg*K)
gamma = 1.4 #ratio of specific heats
Cp = 1005 #J/(kg*K)
SG_kerosene = 0.8 #specific gravity of kerosene
combustionHeatingValue = 43.1e6 #MJ/kg

def calibrateThrust(thrustReading):
	xp = [2,10,18,33,47]
	yp = [0,4,9,18,26]
	return np.interp(thrustReading,xp,yp)

def dataProcess(file):
	data = pd.read_csv(file, sep="\t", header=None, encoding="ISO-8859-1")
	del data[1]
	del data[max(data.columns)]
	data = data.to_numpy()
	headers = data[0]
	units = data[1]
	units[0] = 'sec'
	data = data[2:]
	startTime = datetime.strptime(data[0][0],"%H:%M:%S")
	data = np.vstack(([(datetime.strptime(i,"%H:%M:%S") - startTime).total_seconds() for i in data.T[0]],data.T[1:])).astype(float)
	data = {str(header.replace(' ','')):data[idx] for idx,header in enumerate(headers)}
	units = {str(header.replace(' ','')):units[idx].encode('utf-8') for idx,header in enumerate(headers)}
	data['Thrust'] = calibrateThrust(data['Thrust']) 
	return data,units#[data,headers,units]

fullData, units = dataProcess(fullFile)
midData = dataProcess(midFile)[0]

for testType,data in {'Full Throttle':fullData, 'Medium Throttle':midData}.items():
	fig, ax = plt.subplots(3,1)
	fig.suptitle(testType)
	for idx,yAxis in enumerate(['RPM','Thrust','CompExitPress']):
		ax[idx].plot(data['Time'],data[yAxis])
		ax[idx].set_ylabel("{} [{}]".format(yAxis,units[yAxis]))
		ax[idx].set_xlabel('Time [s]')
		ax[idx].grid()
plt.show()
#medium: 20, full: 25ish
midTimeInstantIdx = np.where(midData['Time']==20)[0][0]
fullTimeInstantIdx = np.where(midData['Time']==24)[0][0]
#embed()
#['Time', 'Comp Inlet Press', 'Comp Exit Press', 'Turb Inlet Press', 'Turb Exit Press', 
#'Noss Exit Press', 'Fuel Flow', 'RPM', 'Thrust', 'Comp Inlet', 'Comp Exit Temp', 'Turb Inlet Temp', 'Turb Exit Temp', 'EGT']

#Analysis of Experimental Data
#----------------------------------------------------------------------
for testType,data in {'Full Throttle':fullData, 'Medium Throttle':midData}.items():
	#Step 1:
	#------------------------------------------------------------
	print "{}\n{}".format(testType,'-'*40)
	if testType == 'Full Throttle':
		tiIdx = fullTimeInstantIdx
	else:
		tiIdx = midTimeInstantIdx
	for key in data.keys():
		data[key] = data[key][tiIdx]
	compSpecificWork = Cp*(data['CompExitTemp']-data['CompInlet'])
	compPressRatio = data['CompExitPress']/data['CompInletPress']
	compTempRatio = data['CompExitTemp']/data['CompInlet']
	compIdealSpecificWork = Cp*data['CompInlet']*((compPressRatio**((gamma-1)/gamma))-1)
	compEffAdibatic = compIdealSpecificWork/compSpecificWork
	print 'Comp. Specific Work: {:.3} J/kg'.format(compSpecificWork)
	print 'Stagnation Temp. Ratio: {:.3}'.format(compTempRatio)
	print 'Stagnation Press. Ratio: {:.3}'.format(compPressRatio)
	print 'Adibatic Efficiency: {:.3}'.format(compEffAdibatic)
	
	#Step 4:
	#------------------------------------------------------------
	combustionPressRatio = data['TurbInletPress']/data['CompExitPress']
	fuelMassFlowFraction = (Cp*(data['TurbInletTemp']-data['TurbExitTemp'])/(combustionHeatingValue-Cp*(data['TurbInletTemp']-data['TurbExitTemp'])))
	fuelHeatTransfer = fuelMassFlowFraction*combustionHeatingValue

	print 'Combustion Chamber pressure ratio, pi_b: {:.3}'.format(combustionPressRatio)
	print 'Fuel Mass Flow Fraction: {:3f}'.format(fuelMassFlowFraction)
	print 'Fuel Heat Transfer: {:3f} J'.format(fuelHeatTransfer)

embed()