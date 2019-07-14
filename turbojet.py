import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt

from IPython import embed

fullFile = "group4data/group4_full.TXT"
midFile = "group4data/group4_mid.TXT"

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
	xp = np.array([2,10,18,33,47])*4.44822
	yp = np.array([0,4,9,18,26])*4.44822
	return np.interp(thrustReading,xp,yp)

def plotCalibration():
	plt.figure()
	plt.plot([2,10,18,33,47],[0,4,9,18,26])
	plt.xlabel("Reading from load cell (lbf)")
	plt.ylabel("Actual Thrust (lbf)")
	plt.title("Thrust Measurement Calibration")
	plt.grid()
	plt.show()

#plotCalibration()

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
	
	data['Thrust'] = calibrateThrust(data['Thrust']*4.44822) 
	units['Thrust'] = 'N'
	
	data['FuelFlow'] /= 951019.388 #m^3/s
	units['FuelFlow'] = 'm^3/s'

	for tempKey in ['TurbExitTemp','CompExitTemp','CompInlet','TurbInletTemp', 'EGT']:
		data[tempKey] += 273.15 #convert to K
		units[tempKey] = 'K'
	for pressKey in ['CompExitPress','TurbInletPress','TurbExitPress','CompInletPress','NossExitPress']:
		data[pressKey] *= 6.89476
		units[pressKey] = 'kPa'


	return data,units#[data,headers,units]

fullData, units = dataProcess(fullFile)
midData = dataProcess(midFile)[0]

def plotTestGraphs():
	for testType,data in {'Full Throttle':fullData, 'Medium Throttle':midData}.items():
		fig, ax = plt.subplots(3,1)
		fig.suptitle(testType)
		for idx,yAxis in enumerate(['RPM','Thrust','CompExitPress']):
			ax[idx].plot(data['Time'],data[yAxis])
			ax[idx].set_ylabel("{} [{}]".format(yAxis,units[yAxis]))
			ax[idx].set_xlabel('Time [s]')
			ax[idx].grid()
	plt.show()

#plotTestGraphs()

#medium: 20, full: 25ish
midDataAvg = {key:np.mean(val) for key,val in midData.items()}
fullDataAvg = {key:np.mean(val) for key,val in fullData.items()}
#embed()
#['Time', 'Comp Inlet Press', 'Comp Exit Press', 'Turb Inlet Press', 'Turb Exit Press', 
#'Noss Exit Press', 'Fuel Flow', 'RPM', 'Thrust', 'Comp Inlet', 'Comp Exit Temp', 'Turb Inlet Temp', 'Turb Exit Temp', 'EGT']

#Analysis of Experimental Data

#----------------------------------------------------------------------
for testType,data in {'Full Throttle':fullDataAvg, 'Medium Throttle':midDataAvg}.items():
	#Step 1:
	#------------------------------------------------------------
	print "\n{}\n{}".format(testType,'-'*40)
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
	print 'Combustion Chamber pressure ratio, pi_b: {:.3}'.format(combustionPressRatio)

	#Step 6:
	#------------------------------------------------------------
	fuelMassFlowFraction = (Cp*(data['TurbInletTemp']-data['TurbExitTemp'])/(combustionHeatingValue-Cp*(data['TurbInletTemp']-data['TurbExitTemp'])))
	fuelHeatTransfer = fuelMassFlowFraction*combustionHeatingValue
	print 'Fuel Mass Flow Fraction: {:.3}'.format(fuelMassFlowFraction)
	print 'Fuel Heat Transfer: {:.3} J'.format(fuelHeatTransfer)

	#Step 8:
	#------------------------------------------------------------
	fuelMassFlowRate = (SG_kerosene*1e3)*data['FuelFlow']
	airMassFlowRate = fuelMassFlowRate/fuelMassFlowFraction
	print 'Measured Fuel Mass Flow Rate: {:.3} kg/s'.format(fuelMassFlowRate)
	print 'Estimated Air Mass Flow Rate: {:.3} kg/s'.format(airMassFlowRate)

	#Step 9:
	#------------------------------------------------------------
	turbPressRatio = data['TurbExitPress']/data['TurbInletPress']
	turbIsentropicEff = (airMassFlowRate*Cp*(data['TurbInletTemp']-data['TurbExitTemp']))/(Cp*data['TurbInletTemp']*(1 - turbPressRatio**((gamma-1)/gamma)))

	print "Turbine Pressure Ratio: {:.3}".format(turbPressRatio)
	print "Turbine Isentropic Efficiency: {:.3}".format(turbIsentropicEff)

	#Step 11:
	#------------------------------------------------------------
	#assume 100% efficient nozzle
	embed()
	nozzleFlowVelocity = np.sqrt(2*Cp*data['EGT']*(1-(data['NossExitPress']/data['TurbExitPress'])**((gamma-1)/gamma)))
	speedOfSound = np.sqrt(gamma*R*data['CompInlet'])
	machNum = nozzleFlowVelocity/speedOfSound

	print "Nozzle Flow Exit Velocty: {:.3} m/s".format(nozzleFlowVelocity)
	print "Mach Number: {:.3}".format(machNum)


embed()