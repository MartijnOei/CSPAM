'''
Martijn Oei, 2017
This program tests the rewritten plotting procedures of the CSPAM STObj Class.
Based on the CASA Calibration Table (CT) / Solution Table (ST) chosen, appropriate plotting procedures are performed.
Use this convenience utility for tests, as it circumvents running the CSPAM pipeline step 'bandpass', which might take a while to complete.
'''

from lib import TableObjects

# Initialize the directory where plots should be stored. This directory could be different from the one where the CSPAM pipeline stores plots.
plotDirectory = "/data1/MartijnOei/martijn_400mugs_pilot/DDTB247-GWB_FULLPOL/DDTB247-GWB_FULLPOL-plot/flux_cal_3C286/"

# Initialize the location of the CT that should be used.
calibrationTablePath = "/data1/MartijnOei/martijn_400mugs_pilot/DDTB247-GWB_FULLPOL/DDTB247-GWB_FULLPOL-cal/flux_cal_3C286/"
calibrationTableStep = "final"
calibrationTableType = 'K'
calibrationTableFile = calibrationTableStep + '.' + calibrationTableType + '/'

# Generate CT plots.
calibrationTable = TableObjects.STObj(calibrationTablePath + calibrationTableFile)
calibrationTable.plot(plotDirectory)