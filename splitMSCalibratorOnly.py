'''
Martijn  Oei, 2017
This program takes a CASA MS, and creates a new MS containing only a calibrator field.
'''

from casat import split
split = split.split


nameFieldCalibrator = "3C286"
pathMSOriginal = "/data1/MartijnOei/martijn_400mugs_pilot/DDTB247-GWB_FULLPOL/DDTB247-GWB_FULLPOL.MS"
pathMSCalibrator = "/data1/MartijnOei/martijn_400mugs_pilot/DDTB247-GWB_FULLPOL/DDTB247-GWB_FULLPOL_" + nameFieldCalibrator + "_ONLY.MS"
dataColumn = "data" # Common choices are "data" and "corrected".

split(vis = pathMSOriginal, outputvis = pathMSCalibrator, field = nameFieldCalibrator, datacolumn = dataColumn, keepflags = True)
print ("Saved calibrator-only MS to: " + pathMSCalibrator)