'''
Martijn Oei, 2016
This program converts (uGMRT) UVFITS data into a CASA MS, using CASANOVA.
'''

from casat import importgmrt
importgmrt = importgmrt.importgmrt

# Set the paths where the UVFITS data can be found, and where the MS data should be stored.
dataPathUVFITS = "/data1/MartijnOei/martijn_400mugs_pilot/DDTB247-GWB_FULLPOL/"
dataPathMS = "/data1/MartijnOei/martijn_400mugs_pilot/DDTB247-GWB_FULLPOL/"

# The original UVFITS file as sent to us by the GMRT staff ('DDTB247-GWB_FULLPOL.UVFITS') could not be converted to a MS using 'importgmrt' -
# the procedure would crash at the end. Huib then imported this UVFITS file in AIPS, and exported it again as 'DDTB247-GWB_FULLPOL.UVFITS2'.
# Using 'importgmrt' then worked without problems.
# Set the names of the UVFITS data and the output MS data.
dataNameUVFITS = "DDTB247-GWB_FULLPOL.UVFITS2" 
dataNameMS = "DDTB247-GWB_FULLPOL.MS"

# Perform the conversion.
importgmrt(fitsfile = dataPathUVFITS + dataNameUVFITS, vis = dataPathMS + dataNameMS)