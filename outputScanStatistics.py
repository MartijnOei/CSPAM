'''
Martijn Oei, 2017
This program loads a (uGMRT logging) text file with information about the observed scans.
It is assumed this file is created beforehand, manually, by taking the uGMRT-provided log file
and removing everything other than the table with scan information.
It produces plots and scan statistics.
'''

import numpy
from matplotlib import pyplot
import datetime

# Set the path and filename to search for.
# In this case, data of the November 26, 2016 observations is loaded - 'ddtb247ScansOnly.log'.
logFilePath = "/data1/MartijnOei/martijn_400mugs_pilot/"
logFileName = "ddtb247ScansOnly.log"

# Set the name of the observing run, to be used in the plots.
observingRunName = "GMRT DDTB247"

# Set visual parameters for the bar plots.
barWidth = 0.9
barAlpha = 0.5
targetNamesRotation = 30 # in degrees



# Load the information.
scanIDs, targetNames, timeStartStrings, timeEndStrings = numpy.loadtxt(logFilePath + logFileName, dtype = { "names": ("scan ID", "target name", "time start", "time end"),
                                                                                                            "formats": ("i2", "a7", "a8", "a8")}, usecols = (1,2,3,5), unpack = True)
scanNumber = len(scanIDs)
targetNumber = len(set(targetNames))

# Generate a list of scan durations from time strings.
scanDurations = [] # in seconds
for i in range(scanNumber):
	timeStart = datetime.datetime.strptime(timeStartStrings[i], "%H:%M:%S")
	timeEnd = datetime.datetime.strptime(timeEndStrings[i], "%H:%M:%S")
	scanDurations.append((timeEnd - timeStart).total_seconds())

# Generate 3 lists, to be used in histograms / bar plots. These contain:
histTargetNames = [] # target names (no duplicates)
histScanNumbers = [] # the number of corresp. scans
histObsTimeCum = [] # the corresp. total observation length

for i in range(scanNumber):
	if (not targetNames[i] in histTargetNames):
		histTargetNames.append(targetNames[i])
		histScanNumbers.append(1)
		histObsTimeCum.append(scanDurations[i])
	else:
		index = histTargetNames.index(targetNames[i])
		histScanNumbers[index] += 1
		histObsTimeCum[index] += scanDurations[i]


# Print the histogram / bar plot statistics.
print (histTargetNames)
print (histScanNumbers)
print (histObsTimeCum)
print ("average duration of scan (s):", numpy.divide(numpy.array(histObsTimeCum), numpy.array(histScanNumbers)))


# Draw histogram bar plot with the total number of scans.
pyplot.bar(numpy.arange(targetNumber), histScanNumbers, barWidth, alpha = barAlpha)
pyplot.ylim(0, numpy.amax(histScanNumbers) + 1)
pyplot.xticks(numpy.arange(targetNumber) + barWidth / 2, histTargetNames)
labels = pyplot.gca().get_xticklabels()
pyplot.setp(labels, rotation = targetNamesRotation)
pyplot.ylabel("number of scans (1)")
pyplot.title(observingRunName + " distribution of scans over targets\ntotal number of scans: " + str(scanNumber))
pyplot.gca().yaxis.grid(True)
pyplot.show()

# Draw bar plot showing the total amount of time each target was observed.
pyplot.bar(numpy.arange(targetNumber), histObsTimeCum, barWidth, alpha = barAlpha)
pyplot.ylim(0, numpy.amax(histObsTimeCum) + 60)
pyplot.xticks(numpy.arange(targetNumber) + barWidth / 2, histTargetNames)
labels = pyplot.gca().get_xticklabels()
pyplot.setp(labels, rotation = targetNamesRotation)
pyplot.ylabel("total observing time (s)")
pyplot.title(observingRunName + " distribution of observing time over targets\ntotal observing time (min): " + str(numpy.round(numpy.sum(histObsTimeCum) / 60, 2)))
pyplot.gca().yaxis.grid(True)
pyplot.show()