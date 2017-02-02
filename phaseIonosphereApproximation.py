'''
Martijn Oei, 2017
This program calculates the quality of a linear fit to the ionospheric plasma opacity phase distortion effect,
which is inversely proportional to the frequency.
'''

from matplotlib import pyplot
import numpy

showGraphPhases = True
showGraphPhaseDifferences = True

deltaTEC = 0.1 # in TECU
frequencyStart = 300 # in MHz (corresponding to the uGMRT)
frequencyEnd = 500 # in MHz (corresponding to the uGMRT)
frequencyNumber = 201

frequencies = numpy.linspace(frequencyStart, frequencyEnd, num = frequencyNumber, endpoint = True)
phaseIonosphere = 8060 * deltaTEC * 60 / frequencies

phaseIonosphereLinearFit = numpy.polyval(numpy.polyfit(frequencies, phaseIonosphere, 1), frequencies)

phaseIonosphereDifference = phaseIonosphere - phaseIonosphereLinearFit


# Plot the inversely proportional ionospheric plasma opacity effect for uGMRT frequencies and a linear fit.
pyplot.figure(figsize = (10, 6))
pyplot.plot(frequencies, phaseIonosphere, label = r"actual effect $\propto\ \nu^{-1}$", color = "blue")
pyplot.plot(frequencies, phaseIonosphereLinearFit, label = r"best fit $\propto\ \nu$", color = "red", ls = "--")
pyplot.xlabel(r"radio frequency $\nu$ (MHz)")
pyplot.ylabel("phase ($\degree$)")
pyplot.grid()
pyplot.title("comparison between ionospheric plasma opacity phase distortion effect and linear fit\n$\Delta$TEC = " + str(deltaTEC))
pyplot.legend()
pyplot.subplots_adjust(left = 0.1, right = 0.98)
pyplot.savefig("./phasesIonosphereLinearFit1.png")

if (showGraphPhases):
    pyplot.show()
else:
    pyplot.close()

# Plot the difference between the two quantitative models of the effect.
pyplot.figure(figsize = (10, 6))
pyplot.fill_between(frequencies, 0, phaseIonosphereDifference, color = "FireBrick", alpha = 0.5)
pyplot.xlabel(r"radio frequency $\nu$ (MHz)")
pyplot.ylabel("phase difference ($\degree$)")
pyplot.grid()
pyplot.title("comparison between ionospheric plasma opacity phase distortion effect and linear fit\nmean absolute phase difference = "
             + str(numpy.round(numpy.mean(numpy.abs(phaseIonosphereDifference)), 1)) + "$\degree$")
pyplot.subplots_adjust(left = 0.1, right = 0.98)
pyplot.savefig("./phasesIonosphereLinearFit2.png")

if (showGraphPhaseDifferences):    
    pyplot.show()
else:
    pyplot.close()