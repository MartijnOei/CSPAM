#
# This file is largely a copy of lines 672-1653 of EVLA_functions.py
# from the EVLA pipeline version 1.3.4 retrieved from
# https://science.nrao.edu/facilities/vla/data-processing/pipeline/scripted-pipeline
# on 16 February 2016 by Kasper van Dam (M.Sc. student Leiden Observatory).
#


# This is needed to run inside CASANOVA
# Kasper van Dam: note that line 239 is changed
# Martijn Oei: use of the string method .lower() in line 246 now depends on whether the EVLA is the telescope used. Initialization argument and public variable 'telescope' added.
# Martijn Oei: line 288 and 289 look superfluous; commented
import casac
casac = casac.casac

from casat import flagdata
flagdata = flagdata.flagdata

def logprint(text, logfileout=None):
	print text


# Description:
# ------------
# This file contains the reference antenna heuristics.

# The present heuristics are geometry and flagging.

# Classes:
# --------
# RefAntHeuristics - This class chooses the reference antenna heuristics.
# RefAntGeometry   - This class contains the geometry heuristics for the
#                    reference antenna.
# RefAntFlagging   - This class contains the flagging heuristics for the
#                    reference antenna.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.
# 2012 Jun 06 - Nick Elias, NRAO
#               Modified to exclude ALMA pipeline classe dependencies.

# ------------------------------------------------------------------------------

# Imports
# -------

import numpy

# ------------------------------------------------------------------------------
# class RefAntHeuristics
# ------------------------------------------------------------------------------

# RefAntHeuristics
# ----------------

# Description:
# ------------
# This class chooses the reference antenna heuristics.

# Public member variables:
# ------------------------
# vis      - This python string contains the MS name.
#
# field    - This python string or list of strings contains the field numbers
#            or IDs.  Presently it is used only for the flagging heuristic.
# spw      - This python string or list of strings contains the spectral
#            window numbers of IDs.  Presently it is used only for the
#            flagging heuristic.
# intent   - This python string or list of strings contains the intent(s).
#            Presently it is used only for the flagging heuristic.
#
# geometry - This python boolean determines whether the geometry heuristic will
#            be used.
# flagging - This python boolean determines whether the flagging heuristic will
#            be used.

# Public member functions:
# ------------------------
# __init__  - This public member function constructs an instance of the
#             RefAntHeuristics() class.
# calculate - This public member function forms the reference antenna list
#             calculated from the selected heuristics.

# Private member functions:
# -------------------------
# _get_names - This private member function gets the antenna names from the MS.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version created with public member variables vis, field,
#               spw, intent, geometry, and flagging; public member functions
#               __init__() and calculate(); and private member function
#               _get_names().
# 2012 Jun 06 - Nick Elias, NRAO
#               api inheritance eliminated.

# ------------------------------------------------------------------------------

class RefAntHeuristics:

# ------------------------------------------------------------------------------

# RefAntHeuristics::__init__

# Description:
# ------------
# This public member function constructs an instance of the RefAntHeuristics()
# class.

# NB: If all of the defaults are chosen, no reference antenna list is returned.

# Inputs:
# -------
# vis        - This python string contains the MS name.
#
# field      - This python string or list of strings contains the field numbers
#              or IDs.  Presently it is used only for the flagging heuristic.
#              The default is ''.
# spw        - This python string or list of strings contains the spectral
#              window numbers of IDs.  Presently it is used only for the
#              flagging heuristic.  The default is ''.
# intent     - This python string or list of strings contains the intent(s).
#              Presently it is used only for the flagging heuristic.  The
#              default is ''.
#
# geometry   - This python boolean determines whether the geometry heuristic
#              will be used in automatic mode.  The default is False.
# flagging   - This python boolean determines whether the flagging heuristic
#              will be used in automatic mode.  The default is False.

# Outputs:
# --------
# None, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.
# 2012 Jun 06 - Nick Eluas, NRAO
#               Input parameter defaults added.

# ------------------------------------------------------------------------------

	def __init__(self, vis, field = '', telescope = "GMRT", spw = '', intent = '', geometry = False,
	    flagging = False):

		# Initialize the public member variables of this class

		self.vis = vis

		self.field = field
		self.telescope = telescope
		self.spw = spw
		self.intent = intent

		self.geometry = geometry
		self.flagging = flagging


		# Returns nothing.
		return None

# ------------------------------------------------------------------------------

# RefAntHeuristics::calculate

# Description:
# ------------
# This public member function forms the reference antenna list calculated from
# the selected heuristics.

# NB: A total score is calculated from all heuristics.  The best antennas have
# the highest scores, so a reverse sort is performed to obtain the final list.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The numpy array of strings containing the ranked reference antenna list,
# returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def calculate( self ):

		# If no heuristics are specified, return no reference antennas

		if not ( self.geometry or self.flagging ): return []


		# Loads (capitalized versions of) the antenna names and initializes the score dictionary.
		names = self._get_names()

		score = dict()
		for n in names: score[n] = 0.0


		# For each selected heuristic, add the score for each antenna

		self.geoScore = 0.0
		self.flagScore = 0.0

		if self.geometry:
			geoClass = RefAntGeometry( self.vis )
			self.geoScore = geoClass.calc_score()
			for n in names: score[n] += self.geoScore[n]

		if self.flagging:
			flagClass = RefAntFlagging( self.vis, self.field,
			    self.spw, self.intent )
			self.flagScore = flagClass.calc_score()
			for n in names:
                            try:
                                score[n] += self.flagScore[n]
                            except KeyError, e:
                                logprint ("WARNING: antenna "+str(e)+", is completely flagged and missing from calibrators.ms", logfileout='logs/refantwarnings.log')


		# Calculate the final score and return the list of ranked
		# reference antennas.  NB: The best antennas have the highest
		# score, so a reverse sort is required.

		keys = numpy.array( score.keys() )
		values = numpy.array( score.values() )
		argSort = numpy.argsort( values )[::-1]
		
		# Martijn: because 'keys[argSort]' = [['C13'],['C13'],...,['C13']], Kasper added '[0]' below.
		refAntUpper = keys[argSort][0] # 'refAntUpper' is now a list with a single element.
		
		# Martijn: Based on the telescope, antenna names are either written uppercase or lowercase.
		# Martijn: EVLA names are lowercase (e.g. 'ea02'); GMRT names are uppercase.
		# Martijn: 'refAntUpper' contains uppercase letters.
		refAnt = list()
		if (self.telescope == "EVLA"):
			for r in refAntUpper: refAnt.append(r.lower())
		else: # It is assumed that all other telescopes (e.g. the GMRT) have uppercase antenna names.
			for r in refAntUpper: refAnt.append(r)

		# Return the list of ranked reference antennas. Martijn: this list has length 1.
		return(refAnt)

# ------------------------------------------------------------------------------

# RefAntHeuristics::_get_names

# Description:
# ------------
# This private member function gets the antenna names from the MS.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The numpy array of strings containing the antenna names, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _get_names( self ):
		# Create the local instance of the table tool and open the MS
		tbLoc = casac.table()
		tbLoc.open( self.vis+'/ANTENNA' )


		# Gets the antenna names and automatically capitalizes them.
		# (Unfortunately, some CASA tools capitalize them and others do not.)
		names = tbLoc.getcol( 'NAME' ).tolist()

		# Martijn: the next two lines do not seem to make sense. I commented them.
		#rNames = range( len(names) )
		#for n in rNames: names[n] = names[n]
		
		# Close the local instance of the table tool and delete it.
		tbLoc.close()
		del tbLoc

		# Return the antenna names
		return names

# ------------------------------------------------------------------------------
# class RefAntGeometry
# ------------------------------------------------------------------------------

# RefAntGeometry
# --------------

# Description:
# ------------
# This class contains the geometry heuristics for the reference antenna.

# Algorithm:
# ----------
# * Calculate the antenna distances from the array center.
# * Normalize the distances by the maximum distance.
# * Calculate the score for each antenna, which is one minus the normalized
#   distance.  The best antennas have the highest score.
# * Sort according to score.

# Public member variables:
# ------------------------
# vis - This python string contains the MS name.

# Public member functions:
# ------------------------
# __init__   - This public member function constructs an instance of the
#              RefAntGeometry() class.
# calc_score - This public member function calculates the geometry score for
#              each antenna.

# Private member functions:
# -------------------------
# _get_info       - This private member function gets the information from the
#                   antenna table of the MS.
# _get_measures   - This private member function gets the measures from the
#                   antenna table of the MS.
# _get_latlongrad - This private member function gets the latitude, longitude
#                   and radius (from the center of the earth) for each antenna.
# _calc_distance  - This private member function calculates the antenna
#                   distances from the array reference from the radii,
#                   longitudes, and latitudes.
# _calc_score     - This private member function calculates the geometry score
#                   for each antenna.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

class RefAntGeometry:

# ------------------------------------------------------------------------------

# RefAntGeometry::__init__

# Description:
# ------------
# This public member function constructs an instance of the RefAntGeometry()
# class.

# Inputs:
# -------
# vis - This python string contains the MS name.

# Outputs:
# --------
# None, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def __init__( self, vis ):

		# Set the public variables

		self.vis = vis


		# Return None

		return None

# ------------------------------------------------------------------------------

# RefAntGeometry::calc_score

# Description:
# ------------
# This public member function calculates the geometry score for each antenna.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def calc_score( self ):

		# Get the antenna information, measures, and locations

		info = self._get_info()
		measures = self._get_measures( info )
		radii, longs, lats = self._get_latlongrad( info, measures )


		# Calculate the antenna distances and scores

		distance = self._calc_distance( radii, longs, lats )
		score = self._calc_score( distance )


		# Return the scores

		return score

# ------------------------------------------------------------------------------

# RefAntGeometry::_get_info

# Description:
# ------------
# This private member function gets the information from the antenna table of
# the MS.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The python dictionary containing the antenna information, returned via the
# function value.  The dictionary format is:
# 'position'          - This numpy array contains the antenna positions.
# 'flag_row'          - This numpy array of booleans contains the flag row
#                       booleans.  NB: This element is of limited use now and
#                       may be eliminated.
# 'name'              - This numpy array of strings contains the antenna names.
# 'position_keywords' - This python dictionary contains the antenna information.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _get_info( self ):

		# Create the local instance of the table tool and open it with
		# the antenna subtable of the MS

		tbLoc = casac.table()
		tbLoc.open( self.vis+'/ANTENNA' )


		# Get the antenna information from the antenna table

		info = dict()

		info['position'] = tbLoc.getcol( 'POSITION' )
		info['flag_row'] = tbLoc.getcol( 'FLAG_ROW' )
		info['name'] = tbLoc.getcol( 'NAME' )
		info['position_keywords'] = tbLoc.getcolkeywords( 'POSITION' )


		# Close the table tool and delete the local instance

		tbLoc.close()
		del tbLoc


		# The flag tool appears to return antenna names as upper case,
		# which seems to be different from the antenna names stored in
		# MSes.  Therefore, these names will be capitalized here.

		rRow = range( len( info['name'] ) )


		# Return the antenna information

		return info

# ------------------------------------------------------------------------------

# RefAntGeometry::_get_measures

# Description:
# ------------
# This private member function gets the measures from the antenna table of the
# MS.

# Inputs:
# -------
# info - This python dictionary contains the antenna information from private
#        member function _get_info().

# Outputs:
# --------
# The python dictionary containing the antenna measures, returned via the
# function value.  The dictionary format is:
# '<antenna name>' - The python dictionary containing the antenna measures.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _get_measures( self, info ):

		# Create the local instances of the measures and quanta tools

		meLoc = casac.measures()
		qaLoc = casac.quanta()


		# Initialize the measures dictionary and the position and
		# position_keywords variables

		measures = dict()

		position = info['position']
		position_keywords = info['position_keywords']

		rf = position_keywords['MEASINFO']['Ref']

		for row,ant in enumerate( info['name'] ):

			if not info['flag_row'][row]:

				p = position[0,row]
				pk = position_keywords['QuantumUnits'][0]
				v0 = qaLoc.quantity( p, pk )

				p = position[1,row]
				pk = position_keywords['QuantumUnits'][1]
				v1 = qaLoc.quantity( p, pk )

				p = position[2,row]
				pk = position_keywords['QuantumUnits'][2]
				v2 = qaLoc.quantity( p, pk )

				measures[ant] = meLoc.position( rf=rf, v0=v0,
				    v1=v1, v2=v2 )


		# Delete the local instances of the measures and quanta tools

		del qaLoc
		del meLoc


		# Return the measures

		return measures

# ------------------------------------------------------------------------------

# RefAntGeometry::_get_latlongrad

# Description:
# ------------
# This private member function gets the latitude, longitude and radius (from the
# center of the earth) for each antenna.

# Inputs:
# -------
# info     - This python dictionary contains the antenna information from
#            private member function _get_info().
# measures - This python dictionary contains the antenna measures from private
#            member function _get_measures().

# Outputs:
# --------
# The python tuple containing containing radius, longitude, and latitude python
# dictionaries, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _get_latlongrad( self, info, measures ):

		# Create the local instance of the quanta tool

		qaLoc = casac.quanta()


		# Get the radii, longitudes, and latitudes

		radii = dict()
		longs = dict()
		lats = dict()

		for ant in info['name']:

			value = measures[ant]['m2']['value']
			unit = measures[ant]['m2']['unit']
			quantity = qaLoc.quantity( value, unit )
			convert = qaLoc.convert( quantity, 'm' )
			radii[ant] = qaLoc.getvalue( convert )

			value = measures[ant]['m0']['value']
			unit = measures[ant]['m0']['unit']
			quantity = qaLoc.quantity( value, unit )
			convert = qaLoc.convert( quantity, 'rad' )
			longs[ant] = qaLoc.getvalue( convert )

			value = measures[ant]['m1']['value']
			unit = measures[ant]['m1']['unit']
			quantity = qaLoc.quantity( value, unit )
			convert = qaLoc.convert( quantity, 'rad' )
			lats[ant] = qaLoc.getvalue( convert )


		# Delete the local instance of the quanta tool

		del qaLoc


		# Return the tuple containing the radius, longitude, and
		# latitude python dictionaries

		return radii, longs, lats

# ------------------------------------------------------------------------------

# RefAntGeometry::_calc_distance

# Description:
# ------------
# This private member function calculates the antenna distances from the array
# reference from the radii, longitudes, and latitudes.

# NB: The array reference is the median location.

# Inputs:
# -------
# radii - This python dictionary contains the radius (from the center of the
#         earth) for each antenna.
# longs - This python dictionary contains the longitude for each antenna.
# lats  - This python dictionary contains the latitude for each antenna.

# Outputs:
# --------
# The python dictionary containing the antenna distances from the array
# reference, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _calc_distance( self, radii, longs, lats ):

		# Convert the dictionaries to numpy float arrays.  The median
		# longitude is subtracted.

		radiusValues = numpy.array( radii.values() )

		longValues = numpy.array( longs.values() )
		longValues -= numpy.median( longValues )

		latValues = numpy.array( lats.values() )


		# Calculate the x and y antenna locations.  The medians are
		# subtracted.

		x = longValues * numpy.cos(latValues) * radiusValues
		x -= numpy.median( x )

		y = latValues * radiusValues
		y -= numpy.median( y )


		# Calculate the antenna distances from the array reference and
		# return them

		distance = dict()
		names = radii.keys()

		for i,ant in enumerate(names):
			distance[ant] = numpy.sqrt( pow(x[i],2) + pow(y[i],2) )

		return distance

# ------------------------------------------------------------------------------

# RefAntGeometry::_calc_score

# Description:
# ------------
# This private member function calculates the geometry score for each antenna.

# Algorithm:
# ----------
# * Calculate the antenna distances from the array center.
# * Normalize the distances by the maximum distance.
# * Calculate the score for each antenna, which is one minus the normalized
#   distance.  The best antennas have the highest score.
# * Sort according to score.

# Inputs:
# -------
# distance - This python dictionary contains the antenna distances from the
#            array reference.  They are calculated in private member function
#            _calc_distance().

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _calc_score( self, distance ):

		# Get the number of good data, calculate the fraction of good
		# data, and calculate the good and bad weights

		far = numpy.array( distance.values(), numpy.float )
		fFar = far / float( numpy.max(far) )

		wFar = fFar * len(far)
		wClose = ( 1.0 - fFar ) * len(far)


		# Calculate the score for each antenna and return them

		score = dict()

		names = distance.keys()
		rName = range( len(wClose) )

		for n in rName: score[names[n]] = wClose[n]

		return score

# ------------------------------------------------------------------------------

# RefAntFlagging
# --------------

# Description:
# ------------
# This class contains the flagging heuristics for the reference antenna.

# Algorithm:
# ----------
# * Get the number of unflagged (good) data for each antenna.
# * Normalize the good data by the maximum good data.
# * Calculate the score for each antenna, which is one minus the normalized
#   number of good data.  The best antennas have the highest score.
# * Sort according to score.

# Public member variables:
# ------------------------
# vis    - This python string contains the MS name.
#
# field  - This python string or list of strings contains the field numbers or
#          or IDs.
# spw    - This python string or list of strings contains the spectral window
#          numbers of IDs.
# intent - This python string or list of strings contains the intent(s).

# Public member functions:
# ------------------------
# __init__   - This public member function constructs an instance of the
#              RefAntFlagging() class.
# calc_score - This public member function calculates the flagging score for
#              each antenna.

# Private member functions:
# -------------------------
# _get_good   - This private member function gets the number of unflagged (good)
#               data from the MS.
# _calc_score - This private member function calculates the flagging score for
#               each antenna.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

class RefAntFlagging:

# ------------------------------------------------------------------------------

# RefAntFlagging::__init__

# Description:
# ------------
# This public member function constructs an instance of the RefAntFlagging()
# class.

# Inputs:
# -------
# vis    - This python string contains the MS name.
#
# field  - This python string or list of strings contains the field numbers or
#          or IDs.
# spw    - This python string or list of strings contains the spectral window
#          numbers of IDs.
# intent - This python string or list of strings contains the intent(s).

# Outputs:
# --------
# None, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def __init__( self, vis, field, spw, intent ):

		# Set the public member functions

		self.vis = vis

		self.field = field
		self.spw = spw
		self.intent = intent


		# Return None

		return None

# ------------------------------------------------------------------------------

# RefAntFlagging::calc_score

# Description:
# ------------
# This public member function calculates the flagging score for each antenna.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def calc_score( self ):

		# Calculate the number of unflagged (good) measurements for each
		# antenna, determine the score, and return them

		good = self._get_good()
		score = self._calc_score( good )

		return( score )

# ------------------------------------------------------------------------------

# RefAntFlagging::_get_good

# Description:
# ------------
# This private member function gets the number of unflagged (good) data from the
# MS.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The dictionary containing the number of unflagged (good) data from the MS,
# returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _get_good( self ):

                #Update April 2015 to use the flagging task
                task_args = {'vis'          : self.vis, 
                           'mode'         : 'summary',
                           'field'        : self.field,
                           'spw'          : self.spw,
                           'intent'       : self.intent,
                           'display'      : '',
                           'flagbackup'   : False,
                           'savepars'     : False}

                d = flagdata(**task_args)                                               



		# Calculate the number of good data for each antenna and return
		# them

		antenna = d['antenna']
		good = dict()

		for a in antenna.keys():
			good[a] = antenna[a]['total'] - antenna[a]['flagged']

		return( good )

# ------------------------------------------------------------------------------

# RefAntFlagging::_calc_score

# Description:
# ------------
# This private member function calculates the flagging score for each antenna.

# Algorithm:
# ----------
# * Get the number of unflagged (good) data for each antenna.
# * Normalize the good data by the maximum good data.
# * Calculate the score for each antenna, which is one minus the normalized
#   number of good data.  The best antennas have the highest score.
# * Sort according to score.

# Inputs:
# -------
# good - This python dictionary contains the number of unflagged (good) data
#        from the MS.  They are obtained in private member function _get_good().

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _calc_score( self, good ):

		# Get the number of good data, calculate the fraction of good
		# data, and calculate the good and bad weights

		nGood = numpy.array( good.values(), numpy.float )
		fGood = nGood / float( numpy.max(nGood) )

		wGood = fGood * len(nGood)
		wBad = ( 1.0 - fGood ) * len(nGood)


		# Calculate the score for each antenna and return them

		score = dict()

		names = good.keys()
		rName = range( len(wGood) )

		for n in rName: score[names[n]] = wGood[n]

		return score

######################################################################
