# Description:
# ------------
# This file contains two classes that deal with CASA tables:
# MSObj (Measurement Set Object) and STObj (Solution Table Object).
# 
# These classes originally come from CSPAM version 0.1 by Francesco
# de Gasperin and Huib Intema. Edited by Kasper van Dam, 2016. Edited by Martijn Oei, 2017.
#
# These objects rely heavily on the CASA environment interfaced using casanova.

import os
import sys
import logging
import numpy
from scipy import interpolate
from matplotlib import pyplot

# CSPAM Modules
import SourceObjects
import utils

# CASA Toolkits
import casac
ms = casac.casac.ms()
tb = casac.casac.table()

# CASA Tasks
from casat import plotcal
plotcal = plotcal.plotcal

class MSObj:
    """
    The MSObj class provides information on CASA Measurement Sets

    Available instance attributes:
    ------------------------------
    name                    type    info

    file_path             - str   - absolute path to the file
    ms_name               - str   - name of the measurement set (without .ms)
    summary               - dict  - summary of the measurement set
    scansummary           - dict  - summary of the main table
    dir_img               - str   - absolute path to the img directory
    dir_cal               - str   - absolute path to the cal directory
    dir_plot              - str   - absolute path to the plot directory
    dir_peel              - str   - absolute path to the peel directory
    minBL_for_cal         - int   - minimum number of baselines needed
    nchan                 - int   - number of channels
    freq                  - float - reference frequency
    telescope             - str   - name of the telescope
    band                  - str   - name of the band
    fluxcalibrator        - user  - user defined class with flux cal info
    phasecalibrator       - user  - user defined class with phase cal info
    leakcalibrator        - user  - user defined class with leakeage cal info
    poscalibrator         - user  - user defined class with position cal info
    anglecalibrator       - user  - user defined class with angle cal info
    targetsources         - list  - list of target sources (user defined class)
    *spw                   - str   - spectral window
    *central_chans         - str   - central channels
    *flag                  - dict  - flag from config in casa style

    * = not used at the moment

    Available public instance methods:
    ----------------------------------
    get_scan_ids_from_field_id(field, user_scan_ids = [])
    
    get_field_name_from_field_id(field):
    
    get_field_name_from_scan_id(scan):
    
    get_field_id_from_scan_id(scan):

    get_direction_from_tgt_field_id(tgt)

    These are pretty much self-explanatory. 
    """

    def __init__(self, path, conf = None):
        # Names and paths
        self.file_path = path
        try:
            occurence_last_slash = path.rfind('/')
            self.ms_name = path[occurence_last_slash+1:]
        except:
            self.ms_name = path
        self.ms_name = self.ms_name.replace('.ms','')
        self.ms_name = self.ms_name.replace('.MS','')
        
        # Measurement set summaries
        ms.open(self.file_path)
        self.summary = ms.summary()
        self.scansummary = ms.getscansummary()
        ms.close()
        
        # Directories where CSPAM data reduction products and human output are stored
        self.dir_img = self.file_path.replace('.ms','')
        self.dir_img = self.dir_img.replace('.MS','')
        self.dir_img = self.dir_img + "-img"
        
        self.dir_cal = self.file_path.replace('.ms','')
        self.dir_cal = self.dir_cal.replace('.MS','')
        self.dir_cal = self.dir_cal + "-cal"
        
        self.dir_plot = self.file_path.replace('.ms','')
        self.dir_plot = self.dir_plot.replace('.MS','')
        self.dir_plot = self.dir_plot + "-plot"

        self.dir_peel = self.file_path.replace('.ms','')
        self.dir_peel = self.dir_peel.replace('.MS','')
        self.dir_peel = self.dir_peel + "-peel"
        
        # Create these directories if they do not exist yet.
        if (not os.path.isdir(self.dir_img)):
            os.makedirs(self.dir_img)
        if (not os.path.isdir(self.dir_plot)):
            os.makedirs(self.dir_plot)
        if (not os.path.isdir(self.dir_cal)):
            os.makedirs(self.dir_cal)
        if (not os.path.isdir(self.dir_peel)):
            os.makedirs(self.dir_peel)
        
        # Minimum number of baselines
        self.minBL_for_cal = self._get_minBL_for_cal()
        
        # Number of channels
        self.nchan = self._getChannelNumber()
        
        # Frequency width of channels
        self.channelWidth = self._getChannelWidth()
        
        # Frequency
        self.freq = self._get_freq()
        
        # Telescope
        self.telescope = self._get_telescope()
        assert self.telescope == "GMRT" or self.telescope == "EVLA"
        
        # Band
        self.band = self._get_band()
        
        # Separate the scans for the calibrator and targets
        # Note that cal_scan_ids is an empty list if there are only targets
        cal_scan_ids = self._determine_cal_scan_ids()
        tgt_scan_ids = list(set(self.scansummary.keys())-set(cal_scan_ids))
        
        cal_field_id = self.get_field_id_from_scan_id(cal_scan_ids)
        
        print ("Calibrator fields found (IDs):", cal_field_id)
        print ("Calibrator scans found (IDs):", cal_scan_ids)
        

        # Start with the assumption that there is one calibrator for all types
        # Override these parameters later if different settings exist in the
        # configuration file.
        
        # Check if there are calibrators
        if len(cal_field_id) > 0:
            # Check if there are multiple calibrators
            if len(cal_field_id) > 1:
                # If there are more calibrators simply take the first and
                # remove the scans associated with the other calibrators
                flux_cal_field_id = cal_field_id[0]
                flux_cal_field_name = self.get_field_name_from_field_id(flux_cal_field_id)
                flux_cal_scan_ids = self.get_scan_ids_from_field_id(flux_cal_field_id, user_scan_ids=cal_scan_ids)
            else:
                # there is only one calibrator
                flux_cal_field_id = cal_field_id[0]
                flux_cal_field_name = self.get_field_name_from_field_id(flux_cal_field_id)
                flux_cal_scan_ids = cal_scan_ids

            # One calibrator assumption so equalize everything
            phase_cal_field_id = flux_cal_field_id
            phase_cal_field_name = self.get_field_name_from_field_id(
                                        phase_cal_field_id)
            phase_cal_scan_ids = flux_cal_scan_ids
            
            leak_cal_field_id = flux_cal_field_id
            leak_cal_field_name = self.get_field_name_from_field_id(
                                       leak_cal_field_id)
            leak_cal_scan_ids = flux_cal_scan_ids
        
            pos_cal_field_id = flux_cal_field_id
            pos_cal_field_name = self.get_field_name_from_field_id(pos_cal_field_id)
            pos_cal_scan_ids = flux_cal_scan_ids
        
            angle_cal_field_id = flux_cal_field_id
            angle_cal_field_name = self.get_field_name_from_field_id(
                                        angle_cal_field_id)
            angle_cal_scan_ids = flux_cal_scan_ids
        
        # Check if there is configuration file. If so, possibly override variable values that were just set.
        # If this mset only contains targets there is no configuration file.
        if conf is not None:
            
            # See if tgt_scans/cal_scans/cal_types are set and override			
            if conf['cal_scans']:
                # Remove the calibrator scans that are unwanted
                cal_scan_ids = list(set(cal_scan_ids).intersection(conf['cal_scans']))
            
            if conf['tgt_scans']:
                # Remove the target scans that are unwanted
                tgt_scan_ids = list(set(tgt_scan_ids).intersection(conf['tgt_scans']))
            
            
            if conf['flux_cal_field']: # This expression will be 'True' if 'conf['flux_cal_field']' is not an empty string.
                flux_cal_field_id = conf['flux_cal_field']
                flux_cal_field_name = self.get_field_name_from_field_id(flux_cal_field_id)
                flux_cal_scan_ids = self.get_scan_ids_from_field_id(flux_cal_field_id, user_scan_ids=cal_scan_ids)
            
            if conf['phase_cal_field']:
                phase_cal_field_id = conf['phase_cal_field']
                phase_cal_field_name = self.get_field_name_from_field_id(phase_cal_field_id)
                phase_cal_scan_ids = self.get_scan_ids_from_field_id(phase_cal_field_id, user_scan_ids=cal_scan_ids)
            
            if conf['leakage_cal_field']:
                leak_cal_field_id = conf['leakage_cal_field']
                leak_cal_field_name = self.get_field_name_from_field_id(leak_cal_field_id)
                leak_cal_scan_ids = self.get_scan_ids_from_field_id(leak_cal_field_id, user_scan_ids=cal_scan_ids)
            
            if conf['position_cal_field']:
                pos_cal_field_id = conf['position_cal_field']
                pos_cal_field_name = self.get_field_name_from_field_id(pos_cal_field_id)
                pos_cal_scan_ids = self.get_scan_ids_from_field_id(pos_cal_field_id, user_scan_ids=cal_scan_ids)
            
            if conf['angle_cal_field']:
                angle_cal_field_id = conf['angle_cal_field']
                angle_cal_field_name = self.get_field_name_from_field_id(angle_cal_field_id)
                angle_cal_scan_ids = self.get_scan_ids_from_field_id(angle_cal_field_id, user_scan_ids=cal_scan_ids)
			
			## The following variables (flag, spw, central_chans) aren't used yet
            # Flags (manual flag to be used with flagdata command)
            self.flag = self._convert_flag(conf['flag']) 
		
            # Spectral windows
            self.spw = conf['spw']
            
            # Central channel
            central_nchan = self.nchan*conf['central_chan_percentage']/100.
            self.central_chans = str(int(round(self.nchan/2.-central_nchan/2.))) \
                                 + '~' + \
                                 str(int(round(self.nchan/2.+central_nchan/2.)))

        # Check if there are calibrators
        if len(cal_field_id) > 0:
            # Calibrators
            self.fluxcalibrator = SourceObjects.CalibratorSource('flux', 
                                  flux_cal_field_name, flux_cal_field_id,
                                  flux_cal_scan_ids)
            self.phasecalibrator = SourceObjects.CalibratorSource('phase',
                                   phase_cal_field_name, phase_cal_field_id,
                                   phase_cal_scan_ids)
            self.leakcalibrator = SourceObjects.CalibratorSource('leakage',
                                  leak_cal_field_name, leak_cal_field_id, 
                                  leak_cal_scan_ids)
            self.poscalibrator = SourceObjects.CalibratorSource('position',
                                 pos_cal_field_name, pos_cal_field_id, 
                                 pos_cal_scan_ids)
            self.anglecalibrator = SourceObjects.CalibratorSource('angle',
                                   angle_cal_field_name, angle_cal_field_id, 
                                   angle_cal_scan_ids)

        # Targets (list of an arbitrary number of target sources)
        self.targetsources = []
        tgt_field_ids = self.get_field_id_from_scan_id(tgt_scan_ids)
        for tgt_field_id in tgt_field_ids:
            tgt_field_name = self.get_field_name_from_field_id(tgt_field_id)
            tgt_scan_ids = self.get_scan_ids_from_field_id(tgt_field_id)
            target = SourceObjects.TargetSource(tgt_field_name, tgt_field_id, 
                                                tgt_scan_ids)
            self.targetsources.append(target)
        
        # List of split off target measurement sets (which are also instances
        # of this class). Note that for a target mset this list is empty.
        self.targetmsets = []
    
    
    # Private Methods
    # (Actually, private methods do not exist in Python. Using the _ is a convention we adhere to rather than a syntactic necessity.)

    def _getChannelNumber(self):
        '''
        Returns: the number of channels (int).
        '''
        tb.open(self.file_path + "/SPECTRAL_WINDOW")
        channelNumber = tb.getcol("NUM_CHAN")[0]
        tb.close()
        return channelNumber
    
    def _getChannelWidth(self):
    		'''
    		Returns: the width of the channels in Hertz (float). It is assumed that all channels have the same width.
    		'''
    		tb.open(self.file_path + "/SPECTRAL_WINDOW")
    		channelWidth = tb.getcol("CHAN_WIDTH")[0]
    		tb.close()
    		return channelWidth
    
    def _get_telescope(self):
        """
        Return: the telescope name
        """
        tb.open(self.file_path+'/OBSERVATION')
        telescope = tb.getcol('TELESCOPE_NAME')[0]
        tb.close()
        return telescope

    def _get_band(self):
        """
        Return telescope band
        Note that if you add more bands later on, you also need to add more
        primary beam attenuations in skymodel.py
        """
        if self.telescope == 'GMRT':
            if self.freq > 650e6: return '1420'
            if self.freq > 550e6 and self.freq < 650e6: return '610'
            if self.freq > 250e6 and self.freq < 350e6: return '325'
            if self.freq > 200e6 and self.freq < 300e6: return '235'
            if self.freq < 200e6: return '151'
        elif self.telescope == 'EVLA':
            if self.freq < 1e9: return 'P'
            if self.freq >= 1e9: return 'L'

    def _get_freq(self):
        """
        Return: the reference frequency
        """
        tb.open(self.file_path+'/SPECTRAL_WINDOW')
        freq = tb.getcol('REF_FREQUENCY')[0]
        tb.close()
        return freq

    def _get_antenna_names(self):
        """
        Return: list of antenna names
        """
        tb.open( '%s/ANTENNA' % self.file_path)
        antenna_names = tb.getcol( 'NAME' )
        tb.close()
        return antenna_names

    def _get_minBL_for_cal(self):
        """
        Return: estimate the minimum BL for calibration steps
        """
        num_antenna = len(self._get_antenna_names())
        return max(3,int(num_antenna/4.0))

    def _convert_flag(self, flag_string=''):
        """
        Convert flag command from the config file to a casa command
        e.g.: C14,E03,E04,S01,W01=; =22:30:00~22:43:00; C03=22:52:30~22:55:30
        """
        flag = {}
        if flag_string != '':
            for flag_group in flag_string.replace(' ','').split(';'):
                ant = flag_group.split('=')[0]
                time = flag_group.split('=')[1]
                flag[ant] = time
            
        return flag

    def _determine_cal_scan_ids(self):
        """
        Save the calibrator scans in a list
        If empy list given, then check for coords
        """
        known_cals = {
        '3C147':{'m0': {'unit': 'rad', 'value': 1.49488177653836},
                'm1': {'unit': 'rad', 'value': 0.87008056907685105},
                'refer': 'J2000',
                'type': 'direction'},
        '3C196':{'m0': {'unit': 'rad', 'value': 2.1537362969610028},
                'm1': {'unit': 'rad', 'value': 0.841554132080366},
                'refer': 'J2000',
                'type': 'direction'},
        '3C286':{'m0': {'unit': 'rad', 'value': 3.5392586514514845},
                'm1': {'unit': 'rad', 'value': 0.53248541037303654},
                'refer': 'J2000',
                'type': 'direction'},
        '3C295':{'m0': {'unit': 'rad', 'value': 3.7146787856873482},
                'm1': {'unit': 'rad', 'value': 0.91111035090915105},
                'refer': 'J2000',
                'type': 'direction'},
        '3C380':{'m0': {'unit': 'rad', 'value': 4.8412379124131713},
                'm1': {'unit': 'rad', 'value': 0.85078013643188044},
                'refer': 'J2000',
                'type': 'direction'},
        '3C48':{'m0': {'unit': 'rad', 'value': 0.42624576436309852},
                'm1': {'unit': 'rad', 'value': 0.57874633182450852},
                'refer': 'J2000',
                'type': 'direction'}
        }
                
        ms.open(self.file_path)
        cal_scan_ids = []
        for cal_scan_id in self.scansummary.keys():
            cal_field_id = self.get_field_id_from_scan_id(cal_scan_id)
            direc = ms.getfielddirmeas(fieldid=int(cal_field_id))
            # if distance with known cal is < than 60" then add it
            for known_cal, known_cal_dir in known_cals.iteritems():
                if utils.angularSeparationOfDirectionsArcsec(
                         direc, known_cal_dir) <= 60:
                    logging.info('Found '+known_cal+' in scan: '
                                 +cal_scan_id)
                    cal_scan_ids.append(cal_scan_id)

                    # update field name for SetJy
                    tb.open('%s/FIELD' % self.file_path, nomodify=False)
                    tb.putcell('NAME', int(cal_field_id), known_cal)
                    source_id = tb.getcell('SOURCE_ID', int(cal_field_id))
                    tb.close()
                    tb.open('%s/SOURCE' % self.file_path, nomodify=False)
                    tb.putcell('NAME', source_id, known_cal)
                    tb.close()

                    break # cal found, useless to keep on iterating
        ms.close()

        if cal_scan_ids == []:
            logging.info('No calibrators found, this ms only contains targets')

        return cal_scan_ids
    
    
    # Public Methods

    def get_scan_ids_from_field_id(self, field_id, user_scan_ids = []):
        """
        This method returns a list of scan ids for a given field_id while
        possibly also excluding scans that are unwanted.
        Note that user_scan_ids is a list of wanted scans.
        """
        scans = []
        for scan_id in self.scansummary.keys():
            test_field_id = self.get_field_id_from_scan_id(scan_id)
            if int(test_field_id) == int(field_id):
                scans.append(scan_id)
        if len(user_scan_ids) > 0:
            scans = list(set(scans).intersection(user_scan_ids))
        return scans

    def get_field_name_from_field_id(self, field):
        field_name = self.summary['field_'+str(field)]['name']
        return field_name
 
    def get_field_name_from_scan_id(self, scan):
        """
        this method returns the field name of a scan id 
        or returns a list of field names of a list of scan ids
        """
        if isinstance(scan, list):
            fieldnames = []
            for i in scan:
                field_name = self.summary['scan_'+str(i)]['0']['FieldName']
                if field_name not in fieldnames:
                    fieldnames.append(field_name)
            return fieldnames
        else:
            field_name = self.summary['scan_'+str(scan)]['0']['FieldName']
            return field_name

    def get_field_id_from_scan_id(self, scan):
        """
        This method returns the field id of a scan id 
        or returns a list of field ids of a list of scan ids
        """
        if isinstance(scan, list):
            fieldids = []
            for i in scan:
                field_id = self.summary['scan_'+str(i)]['0']['FieldId']
                if field_id not in fieldids:
                    fieldids.append(field_id)
            return fieldids
        else:
            field_id = self.summary['scan_'+str(scan)]['0']['FieldId']
            return str(field_id)
        
    def get_direction_from_tgt_field_id(self, tgt):
        direction = self.summary['field_'+str(tgt)]['direction']
        return direction
        
    def get_all_available_timestamps(self):
        tb.open(self.file_path)
        timestamps = tb.getcol('TIME')
        tb.close()
        return timestamps



class STObj:
    '''
    A Solution Table Object (STObj) is used to provide information on (e.g. plot) Solution Tables (i.e. CASA Calibration Tables).
    '''
    def __init__(self, file_path):
        self.file_path = file_path
        self.st_type = self._get_type() # Possibilities are "K Jones", "G Jones" and "B Jones".
        self.calibratorName = self._getCalibratorName() # The name is added whenever a Calibration Table is created. See steps.py --> bandpass_calibration(). Used in plot().
        self.telescopeName = self._getTelescopeName() # E.g. "GMRT"
        self.MSName = self._getMSName() # E.g. "DDTB247-GWB_FULLPOL"
    
    
    # Private Methods
    # (actually private methods don't exist in Python, the _ is a convention)

    def _get_type(self):
        '''
        This method returns the associated Calibration Table type.
        '''
        tb.open(self.file_path)
        st_type = tb.getkeyword("VisCal")
        tb.close()
        return st_type
    
    def _getCalibratorName(self):
        '''
        This method returns the associated Calibration Table calibrator name.
        '''
        tb.open(self.file_path)
        calibratorName = tb.getkeyword("CalibratorName")
        tb.close()
        return calibratorName
    
    def _getTelescopeName(self):
        '''
        This method returns the associated Calibration Table telescope name.
        '''
        tb.open(self.file_path + "/OBSERVATION")
        telescopeName = tb.getcol("TELESCOPE_NAME")[0]
        tb.close()
        return telescopeName
    
    def _getMSName(self):
        '''
        This method returns the associated Calibration Table Measurement Set name.
        '''
        tb.open(self.file_path)
        MSName = tb.getkeyword("MSName")[ : -3] # We choose not to include the extension (".ms" or ".MS").
        tb.close()
        return MSName
    
    
    # Public Methods

    def plot(self, plotDirectory, phase_only = False, amp_only = False):
        '''
        This method plots some of the content of the associated Calibration Table.
        The content plotted depends on the Calibration Table type.
        '''
        
        if (not os.path.isdir(plotDirectory)):
            os.makedirs(plotDirectory)
        
        
        if (self.st_type == "K Jones"):
            # Plot delay
            #tb.open( '%s/ANTENNA' % self.file_path)
            #nameAntenna = tb.getcol( 'NAME' )
            #numAntenna = len(nameAntenna)
            #tb.close()
            
            tb.open(self.file_path)
            ants1 = tb.getcol('ANTENNA1') # Antenna
            nameAntenna = numpy.sort(numpy.unique(ants1))
            numAntenna = len(nameAntenna)
            tb.close()
            
            nplots=int(numAntenna/3)
            for ii in range(nplots):
                filename = plotDirectory + "/delay_" + str(ii) + '.png'
                systemCommand='rm -rf ' + filename
                os.system(systemCommand)
                
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                plotcal(caltable=self.file_path,xaxis='time',yaxis='delay',antenna=antPlot,subplot=311,\
                        overplot=False,clearpanel='All',iteration='antenna',plotrange=[],\
                        plotsymbol='o-',markersize=5.0,fontsize=10.0,showgui=False,\
                        figfile=filename)
        
        if (self.st_type == "G Jones"):
            """
            For G Jones tables plot both amp and phase, if phase_only is True, plot
            only phase.
            """
            #tb.open( '%s/ANTENNA' % self.file_path)
            #nameAntenna = tb.getcol( 'NAME' )
            #numAntenna = len(nameAntenna)
            #tb.close()
            
            tb.open(self.file_path)
            ants1 = tb.getcol('ANTENNA1') # Antenna
            nameAntenna = numpy.sort(numpy.unique(ants1))
            numAntenna = len(nameAntenna)
            tb.close()
            
            nplots=int(numAntenna/3)
            # Plot amp
            if not phase_only:
                tb.open(self.file_path)
                cpar=tb.getcol('CPARAM')
                flgs=tb.getcol('FLAG')
                tb.close()
                amps=numpy.abs(cpar)
                good=numpy.logical_not(flgs)
                plotmax=numpy.max(amps[good])
                BL = False # Not used in this pipeline
                for ii in range(nplots):
                    filename=plotDirectory+'/'+'a_'+str(ii)+'.png'
                    syscommand='rm -rf '+filename
                    os.system(syscommand)
                    antPlot=str(ii*3)+'~'+str(ii*3+2)
                    if BL: xaxis = 'antenna2'
                    else: xaxis = 'time'
                    if BL: plotsymbol = 'o'
                    else: plotsymbol = 'o-'
                    plotcal(caltable=self.file_path,xaxis=xaxis,yaxis='amp',antenna=antPlot,subplot=311,\
                            iteration='antenna',plotrange=[0,0,0,plotmax],plotsymbol=plotsymbol,plotcolor='red',\
                            markersize=5.0,fontsize=10.0,showgui=False,figfile=filename,clearpanel='All')
            # Plot phase
            if not amp_only:
                for ii in range(nplots):
                    filename=plotDirectory+'/'+'p_'+str(ii)+'.png'
                    syscommand='rm -rf '+filename
                    os.system(syscommand)
                    antPlot=str(ii*3)+'~'+str(ii*3+2)
                    BL = False # Not used in this pipeline
                    if BL: xaxis = 'antenna2'
                    else: xaxis = 'time'
                    plotcal(caltable=self.file_path,xaxis=xaxis,yaxis='phase',antenna=antPlot,subplot=311,\
                            iteration='antenna',plotrange=[0,0,-180,180],showflags=False,\
                            plotsymbol='o-',plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,\
                            figfile=filename)
        
        if (self.st_type == "B Jones"):
            '''
            Generate and store bandpass graphs (both amplitude and phase).
            The data is assumed to contain only 2 different polarizations!
            '''
            
            # Open the Calibration Table, extract gains and flags, and close again.
            tb.open(self.file_path)
            gainsComplex = tb.getcol("CPARAM")
            flags = tb.getcol("FLAG")
            tb.close()
            
            # Divide up the data in two different polarizations.
            gainsComplexPol1 = gainsComplex[0]
            gainsComplexPol2 = gainsComplex[1]
            flagsPol1 = flags[0]
            flagsPol2 = flags[1]
            
            # Transpose the data so that 'gainsComplexPoli[j]' is the gain-versus-channel data for polarization i and antenna j.
            # The flags should be transposed accordingly.
            gainsComplexPol1 = numpy.transpose(gainsComplexPol1)
            gainsComplexPol2 = numpy.transpose(gainsComplexPol2)
            flagsPol1 = numpy.transpose(flagsPol1)
            flagsPol2 = numpy.transpose(flagsPol2)
            
            # The first channel of gain-versus-channel data (for some polarization and antenna) is the channel corresponding to the highest frequency.
            # We want the first channel to correspond to the lowest frequency, and thus we mirror the data.
            gainsComplexPol1 = numpy.fliplr(gainsComplexPol1)
            gainsComplexPol2 = numpy.fliplr(gainsComplexPol2)
            flagsPol1 = numpy.fliplr(flagsPol1)
            flagsPol2 = numpy.fliplr(flagsPol2)
            
            # Calculate the gain amplitudes and gain phases.            
            gainAmplitudesPol1 = numpy.absolute(gainsComplexPol1)
            gainAmplitudesPol2 = numpy.absolute(gainsComplexPol2)
            gainPhasesPol1 = numpy.angle(gainsComplexPol1, deg = True)
            gainPhasesPol2 = numpy.angle(gainsComplexPol2, deg = True)
            
            # Create an array containing the central frequency of each channel.
            tb.open(self.file_path + "/SPECTRAL_WINDOW")
            frequencyBand = numpy.flipud(tb.getcol("CHAN_FREQ") / 1e6) # Extract frequencies from table, convert from Hz to MHz, and put in ascending order.
            tb.close()
            
            
            # Set the plot limits. ALERT: Amplitude limits should be made more general!
            plotFrequencyMin = numpy.floor(frequencyBand[0]) # MHz-rounding
            plotFrequencyMax = numpy.ceil(frequencyBand[-1]) # MHz-rounding
            plotGainAmplitudeMin = .5
            plotGainAmplitudeMax = 1.5
            
            # Creates a bandpass amplitude and phase plot for each antenna. Already-existing plots are overwritten.
            antennaNumber = gainsComplex.shape[2]
            for antennaID in range(antennaNumber):
                if (numpy.all(flagsPol1[antennaID]) and numpy.all(flagsPol2[antennaID])):
                    print ("Skipping graph creation for antenna ID " + str(antennaID) + ": all data is flagged.")
                else:
                    print("Starting graph creation for antenna ID " + str(antennaID) + "...")
                    pyplot.figure(figsize = (12, 6))
                    pyplot.scatter(frequencyBand, numpy.ma.masked_array(gainAmplitudesPol1[antennaID], mask = flagsPol1[antennaID]), c = 'b', lw = 0, label = "polarization 1")
                    pyplot.scatter(frequencyBand, numpy.ma.masked_array(gainAmplitudesPol2[antennaID], mask = flagsPol2[antennaID]), c = 'r', lw = 0, label = "polarization 2")
                    pyplot.grid()
                    pyplot.legend()
                    pyplot.xlabel("central channel frequency (MHz)")
                    pyplot.ylabel("gain amplitude (1)")
                    pyplot.xlim(plotFrequencyMin, plotFrequencyMax)
                    pyplot.ylim(plotGainAmplitudeMin, plotGainAmplitudeMax)
                    pyplot.title("bandpass (amplitude)\ndata set: " + self.MSName + " | telescope: " + self.telescopeName + " | antenna ID: " + str(antennaID) + " | calibrator: " + self.calibratorName)
                    pyplot.subplots_adjust(left = .08, right = .98)
                    pyplot.savefig(plotDirectory + "/cal" + self.calibratorName + "BPAmpAnt" + str(antennaID) + ".png")
                    pyplot.close()
                    
                    pyplot.figure(figsize = (12, 6))
                    pyplot.scatter(frequencyBand, numpy.ma.masked_array(gainPhasesPol1[antennaID], mask = flagsPol1[antennaID]), c = 'b', lw = 0, label = "polarization 1")
                    pyplot.scatter(frequencyBand, numpy.ma.masked_array(gainPhasesPol2[antennaID], mask = flagsPol2[antennaID]), c = 'r', lw = 0, label = "polarization 2")
                    pyplot.grid()
                    pyplot.legend()
                    pyplot.xlabel("central channel frequency (MHz)")
                    pyplot.ylabel("gain phase ($\degree$)")
                    pyplot.xlim(plotFrequencyMin, plotFrequencyMax)
                    pyplot.ylim(-180 - 5, 180 + 5)
                    pyplot.yticks(numpy.linspace(-180, 180, num = 9, endpoint = True))
                    pyplot.title("bandpass (phase)\ndata set: " + self.MSName + " | telescope: " + self.telescopeName + " | antenna ID: " + str(antennaID) + " | calibrator: " + self.calibratorName)
                    pyplot.subplots_adjust(left = .08, right = .98)
                    pyplot.savefig(plotDirectory + "/cal" + self.calibratorName + "BPPhaseAnt" + str(antennaID) + ".png")
                    pyplot.close()
    
    
    def invert_table(self):
        '''
        This method inverts a calibration table (by taking the multiplicative inverse of each number)
        and returns the new calibration table path.
        '''
        pathCalTableInverse = self.file_path + "_inv"
        
        # Create a copy of the old calibration table.
        systemCommand = "cp -r " + self.file_path + " " + pathCalTableInverse
        os.system(systemCommand)
        
        # Open the new Calibration Table, and load the (complex) gains.
        tb.open(pathCalTableInverse, nomodify = False)
        gVals = tb.getcol("CPARAM")#, startrow=start, nrow=incr)
        
        # Take the reciprocal of the (complex) gains, and save these to the Table.
        mask = abs(gVals) > 0.0 # only consider non-zero values
        gVals[mask] = 1.0 / gVals[mask] # This is the actual inversion.
        tb.putcol("CPARAM", gVals)#, startrow=start, nrow=incr) # replace the GAIN values with the inverted values
        tb.close()
        
        # Return the path of the inverted calibration table.
        return pathCalTableInverse
    
    
    def re_reference_table(self, refant = 1):
        """
        CASA's gaincal changes the reference antenna if the previous reference
        antenna is flagged. This method sets the varying reference antennas to 
        one antenna.
        """

        # Create a copy of the calibration table
        syscommand = 'rm -rf '+self.file_path+"_reref"
        os.system(syscommand)
        syscommand = "cp -rf "+self.file_path+" "+self.file_path+"_reref"
        os.system(syscommand)
        caltab = self.file_path+"_reref"

        # Open the calibration table and fetch needed parameters
        tb.open(caltab, nomodify=False)
        times = tb.getcol('TIME')
        gains = tb.getcol('CPARAM')
        gainsRR = gains[0,0,:] # RR polarization
        gainsLL = gains[1,0,:] # LL polarization
        ants1 = tb.getcol('ANTENNA1') # Antenna
        ants2 = tb.getcol('ANTENNA2') # Current reference antenna
        flags = tb.getcol('FLAG')
        flagsRR = flags[0,0,:] # RR polarization
        flagsLL = flags[1,0,:] # LL polarization

        # Create empty lists to be filled.
        updatedRefant = []
        updatedgainRR = []
        updatedgainLL = []
        updatedflagRR = []
        updatedflagLL = []
        
        # Create numpy array for easy access (drawback is that all values in
        # this 2D array are complex numbers now)
        table = numpy.column_stack((times, ants1, ants2, gainsRR, gainsLL, flagsRR, flagsLL))
        for time, ant1, ant2, gRR, gLL, fRR, fLL in table:
            if int(numpy.real(ant2)) is not refant:
                # This gain needs to be re-referenced. In order to do this we
                # need the gain corresponding with the new reference antenna.
                # I.e.:  g    = g    . g*   / |g   |
                #         i,n    i,o    n,o     n,o
                # where n is the new reference antenna, o the old one, | |
                # denotes the absolute value and * is the complex conjugate.

                # Find the gain relating the old reference antenna to the new
                # one.
                rerefgainRR = None
                rerefgainLL = None
                rerefflagRR = None
                rerefflagLL = None
                for time_2, ant1_2, ant2_2, gRR_2, gLL_2, fRR_2, fLL_2 in table:
                    if time_2 == time and int(numpy.real(ant1_2)) == refant and ant2_2 == ant2:
                        # This is the gain we need for re-referencing.
                        rerefgainRR = gRR_2
                        rerefgainLL = gLL_2
                        rerefflagRR = fRR_2
                        rerefflagLL = fLL_2
                        break # No need to search any further
                
                # Quit if no re-reference gain was found.
                if rerefgainRR == None:
                    print ("Unable to find correct antenna for re-referencing")
                    sys.exit()
                
                # Calculate new gains with respect to new reference antenna
                newgainRR = gRR * numpy.conj(rerefgainRR) / numpy.absolute(rerefgainRR)
                newgainLL = gLL * numpy.conj(rerefgainLL) / numpy.absolute(rerefgainLL)
                
                # Set flags to false if one of antennas is flagged
                newflagRR = True
                if bool(numpy.real(rerefflagRR)) is False and bool(numpy.real(fRR)) is False:
                    newflagRR = False
                
                newflagLL = True
                if bool(numpy.real(rerefflagLL)) is False and bool(numpy.real(fLL)) is False:
                    newflagLL = False
                
                # Correct for the fact that numpy.column_stack casts all values tp
                # complex numbers.
                updatedRefant.append( refant )
                updatedgainRR.append( newgainRR )
                updatedgainLL.append( newgainLL )
                updatedflagRR.append( newflagRR )
                updatedflagLL.append( newflagLL )
            else:
                # This gain is OK
                # Correct for the fact that numpy.column_stack casts all values tp
                # complex numbers.
                updatedRefant.append( int(numpy.real(ant2)) )
                updatedgainRR.append( gRR )
                updatedgainLL.append( gLL )
                updatedflagRR.append( bool(numpy.real(fRR)) )
                updatedflagLL.append( bool(numpy.real(fLL)) )
                        
        # Add new data to existing calibration table
        gains[0,0,:] = updatedgainRR
        gains[1,0,:] = updatedgainLL
        flags[0,0,:] = updatedflagRR
        flags[1,0,:] = updatedflagLL

        tb.putcol('CPARAM', gains)
        tb.putcol('ANTENNA2', updatedRefant)
        tb.putcol('FLAG', flags)

        tb.close()

        return caltab

    def normalize_reference_antenna(self):
        """
        Somehow, the calibration tables created by CASA's gaincal have a
        changing constant offset with respect to the reference antenna. For
        example, one would expect the reference antenna to have a phase of zero
        everywhere (after all, your measure the phase with respect to the phase
        of the reference antenna), but this is not the case. The phase of the
        reference antenna changes with constant values in time. This method
        corrects for this.
        """

        # Create a copy of the calibration table
        syscommand = 'rm -rf '+self.file_path+"_normref"
        os.system(syscommand)
        syscommand = "cp -rf "+self.file_path+" "+self.file_path+"_normref"
        os.system(syscommand)
        caltab = self.file_path+"_normref"

        # Open the calibration table and fetch needed parameters
        tb.open(caltab, nomodify=False)
        times = tb.getcol('TIME')
        gains = tb.getcol('CPARAM')
        gainsRR = gains[0,0,:] # RR polarization
        gainsLL = gains[1,0,:] # LL polarization
        ants1 = tb.getcol('ANTENNA1') # Antenna
        ants2 = tb.getcol('ANTENNA2') # Current reference antenna

        # Create empty lists to be filled
        updatedgainRR = gainsRR
        updatedgainLL = gainsLL

        # Create numpy array for easy access (drawback is that all values in
        # this 2D array are complex numbers now)
        table = numpy.column_stack((times, ants1, ants2, gainsRR, gainsLL))
        updated_table = table
        for time, ant1, ant2, gRR, gLL in table:
            if ant1 == ant2:
                if gRR != (1+0j) or gLL != (1+0j):
                    # Gains with this timestamp need to be fixed. Find the
                    # the gains with corresponding times.
                    for i, [time_2, ant1_2, ant2_2, gRR_2, gLL_2] in enumerate(table):
                        if time == time_2:
                            # Check which polarizations need to be fixed
                            if gRR != (1+0j):
                                # Fix RR polarization
                                updated_gRR_2 = gRR_2 * numpy.conjugate(gRR) / numpy.absolute(gRR)
                                updatedgainRR[i] = updated_gRR_2
                        
                            if gLL != (1+0j):
                                # Fix LL polarization
                                updated_gLL_2 = gLL_2 * numpy.conjugate(gLL) / numpy.absolute(gLL)
                                updatedgainLL[i] = updated_gLL_2

        # Add new data to existing calibration table
        gains[0,0,:] = updatedgainRR
        gains[1,0,:] = updatedgainLL
        tb.putcol("CPARAM", gains)
        
        tb.close()

        return caltab

    def resample_solutions(self, interpolationtimes, interp_type = 'spline'):
        """
        Calibration tables created by CASA's gaincal can have a smaller number
        of timestamps than the original measurement set (e.g. if the solution
        interval was larger than 'int'). This method resamples the calibration
        table to the original number of timestamps available in the measurement
        set by interpolation.
        
        Note that the calibration table created by this function does not have
        correct data in the following columns: FIELD_ID, INTERVAL, 
        OBSERVATION_ID, PARAMERR, SCAN_NUMBER, SNR, SPECTRAL_WINDOW_ID (these
        data are simply not interpolated). So that means that using these
        columns in applycal is not possible.
        """
        
        # Only use the unique timestamps
        newtimes = numpy.sort(numpy.unique(interpolationtimes))

        # Create a copy of the calibration table
        syscommand = 'rm -rf '+self.file_path+'_resamp'
        os.system(syscommand)
        syscommand = "cp -rf "+self.file_path+" "+self.file_path+"_resamp"
        os.system(syscommand)
        caltab = self.file_path+'_resamp'

        # Open the calibration table and fetch needed parameters
        tb.open(self.file_path, nomodify=False)
        
        tbkeyw = tb.getkeywords()
        tbcoltime = tb.getcolkeywords('TIME')
        tbcolant1 = tb.getcolkeywords('ANTENNA1')
        tbcolant2 = tb.getcolkeywords('ANTENNA2')
        tbcolcpar = tb.getcolkeywords('CPARAM')
        tbcolflag = tb.getcolkeywords('FLAG')
        
        tbcolfieldid = tb.getcolkeywords('FIELD_ID')
        tbcolint = tb.getcolkeywords('INTERVAL')
        tbcolobsid = tb.getcolkeywords('OBSERVATION_ID')
        tbcolpar = tb.getcolkeywords('PARAMERR')
        tbcolscann = tb.getcolkeywords('SCAN_NUMBER')
        tbcolsnr = tb.getcolkeywords('SNR')
        tbcolspwid = tb.getcolkeywords('SPECTRAL_WINDOW_ID')
        tbcolwei = tb.getcolkeywords('WEIGHT')
        
        tbinfo = tb.info()
        times = tb.getcol('TIME')
        gains = tb.getcol('CPARAM')
        gainsRR = gains[0,0,:] # RR polarization
        gainsLL = gains[1,0,:] # LL polarization
        ants1 = tb.getcol('ANTENNA1') # Antenna
        ants2 = tb.getcol('ANTENNA2') # Current reference antenna
        flags = tb.getcol('FLAG')
        flagsRR = flags[0,0,:] # RR polarization
        flagsLL = flags[1,0,:] # LL polarization
        tabdesc = tb.getdesc()  
        dminfo  = tb.getdminfo()
                
        tb.close()
        
        # Check if all reference antennas are the same
        uniqueants = numpy.unique(ants2)
        if len(uniqueants) > 1:
            # Multiple reference antennas in calibration table
            print ("Multiple reference antennas in calibration table")
            sys.exit()
        else:
            refant = uniqueants[0]

        # Separate different antennas: info_per_ant is a list containing
        # for each antenna a list with the time, gains and flags for both
        # polarizations
        uniqueants = numpy.unique(ants1)
        info_per_ant = []
        for unique_ant in uniqueants:
            this_ants_time = []
            this_ants_gRR = []
            this_ants_gLL = []
            this_ants_fRR = []
            this_ants_fLL = []
            for time, gRR, gLL, ant, flagRR, flagLL in zip(times, gainsRR, gainsLL, ants1, flagsRR, flagsLL):
                if ant == unique_ant:
                    this_ants_time.append(time)
                    this_ants_gRR.append(gRR)
                    this_ants_fRR.append(flagRR)
                    this_ants_gLL.append(gLL)
                    this_ants_fLL.append(flagLL)
            
            info_per_ant.append([this_ants_time, this_ants_gRR, this_ants_gLL, this_ants_fRR, this_ants_fLL])
        
        # For each antenna interpolate the data
        updated_ant_info = []
        for ant_info in info_per_ant:
            ant_time = ant_info[0]
            ant_gRR = ant_info[1]
            ant_gLL = ant_info[2]
            ant_fRR = ant_info[3]
            ant_fLL = ant_info[4]
        
            ant_phase_RR = numpy.angle(ant_gRR)
            ant_mag_RR = numpy.abs(ant_gRR)
            ant_phase_LL = numpy.angle(ant_gLL)
            ant_mag_LL = numpy.abs(ant_gLL)

            
            ### REMOVE FLAGGED VALUES ###
                    
            # Although if data is flagged, gains are still stored in the table.
            # This messes up the interpolation so remove these values first.
            updated_ant_timeRR = []
            updated_ant_phaseRR = []
            updated_ant_magRR = []
            for time, phaseRR, magRR, fRR in zip(ant_time, ant_phase_RR, ant_mag_RR, ant_fRR):
                if not fRR:
                    updated_ant_timeRR.append(time)
                    updated_ant_phaseRR.append(phaseRR)
                    updated_ant_magRR.append(magRR)

            updated_ant_timeLL = []
            updated_ant_phaseLL = []
            updated_ant_magLL = []
            for time, phaseLL, magLL, fLL in zip(ant_time, ant_phase_LL, ant_mag_LL, ant_fLL):
                if not fLL:
                    updated_ant_timeLL.append(time)
                    updated_ant_phaseLL.append(phaseLL)
                    updated_ant_magLL.append(magLL)

            ### UNWRAP PHASES ###

            updated_ant_phaseRR = numpy.unwrap(updated_ant_phaseRR)
            updated_ant_phaseLL = numpy.unwrap(updated_ant_phaseLL)


            ### INTERPOLATE VALUES ###
            
            if interp_type == 'spline':
                # Create the interpolation functions and new values
                if len(updated_ant_timeRR) > 4:
                    # Total number of points must be greater than the degree of the
                    # spline
                    splineRepPhaseRR = interpolate.splrep(updated_ant_timeRR, updated_ant_phaseRR, s=0)
                    splineRepMagRR = interpolate.splrep(updated_ant_timeRR, updated_ant_magRR, s=0)
                    newphasesRR = interpolate.splev(newtimes, splineRepPhaseRR, der=0)
                    newmagRR = interpolate.splev(newtimes, splineRepMagRR, der=0)
                else:
                    # (almost) everything is flagged, so just fill it with dummy
                    # values
                    newphasesRR = numpy.linspace(0,0,len(newtimes))
                    newmagRR = numpy.linspace(1,1,len(newtimes))
                
                if len(updated_ant_timeLL) > 4:
                    # Total number of points must be greater than the degree of the
                    # spline
                    splineRepPhaseLL = interpolate.splrep(updated_ant_timeLL, updated_ant_phaseLL, s=0)
                    splineRepMagLL = interpolate.splrep(updated_ant_timeLL, updated_ant_magLL, s=0)
                    newphasesLL = interpolate.splev(newtimes, splineRepPhaseLL, der=0)
                    newmagLL = interpolate.splev(newtimes, splineRepMagLL, der=0)
                else:
                    # (almost) everything is flagged, so just fill it with dummy
                    # values
                    newphasesLL = numpy.linspace(0,0,len(newtimes))
                    newmagLL = numpy.linspace(1,1,len(newtimes))

            elif interp_type == 'linear':
                # Create the interpolation functions and new values
                if len(updated_ant_timeRR) > 2:
                    newphasesRR = numpy.interp(newtimes, updated_ant_timeRR, updated_ant_phaseRR)
                    newmagRR = numpy.interp(newtimes, updated_ant_timeRR, updated_ant_magRR)
                else:
                    # (almost) everything is flagged, so just fill it with dummy
                    # values
                    newphasesRR = numpy.linspace(0,0,len(newtimes))
                    newmagRR = numpy.linspace(1,1,len(newtimes))
                
                if len(updated_ant_timeLL) > 2:
                    newphasesLL = numpy.interp(newtimes, updated_ant_timeLL, updated_ant_phaseLL)
                    newmagLL = numpy.interp(newtimes, updated_ant_timeLL, updated_ant_magLL)
                else:
                    # (almost) everything is flagged, so just fill it with dummy
                    # values
                    newphasesLL = numpy.linspace(0,0,len(newtimes))
                    newmagLL = numpy.linspace(1,1,len(newtimes))

            new_ant_gRR = newmagRR * numpy.exp(1j * newphasesRR)
            new_ant_gLL = newmagLL * numpy.exp(1j * newphasesLL)

            
            ### FIX INTERPOLATION BETWEEN FLAGS ###
                    
            # Fix interpolation between flagged values. This comes basically
            # down to linear interpolation between flag values (zero and one).
            # Values < 1 will next be regarded as flagged.
            ant_fRR_num = map(int, ant_fRR)
            new_ant_fRR = numpy.interp(newtimes, ant_time, ant_fRR_num)
            new_ant_fRR = numpy.ceil(new_ant_fRR)
            new_ant_fRR = map(bool, new_ant_fRR)
            
            ant_fLL_num = map(int, ant_fLL)
            new_ant_fLL = numpy.interp(newtimes, ant_time, ant_fLL_num)
            new_ant_fLL = numpy.ceil(new_ant_fLL)
            new_ant_fLL = map(bool, new_ant_fLL)


            ### FLAG LARGE GAPS ###

            # Sometimes, for larger gaps, timestamps exist, flags are not
            # true but data is not available. This means that for the entire
            # gap the interpolater can do what it wants. So these large gaps
            # need to be flagged as well.
            # First RR
            for i in xrange(len(updated_ant_timeRR)-1):
                # Find out how many interpolated steps are between two data
                # points.
                time = updated_ant_timeRR[i]
                nexttime = updated_ant_timeRR[i+1]
                
                # Since we're dealing with a calibration table with less
                # timestamps, we need to find the timestamps in the measurement
                # set which correspond best with those in the cal table.
                time_matches = []
                nexttime_matches = []
                for match in newtimes:
                    time_matches.append([match-time, match])
                    nexttime_matches.append([match-nexttime, match])
                time_matches.sort(key=lambda x: x[0])
                nexttime_matches.sort(key=lambda x: x[0])
                best_match_with_time = time_matches[0][1]
                best_match_with_nexttime = nexttime_matches[0][1]
                
                # Calculate how many intermediate steps there are
                counter = 0
                start_counter = False
                for interp_time in newtimes:
                    if interp_time == best_match_with_time:
                        start_counter = True
                    if start_counter:
                        counter += 1
                    if interp_time == best_match_with_nexttime:
                        start_counter = False
            
                # If the numbers of steps is too large, flag the gap.
                expected_steps = int(len(newtimes)/numpy.float32(len(updated_ant_timeRR)) + 2)
                if counter > expected_steps:
                    # Flag this gap
                    do_flag = False
                    for i, interp_time in enumerate(newtimes):
						# Only flag values between the two data steps.
                        if interp_time == best_match_with_nexttime:
                            do_flag = False
                        if do_flag:
                            new_ant_fRR[i] = True
                        if interp_time == best_match_with_time:
                            do_flag = True
            # Also do LL
            for i in xrange(len(updated_ant_timeLL)-1):
                # Find out how many interpolated steps are between two data
                # points.
                time = updated_ant_timeLL[i]
                nexttime = updated_ant_timeLL[i+1]
                
                # Since we're dealing with a calibration table with less
                # timestamps, we need to find the timestamps in the measurement
                # set which correspond best with those in the cal table.
                time_matches = []
                nexttime_matches = []
                for match in newtimes:
                    time_matches.append([match-time, match])
                    nexttime_matches.append([match-nexttime, match])
                time_matches.sort(key=lambda x: x[0])
                nexttime_matches.sort(key=lambda x: x[0])
                best_match_with_time = time_matches[0][1]
                best_match_with_nexttime = nexttime_matches[0][1]
                
                # Calculate how many intermediate steps there are
                counter = 0
                start_counter = False
                for interp_time in newtimes:
                    if interp_time == best_match_with_time:
                        start_counter = True
                    if start_counter:
                        counter += 1
                    if interp_time == best_match_with_nexttime:
                        start_counter = False
            
                # If the numbers of steps is too large, flag the gap.
                expected_steps = int(len(newtimes)/numpy.float32(len(updated_ant_timeLL)) + 2)
                if counter > expected_steps:
                    # Flag this gap
                    do_flag = False
                    for i, interp_time in enumerate(newtimes):
						# Only flag values between the two data steps.
                        if interp_time == best_match_with_nexttime:
                            do_flag = False
                        if do_flag:
                            new_ant_fLL[i] = True
                        if interp_time == best_match_with_time:
                            do_flag = True


            ### SET FLAGGED VALUES TO (1+0j) ###

            # This is just the casa convention. Doesn't really matter since
            # the data is flagged.
            for i, flag in enumerate(new_ant_fRR):
                if flag:
                    new_ant_gRR[i] = (1+0j)
            for i, flag in enumerate(new_ant_fLL):
                if flag:
                    new_ant_gLL[i] = (1+0j)

            # Store this new interpolation data
            updated_ant_info.append([newtimes, new_ant_gRR, new_ant_gLL, new_ant_fRR, new_ant_fLL])

        # Put the updated info back in the correct places in a new cal table
        correct_format_times = []
        correct_format_antids = []
        correct_format_antrefids = []
        correct_format_gRR = []
        correct_format_gLL = []
        correct_format_fRR = []
        correct_format_fLL = []
        for time in newtimes:
            for ant_id, updated_info in enumerate(updated_ant_info):
                ant_times = updated_info[0]
                ant_gRR = updated_info[1]
                ant_gLL = updated_info[2]
                ant_fRR = updated_info[3]
                ant_fLL = updated_info[4]
                for i, ant_time in enumerate(ant_times):
                    if ant_time == time:
                        correct_format_times.append(ant_time)
                        correct_format_antids.append(ant_id)
                        correct_format_antrefids.append(refant)
                        correct_format_gRR.append(ant_gRR[i])
                        correct_format_gLL.append(ant_gLL[i])
                        correct_format_fRR.append(ant_fRR[i])
                        correct_format_fLL.append(ant_fLL[i])

        updatedgains = numpy.zeros([2,1,len(correct_format_gRR)], dtype=numpy.complex)
        updatedflags = numpy.zeros([2,1,len(correct_format_gRR)], dtype=bool)
        updatedgains[0,0,:] = correct_format_gRR
        updatedgains[1,0,:] = correct_format_gLL
        updatedflags[0,0,:] = correct_format_fRR
        updatedflags[1,0,:] = correct_format_fLL
        
        # Copy existing table information
        tb.create(caltab, tabdesc, dminfo=dminfo)
        tb.addrows(len(correct_format_times))
        tb.putinfo(tbinfo)
        tb.putkeywords(tbkeyw)
        tb.putcolkeywords('TIME', tbcoltime)
        tb.putcolkeywords('ANTENNA1', tbcolant1)
        tb.putcolkeywords('ANTENNA2', tbcolant2)
        tb.putcolkeywords('CPARAM', tbcolcpar)
        tb.putcolkeywords('FLAG', tbcolflag)
        tb.putcolkeywords('FIELD_ID', tbcolfieldid)
        tb.putcolkeywords('INTERVAL', tbcolint)
        tb.putcolkeywords('OBSERVATION_ID', tbcolobsid)
        tb.putcolkeywords('PARAMERR', tbcolpar)
        tb.putcolkeywords('SCAN_NUMBER', tbcolscann)
        tb.putcolkeywords('SNR', tbcolsnr)
        tb.putcolkeywords('SPECTRAL_WINDOW_ID', tbcolspwid)
        tb.putcolkeywords('WEIGHT', tbcolwei)
        
        # Useful stuff
        tb.putcol('TIME', correct_format_times)
        tb.putcol('CPARAM', updatedgains)
        tb.putcol('ANTENNA1', correct_format_antids)
        tb.putcol('ANTENNA2', correct_format_antrefids)
        tb.putcol('FLAG', updatedflags)
        
        # Needed but not useful
        emptylist = numpy.linspace(0,0,len(correct_format_times))
        secondemptylist = numpy.zeros_like(updatedgains, dtype=numpy.float)
        tb.putcol('FIELD_ID', emptylist)
        tb.putcol('INTERVAL', emptylist)
        tb.putcol('OBSERVATION_ID', emptylist)
        tb.putcol('PARAMERR', secondemptylist) # different format
        tb.putcol('SCAN_NUMBER', emptylist)
        tb.putcol('SNR', secondemptylist) # different format
        tb.putcol('SPECTRAL_WINDOW_ID', emptylist)
        # The WEIGHT COLUMN WAS EMPTY SO LEAVE IT EMPTY
        tb.close()

        return caltab