import sys
import os
import numpy as np
import logging
import ConfigParser # In Python 3, this module was renamed to 'configparser'
import argparse

# CSPAM Modules
from lib import TableObjects
from lib import utils
import steps

def get_conf(config_file = "cspam.config"):
    """
    Prepare a config dictionary with all the user-defined parameters
    """
    if (not os.path.isfile(config_file)):
        logging.critical('Configuration file '+config_file+' not found.')
        sys.exit(1)
    
    confP = ConfigParser.ConfigParser()
    confP.read(config_file)

    conf = {}
    conf['steps'] = confP.get('DEFAULT','steps').replace(' ','').split(',')
    conf['data_dir'] = confP.get('DEFAULT','data_dir')

    # creating MSs from sections
    conf['MSs'] = confP.sections()
    for MS in conf['MSs']:
        conf[MS] = {}
        conf[MS]['file_path'] = conf['data_dir'] + '/' + MS
        conf[MS]['flag'] = confP.get(MS, 'flag')
        conf[MS]['cal_scans'] = confP.get(
				MS, 'cal_scans').replace(' ','').split(',')
        conf[MS]['tgt_scans'] = confP.get(
				MS, 'tgt_scans').replace(' ','').split(',')
        conf[MS]['spw'] = confP.get(MS, 'spw').replace(' ','')
        conf[MS]['central_chan_percentage'] = confP.getint(
					      MS, 'central_chan_percentage')
        conf[MS]['flux_cal_field'] = confP.get(MS, 'flux_cal_field')
        conf[MS]['phase_cal_field'] = confP.get(MS, 'phase_cal_field')
        conf[MS]['leakage_cal_field'] = confP.get(MS, 'leakage_cal_field')
        conf[MS]['position_cal_field'] = confP.get(MS, 'position_cal_field')
        conf[MS]['angle_cal_field'] = confP.get(MS, 'angle_cal_field')
    
    return conf

def parseCommandLineArguments():
    """
    Set up command line parsing.
    """
    parser = argparse.ArgumentParser(description="""CSPAM stands for CASA SPAM: 
             Source Peeling and Atmospheric Modeling""")
    
    parser.add_argument("-c", "--config", help = """Path to the configuration file""", 
                        type = str, required = False)
    parser.add_argument("-v", "--verbose", action = "store_true", 
                        dest = "verbosity", help = """Make PNG plot""")
    
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    # Start writing debug messages to the file 'cspam.log'.
    logging.basicConfig(filename = 'cspam.log', level = logging.DEBUG)
    
    # Load the arguments provided by the user, such as an alternative path to the configuration file to use.
    args = parseCommandLineArguments()    
    config_file = args['config']
    #verbose = args['verbosity'] # This argument is not yet used.
    
    if (config_file):
        conf = get_conf(config_file = config_file)
    else: # Assume use of default configuration file 'cspam.config'.
        conf = get_conf()
    
    print ("Starting CSPAM version: 0.x")
    print ("Parsed configuration file:")
    utils.print_dict(conf)
    
    
    # Create a list of instances of MSObj - one for each MS to act upon.
    MSs = []
    for MSinConf in conf['MSs']:
        print ("MS name: " + MSinConf)
        print ("MS file path: " + conf[MSinConf]['file_path'])
        print ("MS configuration file input: " + str(conf[MSinConf]))
        MSs.append(TableObjects.MSObj(conf[MSinConf]['file_path'], conf = conf[MSinConf])) # Note that the MSObj class might update field names in the MS
    
    
    # Carry out the wanted steps per MS.
    MSIndex = 0 # Only used in print function
    for mset in MSs:
        # Announce which MS to work on.
        MSString = "MS " + str(MSIndex + 1)
        print ("Starting work on " + MSString + " of " + str(len(MSs)) + ".")
        MSIndex += 1
        
        
        # Execute the wanted steps.
        if ('plots' in conf['steps']):
            print (MSString + ") Starting step: plots")
            steps.plots(mset)

        if ('preflag' in conf['steps']):
            print (MSString + ") Starting step: preflag")
            steps.preflag(mset)

        if ('setjy' in conf['steps']):
            print (MSString + ") Starting step: setjy")
            steps.set_flux_density_scale(mset)

        if ('bandpass' in conf['steps']):
            print (MSString + ") Starting step: bandpass")
            steps.bandpass_calibration(mset)

        if ('cal' in conf['steps']):
            print (MSString + ") Starting step: cal")
            steps.calib(mset)

        if ('selfcal' in conf['steps']):
            print (MSString + ") Starting step: selfcal")
            steps.selfcal(mset)

        if ('peeling' in conf['steps']):
            print (MSString + ") Starting step: peeling")
            steps.peeling(mset)

        if ('createimage' in conf['steps']):
            print (MSString + ") Starting step: createimage")
            steps.createimage(mset)