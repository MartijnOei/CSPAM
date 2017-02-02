'''
This file contains two class definitions.
Both 'CalibratorSource' and 'TargetSource' are simple classes, containing only variables.
Written by Kasper van Dam (2016), edited by Martijn Oei (2017).
'''


class CalibratorSource:
    
    def __init__(self, caltype, field_name, field_id, scan_ids):
        
        self.caltype = caltype # String indicating calibration type (e.g. 'flux' or 'phase')
        
        self.field_name = field_name # String indicating field name
        self.field_id = field_id
        self.scan_ids = scan_ids # List of strings indicating scans corresp. to this calibrator
        self.scans = ','.join(self.scan_ids) # String with the scan ids (CASA format)
        
        # The files associated with this calibrator are saved in subdirectories. Example: "/flux_cal_3C147"
        self.extend_dir = '/' + caltype + "_cal_" + self.field_name # Path extension (string used for easy directory creation and access)
        
        self.derived_cal_tables = [] # List of calibration tables derived from this calibrator


class TargetSource:

    def __init__(self, field_name, field_id, scan_ids):
        
        self.field_name = field_name
        self.field_id = field_id
        self.scan_ids = scan_ids # List of strings indicating scans corresp. to this source
        self.scans = ','.join(self.scan_ids) # String with the scan ids (CASA format)
        
        self.extend_dir = '/self_cal_'+self.field_name # Path extension (string used for easy directory creation and access)
        
        self.self_cal_gp_tables = [] # List of self-cal tables derived from this target
        self.self_cal_ga_tables = [] # List of self-cal tables derived from this target