import numpy as np
import pandas as pd
import math
import glob
import sys
import os

sys.path.append('../PMI_analysis/pyext/src/')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################

nproc = 10
top_dir =  sys.argv[1] 
analysis_dir = f'{top_dir}/analysis/'

# Check if analysis dir exists
if not os.path.isdir(analysis_dir):
    os.makedirs(analysis_dir)

# How are the trajectories dir names
dir_head = 'run_'
out_dirs = glob.glob(f'{top_dir}/{dir_head}*/output/')

################################
# Get and organize fields for
# analysis
################################
# Read the total score, plot
# and check for score convengence
XLs_cutoffs = {'DSS':30.0}

# Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name=dir_head,
                          analysis_dir = analysis_dir,
                          nproc=nproc,
                          number_models_out=10000,
                          nskip=5)

# Define restraints to analyze
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs,
                             ambiguous_XLs_restraint=True)
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()

# Read stat files
AT.read_stat_files()
AT.write_models_info()


AT.hdbscan_clustering([ 'XLs_sum', 'EV_sum','CR_sum'],
                      min_cluster_size=50000)
AT.summarize_XLs_info()
exit()


