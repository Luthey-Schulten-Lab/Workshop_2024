from datetime import datetime
import numpy as np
import os
from Bio import SeqIO
import importlib
import sys
# add plotting scripts' folder as system path
sys.path.append('/home/enguang/CMEODE/CMEODE_Hook/analysis/wcm_analyses')


import WCM_analysis
importlib.reload(WCM_analysis)

import WCM_gene as gene

log = open('./analyse_initiation.log', 'w')
origin_stdout = sys.stdout
origin_stderr = sys.stderr

sys.stdout = log
sys.stderr = log

merging_required = True

in_dirs = ['./WithVolumeChange/', './WithoutVolumeChange/']

in_label = 'initiation'

reps = np.arange(1, 50+1)

pkl_dirs = in_dirs

pkl_label = 'intiation'

if merging_required == True:
        
        for i_change, in_dir in enumerate(in_dirs):
                 

            w = WCM_analysis.WCM_ensemble()

            # this assumes trajectories are in the below format
            # (in_dir)/(in_label).(replicate).csv
            w.set_traj_files(in_dir,in_label,reps)

            w.load_trajs()

            w.merge_trajs()

            pkl_dir = pkl_dirs[i_change]

            w.write_merged_ensemble(pkl_dir,pkl_label)

            print('New ensemble of trajectories Created in folder {0}'.format(pkl_dir))


sys.stdout = origin_stdout
sys.stderr = origin_stderr

