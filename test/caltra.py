import datetime
import numpy as np
from lagranto import caltra

mapping = {datetime.datetime(2009, 11, 30, 0): 'datadir/xjjhq/xjjhqa_036.pp',
           datetime.datetime(2009, 11, 30, 1): 'datadir/xjjhq/xjjhqa_036.pp'}

trainp = np.zeros([2, 3])

trainp[:, 0] = 360
trainp[0, 2] = 5000
trainp[1, 2] = 10000

traout = caltra.caltra(trainp, mapping, tracers=['advection_only_pv',
                                                 'ertel_potential_vorticity'])

for traj in traout:
    print traj
