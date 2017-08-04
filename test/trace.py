import datetime
import numpy as np
from lagranto import trajectory, trace
from scripts import case_studies

# IOP5
forecast = case_studies.iop5b.copy()
datadir = '/home/lsaffin/Documents/meteorology/data/iop5/'
start_time = datetime.datetime(2011, 11, 28, 12)
mapping = {start_time: datadir + '20111128_analysis12.nc'}


# IOP8
"""
forecast = case_studies.iop8.copy()
datadir = '/projects/diamet/lsaffi/iop8/'
start_time = datetime.datetime(2011, 12, 7, 12)
mapping = {start_time: datadir + '20111207_analysis12.nc'}
"""

# Create mapping of times to files
for n in range(6, 25, 6):
    time = start_time + datetime.timedelta(hours=n)
    mapping[time] = datadir + '*_' + str(n).zfill(3) + '.nc'

# Define variables to calculate along trajectories
tracers = ['long_wave_radiation_pv', 'microphysics_pv']

# Load trajectories
filename = datadir + 'backward_trajectories_from_circle_320K.pkl'
trainp = trajectory.load(filename)

# Calculate the trajectories
print len(trainp)
traout = trace.trace(trainp, tracers, mapping)

# Save the trajectories
traout.save(filename)
