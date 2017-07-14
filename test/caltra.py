import datetime
from lagranto import caltra, geometry
from scripts import case_studies

# Set up the mapping of times to files:
datadir = '/home/lsaffin/Documents/meteorology/data/iop5/'
start_time = datetime.datetime(2011, 11, 28, 12)
mapping = {start_time: datadir + '20111128_analysis12.nc'}
for n in range(6, 25, 6):
    time = start_time + datetime.timedelta(hours=n)
    mapping[time] = datadir + 'prognostics_' + str(n).zfill(3) + '.nc'

# Define the start points
forecast = case_studies.iop5b.copy()
cubes = forecast.set_lead_time(hours=24)
levels = ('air_potential_temperature', [315])
trainp = geometry.select(cubes, 'total_minus_advection_only_theta', '>', 10,
                         levels=levels)
#trainp = geometry.select_inflow_region(datadir + 'backward_trajectories.pkl')
#trainp = geometry.circle((366, 6.5), 6, [320], 0.5)

# Calculate the trajectories
traout = caltra.caltra(trainp, mapping, fbflag=-1,
                       tracers=['air_potential_temperature', 'air_pressure'])

# Save the trajectories
traout.to_pickle(
    datadir + 'backward_trajectories_from_heated_region_315K.pkl')
