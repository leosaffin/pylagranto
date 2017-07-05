import datetime
import numpy as np
import matplotlib.pyplot as plt
import iris
from mymodule import convert, grid, plot
from lagranto import caltra
from scripts import case_studies

datadir = '/home/lsaffin/Documents/meteorology/data/iop5/'
mapping = {
    datetime.datetime(2011, 11, 28, 12): datadir + '20111128_analysis12.nc',
    datetime.datetime(2011, 11, 28, 18): datadir + 'prognostics_006.nc',
    datetime.datetime(2011, 11, 29, 0): datadir + 'prognostics_012.nc',
    datetime.datetime(2011, 11, 29, 6): datadir + 'prognostics_018.nc',
    datetime.datetime(2011, 11, 29, 12): datadir + 'prognostics_024.nc',
    # datetime.datetime(2011, 11, 29, 18): datadir + 'prognostics_030.nc',
    # datetime.datetime(2011, 11, 30, 0): datadir + 'prognostics_036.nc',
}


# Select all start positions that have experiences large latent heating
forecast = case_studies.iop5b.copy()
cubes = forecast.set_lead_time(hours=24)
theta_adv = convert.calc('advection_only_theta', cubes,
                         levels=('air_potential_temperature', [320]))[0]
z = convert.calc('altitude', cubes,
                 levels=('air_potential_temperature', [320]))[0]

ny, nx = theta_adv.shape
lon, lat = grid.get_xy_grids(theta_adv)

# Set up an input array of trajectories
trainp = []
for i in range(nx):
    for j in range(ny):
        if theta_adv.data[j, i] < 310:
            trainp.append([lon[j, i], lat[j, i], z[j, i].data])
trainp = np.array(trainp)


# Calculate the trajectories
traout = caltra.caltra(trainp, mapping, fbflag=-1,
                       tracers=['air_potential_temperature',
                                'air_pressure'])

traout.to_pickle(datadir + 'test_trajectories.pkl')

# Plot the lat/lon positions of trajectories
plt.figure()
for n in range(len(traout)):
    lc = plot.colored_line_plot(
        traout[n].grid_longitude, traout[n].grid_latitude,
        traout[n].air_potential_temperature,
        vmin=300, vmax=350, cmap='viridis')

plt.xlim(330, 370)
plt.ylim(-15, 15)
plt.colorbar(lc)
plt.show()
