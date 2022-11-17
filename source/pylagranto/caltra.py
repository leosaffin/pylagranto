import numpy as np
from pylagranto import trajectory
import pylagranto.fortran


def caltra(trainp, times, datasource,
           imethod=1, numit=3, nsubs=4, fbflag=1, jflag=False,
           tracers=[]):
    """Calculate a set of trajectories from the input points

    args:
        trainp (np.array): An array of start positions of the trajectories.
            Must have shape (number of trajectories x 3) with the 3 being the
            (x,y,z) co-ordinates

        times (list): List of times in the wind data used for the trajectory
            calculations

        datasource (pylagranto.datasets.DataSource):

        imethod (int, optional): Numerical method of integration. 1=Euler,
            2=Runge-Kutta. Default is 1.

        numit (int, optional): Number of iterations for the Euler method.
            Default is 3

        nsubs (int, optional): Number of substeps to take between files.
            Default is 4.

        fbflag (int, optional): Forward trajectories (1) or reverse
            trajectories (-1). Default is 1.

        jflag (logical, optional): Flag for whether trajectories re-enter the
            atmosphere on hitting the ground. Default is False

        tracers (list): A list of variable names to trace at each point along
            the trajectory

    returns:
        traout (Lagranto.trajectory.TrajectoryEnsemble):
    """
    # Initialise output
    ntra = len(trainp)   # Number of trajectories
    ntim = len(times)  # Number of files at different times
    traout = np.zeros([ntra, ntim, 3 + len(tracers)])

    # Extract starting trajectory positions
    x = trainp[:, 0]
    y = trainp[:, 1]
    z = trainp[:, 2]

    # Initialise the flag and the counter for trajectories leaving the domain
    leftflag = np.zeros(ntra)

    # Calulate the timestep in seconds between input files and divide by
    # number of substeps
    ts = abs((times[1] - times[0]).total_seconds()) / nsubs

    # Reverse file load order for backward trajectories
    times.sort()
    if fbflag == -1:
        times.reverse()

    # Loop over all input files
    for n, time in enumerate(times):
        print(time)

        if n > 0:
            # Copy old velocities and pressure fields to new ones
            uut0 = uut1.copy()
            vvt0 = vvt1.copy()
            wwt0 = wwt1.copy()
            p3t0 = p3t1.copy()
            spt0 = spt1.copy()

        # Read wind fields and surface pressure at next time
        datasource.set_time(time)
        spt1, uut1, vvt1, wwt1, p3t1 = datasource.winds()

        if n == 0:
            # Load grid parameters at first timestep
            nx, ny, nz, xmin, ymin, dx, dy, hem, per, names = \
                datasource.grid_parameters()

            # Add the list of traced variables to the names in the output
            names += tracers

        elif n > 0:
            # Call fortran routine to update trajectory positions
            x, y, z, leftflag = pylagranto.fortran.caltra.main(
                x, y, z, leftflag, ts, nsubs, imethod, numit, jflag, fbflag,
                spt0, spt1, p3t0, p3t1, uut0, uut1, vvt0, vvt1, wwt0, wwt1,
                xmin, ymin, dx, dy, per, hem, nx, ny, nz, ntra)

        # Save positions
        traout[:, n, 0] = x
        traout[:, n, 1] = y
        traout[:, n, 2] = z

        # Trace additional fields
        for m, tracer in enumerate(tracers):
            try:
                array = datasource.get_variable(tracer)
            except ValueError:
                # If variable can't be loaded print a warning and put zero
                print("Variable {} not available at this time."
                      "Replacing with zeros".format(tracer))
                array = np.zeros_like(uut1)
            traout[:, n, m + 3] = pylagranto.fortran.trace.interp_to(
                array, x, y, z, leftflag, p3t1, spt1, xmin, ymin,
                dx, dy, nx, ny, nz, ntra)

    return trajectory.TrajectoryEnsemble(traout, times, names)

