import numpy as np
import iris
from iris.analysis import Linear
from mymodule import convert, grid
from lagranto import pyLagranto
from lagranto.trajectory import TrajectoryEnsemble


def caltra(trainp, mapping, imethod=1, numit=3, nsubs=4, fbflag=1, jflag=False,
           tracers=[], levels=None):
    """Calculate a set of trajectories from the input points

    args:
        trainp (np.array): An array of start positions of the trajectories.
            Must have shape (number of trajectories x 3) with the 3 being the
            (x,y,z) co-ordinates

        mapping (dict): A mapping between datetime objects and filenames to be
            loaded for the trajectory calculations

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

        levels (tuple or None): The name of a vertical coordinate and list of
            values to interpolate fields to prior to performing the trajectory
            calculations

    returns:
        traout (Lagranto.trajectory.TrajectoryEnsemble):
    """
    # Initialise output
    ntra = len(trainp)   # Number of trajectories
    ntim = len(mapping)  # Number of files at different times
    traout = np.zeros([ntra, ntim, 3 + len(tracers)])

    # Extract starting trajectory positions
    x = trainp[:, 0]
    y = trainp[:, 1]
    z = trainp[:, 2]

    # Initialise the flag and the counter for trajectories leaving the domain
    leftflag = np.zeros(ntra)

    # Extract times relating to filenames
    times = sorted(list(mapping))

    # Calulate the timestep in seconds between input files and divide by
    # number of substeps
    ts = (times[1] - times[0]).total_seconds() / nsubs

    # Reverse file load order for reverse trajectories
    if (fbflag == -1):
        times.reverse()

    # Loop over all input files
    for n, time in enumerate(times):
        print time

        if n > 0:
            # Copy old velocities and pressure fields to new ones
            uut0 = uut1.copy()
            vvt0 = vvt1.copy()
            wwt0 = wwt1.copy()
            p3t0 = p3t1.copy()
            spt0 = spt1.copy()

        # Read wind fields and surface pressure at next time
        cubes = iris.load(mapping[time])
        spt1, uut1, vvt1, wwt1, p3t1 = load_winds(cubes, levels)

        if n == 0:
            # Load grid parameters at first timestep
            example_cube = convert.calc('upward_air_velocity', cubes,
                                        levels=levels)
            nx, ny, nz, xmin, ymin, dx, dy, hem, per, names = \
                grid_parameters(example_cube, levels)

            # Add the list of traced variables to the names in the output
            names += tracers

        elif n > 0:
            # Call fortran routine to update trajectory positions
            x, y, z, leftflag = pyLagranto.caltra.main(
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
                cube = convert.calc(tracer, cubes, levels=levels)
                array = cube.data.transpose().flatten(order='F')
            except ValueError:
                # If variable can't be loaded print a warning and put zero
                print ('Variable ' + tracer + ' not available at this time. ' +
                       'Replacing with zeros')
                array = np.zeros_like(uut1)
            traout[:, n, m + 3] = pyLagranto.trace.interp_to(
                array, x, y, z, leftflag, p3t1, spt1, xmin, ymin,
                dx, dy, nx, ny, nz, ntra)

    return TrajectoryEnsemble(traout, times, names)


def grid_parameters(cube, levels):
    """Extract grid parameters for calculations from cube

    args:
        cube (iris.cube.Cube):

    returns:
        nz, ny, nx (int): Grid dimensions.

        xmin, ymin (float): Minimum longitude and latitude.

        dx, dy (float): Grid spacing in degrees.

        hem, per (int): Flag for whether the domain is hemispheric and/or
            periodic.

        names (list of str): The names of the dimensional coordinates
    """
    # Extract grid dimesions
    nz, ny, nx = cube.shape
    x = grid.extract_dim_coord(cube, 'x')
    y = grid.extract_dim_coord(cube, 'y')

    if levels is None:
        names = [x.name(), y.name(), 'altitude']
    else:
        names = [x.name(), y.name(), levels[0]]

    # Find minimum latitude and longitude
    xmin = x.points.min()
    xmax = x.points.max()
    ymin = y.points.min()
    ymax = y.points.max()

    # Grid spacing
    dx = (x.points[1:] - x.points[:-1]).mean()
    dy = (y.points[1:] - y.points[:-1]).mean()

    # Set logical flag for periodic data set (hemispheric or not)
    if abs(xmax - xmin - dx - 360) < dx:
        hem = 1
        per = 360
    else:
        hem = 0
        per = 0

    return nx, ny, nz, xmin, ymin, dx, dy, hem, per, names


def load_winds(cubes, levels):
    """Load the wind fields from a cubelist
    """
    # Extract fields needed as cubes
    u = convert.calc('x_wind', cubes, levels=levels)
    v = convert.calc('y_wind', cubes, levels=levels)

    if levels is None:
        # Default is height based coordinate
        w = convert.calc('upward_air_velocity', cubes, levels=levels)
        z = grid.make_cube(w, 'altitude')
        surface = grid.make_cube(w[0], 'surface_altitude')
    else:

        z = np.zeros_like(v.data)
        # Create a uniform height coordinate for the levels interpolated to
        for n, level in enumerate(levels[1]):
            z[n, :, :] = level
        z = v.copy(data=z)

        # Set vertical velocity and surface height to zero
        w = v.copy(data=np.zeros_like(v.data))
        surface = w[0].copy(data=np.zeros_like(w[0].data))

    # Return fields as 1d arrays with size nx*ny*nz
    return [x.data.transpose().flatten(order='F') for x in surface, u, v, w, z]
