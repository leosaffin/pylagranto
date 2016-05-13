import numpy as np
from iris.analysis import Linear
from mymodule import convert, files, grid
from lagranto import pyLagranto


def caltra(trainp, mapping, imethod=1, numit=3, nsubs=4, fbflag=1, jflag=False,
           wfactor=100, tracers=[]):
    """
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

        wfactor (float, optional): Factor for difference in units for vertical
            velocity. Default is 100 because Lagranto is retarded and divides
            by 100.

        tracers (list): A list of variable names to trace at each point along
            the trajectory
    returns:
    """
    # Initialise output
    ntra = len(trainp)   # Number of trajectories
    ntim = len(mapping)  # Number of files at different times
    traout = np.zeros([ntra, ntim, 4 + len(tracers)])

    # Extract trajectory positions
    xx0 = trainp[:, 0]
    yy0 = trainp[:, 1]
    pp0 = trainp[:, 2]

    # Initialise the flag and the counter for trajectories leaving the domain
    leftflag = np.zeros(ntra)

    # Save starting positions
    traout[:, 0, 0] = 0.
    traout[:, 0, 1] = xx0
    traout[:, 0, 2] = yy0
    traout[:, 0, 3] = pp0

    # Extract times relating to filenames
    times = mapping.keys()
    times.sort()

    # Calulate the timestep in seconds between input files
    ts = (times[1] - times[0]).total_seconds()

    # Reverse file load order for reverse trajectories
    if (fbflag == -1):
        times.reverse()

    # Read wind field and grid from first file
    cubes = files.load(mapping[times[0]])
    spt1, uut1, vvt1, wwt1, p3t1 = load_winds(cubes)

    example_cube = convert.calc('upward_air_velocity', cubes)
    nx, ny, nz, xmin, ymin, dx, dy, hem, per = grid_parameters(example_cube)

    # Loop over all input files
    for n, time in enumerate(times[1:], start=1):
        print time
        # Copy old velocities and pressure fields to new ones
        uut0 = uut1.copy()
        vvt0 = vvt1.copy()
        wwt0 = wwt1.copy()
        p3t0 = p3t1.copy()
        spt0 = spt1.copy()

        # Read wind fields and surface pressure at next time
        cubes = files.load(mapping[time])
        spt1, uut1, vvt1, wwt1, p3t1 = load_winds(cubes)

        # Call fortran routine
        xx0, yy0, pp0, leftflag = pyLagranto.caltra.main(
            xx0, yy0, pp0, leftflag, ts, nsubs, imethod, numit, jflag, wfactor,
            fbflag, spt0, spt1, p3t0, p3t1, uut0, uut1, vvt0, vvt1, wwt0, wwt1,
            xmin, ymin, dx, dy, per, hem, nx, ny, nz, ntra)

        # Save positions
        traout[:, n, 0] = (time - times[0]).total_seconds()
        traout[:, n, 1] = xx0
        traout[:, n, 2] = yy0
        traout[:, n, 3] = pp0

        # Trace additional fields
        for m, tracer in enumerate(tracers):
            cube = convert.calc(tracer, cubes)
            cube = remap_3d(cube, example_cube)
            array = cube.data.transpose().flatten(order='F')
            traout[:, n, m + 4] = pyLagranto.trace.interp_to(
                array, xx0, yy0, pp0, leftflag, p3t1, spt1, xmin, ymin,
                dx, dy, nx, ny, nz, ntra)

    return traout


def grid_parameters(cube):
    """Extract grid parameters for calculations from cube

    args:
        cube (iris.cube.Cube):

    returns:
        nz, ny, nx (int): Grid dimensions.

        xmin, ymin (float): Minimum longitude and latitude.

        dx, dy (float): Grid spacing in degrees.

        hem, per (int): Flag for whether the domain is hemispheric and/or
            periodic.
    """
    # Extract grid dimesions
    nz, ny, nx = cube.shape
    x = grid.extract_dim_coord(cube, 'x').points
    y = grid.extract_dim_coord(cube, 'y').points

    # Find minimum latitude and longitude
    xmin = x.min()
    ymin = y.min()

    # Grid spacing
    dx = (x[1:] - x[:-1]).mean()
    dy = (y[1:] - y[:-1]).mean()

    # Set logical flag for periodic data set (hemispheric or not)
    hem = 0
    per = 0
    """
    if (per == 0):
        delta = xmax - xmin - 360.
        if (abs(delta + dx).lt.eps):
            raise Exception('arrays must be closed')
        elif (np.abs(delta) < eps):
            # Periodic and hemispheric
            hem = 1
            per = 360

    else:
        # Periodic and hemispheric
        hem = 1
    """

    return nx, ny, nz, xmin, ymin, dx, dy, hem, per


def load_winds(cubes):
    """Load the wind fields from a cubelist
    """
    # Extract fields needed as cubes
    u = convert.calc('x_wind', cubes)
    v = convert.calc('y_wind', cubes)
    w = convert.calc('upward_air_velocity', cubes)
    z = grid.make_cube(w, 'altitude')
    surface = convert.calc('surface_altitude', cubes)

    # Remap staggered grid cells
    u = remap_3d(u, w)
    v = remap_3d(v, w)

    # Return fields as 1d arrays with size nx*ny*nz
    return [x.data.transpose().flatten(order='F') for x in surface, u, v, w, z]


def remap_3d(cube, target):
    cube = cube.regrid(target, Linear())
    cube = cube.interpolate(
        [('level_height', target.coord('level_height').points)], Linear())

    return cube
