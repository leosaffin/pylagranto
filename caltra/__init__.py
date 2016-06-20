import numpy as np
from iris.analysis import Linear
from mymodule import convert, files, grid
from lagranto import pyLagranto
from lagranto.trajectory import TrajectoryEnsemble


def caltra(trainp, mapping, imethod=1, numit=3, nsubs=4, fbflag=1, jflag=False,
           tracers=[]):
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

        tracers (list): A list of variable names to trace at each point along
            the trajectory
    returns:
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

    # Save starting positions
    traout[:, 0, 0] = x
    traout[:, 0, 1] = y
    traout[:, 0, 2] = z

    # Extract times relating to filenames
    times = sorted(list(mapping))

    # Calulate the timestep in seconds between input files and divide by
    # number of substeps
    ts = (times[1] - times[0]).total_seconds() / nsubs

    # Reverse file load order for reverse trajectories
    if (fbflag == -1):
        times.reverse()

    # Read wind field and grid from first file
    cubes = files.load(mapping[times[0]])
    spt1, uut1, vvt1, wwt1, p3t1 = load_winds(cubes)

    example_cube = convert.calc('upward_air_velocity', cubes)
    nx, ny, nz, xmin, ymin, dx, dy, hem, per, names = \
        grid_parameters(example_cube)

    # Add the list of traced variables to the names in the output
    names += tracers

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

        # Call fortran routine to update trajectory positions
        x, y, z, leftflag = pyLagranto.caltra.main(
            x, y, z, leftflag, ts, nsubs, imethod, numit, jflag,
            fbflag, spt0, spt1, p3t0, p3t1, uut0, uut1, vvt0, vvt1, wwt0, wwt1,
            xmin, ymin, dx, dy, per, hem, nx, ny, nz, ntra)

        # Save positions
        traout[:, n, 0] = x
        traout[:, n, 1] = y
        traout[:, n, 2] = z

        # Trace additional fields
        for m, tracer in enumerate(tracers):
            cube = convert.calc(tracer, cubes)
            cube = remap_3d(cube, example_cube)
            array = cube.data.transpose().flatten(order='F')
            traout[:, n, m + 3] = pyLagranto.trace.interp_to(
                array, x, y, z, leftflag, p3t1, spt1, xmin, ymin,
                dx, dy, nx, ny, nz, ntra)

    return TrajectoryEnsemble(traout, times, names)


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

        names (list of str): The names of the dimensional coordinates
    """
    # Extract grid dimesions
    nz, ny, nx = cube.shape
    x = grid.extract_dim_coord(cube, 'x')
    y = grid.extract_dim_coord(cube, 'y')
    names = [x.name(), y.name(), 'altitude']

    # Find minimum latitude and longitude
    xmin = x.points.min()
    ymin = y.points.min()

    # Grid spacing
    dx = (x.points[1:] - x.points[:-1]).mean()
    dy = (y.points[1:] - y.points[:-1]).mean()

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

    return nx, ny, nz, xmin, ymin, dx, dy, hem, per, names


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
