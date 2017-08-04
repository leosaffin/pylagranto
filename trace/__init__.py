import numpy as np
import iris
from mymodule import convert
from lagranto import caltra, trajectory, pyLagranto


def trace(trajectories, tracers, mapping, levels=None):
    """Trace fields using existing trajectories

    Args:
        trajectories (trajectory.TrajectoryEnsemble): Output object of
            trajectory calculations

        tracers (list): Names of new variables to be traced.

        mapping (dict): Mapping of times to filenames with data.

        levels (tuple or None): The name of a vertical coordinate used in the
            trajectory calculations and list of values to interpolate fields to
            prior to performing interpolating to tracer positions.

    Returns:
        traout (trajectory.TrajectoryEnsemble): Same as trajectories input but
            with tracers appended to the data
    """

    # Extract trajectory data
    times = trajectories.times
    names = trajectories.names

    # Only trace new variables
    new_tracers = []
    for tracer in tracers:
        if tracer not in names:
            new_tracers.append(tracer)
            names.append(tracer)

    tracers = new_tracers
    if len(tracers == 0):
        raise ValueError('No new variables to tracer')
    else:
        print('Tracing', tracers)

    # Increase the size of the array holding trajectory data
    ntra, ntim, nvar = trajectories.data.shape
    traout = np.zeros([ntra, ntim, nvar + len(tracers)])

    # Copy old data
    traout[:, :, 0:nvar] = trajectories.data

    for n, time in enumerate(times):
        print time

        # Extract trajectory positions
        x = traout[:, n, 0]
        y = traout[:, n, 1]
        z = traout[:, n, 2]

        # Check which trajectories have left the domain
        leftflag = (traout[:, n, 2] == -1000.).astype(int)

        # Read wind fields and surface pressure at next time
        cubes = iris.load(mapping[time])
        spt1, uut1, vvt1, wwt1, p3t1 = caltra.load_winds(cubes, levels)

        if n == 0:
            # Load grid parameters
            example_cube = convert.calc('upward_air_velocity', cubes,
                                        levels=levels)
            nx, ny, nz, xmin, ymin, dx, dy, hem, per, varnames = \
                caltra.grid_parameters(example_cube, levels)

        # Add each tracer to the trajectory data
        for m, tracer in enumerate(tracers, start=nvar):
            try:
                cube = convert.calc(tracer, cubes, levels=levels)
                array = cube.data.transpose().flatten(order='F')
            except ValueError:
                # If variable can't be loaded print a warning and put zero
                print ('Variable ' + tracer + ' not available at this time. ' +
                       'Replacing with zeros')
                array = np.zeros_like(uut1)
            traout[:, n, m] = pyLagranto.trace.interp_to(
                array, x, y, z, leftflag, p3t1, spt1, xmin, ymin,
                dx, dy, nx, ny, nz, ntra)

    return trajectory.TrajectoryEnsemble(traout, times, names)
