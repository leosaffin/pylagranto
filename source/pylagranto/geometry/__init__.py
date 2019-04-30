"""Functions for defining start points of trajectories
"""

import numpy as np
import pandas as pd
import iris.plot as iplt
from irise import convert, grid
from mymodule.detection.rossby_waves import tropopause_contour as trop
from pylagranto import caltra, operator_dict
import pylagranto.fortran


def circle(centre, radius, vertical_levels, resolution):
    """Set up an array of points bounded by a horizontal circle

    Args:
        centre (tuple of floats): The lon/lat positions of the circle centre.

        radius (float): The radius of the circle in degrees.

        vertical_levels (list): A list of vertical level values to insert in
            the z position for each point.

        resolution (float): The spacing between points in degrees.

    Returns:
        trainp (np.array): Array of shape (ntra, 3) with all the points in the
            circle.
    """
    # Set up a square a points
    points = np.arange(resolution, radius, resolution)

    # Select only points within a circle
    trainp = []
    cos_phi = np.cos(centre[1] * np.pi / 180)
    for zp in vertical_levels:
        trainp.append([centre[0], centre[1], zp])
        for xp in points:
            # Take all points along zero line
            trainp.append([centre[0] + xp, centre[1], zp])
            trainp.append([centre[0] - xp, centre[1], zp])
            trainp.append([centre[0], centre[1] + xp, zp])
            trainp.append([centre[0], centre[1] - xp, zp])
            for yp in points:
                if (cos_phi * xp)**2 + yp**2 < radius**2:
                    # Add a point for each quarter circle
                    trainp.append([centre[0] + xp, centre[1] + yp, zp])
                    trainp.append([centre[0] + xp, centre[1] - yp, zp])
                    trainp.append([centre[0] - xp, centre[1] + yp, zp])
                    trainp.append([centre[0] - xp, centre[1] - yp, zp])

    return np.array(trainp)


def ring(centre, radius, vertical_levels, resolution):
    """Set up an array of points along a horizontal circle

    Args:
        centre (tuple of floats): The lon/lat positions of the circle centre.

        radius (float): The radius of the circle in degrees.

        vertical_levels (list): A list of vertical level values to insert in
            the z position for each point.

        resolution (float): The angular spacing between points in degrees.

    Returns:
        trainp (np.array): Array of shape (ntra, 3) with all the points in the
            circle.
    """
    trainp = []
    for zp in vertical_levels:
        # Select all points around the edge of a circle
        for angle in np.arange(0, 360, resolution):
            xp = centre[0] + radius * np.cos(angle * np.pi / 180)
            yp = centre[1] + radius * np.sin(angle * np.pi / 180)
            trainp.append([xp, yp, zp])

    return np.array(trainp)


def select(cubes, variable, criteria, value, levels=None):
    """Select start points where the variable satisfies a given criteria

    e.g.
    >>> select(cubes, 'altitude', '>', 10000,
               levels=('air_potential_temperature', [320]))

    Will select all gridpoints on the 320 K isentrope above 10 km.

    Args:
        cubes (iris.cube.CubeList): A cubelist containing the data required
            to calculate start positions.

        variable (str): The variable to calculate the criteria on.

        criteria (str): A string defining a comparison operator, e.g. 'lt' and
            '<' are equivalent. Can also define the operator using the built
            in python definitions.

        value (float): The value to use with the comparison criteria.

        levels (optional tuple): A set of vertical levels to calculate the
            start points at.

    Returns:
        trainp (np.array): Array of shape (ntra, 3) with all the grid points
            satisfying the give criteria.
    """
    # Convert the criteria to the relevant operator
    if type(criteria) is str:
        criteria = operator_dict[criteria]

    # Extract the variable
    cube = convert.calc(variable, cubes, levels=levels)

    # Get the corresponding grid data
    nz, ny, nx = cube.shape
    lon, lat = grid.get_xy_grids(cube)
    z = convert.calc('altitude', cubes, levels=levels).data

    # Set up an input array of trajectories
    trainp = []

    # Select all grid points that satisfy the criteria
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if criteria(cube.data[k, j, i], value):
                    trainp.append([lon[j, i], lat[j, i], z[k, j, i]])

    return np.array(trainp)


def contour(cubes, varname, value, levels=None):
    cube = convert.calc(varname, cubes, levels=levels)

    trainp = []
    for n, level in enumerate(levels[1]):
        cs = iplt.contour(cube[n], [value])
        contours = trop.get_contour_verts(cs)[0]
        path = trop.get_tropopause_contour(contours)

        for x, y in path.vertices():
            trainp.append([x, y, level])

    return np.array(trainp)


def select_2d(cubes, variable, criteria, value, levels):
    """Select start points where the variable satisfies a given criteria

    Same as select but for point on a2d cube such as sea-level pressure

    Args:
        cubes (iris.cube.CubeList): A cubelist containing the data required
            to calculate start positions.

        variable (str): The variable to calculate the criteria on.

        criteria (str): A string defining a comparison operator, e.g. 'lt' and
            '<' are equivalent. Can also define the operator using the built
            in python definitions.

        value (float): The value to use with the comparison criteria.

        levels (np.array): Vertical levels for the start points.

    Returns:
        trainp (np.array): Array of shape (ntra, 3) with all the grid points
            satisfying the give criteria.
    """
    # Convert the criteria to the relevant operator
    if type(criteria) is str:
        criteria = operator_dict[criteria]

    # Extract the variable
    cube = convert.calc(variable, cubes)

    # Get the corresponding grid data
    nz = len(levels)
    ny, nx = cube.shape
    lon, lat = grid.get_xy_grids(cube)

    # Set up an input array of trajectories
    trainp = []

    # Select all grid points that satisfy the criteria
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if criteria(cube.data[j, i], value):
                    trainp.append([lon[j, i], lat[j, i], levels[k]])

    return np.array(trainp)


def select_from_trajectories(filename):
    """Create start points from the end of existing trajectories

    Args:
        filename (str): The '.pkl' filename of the saved trajectories

    Returns:
        trainp (np.array): Array of shape (ntra, 3) with the end points of the
            trajectories.
    """
    # Load the trajectories
    trajectories = pd.read_pickle(filename)

    trainp = []
    for T in trajectories:
        # Only take trajectories that stay in the domain
        if (T.altitude.values > 0).all():
            trainp.append([T.grid_longitude[-1], T.grid_latitude[-1],
                           T.altitude[-1]])

    return np.array(trainp)


def replace_z(trainp, cubes, levels, name='altitude'):
    """Change the z coordinate on the trajectory start positions

    Example:

    >>> trainp = geometry.circle([0,0], 1, [310, 320], [0.1])

    Gives a set of points inside a 1 degree radius circle with the vertical
    coordinates set to 310 and 320. To change these vertical points from
    isentropic to altitude based,

    >>> levels = ('air_potential_temperature', range(300, 330, 2))
    >>> trainp = replace_z(trainp, cubes, levels, name='altitude')

    Or pressure based

    >>> levels = ('air_potential_temperature', range(300, 330, 2))
    >>> trainp = replace_z(trainp, cubes, levels, name='air_pressure')

    etc.


    Args:
        trainp (np.array): Array of shape (ntra, 3) with the initial trajectory
            start points.

        cubes (iris.cube.CubeList): A cubelist containing the data required
            to calculate start positions.

        levels (tuple): A set of vertical levels of the original
            coordinate to calculate the new coordinate from.

        name (str): The name of the new coordinate

    Returns
        new_trainp (np.array): The original trainp array with the z points
            replaced.
    """
    # Load new vertical coordinate
    z = convert.calc(name, cubes, levels=levels)
    array = z.data.transpose().flatten(order='F')

    # Extract grid parameters
    nx, ny, nz, xmin, ymin, dx, dy, hem, per, names = caltra.grid_parameters(
        z, levels)

    # Load surface fields
    spt1, uut1, vvt1, wwt1, p3t1 = caltra.load_winds(cubes, levels)
    zp = pylagranto.fortan.trace.interp_to(
        array, trainp[:, 0], trainp[:, 1], trainp[:, 2],
        np.zeros_like(trainp[:, 0]), p3t1, spt1, xmin, ymin, dx, dy,
        nx, ny, nz, len(trainp))

    new_trainp = trainp.copy()
    new_trainp[:, 2] = zp

    return new_trainp
