"""Functions for defining start points of trajectories
"""

import numpy as np
import pandas as pd
from mymodule import convert, grid
from lagranto import operator_dict


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
                if ((cos_phi * xp)**2 + yp**2 < radius**2):
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
        for angle in range(0, 360, resolution):
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
