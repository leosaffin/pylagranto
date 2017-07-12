"""Functions for defining start points of trajectories
"""

import numpy as np


def circle(centre, radius, vertical_levels, resolution):
    """Set up an array of points bounded by a horizontal circle
    """
    # Set up a square a points
    points = np.arange(resolution, radius, resolution)

    # Select only points within a circle
    trainp = []
    cos_phi = np.cos(centre[1] * np.pi / 180)
    for zp in vertical_levels:
        trainp.append(centre[0], centre[1], zp)
        for xp in points:
            for yp in points:
                if ((cos_phi * xp)**2 + yp**2 < radius**2):
                    # Add a point for each quarter circle
                    trainp.append([centre[0] + xp, centre[1] + yp, zp])
                    trainp.append([centre[0] + xp, centre[1] - yp, zp])
                    trainp.append([centre[0] - xp, centre[1] + yp, zp])
                    trainp.append([centre[0] - xp, centre[1] - yp, zp])

    return np.array(trainp)
