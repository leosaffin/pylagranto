from pandas import DataFrame, Panel
from datetime import timedelta as dt


class Trajectory(DataFrame):
    """A class for holding data for a trajectory

    Example:

    time                lon       lat       p
    -------------------------------------------
    2011-11-29 22:00    37.249    24.455    500
    2011-11-29 23:00    37.249    24.455    500

    This would represent a numpy.ndarray of shape (2,4).
    The times are stored as differences from a start time.
    """

    def __init__(self, data, varnames, times):
        DataFrame.__init__(self, data, index=times, columns=varnames)


class TrajectoryEnsemble(Panel):
    """A class for holding multiple trajectory objects

    Args:
        trajectories (list): A list of Trajectory objects
    """

    def __init__(self, trajectories):
        data = {n: trajectory for n, trajectory in enumerate(trajectories)}
        Panel.__init__(self, data)
