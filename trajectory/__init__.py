import numpy as np
from pandas import DataFrame, Panel
from datetime import timedelta as dt
from lagranto import operator_dict


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
        data (np.array): A 3d array with the data from multiple trajectories

        times (list of datetime.datetime):

        names (list of strings):
    """

    def __init__(self, data, times, names):
        Panel.__init__(self, data=data, items=range(len(data)),
                       major_axis=times, minor_axis=names)

    @property
    def times(self):
        return list(self.major_axis)

    @property
    def relative_times(self):
        return [T - self.times[0] for T in self.times]

    @property
    def names(self):
        return list(self.minor_axis)

    def __iter__(self):
        for index in self.items:
            yield self[index]

    def select(self, variable, criteria, value, time):
        """Select all trajectories where the variable matches the criteria
        """
        # Convert the criteria to the relevant operator
        if type(criteria) is str:
            criteria = operator_dict[criteria]

        # Extract the indices for the data to be checked
        var_index = self.names.index(variable)
        time_index = self.relative_times.index(time)

        # Get the indices of the trajectories that match the criteria
        indices = np.where(
            criteria(self.values[:, time_index, var_index], value))

        # Create a new trajectory ensemble with the subset of trajectories
        subset = TrajectoryEnsemble(
            self.values[indices], self.times, self.names)

        return subset
