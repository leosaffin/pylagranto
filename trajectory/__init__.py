import cPickle as pickle
import numpy as np
from datetime import timedelta as dt
from lagranto import operator_dict


def load(filename):
    with open(filename, 'rb') as savefile:
        data = pickle.load(savefile)

    return data


class Trajectory(object):
    """A class for holding data for a trajectory

    Example:

    time                lon       lat       p
    -------------------------------------------
    2011-11-29 22:00    37.249    24.455    500
    2011-11-29 23:00    37.249    24.455    500

    This would represent a numpy.ndarray of shape (2,4).
    The times are stored as differences from a start time.
    """

    def __init__(self, data, times, names):
        self.data = data
        self.times = times
        self.names = names

    def __getitem__(self, key):
        if type(key) is str:
            index = self.names.index(key)

        return self.data[:, index]


class TrajectoryEnsemble(object):
    """A class for holding multiple trajectory objects

    Args:
        data (np.array): A 3d array with the data from multiple trajectories

        times (list of datetime.datetime):

        names (list of strings):
    """

    def __init__(self, data, times, names):
        self.data = data
        self.times = times
        self.names = names

    @property
    def relative_times(self):
        return [T - self.times[0] for T in self.times]

    def __len__(self):
        return len(self.data)

    def __getitem__(self, key):
        if type(key) is int:
            return Trajectory(self.data[key], self.times, self.names)
        elif type(key) is str:
            index = self.names.index(key)
            return self.data[:, :, index]
        else:
            raise TypeError

    def __iter__(self):
        for n in range(len(self)):
            yield self[n]

    def save(self, filename):
        with open(filename, 'wb') as savefile:
            pickle.dump(self, savefile)

        return

    def select(self, variable, criteria, value, time=[]):
        """Select all trajectories where the variable matches the criteria
        """
        # Convert the criteria to the relevant operator
        if type(criteria) is str:
            criteria = operator_dict[criteria]

        # Extract the indices for the data to be checked
        var_index = self.names.index(variable)

        if len(time) == 0:
            # Take trajectories that always match the criteria
            indices = np.where(
                criteria(self.data[:, :, var_index], value).all(axis=1))

        elif len(time) == 1:
            # Take trajectories that match the criteria at the given time
            time_index = self.relative_times.index(time[0])

            # Get the indices of the trajectories that match the criteria
            indices = np.where(
                criteria(self.data[:, time_index, var_index], value))

        elif len(time) == 2:
            # Take trajectories where the difference between the two times
            # matches the criteria
            time_index_1 = self.relative_times.index(time[0])
            time_index_2 = self.relative_times.index(time[1])

            diff = (self.data[:, time_index_1, var_index] -
                    self.data[:, time_index_2, var_index])

            indices = criteria(diff, value)

        # Create a new trajectory ensemble with the subset of trajectories
        subset = TrajectoryEnsemble(
            self.data[indices], self.times, self.names)

        return subset
