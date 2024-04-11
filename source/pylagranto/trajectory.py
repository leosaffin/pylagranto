import datetime
import pickle
import numpy as np
import xarray as xr
from pylagranto import operator_dict


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

    @property
    def x(self):
        return self.data[:, 0]

    @property
    def y(self):
        return self.data[:, 1]

    @property
    def z(self):
        return self.data[:, 2]

    def __getitem__(self, key):
        if type(key) is str:
            index = self.names.index(key)

        return self.data[:, index]


class TrajectoryEnsemble(object):
    """A class for holding multiple trajectory objects

    Example:

    176 Trajectories

    Times
    ---------------------------
    2011-11-29 00:00    0:00:00
    2011-11-29 01:00    1:00:00
    2011-11-29 02:00    2:00:00
    2011-11-29 03:00    3:00:00

    Variables
    ---------------------------
    Grid Latitude
    Grid Longitude
    Altitude
    Air potential temperature


    Args:
        data (np.array): A 3d array with the data from multiple trajectories
            The dimensions are (Trajectory Number, Time, Variable).

        times (list of datetime.datetime):

        names (list of strings):
    """

    def __init__(self, data, times, names, misdat=-1000):
        self.data = np.ma.masked_where(data==misdat, data)
        self.times = times
        self.names = names

    @property
    def x(self):
        return self.data[:, :, 0]

    @property
    def y(self):
        return self.data[:, :, 1]

    @property
    def z(self):
        return self.data[:, :, 2]

    @property
    def relative_times(self):
        return [T - self.times[0] for T in self.times]

    def __len__(self):
        return len(self.data)

    def __str__(self, *args, **kwargs):
        string = ' ' + str(len(self)) + ' Trajectories \n \n '
        string += 'Times \n ' + '-' * 50 + ' \n '
        for time, dt in zip(self.times, self.relative_times):
            string += str(time) + '    ' + str(dt) + ' \n '
        string += '\n Variables \n ' + '-' * 50 + ' \n '
        for name in self.names:
            string += str(name).replace('_', ' ') + ' \n '
        return string

    def __getitem__(self, key):
        if type(key) is int:
            return Trajectory(self.data[key], self.times, self.names)
        elif type(key) is np.ndarray:
            return TrajectoryEnsemble(self.data[key], self.times, self.names)
        elif type(key) is str:
            index = self.names.index(key)
            return self.data[:, :, index]
        elif type(key) is datetime.datetime:
            index = self.times.index(key)
            return self.data[:, index, :]
        else:
            raise TypeError

    def __iter__(self):
        for n in range(len(self)):
            yield self[n]

    def __add__(self, other):
        if type(other) is TrajectoryEnsemble:
            # The addition of two trajectory ensembles is useful for stitching
            # together a set of forward and backward trajectories initialised
            # at the same point. In this case the zeroth time entry should be
            # identical but following times go in opposite directions. The
            # general use of 'add' should be to put together two
            # TrajectoryEnsembles with 0-N overlapping times at which all the
            # points are identical

            # No point adding a TrajectoryEnsemble to itself
            if self is other:
                raise ValueError('Cannot add a TrajectoryEnsemble to itself')

            # First check that the number of trajectories and variables along
            # each trajectory match.
            if len(self) != len(other) or self.names != other.names:
                raise ValueError('TrajectoryEnsemble shapes do not match')

            # Check that matching times contain identical points
            for time in self.times:
                if time in other.times:
                    if (self[time] != other[time]).all():
                        raise ValueError('Trajectories at time ' + str(time) +
                                         'do not match')

            # Determine the set of unique times in both trajectories
            times = list(set(self.times + other.times))
            times.sort()

            # Create and populate a data array for the combined
            # TrajectoryEnsemble
            new_array = np.zeros([len(self), len(times), len(self.names)])

            for n, time in enumerate(times):
                if time in self.times:
                    new_array[:, n, :] = self[time]
                else:
                    new_array[:, n, :] = other[time]

            return TrajectoryEnsemble(new_array, times, self.names)
        else:
            raise NotImplementedError('Cannot add ' + str(type(other)) +
                                      'to TrajectoryEnsemble')

    def __radd__(self, other):
        self.__add__(other)

    def save(self, filename):
        with open(filename, 'wb') as savefile:
            pickle.dump(self, savefile)

        return

    def append(self, name, data):
        """Add a new variable to the trajectory

        Args:
            name (str): Name of the variable

            data (np.array): Trajectory variable vs time with shape (ntraj, nt)
        """
        ntraj, nt, nvar = self.data.shape

        if data.shape != (ntraj, nt):
            raise ValueError('Shape of new data does not match trajectories')

        # Create a new data array with space for an extra variable
        newdata = np.zeros([ntraj, nt, nvar + 1])
        newdata[:, :, :-1] = self.data
        newdata[:, :, -1] = data

        # Overwrite existing values
        self.names.append(name)
        self.data = newdata

        return

    def select(self, variable, criteria, value, time=[]):
        """Select all trajectories where the variable matches the criteria

        e.g.

        >>> t0 = datetime.timedelta(hours=0)
        >>> t48 = datetime.timedelta(hours=48)
        >>> wcb_trajectories = trajectories.select('air_pressure', '>', 600,
        >>>                                        time=[t48, t0])

        Would select all trajectories that have decrease in pressure of more
        than 600 between 0-48 hours.

        Args:
            variable (str): The name of the variable in the trajectory to check
                the criteria against.

            criteria (str): A string representation of a python comparison
                operator (<, >, ==, !=, <=, >=).

            value (scalar): The value to compare the variable to.

            time (list): A list of zero, one or two datetime.timedelta objects
                representing the relative times to apply the criteria at. For
                zero times the criteria must be satisfied at all times. For one
                time the criteria must be satisfied at that time. For two times
                the criteria refers to the difference between those two times.
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

    def to_xarray(self):
        data = {
            name: (("trajectory", "time"), self.data[:, :, n])
            for n, name in enumerate(self.names)
        }

        return xr.Dataset(
            data, coords=dict(time=self.times, trajectory=range(len(self)))
        )


