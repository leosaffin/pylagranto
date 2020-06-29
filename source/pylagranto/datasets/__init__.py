from abc import ABC

import numpy as np

import xarray as xr
import iris

from irise import convert, grid

from . import era5_heights_pressures as lagtraj


class DataSource():
    """Base class for loading data into the trajectory calculations

    Caltra requires two functions for getting informations from a data source.

    grid_parameters() -> nx, ny, nz, xmin, ymin, dx, dy, hem, per, names

    load_winds(time) -> z0, u, v, w, z
        Here z0 and z can represent any arbitrary vertical coordinate with w
        representing the vertical velocity in that coordinate system,

    """

    x_name = "longitude"
    y_name = "latitude"
    z_name = "altitude"
    surface_name = "surface_altitude"

    u_name = "eastward_wind"
    v_name = "northward_wind"
    w_name = "upward_air_velocity"

    def __init__(self, mapping):
        pass

    @property
    def shape(self):
        raise NotImplementedError

    def set_time(self, time):
        """Load the underlying data at the requested time

        Args:
            time (datetime.datetime):
        """
        raise NotImplementedError

    def get_variable(self, name):
        """

        Args:
            name (str):

        Returns:
            np.Array:
        """
        raise NotImplementedError

    def grid_parameters(self):
        """Extract grid parameters for calculations from cube

        Returns
            nz, ny, nx (int): Grid dimensions.

            xmin, ymin (float): Minimum longitude and latitude.

            dx, dy (float): Grid spacing in degrees.

            hem, per (int): Flag for whether the domain is hemispheric and/or
                periodic.

            names (list of str): The names of the dimensional coordinates
        """
        # Extract grid dimesions
        _, nz, ny, nx = self.shape

        # Find minimum latitude and longitude
        x = self.get_variable(self.x_name)
        y = self.get_variable(self.y_name)
        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()

        # Grid spacing
        # TODO check that the grid spacing is uniform
        dx = abs((x[1:] - x[:-1]).mean())
        dy = abs((y[1:] - y[:-1]).mean())

        # Set logical flag for periodic data set (hemispheric or not)
        if abs(xmax + dx - xmin - 360) < dx:
            hem = 1
            per = 360
        else:
            hem = 0
            per = 0

        names = [self.x_name, self.y_name, self.z_name]

        return nx, ny, nz, xmin, ymin, dx, dy, hem, per, names

    def winds(self):
        """Load the wind fields
        """
        u = self.get_variable(self.u_name)
        v = self.get_variable(self.v_name)
        w = self.get_variable(self.w_name)
        z = self.get_variable(self.z_name)
        surface = self.get_variable(self.surface_name)

        # Return fields as 1d arrays with size nx*ny*nz
        return [x.transpose().flatten(order='F') for x in (surface, u, v, w, z)]


class LagTraj(DataSource):
    x_name = "longitude"
    y_name = "latitude"
    z_name = "p_f"
    surface_name = "sp"

    u_name = "u"
    v_name = "v"
    w_name = "w"

    def __init__(self, filenames):
        self.filenames = filenames

    @property
    def shape(self):
        return self.ds.u.shape

    def set_time(self, time):
        ds = xr.open_dataset(self.filenames[time])
        ds = ds.drop_vars("z")
        ds = ds.sel(time=[time])

        ds0 = xr.open_dataset(self.filenames[time].replace("an_model", "an_single"))
        ds0 = ds0.sel(time=[time])
        ds0["longitude"] = ds["longitude"]

        ds = ds.merge(ds0)
        lagtraj.add_heights_and_pressures(ds)

        self.ds = lagtraj.era5_on_pressure_levels(ds, np.array([95000.]))
        self.ds = self.ds.reindex(latitude=list(reversed(ds.latitude)))

    def get_variable(self, name):
        if name == "p_f":
            return np.ones_like(self.ds.u.values)*95000
        elif name == "w":
            return np.zeros_like(self.ds.u.values)
        else:
            return self.ds[name].values
