import numpy as np
import iris

from irise import convert, grid, interpolate, variable


class DataSource:
    """Base class for loading data into the trajectory calculations

    Caltra requires two functions for getting information from a data source.

    grid_parameters() -> nx, ny, nz, xmin, ymin, dx, dy, hem, per, names

    load_winds(time) -> z0, u, v, w, z
        Here z0 and z can represent any arbitrary vertical coordinate with w
        representing the vertical velocity in that coordinate system,

    """

    example_cube = None
    x_name = "longitude"
    y_name = "latitude"
    z_name = "altitude"
    surface_name = "surface_altitude"

    u_name = "eastward_wind"
    v_name = "northward_wind"
    w_name = "upward_air_velocity"

    def __init__(self, mapping, levels=None):
        self.mapping = mapping
        self.data = None
        self.levels = levels

    @property
    def shape(self):
        return self.example_cube.shape

    def set_time(self, time):
        self.data = iris.load(
            self.mapping[time], iris.Constraint(time=lambda x: x.point == time)
        )
        self.example_cube = convert.calc(self.w_name, self.data, levels=self.levels)

    def get_variable(self, name):
        """
        Args:
            name (str):

        Returns:
            np.Array:
        """
        return convert.calc(name, self.data).data

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
        nz, ny, nx = self.shape
        x = self.example_cube.coord(axis="x", dim_coords=True)
        y = self.example_cube.coord(axis="y", dim_coords=True)

        if self.levels is None:
            names = [x.name(), y.name(), self.z_name]
        else:
            names = [x.name(), y.name(), self.levels[0]]

        # Find minimum latitude and longitude
        xmin = x.points.min()
        xmax = x.points.max()
        ymin = y.points.min()
        ymax = y.points.max()

        # Grid spacing
        dx = (x.points[1:] - x.points[:-1]).mean()
        dy = (y.points[1:] - y.points[:-1]).mean()

        # Set logical flag for periodic data set (hemispheric or not)
        if abs(xmax + dx - xmin - 360) < dx:
            hem = 1
            per = 360
        else:
            hem = 0
            per = 0

        return nx, ny, nz, xmin, ymin, dx, dy, hem, per, names

    def winds(self):
        """Load the wind fields
        """
        u = self.get_variable(self.u_name)
        v = self.get_variable(self.v_name)
        w = self.get_variable(self.w_name)
        z = self.get_variable(self.z_name)
        surface = self.get_variable(self.surface_name)

        return [surface, u, v, w, z]


class MetUMStaggeredGrid(DataSource):

    vert_coord = "atmosphere_hybrid_height_coordinate"

    u_name = "x_wind"
    v_name = "y_wind"

    def get_variable(self, name):
        # Return fields as 1d arrays with size nx*ny*nz
        return self._get_variable(name).transpose()

    def _get_variable(self, name):
        if name == self.surface_name:
            if self.levels is None:
                # Extract the surface coordinate
                surface = self.example_cube.coord(self.surface_name).points
            else:
                # Set the surface coordinate to zeros
                surface = np.zeros_like(self.example_cube[0].data)
            return surface

        if name == self.z_name:
            if self.levels is None:
                # Default is height based coordinate
                z = variable.height(self.example_cube).data
            else:
                z = np.zeros_like(self.example_cube.data)
                # Create a uniform height coordinate for the levels interpolated to
                for n, level in enumerate(self.levels[1]):
                    z[n, :, :] = level

            return z

        if name == self.w_name and self.levels is not None:
            # Set vertical velocity to zero
            return np.zeros_like(self.example_cube.data)

        cube = convert.calc(name, self.data, levels=self.levels)

        if name in [self.u_name, self.v_name]:
            if self.levels is None:
                cube = interpolate.remap_3d(cube, self.example_cube, vert_coord=self.vert_coord)
            else:
                cube = cube.regrid(self.example_cube, iris.analysis.Linear())

        return cube.data


class ERA5(DataSource):
    z_name = "air_pressure"
    surface_name = "surface_air_pressure"

    w_name = "lagrangian_tendency_of_air_pressure"

    def set_time(self, time):
        super().set_time(time)
        # Pressure is stored as a 1d coordinate "pressure_level" and in hPa while
        # other variables (e.g. omega velocity) use Pa
        p = grid.make_cube(self.example_cube, "pressure_level")
        p.convert_units("Pa")
        p = grid.broadcast_to_cube(p, self.example_cube)
        p.rename("air_pressure")
        self.data.append(p)

        # Vertical coordinate and latitude are reversed compared to what Lagranto
        # expects
        self.example_cube = self.example_cube[::-1, ::-1, :]

    def get_variable(self, name):
        if name == self.surface_name:
            return super().get_variable(name)[::-1, :].transpose()
        else:
            return super().get_variable(name)[::-1, ::-1, :].transpose()
