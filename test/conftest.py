import pathlib
import datetime

import pytest
import numpy as np

import iris
import iris.cube
import iris.coords
import iris.aux_factory


# A testdata folder in this directory
testdata_dir = pathlib.Path(__file__).parent / "testdata"
t0 = datetime.datetime(2000, 1, 1)
dt = datetime.timedelta(hours=1)


def generate_test_netcdf():

    nz, ny, nx = 360, 181, 30

    # Create a 1 degree grid
    longitude = iris.coords.DimCoord(
        np.linspace(0, 359, nx),
        standard_name="longitude",
        units="degrees",
        circular=True,
    )
    latitude = iris.coords.DimCoord(
        np.linspace(-90, 90, ny),
        standard_name="latitude",
        units="degrees",
    )
    level_height = iris.coords.DimCoord(
        np.linspace(500, 15000, nz),
        long_name="level_height",
        units="m",
    )

    # Not sure how this is typically defined but that's not too important. It just needs
    # to be a function between 1 and 0 that tapers off with height
    sigma = iris.coords.AuxCoord(
        1 - (level_height.points / level_height.points[-1])**2,
        standard_name="atmosphere_sigma_coordinate",
        units="",

    )
    surface_altitude = iris.coords.AuxCoord(
        np.zeros([ny, nx]),
        standard_name="surface_altitude",
        units="m",
    )

    altitude = iris.aux_factory.HybridHeightFactory(
        delta=level_height,
        sigma=sigma,
        orography=surface_altitude
    )

    # Generate idealised wind fields
    u = np.ones([nz, ny, nx])
    v = np.ones([nz, ny, nx])
    w = np.ones([nz, ny, nx]) * 0.01

    for n in range(2):
        time = iris.coords.AuxCoord(n, long_name="time", units="hours since 2000-01-01",)
        cubes = iris.cube.CubeList()
        for data, name in [(u, "x_wind"), (v, "y_wind"), (w, "upward_air_velocity")]:
            cube = iris.cube.Cube(
                data,
                standard_name=name,
                units="m s-1",
                dim_coords_and_dims=[
                    (level_height, 0),
                    (latitude, 1),
                    (longitude, 2),
                ],
                aux_coords_and_dims=[
                    (time, None),
                    (sigma, 0),
                    (surface_altitude, [1, 2]),
                ],
                aux_factories=[altitude],
            )
            cubes.append(cube)

        iris.save(cubes, str(testdata_dir/"testdata_{}.nc".format(n)))


@pytest.fixture
def testdata(scope="session"):
    if not testdata_dir.exists():
        testdata_dir.mkdir()
        generate_test_netcdf()

    return {t0 + n*dt: testdata_dir/"testdata_{}.nc".format(n) for n in range(2)}
