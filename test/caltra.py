import pytest

import numpy as np

import iris

from pylagranto.fortran import caltra


def test_caltra():
    x, y, z = 302.33, 13.33, 95000.
    ntra = 1
    leftflag = np.zeros(1)
    ts = 900.
    nsubs = 4
    imethod = 1
    numit = 3
    jflag = False
    fbflag = -1

    hem = 0
    per = 0

    cubes0 = iris.load("ERA5_950hPa_2020-02-05T15:00.nc")
    cubes1 = iris.load("ERA5_950hPa_2020-02-05T14:00.nc")

    xmin, ymin, dx, dy, nx, ny, nz = grid_parameters(cubes0)
    spt0, p3t0, uut0, vvt0, wwt0 = load_winds(cubes0)
    spt1, p3t1, uut1, vvt1, wwt1 = load_winds(cubes1)

    x, y, z, leftflag = caltra.main(
                    x, y, z, leftflag, ts, nsubs, imethod, numit, jflag, fbflag,
                    spt0, spt1, p3t0, p3t1, uut0, uut1, vvt0, vvt1, wwt0, wwt1,
                    xmin, ymin, dx, dy, per, hem, nx, ny, nz, ntra)

    assert x == 302.69702148
    assert y == 13.37627983


def grid_parameters(cubes):
    cube = cubes.extract_strict("northward_wind")

    x = cube.coord("longitude").points
    y = cube.coord("latitude").points

    # Find minimum latitude and longitude
    xmin = x.min()
    ymin = y.min()

    # Grid spacing
    dx = abs((x[1:] - x[:-1]).mean())
    dy = abs((y[1:] - y[:-1]).mean())

    # Extract grid dimesions
    _, nz, ny, nx = cube.shape

    return xmin, ymin, dx, dy, nx, ny, nz


def load_winds(cubes):
    sp = cubes.extract_strict("surface_air_pressure")

    u = cubes.extract_strict("eastward_wind")
    v = cubes.extract_strict("northward_wind")
    w = v.copy(data=np.zeros_like(v.data))
    p = v.copy(data=np.ones_like(v.data)*95000.)

    return [cube.data.transpose().flatten(order="F").data.astype(np.double) for cube in [sp, p, u, v, w]]


if __name__=="__main__":
    test_caltra()