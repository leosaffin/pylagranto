import pytest
import iris
from irise import convert

from pylagranto.caltra import load_winds, grid_parameters
from pylagranto.fortran import caltra


@pytest.mark.parametrize(
    "z,plevs,leftflag",
    [
        # Height coordinate in bounds
        (1000, False, 0),
        # Height coordinate below surface
        (-1, False, 1),
        # Height coordinate above top
        (1e10, False, 1),
        # Pressure coordinate in bounds
        (1000, True, 0),
        # Pressure coordinate below surface
        (1e6, True, 1),
        # Pressure coordinate above top
        (1, True, 1),
    ]
)
def test_check_boundaries(testdata, z, plevs, leftflag):
    levels = None
    cubes = iris.load(testdata[list(testdata.keys())[0]])
    spt1, uut1, vvt1, wwt1, p3t1 = load_winds(cubes, levels)
    example_cube = convert.calc('upward_air_velocity', cubes, levels=levels)
    nx, ny, nz, xmin, ymin, dx, dy, hem, per, names = grid_parameters(example_cube, levels)
    reltpos = 0.0
    jump = 0

    # z-coordinate decreasing with height (like pressure)
    if plevs:
        p3t1 = p3t1[:, :, ::-1]
        spt1[:] = 1e5

    left = caltra.check_boundaries(
        45, 45, z, reltpos, jump, spt1, spt1, p3t1, p3t1, xmin, ymin, dx, dy, per, hem
    )
    assert left == leftflag
