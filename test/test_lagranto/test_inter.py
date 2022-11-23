import pytest
import numpy as np

from pylagranto.fortran import inter


@pytest.mark.parametrize("mode", [1, 2, 3])
@pytest.mark.parametrize(
    "position,indices",
    [
        ([1.5, 20.0, 900.0], [2, 3, 4]),
        ([1.5, 20.0, 300.0], [2, 3, 0.5]),
        ([1.5, 20.0, 1000.0], [2, 3, 0]),
    ]
)
def test_get_index3(mode, position, indices):
    nx, ny, nz = 3, 4, 5

    lats = 15.0
    lonw = 0.0
    dlon = 1.5
    dlat = 2.5

    xpo, ypo, ppo = position

    surf = np.zeros([nx, ny], dtype=float)
    vert = np.linspace(600, 1000, 5, dtype=float)
    vert = np.broadcast_to(vert, [nx, ny, nz])

    rid, rjd, rkd = inter.get_index3(
        xpo, ypo, ppo, mode, vert, surf, lonw, lats, dlon, dlat
    )

    assert rid == indices[0]
    assert rjd == indices[1]
    assert rkd == indices[2]


