import pytest
import numpy as np

from pylagranto.fortran import inter


@pytest.mark.parametrize("mode", [1, 2, 3])
@pytest.mark.parametrize(
    "position,indices",
    [
        ([1.5, 20.0, 900.0], [2, 3, 4.5]),
        ([1.5, 20.0, 300.0], [2, 3, 0.5]),
        ([1.5, 20.0, 1000.0], [2, 3, 5]),
    ]
)
@pytest.mark.parametrize("inverted", [True, False])
def test_get_index3(mode, position, indices, inverted):
    nx, ny, nz = 3, 4, 5

    lats = 15.0
    lonw = 0.0
    dlon = 1.5
    dlat = 2.5

    xpo, ypo, ppo = position

    surf = np.zeros([nx, ny], dtype=float)
    vert = np.array([600, 650, 700, 800, 1000], dtype=float)

    if inverted:
        surf[:, :] = 1050.
        vert = vert[::-1]
        indices = indices.copy()
        indices[-1] = 6 - indices[-1]

    vert = np.broadcast_to(vert, [nx, ny, nz])

    rid, rjd, rkd = inter.get_index3(
        xpo, ypo, ppo, mode, vert, surf, lonw, lats, dlon, dlat
    )

    assert rid == indices[0]
    assert rjd == indices[1]
    assert rkd == indices[2]


