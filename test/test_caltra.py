import pytest

import numpy as np

from pylagranto import caltra


@pytest.mark.parametrize("imethod", [1, 2])
def test_caltra(testdata, imethod):
    trainp = np.array([[180, 0, 5000]])

    # Calculate the trajectories
    traout = caltra.caltra(trainp, testdata, imethod=imethod)

    assert traout.x[0, 1] == pytest.approx(180.0323, abs=0.001)
    assert traout.y[0, 1] == pytest.approx(0.0323, abs=0.001)
    assert traout.z[0, 1] == pytest.approx(5036.0, abs=0.001)
