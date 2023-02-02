import pytest

import datetime

import numpy as np

from pylagranto import caltra


@pytest.mark.parametrize("imethod", [1, 2])
def test_caltra(testdata, imethod):
    # Set up the mapping of times to files:
    start_time = datetime.datetime(2020, 1, 1)
    dt = datetime.timedelta(hours=1)
    mapping = {start_time + n*dt: str(testdata) for n in range(2)}

    trainp = np.array([[180, 0, 5000]])

    print(trainp)

    # Calculate the trajectories
    traout = caltra.caltra(trainp, mapping, imethod=imethod)

    assert traout.x[0, 1] == pytest.approx(180.0323, abs=0.001)
    assert traout.y[0, 1] == pytest.approx(0.0323, abs=0.001)
    assert traout.z[0, 1] == pytest.approx(5036.0, abs=0.001)
