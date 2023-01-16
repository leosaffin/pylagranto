import pytest

import numpy as np

from pylagranto import caltra, datasets


def test_caltra(testdata):
    # Set up the mapping of times to files:
    datasource = datasets.MetUMStaggeredGrid(testdata)

    trainp = np.array([[180, 0, 5000]])

    print(trainp)

    # Calculate the trajectories
    traout = caltra.caltra(trainp, list(testdata), datasource)

    assert traout.x[0, 1] == pytest.approx(180.0323, abs=0.001)
    assert traout.y[0, 1] == pytest.approx(0.0323, abs=0.001)
    assert traout.z[0, 1] == pytest.approx(5036.0, abs=0.001)
