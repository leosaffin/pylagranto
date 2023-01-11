import unittest
import datetime
import numpy as np
from pylagranto import trajectory


class TestTrajectoryEnsemble(unittest.TestCase):
    def setUp(self):
        data = np.array(
            [
                [[350., 45., 500.], [354., 45.5, 502.], [361., 46.5, 510]],
                [[355., 48., 550.], [359., 47.5, 511.], [363., 46., 487.]],
            ])
        t0 = datetime.datetime(2000, 1, 1)
        times = [t0 + datetime.timedelta(hours=n)
                 for n in np.linspace(0, 12, 3)]
        names = ['longitude', 'latitude', 'pressure']
        self.TrajectoryEnsemble1 = trajectory.TrajectoryEnsemble(
            data, times, names)

        times = [t0 - datetime.timedelta(hours=n)
                 for n in np.linspace(0, 12, 3)]
        self.TrajectoryEnsemble2 = trajectory.TrajectoryEnsemble(
            data, times, names)
        return

    def tearDown(self):
        return

    def test_add_same_fails(self):
        with self.assertRaises(ValueError):
            self.TrajectoryEnsemble1 + self.TrajectoryEnsemble1

    def test_add(self):
        result = self.TrajectoryEnsemble1 + self.TrajectoryEnsemble2
        print(result)


if __name__ == '__main__':
    unittest.main()
