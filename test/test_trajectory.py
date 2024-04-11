import pytest
import xarray


def test_add(trajectory_ensemble_forward, trajectory_ensemble_backward):
    result = trajectory_ensemble_forward + trajectory_ensemble_backward

    assert len(result) == len(trajectory_ensemble_forward)
    assert (
        len(result.times)
        == len(trajectory_ensemble_forward.times)
        + len(trajectory_ensemble_backward.times)
        - 1
    )


def test_add_same_fails(trajectory_ensemble_forward):
    with pytest.raises(ValueError):
        trajectory_ensemble_forward + trajectory_ensemble_forward


def test_to_xarray(trajectory_ensemble_forward):
    ds = trajectory_ensemble_forward.to_xarray()
    assert isinstance(ds, xarray.Dataset)
    assert len(ds) == len(trajectory_ensemble_forward.names)
    assert len(ds.time) == len(trajectory_ensemble_forward.times)
    assert len(ds.trajectory) == len(trajectory_ensemble_forward)
