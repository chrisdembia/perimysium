import numpy as np
from numpy import testing

from epimysium import postprocessing as pproc

def test_shift_data_to_cycle():

    ordinate = np.array([2, 1, 2, 3, 4, 5, 6])
    time = np.array([0.5, 1.0, 1.2, 1.35, 1.4, 1.5, 1.8])
    arbitrary_cycle_start_time = 1.0
    arbitrary_cycle_end_time = 1.5
    new_cycle_start_time = 1.35
    shifted_time, shifted_ordinate = pproc.shift_data_to_cycle(
        arbitrary_cycle_start_time, arbitrary_cycle_end_time,
        new_cycle_start_time,
        time, ordinate)

    testing.assert_allclose(shifted_ordinate, np.array([3, 4, 5, 1, 2]))
    testing.assert_allclose(shifted_time,
            np.array([0.0, 0.05, 0.15, 0.30, 0.5]))

    testing.assert_equal(shifted_time.size, 5)
    testing.assert_equal(shifted_ordinate.size, 5)

    assert shifted_time[0] == 0.0
    assert shifted_time[-1] == (arbitrary_cycle_end_time -
            arbitrary_cycle_start_time)

    # TODO add tests where the 3 specified times do not exactly match those
    # in the time array.
