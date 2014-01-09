import os

import numpy as np
from numpy import testing

from epimysium import postprocessing as pproc

parentdir = os.path.abspath(os.path.dirname(__file__))

def test_filter_emg():
    dat = np.loadtxt(os.path.join(parentdir, 'emg_raw_and_filtered.txt'))
    raw = dat[:, 0]
    filtered_des = dat[:, 1]
    sampling_rate = 2000
    filtered = pproc.filter_emg(raw, sampling_rate)
    testing.assert_allclose(filtered, filtered_des, atol=0.001)

def test_shift_data_to_cycle():

    ordinate = np.array([2.0, 1., 2., 3., 4., 5., 6.])
    time = np.array([0.5, 1.0, 1.2, 1.35, 1.4, 1.5, 1.8])
    arbitrary_cycle_start_time = 1.0
    arbitrary_cycle_end_time = 1.5
    new_cycle_start_time = 1.35
    shifted_time, shifted_ordinate = pproc.shift_data_to_cycle(
        arbitrary_cycle_start_time, arbitrary_cycle_end_time,
        new_cycle_start_time,
        time, ordinate)

    testing.assert_allclose(shifted_ordinate, np.array([3, 4, np.nan, 1, 2]))
    testing.assert_allclose(shifted_time,
            np.array([0.0, 0.05, 0.15, 0.30, 0.5]))

    testing.assert_equal(shifted_time.size, 5)
    testing.assert_equal(shifted_ordinate.size, 5)

    assert shifted_time[0] == 0.0
    assert shifted_time[-1] == (arbitrary_cycle_end_time -
            arbitrary_cycle_start_time)

    # TODO add tests where the 3 specified times do not exactly match those
    # in the time array.
    arbitrary_cycle_start_time = 0.9
    arbitrary_cycle_end_time = 1.6
    new_cycle_start_time = 1.35
    shifted_time, shifted_ordinate = pproc.shift_data_to_cycle(
        arbitrary_cycle_start_time, arbitrary_cycle_end_time,
        new_cycle_start_time,
        time, ordinate)
    testing.assert_allclose(shifted_ordinate, np.array([3, 4, np.nan, 1, 2]))
    testing.assert_allclose(shifted_time,
            np.array([0.0, 0.05, 0.25, 0.40, 0.7]))

    testing.assert_equal(shifted_time.size, 5)
    testing.assert_equal(shifted_ordinate.size, 5)

    # Don't actually want this to be true.
    testing.assert_equal(shifted_time[0], 0.0)
    testing.assert_equal(shifted_time[-1] - shifted_time[0],
            (arbitrary_cycle_end_time - arbitrary_cycle_start_time))

def test_shift_data_to_cycle_for_less_than_full_cycle():
    # When [arbitrary_cycle_start_time, arbitrary_cycle_end_time] is NOT a
    # subset of the time interval of the given data. That is, we don't have
    # data for the complete gait cycle.
    import pylab as pl
    arbitrary_cycle_start_time = 0.2
    arbitrary_cycle_end_time = 1.3
    duration = arbitrary_cycle_end_time - arbitrary_cycle_start_time

    # Case 1: data ends early.
    data_starts_at = 0.1
    data_ends_at = 0.8

    x1 = np.linspace(data_starts_at, data_ends_at, 81)
    x10 = x1[pproc.nearest_index(x1, arbitrary_cycle_start_time)]
    y1 = 2 * (x1 - x10)
    pl.subplot(2, 2, 1)
    pl.plot(x1, y1) 
    pl.xlim(xmin=arbitrary_cycle_start_time, xmax=arbitrary_cycle_end_time)
    pl.ylim(0, y1[-1])

    # Case 1a: don't shift, just represent as percent gait cycle.
    xs1a, ys1a = pproc.shift_data_to_cycle(arbitrary_cycle_start_time,
            arbitrary_cycle_end_time, arbitrary_cycle_start_time, x1, y1)
    # Don't expect to satisfy this in this case; see docstring:
    # testing.assert_equal(xs1a[-1], duration)
    pl.plot(xs1a, ys1a)
    testing.assert_almost_equal(pproc.percent_duration(xs1a, 0, duration)[-1],
            54.545454545454)
    pl.subplot(2, 2, 2)
    pl.plot(pproc.percent_duration(xs1a, 0, duration), ys1a, label='1a')
    pl.xlim(0, 100)
    #pl.ylim(0, y1[-1])

    # Case 1b: shift to a new cycle start time.
    xs1b, ys1b = pproc.shift_data_to_cycle(arbitrary_cycle_start_time,
            arbitrary_cycle_end_time, 0.4, x1, y1)
    # When the gap is in the middle, we DO want to satisfy this:
    testing.assert_almost_equal(xs1b[-1], duration)
    pl.plot(pproc.percent_duration(xs1b, 0, duration), ys1b, label='1b')

    # Case 1c: the new cycle start time is in the missing-data portion.
    xs1c, ys1c = pproc.shift_data_to_cycle(arbitrary_cycle_start_time,
            arbitrary_cycle_end_time, 1.0, x1, y1)
    testing.assert_almost_equal(xs1c[-1] + 0.2, duration) 
    testing.assert_almost_equal(xs1c[0], arbitrary_cycle_end_time - 1.0)
    # Increasing (last time value is the maximum value).
    assert np.all(np.diff(xs1c) > 0)
    pl.plot(pproc.percent_duration(xs1c, 0, duration), ys1c, label='1c')

    pl.legend()

    # Case 2: data starts late.
    data_starts_at = 0.4
    data_ends_at = 1.5

    x2 = np.linspace(data_starts_at, data_ends_at, 111)
    x20 = x2[pproc.nearest_index(x2, arbitrary_cycle_start_time)]

    pl.subplot(2, 2, 3)
    y2 = 2 * (x2 - x20)
    pl.plot(x2, y2)
    pl.xlim(xmin=arbitrary_cycle_start_time, xmax=arbitrary_cycle_end_time)
    pl.ylim(ymin=0, ymax=1.8)

    # Case 2a: don't shift.
    xs2a, ys2a = pproc.shift_data_to_cycle(arbitrary_cycle_start_time,
            arbitrary_cycle_end_time, arbitrary_cycle_start_time, x2, y2)
    pl.plot(xs2a, ys2a)
    testing.assert_equal(xs2a[-1], duration)
    pl.subplot(2, 2, 4)
    pl.plot(pproc.percent_duration(xs2a, 0, duration), ys2a, label='2a')
    pl.xlim(0, 100)
    pl.ylim(ymin=0, ymax=1.8)

    # Case 2b:  shift to a new cycle start time.
    xs2b, ys2b = pproc.shift_data_to_cycle(arbitrary_cycle_start_time,
            arbitrary_cycle_end_time, 0.8, x2, y2)
    testing.assert_equal(xs2b[-1], duration)
    pl.plot(pproc.percent_duration(xs2b, 0, duration), ys2b, label='2b')

    # Case 2c: the new cycle start time is in the missing-data portion.
    xs2c, ys2c = pproc.shift_data_to_cycle(arbitrary_cycle_start_time,
            arbitrary_cycle_end_time, 0.3, x2, y2)
    # Can't satisfy this when the gap is not in the middle:
    #testing.assert_equal(xs2c[-1] - xs2c[0], duration)
    pl.plot(pproc.percent_duration(xs2c, 0, duration), ys2c, label='2c')

    pl.legend()
    # Case 3: data starts late AND ends early.
    # TODO

if __name__ == '__main__':
    #import pylab as pl
    #test_shift_data_to_cycle_for_less_than_full_cycle()
    #pl.show()
    pass
