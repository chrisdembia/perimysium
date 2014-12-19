"""Contains methods for postprocessing simulations, or computing derived
quantities, such as sum of squared activations, from simulation (e.g., CMC)
results.

"""

import collections
import copy
import os
import re

import numpy as np
import matplotlib
import pylab as pl
import tables
from scipy.signal import butter, filtfilt
from scipy.stats import nanmean, nanstd

from perimysium import dataman

def savefigtolog(figname, *args, **kwargs):
    pl.savefig(os.path.join(os.environ['LOGFIGS'], figname), *args, **kwargs)


def nearest_index(array, val):
    return np.abs(array - val).argmin()

def filter_emg(raw_signal, sampling_rate, bandpass_order=6,
        bandpass_lower_frequency=50, bandpass_upper_frequency=500,
        lowpass_order=4,
        lowpass_frequency=7.5,
        cd_lowpass_frequency=15.0):
    """Filters a raw EMG signal. The signal must have been sampled at a
    constant rate. We perform the following steps:

    1. Butterworth bandpass filter.
    2. Butterworth lowpass filter.
    3. Critically damped lowpass filter.

    The signals are applied forward and backward (`filtfilt`), which should
    prevent a time delay.

    Parameters
    ----------
    raw_signal : array_like
        Raw EMG signal.
    sampling_rate : float
        In Hertz.
    bandpass_order : int, optional
    bandpass_lower_frequency : float, optional
        In the bandpass filter, what is the lower cutoff frequency? In Hertz.
    bandpass_upper_frequency : float, optional
        In the bandpass filter, what is the upper cutoff frequency? In Hertz.
    lowpass_order : int, optional
    lowpass_frequency : float, optional
        In the lowpass filter, what is the cutoff frequency? In Hertz.
    cd_lowpass_frequency : float, optional
        In the Critically damped lowpass filter, what is the cutoff frequency?
        In Hertz.

    Returns
    -------
    filtered_signal : array_like

    """
    nyquist_frequency = 0.5 * sampling_rate

    # Bandpass.
    # ---------
    normalized_bandpass_lower = bandpass_lower_frequency / nyquist_frequency
    normalized_bandpass_upper = bandpass_upper_frequency / nyquist_frequency
    bandpass_cutoffs = [normalized_bandpass_lower, normalized_bandpass_upper]
    bandpass_b, bandpass_a = butter(bandpass_order, bandpass_cutoffs,
            btype='bandpass')

    bandpassed = filtfilt(bandpass_b, bandpass_a, raw_signal)

    # Rectify.
    # --------
    rectified = np.abs(bandpassed)

    # Lowpass.
    # --------
    lowpass_cutoff = lowpass_frequency / nyquist_frequency
    lowpass_b, lowpass_a = butter(lowpass_order, lowpass_cutoff)

    lowpassed = filtfilt(lowpass_b, lowpass_a, rectified)

    # Critically damped filter.
    # -------------------------
    cd_order = 4
    cdfed = filter_critically_damped(lowpassed, sampling_rate,
            cd_lowpass_frequency, order=4)

    return cdfed

def filter_critically_damped(data, sampling_rate, lowpass_cutoff_frequency,
        order=4):
    """See Robertson, 2003. This code is transcribed from some MATLAB code that
    Amy Silder gave me. This implementation is slightly different from that
    appearing in Robertson, 2003. We only allow lowpass filtering.

    Parameters
    ----------
    data : array_like
        The signal to filter.
    sampling_rate : float
    lowpass_cutoff_frequency : float
        In Hertz (not normalized).
    order : int, optional
        Number of filter passes.

    Returns
    -------
    data : array_like
        Filtered data.

    """
    # 3 dB cutoff correction.
    Clp = (2.0 ** (1.0 / (2.0 * order)) - 1.0) ** (-0.5)

    # Corrected cutoff frequency.
    flp = Clp * lowpass_cutoff_frequency / sampling_rate

    # Warp cutoff frequency from analog to digital domain.
    wolp = np.tan(np.pi * flp)

    # Filter coefficients, K1 and K2.
    # lowpass: a0 = A0, a1 = A1, a2 = A2, b1 = B2, b2 = B2
    K1lp = 2.0 * wolp
    K2lp = wolp ** 2

    # Filter coefficients.
    a0lp = K2lp / (1.0 + K1lp + K2lp)
    a1lp = 2.0 * a0lp
    a2lp = a0lp
    b1lp = 2.0 * a0lp  * (1.0 / K2lp - 1.0)
    b2lp = 1.0 - (a0lp + a1lp + a2lp + b1lp)

    num_rows = len(data)
    temp_filtered = np.zeros(num_rows)
    # For order = 4, we go forward, backward, forward, backward.
    for n_pass in range(order):
        for i in range(2, num_rows):
            temp_filtered[i] = (a0lp * data[i] +
                    a1lp * data[i - 1] +
                    a2lp * data[i - 2] +
                    b1lp * temp_filtered[i - 1] +
                    b2lp * temp_filtered[i - 2])
        # Perform the filter backwards.
        data = np.flipud(temp_filtered)
        temp_filtered = np.zeros(num_rows)

    return data


def metabolic_expenditure_const_eff(power, exclude=None, concentric_eff=0.25,
        eccentric_eff=-1.20):
    """Computes metabolic muscle power expenditure using constant muscle
    efficiencies for concentric and eccentric contraction. Defaut efficiencies
    are accepted in the literature (Voinescu, 2012). Assumes muscles are the
    actuators with names ending in '_l' or '_r'.

    Parameters
    ----------
    power : pytables.Table of an Actuation_power.sto Storage file.
    exclude : list of str's, optional
        Names of columns in `power` to exclude in the computation, in case
        there are columns for actuators that are not muscles but whose name
        still ends in '_l' or '_r'.
    concentric_eff : float, optional
        Efficiency of concentric contraction (shortening).
    eccentric_eff : float, optional
        Efficiency of eccentric contraction (lengthening). This is when a
        muscle is generating force even though it is lengthening. That is, the
        muscle is "braking". This should be negative, because it still costs
        energy to "brake" a muscle.

    Returns
    -------
    met_expenditure_rate : np.array
        Metabolic rate throughout the simulation.
    met_expenditure_avg_rate : float
        Average metabolic rate over the time interval.

    """
    met_expenditure_rate = np.zeros(len(power.cols.time))
    for actu_name in power.colnames:
        if (actu_name.endswith('_l') or actu_name.endswith('_r') and (
            (actu_name not in exclude) if exclude else None)):
            # Boolean expressions to sort out positive and negative work.
            met_expenditure_rate += (
                    ((power.col(actu_name) > 0) / concentric_eff +
                            (power.col(actu_name) < 0) / eccentric_eff) *
                    power.col(actu_name))
    return met_expenditure_rate, avg(power.cols.time, met_expenditure_rate)


def sum_of_squared_activations(states_table, weight_map=None):
    """Computes the sum, across all muscles, of the square of muscle
    activations.

    Parameters
    ----------
    states_table : tables.Table
        A pyTables table containing muscle activation time series data.
    weight_map : dict, optional
        If you want to weight muscles by, e.g. muscle volume. Keys are muscle
        names, values are the weights as floats.

    Returns
    -------
    SSA : numpy.ndarray
        Time series of the sum of squared muscle activations.

    """
    SSA = np.zeros(states_table.col('time').shape)
    for col_name in states_table.colnames:
        if col_name.endswith('activation'):
            if weight_map == None:
                weight = 1
            else:
                weight = weight_map[col_name.replace('_activation', '')]
            SSA += weight * states_table.col(col_name)**2
    return SSA

def avg_sum_of_squared_activations(states_table, cycle_duration=None,
        cycle_start=None, weight_map=None):
    """Computes the average value of the sum of squared activations.

    Parameters
    ----------
    states_table : tables.Table
        A pyTables table containing muscle activation time series data.
    cycle_duration : optional, float
        If provided, the average will only be over this duration, starting
        either at the initial time, or at `cycle_start` if provided.
    cycle_start : optional, float
        Only used if cycle_duration is not None. The start time of the
        interval over which to conduct the average.
    weight_map : dict, optional
        If you want to weight muscles by, e.g. muscle volume. Keys are muscle
        names, values are the weights as floats.

    Returns
    -------
    ASSA : float
        Average value, in time, of the sum of squared muscle activations.

    """
    time = states_table.cols.time
    ssa = sum_of_squared_activations(states_table, weight_map=weight_map)
    if cycle_duration == None:
        duration = time[-1] - time[0]
        integral = np.trapz(ssa, x=time)
        return integral / duration
    else:
        return avg_over_gait_cycle(states_table.cols.time, ssa,
                cycle_duration, cycle_start=cycle_start)


def avg(time, value, init_time=None, final_time=None, interval=None):
    """Finds the average value of `value` in time, using np.trapz.

    Parameters
    ----------
    time : numpy.array
        Array of time values, probably in seconds.
    value : numpy.array
        Array of the quantity to be averaged.
    init_time : float, optional
        The lower bound on the integral, in units of time (probably seconds).
        By default, the first time in the series is used.
    final_time : float, optional
        The upper bound on the integral, in units of time (probably seconds).
        By default, the last time in the series is used.
    interval : int, optional
        Interval of data points to skip/include in the plot (states files have
        way too many data points for plots). By default, no points are skipped.

    Returns
    -------
    avg : float
        Time-averaged value of the data between the times specified.

    """
    if init_time == None:
        init_idx = 0
    else:
        init_idx = np.abs(time - init_time).argmin()
    if final_time == None:
        final_idx = len(time)
    else:
        final_idx = np.abs(time - final_time).argmin()

    duration = time[final_idx-1] - time[init_idx]
    return np.trapz(value[init_idx:final_idx:interval],
            x=time[init_idx:final_idx:interval]) / duration

def avg_over_gait_cycle(time, value, cycle_duration, cycle_start=None,
        avg_over_half_if_partial_data=True):
    """Average a quantity over a gait cycle with duration/period
    `cycle_duration`. The cycle starts at time[0], unless you specify a
    `cycle_start` time. If `time` does not contain a full gait cycle's worth of
    data (time[-1] - time[0] < cycle_duration), we only average over half the
    gait cycle.

    Parameters
    ----------
    time : array_like
    value : array_like
    cycle_duration : float
    cycle_start : float, optional
    avg_over_gait_cycle : bool, optional

    """
    avail_duration = time[-1] - time[0]

    if cycle_start == None:
        cycle_start = time[0]

    if avail_duration >= cycle_duration:
        if cycle_start + cycle_duration > time[-1]:
            raise Exception('Requested time for integration is unavailable '
                    'in the data.')
        return avg(time, value, init_time=cycle_start,
                final_time=cycle_start + cycle_duration)
    else:
        if cycle_start + 0.5 * cycle_duration > time[-1]:
            raise Exception('Requested time for integration is unavailable '
                    'in the data.')
        return avg(time, value, init_time=cycle_start,
                final_time=cycle_start + 0.5 * cycle_duration)


def specific_metabolic_cost(subject_mass,
        time, value, cycle_duration, cycle_start=None):
    return 1.0 / subject_mass * avg_over_gait_cycle(time, value,
            cycle_duration, cycle_start=cycle_start)

def cost_of_transport(subject_mass,
        forward_speed,
        time, value, cycle_duration, cycle_start=None,
        accel_due_to_gravity=9.81):
    return avg_over_gait_cycle(time, value, cycle_duration,
            cycle_start=cycle_start) / (
                    subject_mass * accel_due_to_gravity * forward_speed)

def average_whole_body_power(actu_power, cycle_duration, cycle_start=None,
        ignore='reserve|^F|^M'):
    total_power = np.zeros(len(actu_power.cols.time[:]))
    for coln in actu_power.colnames:
        if coln != 'time':
            if ignore == None or re.search(ignore, coln) == None:
                total_power += actu_power.col(coln)
    return avg_over_gait_cycle(actu_power.cols.time[:], total_power,
            cycle_duration, cycle_start=cycle_start)


def sorted_maxabs(table, init_time=None, final_time=None, exclude=None,
        include_only=None):
    """Returns sort and argsort for all columns given. The quantity
    bieng sorted is the max absolute value of each column either over all
    times, or over the time interval specified.

    Parameters
    ----------
    table : tables.Table
        A pyTables table with, probably, more than 1 column. One column must be
        'time'. Also, the 'time' column is skipped in the sorting.
    init_time : float, optional
        The lower bound of times at which to check max abs values, in units of
        time (probably seconds).  By default, the first time in the series is
        used.
    final_time : float, optional
        The upper bound of times at which to check max abs values, in units of
        time (probably seconds).  By default, the last time in the series is
        used.
    exclude : str or list of str's, optional
        If str, it's a regular expression. Otherwise, it's the exact
        names of columns to exclude from the sorting (e.g.,
        ['umb_val_wholebody']). Automatically excludes 'time' column.
    include_only : str, optional
        A regular expression that a column name must match to be included in
        the output.

    Returns
    -------
    sorted_vals : numpy.array
        The max absolute value for all non-excluded columns, sorted in
        ascending order.
    sorted_args : list of str's
        The column names corresponding to the sorted values, in ascending
        order.

    """
    if exclude != None and not (type(exclude) == list or type(exclude) == str):
        raise Exception("'exclude' must be None, a str, or a list.")

    time = table.cols.time
    if init_time == None:
        init_idx = 0
    else:
        init_idx = np.abs(time - init_time).argmin()
    if final_time == None:
        final_idx = len(time)
    else:
        final_idx = np.abs(time - final_time).argmin()

    maxs = dict()

    def do_exclude(coln):
        if coln == 'time':
            return True
        elif exclude == None:
            return False
        elif type(exclude) == list:
            return coln in exclude
        elif type(exclude) == str:
            return re.search(exclude, coln)

    def do_include(coln):
        if include_only == None:
            return True
        else:
            return re.search(include_only, coln)

    for coln in table.colnames:
        if do_include(coln) and not do_exclude(coln):
            maxs[coln] = np.max(np.abs(table.col(coln)[init_idx:final_idx]))

    sorted_vals = np.sort(maxs.values())
    sorted_args = [maxs.keys()[idx] for idx in np.argsort(maxs.values())]

    return sorted_vals, sorted_args


def sorted_avg(table, init_time=None, final_time=None, exclude=None,
        include_only=None):
    """Returns sort and argsort for all columns given. The quantity
    bieng sorted is the average of each column (i.e., probe value)
    either over all times, or over the time interval specified.

    Parameters
    ----------
    table : tables.Table
        A pyTables table with, probably, more than 1 column. One column must be
        'time'. Also, the 'time' column is skipped in the sorting.
    init_time : float, optional
        The lower bound on the integral, in units of time (probably seconds).
        By default, the first time in the series is used.
    final_time : float, optional
        The upper bound on the integral, in units of time (probably seconds).
        By default, the last time in the series is used.
    exclude : str or list of str's, optional
        If str, it's a regular expression. Otherwise, it's the exact
        names of columns to exclude from the sorting (e.g.,
        ['umb_val_wholebody']). Automatically excludes 'time' column.
    include_only : str, optional
        A regular expression that a column name must match to be included in
        the output.

    Returns
    -------
    sorted_vals : numpy.array
        The average value for all non-excluded columns, sorted in ascending
        order.
    sorted_args : list of str's
        The column names corresponding to the sorted values, in ascending
        order.

    """
    if exclude != None and not (type(exclude) == list or type(exclude) == str):
        raise Exception("'exclude' must be None, a str, or a list.")

    time = table.cols.time

    avgs = dict()

    def do_exclude(coln):
        if coln == 'time':
            return True
        elif exclude == None:
            return False
        elif type(exclude) == list:
            return coln in exclude
        elif type(exclude) == str:
            return re.search(exclude, coln)

    def do_include(coln):
        if include_only == None:
            return True
        else:
            return re.search(include_only, coln)

    for coln in table.colnames:
        if do_include(coln) and not do_exclude(coln):
            avgs[coln] = avg(time, table.col(coln), init_time, final_time)

    sorted_vals = np.sort(avgs.values())
    sorted_args = [avgs.keys()[idx] for idx in np.argsort(avgs.values())]

    return sorted_vals, sorted_args


def sorted_avg_difference(table1, table2,
        init_time=None, final_time=None, exclude=None, include_only=None,
        interval=None):
    """Returns sort and argsort for a difference in average value of each
    column between two tables (table1 minus table2). Both tables must have the
    same columns.

    Parameters
    ----------
    table1 : tables.Table
        A pyTables table with, probably, more than 1 column. One column must be
        'time'. Also, the 'time' column is skipped in the sorting.
    table2 : tables.Table
        Like table2. Average values for the columns of this table are
        subtracted from average values for the columns of table1.
    init_time : float, optional
        The lower bound on the integral, in units of time (probably seconds).
        By default, the first time in the series is used.
    final_time : float, optional
        The upper bound on the integral, in units of time (probably seconds).
        By default, the last time in the series is used.
    exclude : list of str's, optional
        The names of columns to exclude from the sorting (e.g.,
        ['umb_val_wholebody']). Automatically excludes 'time' column.
    include_only : str, optional
        A substring that a column name must contain to be included in the
        output.
    interval : int, optional
        Interval of data points to skip/include in the plot (states files have
        way too many data points for plots). By default, no points are skipped.

    Returns
    -------
    sorted_vals : numpy.ndarray
        The average value for all non-excluded columns, sorted in ascending
        order.
    sorted_args : list of str's
        The column names corresponding to the sorted values, in ascending
        order.

    """
    if exclude != None and not type(exclude) == list:
        raise Exception("'exclude' must be a list.")

    diff_avgs = dict()

    for coln in table1.colnames:
        if (coln != 'time' and (exclude == None or not coln in exclude) and
                (include_only == None or coln.count(include_only) != 0)):
            diff_avgs[coln] = (
                    avg(table1.cols.time[::interval],
                        table1.col(coln)[::interval],
                        init_time, final_time) -
                    avg(table2.cols.time[::interval],
                        table2.col(coln)[::interval],
                        init_time, final_time))

    sorted_vals = np.sort(diff_avgs.values())
    sorted_args = [diff_avgs.keys()[idx] for idx in
            np.argsort(diff_avgs.values())]

    return sorted_vals, sorted_args


def plot_simulation_verification(sim_group, **kwargs):
    """Plots kinematics, residual, and reserve verification information, using
    matplotlib.

    Parameters
    ----------
    sim_group : tables.Group
        A pyTables group holding tables docked from a simulation (e.g., CMC).
    n_max : int, optional
        Plots the n_max maximum columns. By default, all columns are plotted.
        If n_max = 5, only the 5 most errorful columns are plotted.
    violators_only : bool, optional
        Only the columns violating the "good" threshold are plotted. This
        overrides the n_max setting. By default, this is False.
    show_legend : bool, optional
        Show a legend for all columns plotted. By default, False.

    Returns
    -------
    fig : pl.figure
        Figure handle, or list of handles.

    """
    if 'show_legend' in kwargs and kwargs['show_legend']:
        figk = plot_kinematics_verification(sim_group.pErr, **kwargs)
        figrd = plot_residuals_verification(sim_group.Actuation_force, **kwargs)
        figrv = plot_reserves_verification(sim_group.Actuation_force, **kwargs)
        return [figk, figrd, figrv]
    else:
        fig = pl.figure(figsize=(15, 8))
        plot_kinematics_verification(sim_group.pErr, big_picture=(3, 0),
                **kwargs)
        plot_residuals_verification(sim_group.Actuation_force,
                big_picture=(3, 1), **kwargs)
        plot_reserves_verification(sim_group.Actuation_force,
                big_picture=(3, 2), **kwargs)
        return fig


def rms(array):
    """Root-mean-square of the data. Behavior for multidimensional arrays is
    unknown.

    Parameters
    ----------
    array : numpy.array
        Data to be root-mean-squared.

    Returns
    -------
    rms : float
        Root-mean-square of data.

    """
    return np.sqrt(np.mean(array**2))


def plot_rra_gait_info(rra_results_dir):
    """Plots the information useful to tweaking RRA tasks and actuators:
    residual forces and moments, as well as the kinematics errors that are
    important for gait.

    Parameters
    ----------
    rra_results_dir : str
        Path to folder containing results from RRA. We'll grab the
        Actuation_force and pErr Storage files.

    Returns
    -------
    fig : pylab.figure

    """
    for fname in os.listdir(rra_results_dir):
        if fname.endswith('pErr.sto'):
            pErr_fpath = os.path.join(rra_results_dir, fname)
        elif fname.endswith('Actuation_force.sto'):
            actu_fpath = os.path.join(rra_results_dir, fname)

    def plot_thresholds(data, val):
        pl.plot([data['time'][0], data['time'][-1]], [val, val],
                c=[0.7, 0.7, 0.7])
        pl.plot([data['time'][0], data['time'][-1]], [-val, -val],
                c=[0.7, 0.7, 0.7])

    def make_pretty_perr(ylabel='rotation error (deg)'):
        pl.axhline(0, c='k')
        pl.ylabel(ylabel)
        pl.xlabel('time (s)')
        pl.xlim(xmin=pErr['time'][0], xmax=pErr['time'][-1])
        pl.ylim((-2, 2))
        pl.legend(**legend_kwargs)
        pl.grid(axis='y')
        plot_thresholds(pErr, 1)

    legend_kwargs = {'loc': 'best', 'prop': {'size': 12}, 'frameon': False}

    pErr = dataman.storage2numpy(pErr_fpath)
    actu = dataman.storage2numpy(actu_fpath)

    fig = pl.figure(figsize=(12, 24))
    pl.subplot(621)
    pl.title('residual forces')
    for coln in ['FX', 'FY', 'FZ']:
        pl.plot(actu['time'], actu[coln], label=coln)
    pl.axhline(0, c='k')
    pl.ylabel('force (N)')
    pl.legend(**legend_kwargs)
    pl.ylim((-40, 40))
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.grid(axis='y')
    plot_thresholds(actu, 10)
    plot_thresholds(actu, 25)

    pl.subplot(622)
    pl.title('residual moments')
    for coln in ['MX', 'MY', 'MZ']:
        pl.plot(actu['time'], actu[coln], label=coln)
    pl.axhline(0, c='k')
    pl.ylabel('torque (N-m)')
    pl.legend(**legend_kwargs)
    pl.ylim((-40, 40))
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.grid(axis='y')
    plot_thresholds(actu, 30)

    pl.subplot(623)
    m2cm = 100
    pl.title('pelvis translation')
    for coln in ['pelvis_tx', 'pelvis_ty', 'pelvis_tz']:
        pl.plot(pErr['time'], pErr[coln] * m2cm, label=coln[-1])
    make_pretty_perr('translation error (cm)')

    pl.subplot(624)
    rad2deg = np.rad2deg(1.0)
    pl.title('pelvis rotations')
    for coln in ['pelvis_tilt', 'pelvis_list', 'pelvis_rotation']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[-1])
    make_pretty_perr()

    pl.subplot(625)
    pl.title('left lower limb')
    for coln in ['hip_flexion_l', 'knee_angle_l', 'ankle_angle_l']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[0])
    make_pretty_perr()

    pl.subplot(626)
    pl.title('right lower limb')
    for coln in ['hip_flexion_r', 'knee_angle_r', 'ankle_angle_r']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[0])
    make_pretty_perr()

    pl.subplot(627)
    pl.title('lumbar rotations')
    for coln in ['lumbar_bending', 'lumbar_extension', 'lumbar_rotation']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[-1])
    make_pretty_perr()

    pl.subplot(628)
    pl.title('hips')
    labels = ['rot_r', 'rot_l', 'add_r', 'add_l']
    for i, coln in enumerate(['hip_rotation_r', 'hip_rotation_l',
        'hip_adduction_r', 'hip_adduction_l']):
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=labels[i])
    make_pretty_perr()

    pl.subplot(629)
    pl.title('left upper arm')
    for coln in ['arm_flex_l', 'arm_add_l', 'arm_rot_l']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[1])
    make_pretty_perr()

    pl.subplot(6, 2, 10)
    pl.title('right upper arm')
    for coln in ['arm_flex_r', 'arm_add_r', 'arm_rot_r']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[1])
    make_pretty_perr()

    pl.subplot(6, 2, 11)
    pl.title('left lower arm')
    for coln in ['elbow_flex_l', 'pro_sup_l']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[1])
    make_pretty_perr()

    pl.subplot(6, 2, 12)
    pl.title('right lower arm')
    for coln in ['elbow_flex_r', 'pro_sup_r']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[1])
    make_pretty_perr()

    pl.tight_layout()
    return fig


def plot_so_gait_info(so_results_dir):
    """Plots residuals and reserves from static optimization results for a
    model like gait2392. Reserves must be named as 'reserve_<coord-name>'.

    Parameters
    ----------
    so_results_dir : str
        Path to folder containing results from Static Optimization. It must
        have the following file:
            - <name>_StaticOptimization_force.sto

    Returns
    -------
    fig : pylab.figure

    """
    for fname in os.listdir(so_results_dir):
        if fname.endswith('_force.sto'):
            fpath = os.path.join(so_results_dir, fname)
            break

    def plot_thresholds(data, val):
        pl.plot([data['time'][0], data['time'][-1]], [val, val],
                c=[0.7, 0.7, 0.7])
        pl.plot([data['time'][0], data['time'][-1]], [-val, -val],
                c=[0.7, 0.7, 0.7])

    legend_kwargs = {'loc': 'best', 'prop': {'size': 12}}

    actu = dataman.storage2numpy(fpath, 2)

    fig = pl.figure(figsize=(10, 12))
    grid_size = (3, 2)
    pl.subplot2grid(grid_size, (0, 0))
    pl.title('residual forces')
    for coln in ['FX', 'FY', 'FZ']:
        pl.plot(actu['time'], actu[coln], label=coln)
    pl.axhline(0, c='k')
    pl.ylabel('force (N)')
    pl.legend(**legend_kwargs)
    pl.ylim((-40, 40))
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.grid(axis='y')
    plot_thresholds(actu, 10)
    plot_thresholds(actu, 25)

    pl.subplot2grid(grid_size, (0, 1))
    pl.title('residual moments')
    for coln in ['MX', 'MY', 'MZ']:
        pl.plot(actu['time'], actu[coln], label=coln)
    pl.axhline(0, c='k')
    pl.ylabel('torque (N-m)')
    pl.legend(**legend_kwargs)
    pl.ylim((-40, 40))
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.grid(axis='y')
    plot_thresholds(actu, 30)

    pl.subplot2grid(grid_size, (1, 0))
    pl.title('left lower limb reserves')
    for coln in ['hip_flexion_l', 'knee_angle_l', 'ankle_angle_l']:
        pl.plot(actu['time'], actu['reserve_%s' % coln],
                label=coln.split('_')[0])
    pl.axhline(0, c='k')
    pl.ylabel('torque (N-m)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.ylim((-25, 25))
    pl.grid(axis='y')
    plot_thresholds(actu, 10)

    pl.subplot2grid(grid_size, (1, 1))
    pl.title('right lower limb reserves')
    for coln in ['hip_flexion_r', 'knee_angle_r', 'ankle_angle_r']:
        pl.plot(actu['time'], actu['reserve_%s' % coln],
                label=coln.split('_')[0])
    pl.axhline(0, c='k')
    pl.ylabel('torque (N-m)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.ylim((-25, 25))
    pl.grid(axis='y')
    plot_thresholds(actu, 10)

    pl.subplot2grid(grid_size, (2, 0))
    pl.title('lumbar rotations reserves')
    for coln in ['lumbar_bending', 'lumbar_extension', 'lumbar_rotation']:
        actu_name = 'reserve_%s' % coln
        if actu_name in actu.dtype.names:
            pl.plot(actu['time'], actu[actu_name],
                    label=coln.split('_')[-1])
        else:
            # Apoorva's model does not have lumbar muscles.
            pl.text(0, 0, 'no lumbar reserves', ha='left', va='bottom')
    pl.axhline(0, c='k')
    pl.ylabel('moment (N-m)')
    pl.xlabel('time (s)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.ylim((-25, 25))
    pl.grid(axis='y')
    plot_thresholds(actu, 10)

    pl.subplot2grid(grid_size, (2, 1))
    pl.title('hip reserves')
    labels = ['rot_r',' rot_l', 'add_r',' add_l']
    for i, coln in enumerate(['hip_rotation_r', 'hip_rotation_l',
        'hip_adduction_r', 'hip_adduction_l']):
        pl.plot(actu['time'], actu['reserve_%s' % coln], label=labels[i])
    pl.axhline(0, c='k')
    pl.ylabel('moment (N-m)')
    pl.xlabel('time (s)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.ylim((-25, 25))
    pl.grid(axis='y')
    plot_thresholds(actu, 10)

    pl.tight_layout()
    return fig


def plot_cmc_gait_info(cmc_results_dir):
    """

    Parameters
    ----------
    cmc_results_dir : str
        Path to folder containing results from CMC. We'll grab the
        Actuation_force and pErr Storage files.

    Returns
    -------
    fig : pylab.figure

    """
    for fname in os.listdir(cmc_results_dir):
        if fname.endswith('pErr.sto'):
            pErr_fpath = os.path.join(cmc_results_dir, fname)
        elif fname.endswith('Actuation_force.sto'):
            actu_fpath = os.path.join(cmc_results_dir, fname)

    def plot_thresholds(data, val):
        pl.plot([data['time'][0], data['time'][-1]], [val, val],
                c=[0.7, 0.7, 0.7])
        pl.plot([data['time'][0], data['time'][-1]], [-val, -val],
                c=[0.7, 0.7, 0.7])

    legend_kwargs = {'loc': 'best', 'prop': {'size': 12}}

    pErr = dataman.storage2numpy(pErr_fpath)
    actu = dataman.storage2numpy(actu_fpath)

    fig = pl.figure(figsize=(16, 12))
    pl.subplot2grid((3, 4), (0, 2))
    pl.title('residual forces')
    for coln in ['FX', 'FY', 'FZ']:
        pl.plot(actu['time'], actu[coln], label=coln)
    pl.axhline(0, c='k')
    pl.ylabel('force (N)')
    pl.legend(**legend_kwargs)
    pl.ylim((-40, 40))
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.grid(axis='y')
    plot_thresholds(actu, 10)
    plot_thresholds(actu, 25)

    pl.subplot2grid((3, 4), (0, 3))
    pl.title('residual moments')
    for coln in ['MX', 'MY', 'MZ']:
        pl.plot(actu['time'], actu[coln], label=coln)
    pl.axhline(0, c='k')
    pl.ylabel('torque (N-m)')
    pl.legend(**legend_kwargs)
    pl.ylim((-40, 40))
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.grid(axis='y')
    plot_thresholds(actu, 30)

    pl.subplot2grid((3, 4), (0, 0))
    m2cm = 100
    pl.title('pelvis translation')
    for coln in ['pelvis_tx', 'pelvis_ty', 'pelvis_tz']:
        pl.plot(pErr['time'], pErr[coln] * m2cm, label=coln[-1])
    pl.axhline(0, c='k')
    pl.ylabel('translation error (cm)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=pErr['time'][0], xmax=pErr['time'][-1])
    pl.ylim((-2, 2))
    pl.grid(axis='y')
    plot_thresholds(pErr, 1)

    pl.subplot2grid((3, 4), (0, 1))
    rad2deg = np.rad2deg(1.0)
    pl.title('pelvis rotations')
    for coln in ['pelvis_tilt', 'pelvis_list', 'pelvis_rotation']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[-1])
    pl.axhline(0, c='k')
    pl.ylabel('rotation error (deg)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=pErr['time'][0], xmax=pErr['time'][-1])
    pl.ylim((-2, 2))
    pl.grid(axis='y')
    plot_thresholds(pErr, 1)

    pl.subplot2grid((3, 4), (1, 0))
    pl.title('left lower limb')
    for coln in ['hip_flexion_l', 'knee_angle_l', 'ankle_angle_l']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[0])
    pl.axhline(0, c='k')
    pl.ylabel('rotation error (deg)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=pErr['time'][0], xmax=pErr['time'][-1])
    pl.ylim((-2, 2))
    pl.grid(axis='y')
    plot_thresholds(pErr, 1)

    pl.subplot2grid((3, 4), (1, 2))
    pl.title('left lower limb reserves')
    for coln in ['hip_flexion_l', 'knee_angle_l', 'ankle_angle_l']:
        pl.plot(actu['time'], actu['reserve_%s' % coln],
                label=coln.split('_')[0])
    pl.axhline(0, c='k')
    pl.ylabel('torque (N-m)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.ylim((-25, 25))
    pl.grid(axis='y')
    plot_thresholds(actu, 10)

    pl.subplot2grid((3, 4), (1, 1))
    pl.title('right lower limb')
    for coln in ['hip_flexion_r', 'knee_angle_r', 'ankle_angle_r']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[0])
    pl.axhline(0, c='k')
    pl.ylabel('rotation error (deg)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=pErr['time'][0], xmax=pErr['time'][-1])
    pl.ylim((-2, 2))
    pl.grid(axis='y')
    plot_thresholds(pErr, 1)

    pl.subplot2grid((3, 4), (1, 3))
    pl.title('right lower limb reserves')
    for coln in ['hip_flexion_r', 'knee_angle_r', 'ankle_angle_r']:
        pl.plot(actu['time'], actu['reserve_%s' % coln],
                label=coln.split('_')[0])
    pl.axhline(0, c='k')
    pl.ylabel('torque (N-m)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.ylim((-25, 25))
    pl.grid(axis='y')
    plot_thresholds(actu, 10)

    pl.subplot2grid((3, 4), (2, 0))
    pl.title('lumbar rotations')
    for coln in ['lumbar_bending', 'lumbar_extension', 'lumbar_rotation']:
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=coln.split('_')[-1])
    pl.axhline(0, c='k')
    pl.ylabel('rotation error (deg)')
    pl.xlabel('time (s)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=pErr['time'][0], xmax=pErr['time'][-1])
    pl.ylim((-2, 2))
    pl.grid(axis='y')
    plot_thresholds(pErr, 1)

    pl.subplot2grid((3, 4), (2, 2))
    pl.title('lumbar rotations reserves')
    for coln in ['lumbar_bending', 'lumbar_extension', 'lumbar_rotation']:
        actu_name = 'reserve_%s' % coln
        if actu_name in actu.dtype.names:
            pl.plot(actu['time'], actu[actu_name],
                    label=coln.split('_')[-1])
        else:
            # Apoorva's model does not have lumbar muscles.
            pl.text(0, 0, 'no lumbar reserves', ha='left', va='bottom')
    pl.axhline(0, c='k')
    pl.ylabel('torque (N-m)')
    pl.xlabel('time (s)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.ylim((-25, 25))
    pl.grid(axis='y')
    plot_thresholds(actu, 10)

    pl.subplot2grid((3, 4), (2, 1))
    pl.title('hips')
    labels = ['rot_r',' rot_l', 'add_r',' add_l']
    for i, coln in enumerate(['hip_rotation_r', 'hip_rotation_l',
        'hip_adduction_r', 'hip_adduction_l']):
        pl.plot(pErr['time'], pErr[coln] * rad2deg, label=labels[i])
    pl.axhline(0, c='k')
    pl.ylabel('rotation error (deg)')
    pl.xlabel('time (s)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=pErr['time'][0], xmax=pErr['time'][-1])
    pl.ylim((-2, 2))
    pl.grid(axis='y')
    plot_thresholds(pErr, 1)

    pl.subplot2grid((3, 4), (2, 3))
    pl.title('hip reserves')
    labels = ['rot_r',' rot_l', 'add_r',' add_l']
    for i, coln in enumerate(['hip_rotation_r', 'hip_rotation_l',
        'hip_adduction_r', 'hip_adduction_l']):
        pl.plot(actu['time'], actu['reserve_%s' % coln], label=labels[i])
    pl.axhline(0, c='k')
    pl.ylabel('torque (N-m)')
    pl.xlabel('time (s)')
    pl.legend(**legend_kwargs)
    pl.xlim(xmin=actu['time'][0], xmax=actu['time'][-1])
    pl.ylim((-25, 25))
    pl.grid(axis='y')
    plot_thresholds(actu, 10)

    pl.tight_layout()
    return fig

def plot_lower_limb_kinematics(kinematics_q_fpath, gl=None):
    """Plots pelvis tilt, pelvis list, pelvis rotation, hip adduction, hip
    flexion, knee angle, and ankle angle for both limbs.

    Parameters
    ----------
    kinematics_q_fpath : str
        Path to a Kinematics_q.sto file.
    gl : dataman.GaitLandmarks, optional
        If provided, the plots are for a single gait cycle.

    Returns
    -------
    fig : pylab.figure

    """
    fig = pl.figure(figsize=(7, 10))
    dims = (4, 2)

    sto = dataman.storage2numpy(kinematics_q_fpath)

    def plot(time, y, label, side):
        if gl != None:
            plot_pgc(time, y, gl, side=side, plot_toeoff=True, label=label)
        else:
            pl.plot(time, y, label=label)

    def plot_coord(coord, side='right'):
        plot(sto['time'], sto[coord], side, side)
    def plot_one(loc, coord, ylim):
        ax = pl.subplot2grid(dims, loc)
        plot_coord(coord)
        pl.ylim(ylim)
        pl.axhline(0, color='gray', zorder=0)
        pl.title(coord)
    def plot_both_sides(loc, coord_pre, ylim):
        ax = pl.subplot2grid(dims, loc)
        for side in ['left', 'right']:
            coord = '%s_%s' % (coord_pre, side[0])
            plot_coord(coord, side)
        pl.legend(frameon=False)
        pl.ylim(ylim)
        pl.axhline(0, color='gray', zorder=0)
        pl.title(coord_pre)

    plot_one((0, 0), 'pelvis_tilt', [-10, 10])
    plot_one((1, 0), 'pelvis_list', [-15, 15])
    plot_one((2, 0), 'pelvis_rotation', [-10, 10])
    pl.xlabel('time (% gait cycle)')

    plot_both_sides((0, 1), 'hip_adduction', [-20, 20])
    plot_both_sides((1, 1), 'hip_flexion', [-30, 50])
    plot_both_sides((2, 1), 'knee_angle', [-10, 90])
    plot_both_sides((3, 1), 'ankle_angle', [-40, 25])
    pl.xlabel('time (% gait cycle)')

    pl.tight_layout() #fig) #, rect=[0, 0, 1, 0.95])
    return fig

def plot_kinematics_verification(pErr_table,
        n_max=None, violators_only=False, show_legend=False, big_picture=None):
    """Plots kinematics verification information, with comparisons to "good"
    and "okay", and "bad" thresholds, using matplotlib. Assumes that all
    translation-related degrees of freedom have names ending with either 'tx',
    'ty', or 'tz'; all other degrees of freedom are treated as rotational.

    Regions on the plot are colored as follows:
    green: good
    yellow: okay
    red: bad (only shown if necessary)
    The thresholds come from OpenSim's user manual online.

    Parameters
    ----------
    pErr_table : tables.Table
        A pyTables table containing pErr data; probably translational and
        rotational.
    n_max : int, optional
        Plots the n_max maximum columns. By default, all columns are plotted.
        If n_max = 5, only the 5 most errorful columns are plotted.
    violators_only : bool, optional
        Only the columns violating the "good" threshold are plotted. This
        overrides the n_max setting. By default, this is False.
    show_legend : bool, optional
        Show a legend for all columns plotted. By default, False.
    big_picture : tuple of 2 int's, optional (default: None)
        If these two plots are to make up part of a plot grid, specify a tuple
        with the first element as the number of columns in the figure, and the
        second element as the column index (index starts at 0) for these plots.
        Example is: (3, 0).

    Returns
    -------
    fig : pl.figure
        Figure handle if this method creates a new figure.

    """
    legend_kwargs = {'loc': 'upper left', 'bbox_to_anchor': (1, 1)}
    MtoCM = 100.0

    trans_good_thresh = 1 # centimeters
    trans_okay_thresh = 2
    rot_good_thresh = 2 # degrees
    rot_okay_thresh = 5

    time = pErr_table.cols.time
    duration = time[-1] - time[0]

    sorted_pErr, sorted_colns = sorted_maxabs(pErr_table)

    max_trans = -np.inf
    for i, coln in enumerate(sorted_colns):
        if coln.endswith('tx') or coln.endswith('ty') or coln.endswith('tz'):
            max_trans = max(max_trans, sorted_pErr[i] * MtoCM)
    ylim_trans = max(trans_okay_thresh, 1.1 * max_trans)

    fig = None
    if big_picture != None:
        n_col = big_picture[0]
        this_col = big_picture[1]
    else:
        n_col = 1
        this_col = 0
        fig = pl.figure(figsize=(5, 8))

    # Translation
    # -----------
    pl.subplot2grid((2, n_col), (0, this_col))

    # -- Reference lines.
    pl.plot([time[0], time[-1]], [0, 0], color=(0.5, 0.5, 0.5))
    # Good.
    pl.gca().add_patch(pl.Rectangle((time[0], -trans_good_thresh),
        duration, 2 * trans_good_thresh,
        color='g', alpha=0.1))
    for sgn in [-1, 1]:
        # Okay.
        pl.gca().add_patch(pl.Rectangle((time[0], sgn * trans_okay_thresh),
            duration, -sgn * (trans_okay_thresh - trans_good_thresh),
            color='y', alpha=0.1))
        # Bad.
        pl.gca().add_patch(pl.Rectangle((time[0], sgn * ylim_trans),
            duration, -sgn * (ylim_trans - trans_okay_thresh),
            color='r', alpha=0.1))

    # -- Actual data.
    def plot(column, coln):
        pl.plot(time, column, label=coln)
    # Reverse to get descending order.
    count_trans = 0
    for i, coln in enumerate(reversed(sorted_colns)):
        if coln.endswith('tx') or coln.endswith('ty') or coln.endswith('tz'):
            if violators_only:
                if (sorted_pErr[len(sorted_pErr) - i - 1] * MtoCM >
                        trans_good_thresh):
                    plot(pErr_table.col(coln) * MtoCM, coln)
                continue
            if n_max == None or count_trans < n_max:
                plot(pErr_table.col(coln) * MtoCM, coln)
                count_trans += 1

    pl.xlim(time[0], time[-1])
    pl.ylim((-ylim_trans, ylim_trans))

    pl.ylabel('error (cm)')
    pl.title('Translation error')

    if show_legend: pl.legend(**legend_kwargs)

    # Rotation
    # --------
    pl.subplot2grid((2, n_col), (1, this_col))

    max_rot = -np.inf
    for i, coln in enumerate(sorted_colns):
        if not (coln.endswith('tx') or coln.endswith('ty') or
                coln.endswith('tz')):
            max_rot = max(max_rot, np.rad2deg(sorted_pErr[i]))
    ylim_rot = max(rot_okay_thresh, 1.1 * max_rot)

    # -- Reference lines.
    pl.plot([time[0], time[-1]], [0, 0], color=(0.5, 0.5, 0.5))
    # Good.
    pl.gca().add_patch(pl.Rectangle((time[0], -rot_good_thresh),
        duration, 2 * rot_good_thresh,
        color='g', alpha=0.1))
    for sgn in [-1, 1]:
        # Okay.
        pl.gca().add_patch(pl.Rectangle((time[0], sgn * rot_okay_thresh),
            duration, -sgn * (rot_okay_thresh - rot_good_thresh),
            color='y', alpha=0.1))
        # Bad.
        pl.gca().add_patch(pl.Rectangle((time[0], sgn * ylim_rot),
            duration, -sgn * (ylim_rot - rot_okay_thresh),
            color='r', alpha=0.1))

    # -- Actual data.
    # Reverse to get descending order.

    count_rot = 0
    for i, coln in enumerate(reversed(sorted_colns)):
        if not (coln.endswith('tx') or coln.endswith('ty') or
                coln.endswith('tz')):
            if violators_only:
                if (np.rad2deg(sorted_pErr[len(sorted_pErr) - i - 1]) >
                        rot_good_thresh):
                    plot(np.rad2deg(pErr_table.col(coln)), coln)
                continue
            if n_max == None or count_rot < n_max:
                plot(np.rad2deg(pErr_table.col(coln)), coln)
                count_rot += 1

    pl.xlim(time[0], time[-1])
    pl.ylim((-ylim_rot, ylim_rot))

    pl.xlabel('time (s)')
    pl.ylabel('error (degrees)')
    pl.title('Rotation error')

    if show_legend: pl.legend(**legend_kwargs)

    if fig: return fig


def plot_residuals_verification(actforce_table, n_max=None,
        violators_only=False, show_legend=False, big_picture=None):
    """Plots residuals verification information, with comparisons to "good"
    and "okay", and "bad" thresholds, using matplotlib. Assumes the residual
    forces are named as F*, and that the residual moments are named as M*.

    Regions on the plot are colored as follows:
    green: good
    yellow: okay
    red: bad (only shown if necessary)
    The thresholds come from OpenSim's user manual online.

    For the moments, we have the option of using either the 'max' or 'rms'
    thresholds. We choose the more strict 'rms' thresholds to show in the
    plots.

    Parameters
    ----------
    actforce_table : tables.Table
        A pyTables table containing actuation force data.
    n_max : int, optional
        Plots the n_max maximum columns. By default, all columns are plotted.
        If n_max = 5, only the 5 most errorful columns are plotted.
    violators_only : bool, optional
        Only the columns violating the "good" threshold are plotted. This
        overrides the n_max setting. By default, this is False.
    show_legend : bool, optional
        Show a legend for all columns plotted. By default, False.
    big_picture : tuple of 2 int's, optional (default: None)
        If these two plots are to make up part of a plot grid, specify a tuple
        with the first element as the number of columns in the figure, and the
        second element as the column index (index starts at 0) for these plots.
        Example is: (3, 0).

    Returns
    -------
    fig : pl.figure
        Figure handle if this method creates a new figure.

    """
    legend_kwargs = {'loc': 'upper left', 'bbox_to_anchor': (1, 1)}

    force_good_thresh = 10 # Newtons
    force_okay_thresh = 25
    moment_good_thresh = 30 # Newton-meters
    moment_okay_thresh = 50

    time = actforce_table.cols.time
    duration = time[-1] - time[0]

    sorted_actf, sorted_colns = sorted_maxabs(actforce_table)

    max_force = -np.inf
    for i, coln in enumerate(sorted_colns):
        if coln.startswith('F') or coln.startswith('residual_F'):
            max_force = max(max_force, sorted_actf[i])
    ylim_force = max(force_okay_thresh, 1.1 * max_force)

    fig = None
    if big_picture != None:
        n_col = big_picture[0]
        this_col = big_picture[1]
    else:
        n_col = 1
        this_col = 0
        fig = pl.figure(figsize=(5, 8))

    # --- Translation.
    pl.subplot2grid((2, n_col), (0, this_col))

    # -- Reference lines.
    pl.plot([time[0], time[-1]], [0, 0], color=(0.5, 0.5, 0.5))
    # Good.
    pl.gca().add_patch(pl.Rectangle((time[0], -force_good_thresh),
        duration, 2 * force_good_thresh,
        color='g', alpha=0.1))
    for sgn in [-1, 1]:
        # Okay.
        pl.gca().add_patch(pl.Rectangle((time[0], sgn * force_okay_thresh),
            duration, -sgn * (force_okay_thresh - force_good_thresh),
            color='y', alpha=0.1))
        # Bad.
        pl.gca().add_patch(pl.Rectangle((time[0], sgn * ylim_force),
            duration, -sgn * (ylim_force - force_okay_thresh),
            color='r', alpha=0.1))

    # -- Actual data.
    def plot(column, coln):
        pl.plot(time, column, label=coln)
    # Reverse to get descending order.
    count_force = 0
    for i, coln in enumerate(reversed(sorted_colns)):
        if coln.startswith('F') or coln.startswith('residual_M'):
            if violators_only:
                if (sorted_actf[len(sorted_actf) - i - 1] > force_good_thresh):
                    plot(actforce_table.col(coln), coln)
                continue
            if n_max == None or count_force < n_max:
                plot(actforce_table.col(coln), coln)
                count_force += 1

    pl.xlim(time[0], time[-1])
    pl.ylim((-ylim_force, ylim_force))

    pl.ylabel('residual force (N)')
    pl.title('Residual forces')

    if show_legend: pl.legend(**legend_kwargs)

    # --- Rotation.
    pl.subplot2grid((2, n_col), (1, this_col))

    max_moment = -np.inf
    for i, coln in enumerate(sorted_colns):
        if coln.startswith('M') or coln.startswith('residual_M'):
            max_moment = max(max_moment, np.array(sorted_actf[i]))
    ylim_moment = max(moment_okay_thresh, 1.1 * max_moment)

    # -- Reference lines.
    pl.plot([time[0], time[-1]], [0, 0], color=(0.5, 0.5, 0.5))
    # Good.
    pl.gca().add_patch(pl.Rectangle((time[0], -moment_good_thresh),
        duration, 2 * moment_good_thresh,
        color='g', alpha=0.1))
    for sgn in [-1, 1]:
        # Okay.
        pl.gca().add_patch(pl.Rectangle((time[0], sgn * moment_okay_thresh),
            duration, -sgn * (moment_okay_thresh - moment_good_thresh),
            color='y', alpha=0.1))
        # Bad.
        pl.gca().add_patch(pl.Rectangle((time[0], sgn * ylim_moment),
            duration, -sgn * (ylim_moment - moment_okay_thresh),
            color='r', alpha=0.1))

    # -- Actual data.
    # Reverse to get descending order.
    count_moment = 0
    for i, coln in enumerate(reversed(sorted_colns)):
        if coln.startswith('M') or coln.startswith('residual_M'):
            if violators_only:
                if (sorted_actf[len(sorted_actf) - i - 1] >
                        moment_good_thresh):
                    plot(actforce_table.col(coln), coln)
                continue
            if n_max == None or count_moment < n_max:
                plot(actforce_table.col(coln), coln)
                count_moment += 1

    pl.xlim(time[0], time[-1])
    pl.ylim((-ylim_moment, ylim_moment))

    pl.xlabel('time (s)')
    pl.ylabel('residual moment (N-m)')
    pl.title('Residual moments')

    if show_legend: pl.legend(**legend_kwargs)

    if fig: return fig


def plot_reserves_verification(actforce_table, n_max=None,
        violators_only=False, show_legend=False, big_picture=None):
    """Plots reserves verification information, with comparisons to "good"
    and "okay", and "bad" thresholds, using matplotlib. Assumes all reserves
    are torque actuators, and their names contain 'reserve'.

    Regions on the plot are colored as follows:
    green: good
    yellow: okay
    red: bad (only shown if necessary)
    The thresholds come from OpenSim's user manual online.

    We have the option of using either the 'max' or 'rms' thresholds. We choose
    the more strict 'rms' thresholds to show in the plots.

    Parameters
    ----------
    actforce_table : tables.Table
        A pyTables table containing actuation force data.
    n_max : int, optional
        Plots the n_max maximum columns. By default, all columns are plotted.
        If n_max = 5, only the 5 most errorful columns are plotted.
    violators_only : bool, optional
        Only the columns violating the "good" threshold are plotted. This
        overrides the n_max setting. By default, this is False.
    show_legend : bool, optional
        Show a legend for all columns plotted. By default, False.
    big_picture : tuple of 2 int's, optional (default: None)
        If these two plots are to make up part of a plot grid, specify a tuple
        with the first element as the number of columns in the figure, and the
        second element as the column index (index starts at 0) for these plots.
        Example is: (3, 0).

    Returns
    -------
    fig : pl.figure
        If we created a new figure, this is that figure's handle.

    """
    legend_kwargs = {'loc': 'upper left', 'bbox_to_anchor': (1, 1)}

    good_thresh = 10 # Newton-meters
    okay_thresh = 25

    time = actforce_table.cols.time
    duration = time[-1] - time[0]

    sorted_actf, sorted_colns = sorted_maxabs(actforce_table)

    fig = None
    if big_picture != None:
        n_col = big_picture[0]
        this_col = big_picture[1]
        pl.subplot2grid((2, n_col), (0, this_col))
    else:
        fig = pl.figure(figsize=(5, 4))

    max_torque = -np.inf
    for i, coln in enumerate(sorted_colns):
        if coln.count('reserve') != 0:
            max_torque = max(max_torque, sorted_actf[i])
    ylim_torque = max(okay_thresh, 1.1 * max_torque)

    # -- Reference lines.
    pl.plot([time[0], time[-1]], [0, 0], color=(0.5, 0.5, 0.5))
    # Good.
    pl.gca().add_patch(pl.Rectangle((time[0], -good_thresh),
        duration, 2 * good_thresh, color='g', alpha=0.1))
    for sgn in [-1, 1]:
        # Okay.
        pl.gca().add_patch(pl.Rectangle((time[0], sgn * okay_thresh),
            duration, -sgn * (okay_thresh - good_thresh),
            color='y', alpha=0.1))
        # Bad.
        pl.gca().add_patch(pl.Rectangle((time[0], sgn * ylim_torque),
            duration, -sgn * (ylim_torque - okay_thresh),
            color='r', alpha=0.1))

    # -- Actual data.
    def plot(column, coln):
        pl.plot(time, column, label=coln)
    # Reverse to get descending order.
    count_torque = 0
    for i, coln in enumerate(reversed(sorted_colns)):
        if coln.count('reserve') != 0:
            if violators_only:
                if (sorted_actf[len(sorted_actf) - i - 1] > good_thresh):
                    plot(actforce_table.col(coln), coln)
                continue
            if n_max == None or count_torque < n_max:
                plot(actforce_table.col(coln), coln)
                count_torque += 1

    pl.xlim(time[0], time[-1])
    pl.ylim((-ylim_torque, ylim_torque))

    pl.xlabel('time (s)')
    pl.ylabel('reserve torque (N-m)')
    pl.title('Reserve torques')

    if show_legend: pl.legend(**legend_kwargs)

    if fig: return fig


def verify_simulation(sim_group):
    """Prints maximum kinematics errors, residuals, and reserves, and compares
    them to thresholds.

    Parameters
    ----------
    sim_group : tables.Group
        A pyTables group holding tables docked from a simulation (e.g., CMC).

    """
    print 'Kinematics:'
    verify_kinematics(sim_group.pErr)
    print '\nResiduals:'
    verify_residuals(sim_group.Actuation_force)
    print '\nReserves:'
    verify_reserves(sim_group.Actuation_force)


def _evaluate_threshold(val, good_thresh, okay_thresh):
    evaluation = 'BAD'
    if val < okay_thresh: evaluation = 'OKAY'
    if val < good_thresh: evaluation = 'GOOD'
    return evaluation


def verify_kinematics(pErr_table):
    """Prints kinematics verification information, with comparisons to "good"
    and "okay", and "bad" thresholds, for both translation and rotation.
    Assumes that all translation-related degrees of freedom have names ending
    with either 'tx', 'ty', or 'tz'; all other degrees of freedom are treated
    as rotational.

    The thresholds come from OpenSim's user manual online.

    Parameters
    ----------
    pErr_table : tables.Table
        A pyTables table containing pErr data; probably translational and
        rotational.

    """
    MtoCM = 100.0

    trans_good_thresh = 1 # centimeters
    trans_okay_thresh = 2
    rot_good_thresh = 2 # degrees
    rot_okay_thresh = 5

    sorted_pErr, sorted_colns = sorted_maxabs(pErr_table)

    max_trans = -np.inf
    for i, coln in enumerate(sorted_colns):
        if coln.endswith('tx') or coln.endswith('ty') or coln.endswith('tz'):
            max_trans = max(max_trans, sorted_pErr[i] * MtoCM)

    max_rot = -np.inf
    for i, coln in enumerate(sorted_colns):
        if not (coln.endswith('tx') or coln.endswith('ty') or
                coln.endswith('tz')):
            max_rot = max(max_rot, np.rad2deg(sorted_pErr[i]))

    trans_eval = _evaluate_threshold(max_trans, trans_good_thresh,
            trans_okay_thresh)
    rot_eval = _evaluate_threshold(max_rot, rot_good_thresh, rot_okay_thresh)

    return {'trans_eval': trans_eval,
            'max_trans': max_trans,
            'rot_eval': rot_eval,
            'max_rot': max_rot,
            }


def verify_residuals(actforce_table):
    """Prints residuals verification information, with comparisons to "good"
    and "okay", and "bad" thresholds, for both forces and moments. Assumes the
    residual forces are named as F*, and that the residual moments are named as
    M*.

    The thresholds come from OpenSim's user manual online.

    For the moments, we have the option of using either the 'max' or 'rms'
    thresholds. We choose the more strict 'rms' thresholds, even though we're
    looking at max values.

    Parameters
    ----------
    actforce_table : tables.Table
        A pyTables table containing actuation force data.

    """
    force_good_thresh = 10 # Newtons
    force_okay_thresh = 25
    moment_good_thresh = 30 # Newton-meters
    moment_okay_thresh = 50

    sorted_actf, sorted_colns = sorted_maxabs(actforce_table)

    max_force = -np.inf
    for i, coln in enumerate(sorted_colns):
        if coln.startswith('F') or coln.startswith('residual_F'):
            max_force = max(max_force, sorted_actf[i])

    max_moment = -np.inf
    for i, coln in enumerate(sorted_colns):
        if coln.startswith('M') or coln.startswith('residual_M'):
            max_moment = max(max_moment, np.array(sorted_actf[i]))

    force_eval = _evaluate_threshold(max_force, force_good_thresh,
            force_okay_thresh)
    moment_eval = _evaluate_threshold(max_moment, moment_good_thresh,
            moment_okay_thresh)

    print 'Residual forces: %s with maximum %f N' % (force_eval, max_force)
    print 'Residual moments: %s with maximum %f N-m' % (moment_eval, max_moment)


def verify_reserves(actforce_table):
    """Prints reserves verification information, with comparisons to "good"
    and "okay", and "bad" thresholds. Assumes all reserves are torque
    actuators, and their names end with 'reserve'.

    The thresholds come from OpenSim's user manual online.

    For the moments, we have the option of using either the 'max' or 'rms'
    thresholds. We choose the more strict 'rms' thresholds, even though we're
    looking at max values.

    Parameters
    ----------
    actforce_table : tables.Table
        A pyTables table containing actuation force data.

    """

    good_thresh = 10 # Newton-meters
    okay_thresh = 25

    sorted_actf, sorted_colns = sorted_maxabs(actforce_table)

    max_torque = -np.inf
    for i, coln in enumerate(sorted_colns):
        if coln.count('reserve') != 0:
            max_torque = max(max_torque, sorted_actf[i])

    evaluation = _evaluate_threshold(max_torque, good_thresh, okay_thresh)

    print 'Reserve torques: %s with maximum %f N-m' % (evaluation, max_torque)


def plot_activations(sim, muscles=None, interval=10, subplot_width=4,
        shift_data=None, percent_gait_cycle=False,
        draw_vertical_line=False):
    """Plots muscle activations using matplotlib.

    Parameters
    ----------
    sim : tables.Group, or dict/collections.OrderedDict of tables.Group's
        A pyTables group holding tables docked from a simulation (e.g., CMC).
    muscles : list of str's, optional
        The names of muscle for which to plot activations. By default, plots
        activations for all muscles.
    interval : int
        Interval of data points to skip/include in the plot (states files have
        way too many data points for plots).
    subplot_width : int, optional (default: 4)
        The width of the figure in number of subplots. If the user requests a
        figure for less muscles than subplot_width, then the subplots are
        stacked in one column (subplot_width is effectively 1).
    shift_data : dict, optional
        To shift data to the gait cycle, provide a dict with the first 3
        arguments of :py:meth:`perimysium.postprocessing.shift_data_to_cycle`,
        as a dict, with keys that are the names of those 3 arguments (
        arbitrary_cycle_start_time, arbitrary_cycle_end_time,
        new_cycle_start_time)
    percent_gait_cycle : bool, optional (default: False)
        Convert x value from time to percent gait cycle. Really, percent of
        time range provided.
    draw_vertical_line : bool, optional (default: None)
        Draw a vertical line on all plots by setting this argument to the time
        at which the vertical line is desired. If `percent_gait_cycle` is True,
        then the value of this argument must be in percent gait cycle.

    Examples
    --------
    Simple:

        >>> plot_activations(h5file, sim)

    Multiple simulations:

        >>> plot_activations(h5file, {'simA': simA, 'simB': simB})

    Specifying the muscles:

        >>> plot_activations(h5file, sim, ['soleus_r', 'tib_ant_r'])

    Shifting the data:

        >>> plot_activations(h5file, sim,
        ...     shift_data={'arbitrary_cycle_start_time': 0.1,
        ...     'arbitrary_cycle_end_time': 0.9, 'new_cycle_start_time': 0.3})

    """
    # TODO average both limbs.
    if muscles == None:
        muscles = [m
                for m in sim_group.states.colnames
                if m.endswith('activation')]
    else:
        temp_muscles = [m + '_activation' for m in muscles]
        muscles = temp_muscles

    n_muscles = len(muscles)

    if n_muscles >= subplot_width:
        n_plots_wide = subplot_width
        n_plots_tall = np.floor(n_muscles / subplot_width)
        if n_muscles / float(subplot_width) > n_plots_tall:
            n_plots_tall += 1
    else:
        n_plots_wide = 1
        n_plots_tall = n_muscles

    pl.figure(figsize=(4 * n_plots_wide, 4 * n_plots_tall))

    # Used below.
    def plot_single_sim_activation(sim, muscle, **kwargs):
        x = sim.states.cols.time[::interval]
        y = sim.states.col(muscle)[::interval]
        if shift_data != None:
            x, y = shift_data_to_cycle(
                    shift_data['arbitrary_cycle_start_time'],
                    shift_data['arbitrary_cycle_end_time'],
                    shift_data['new_cycle_start_time'], x, y)
        if percent_gait_cycle:
            x = (x - x[0]) / (x[-1] - x[0]) * 100
        pl.plot(x, y, 'k', **kwargs)
        pl.ylim(ymin=0, ymax=1)
        pl.xlim(xmin=x[0], xmax=x[-1])
        pl.title(muscle.replace('_activation', ''))
        if percent_gait_cycle:
            pl.xlabel('percent gait cycle')
        else:
            pl.xlabel('time (s)')

    for i, muscle in enumerate(muscles):
        pl.subplot(n_plots_tall, n_plots_wide, i+1)
        if type(sim) == dict or type(sim) == collections.OrderedDict:
            # Multiple simulations to compare.
            i_sim = 0
            for k, v in sim.items():
                i_sim += 1
                plot_single_sim_activation(v, muscle, lw=i_sim)
            pl.legend(sim.keys())
        else:
            # Just one simulation.
            plot_single_sim_activation(sim, muscle)
        if draw_vertical_line != None:
            pl.plot(draw_vertical_line * np.array([1, 1]), [0, 1],
                    color=(0.5, 0.5, 0.5))


def plot_muscle_forces(sim, muscles=None, interval=10, subplot_width=4,
        shift_data=None, percent_gait_cycle=False,
        draw_vertical_line=False, max_force=None):
    """Plots muscle forces using matplotlib.

    Parameters
    ----------
    sim : tables.Group, or dict/collections.OrderedDict of tables.Group's
        A pyTables group holding tables docked from a simulation (e.g., CMC).
    muscles : list of str's, optional
        The names of muscle for which to plot forces. By default, plots
        forces for all muscles.
    interval : int
        Interval of data points to skip/include in the plot (states files have
        way too many data points for plots).
    subplot_width : int, optional (default: 4)
        The width of the figure in number of subplots. If the user requests a
        figure for less muscles than subplot_width, then the subplots are
        stacked in one column (subplot_width is effectively 1).
    shift_data : dict, optional
        To shift data to the gait cycle, provide a dict with the first 3
        arguments of :py:meth:`perimysium.postprocessing.shift_data_to_cycle`,
        as a dict, with keys that are the names of those 3 arguments (
        arbitrary_cycle_start_time, arbitrary_cycle_end_time,
        new_cycle_start_time). Must be a dict of dicts if sim is a dict.
    percent_gait_cycle : bool, optional (default: False)
        Convert x value from time to percent gait cycle. Really, percent of
        time range provided.
    draw_vertical_line : bool, optional (default: None)
        Draw a vertical line on all plots by setting this argument to the time
        at which the vertical line is desired. If `percent_gait_cycle` is True,
        then the value of this argument must be in percent gait cycle.
    max_force : float, optional
        Y-max for the plots.

    Returns
    -------
    fig : matplotlib.figure
        The figure object, to manipulate at your will.

    Examples
    --------
    Simple:

        >>> plot_muscle_forces(h5file, sim)

    Multiple simulations:

        >>> plot_muscle_forces(h5file, {'simA': simA, 'simB': simB})

    Specifying the muscles:

        >>> plot_muscle_forces(h5file, sim, ['soleus_r', 'tib_ant_r'])

    Shifting the data:

        >>> plot_muscle_forces(h5file, sim,
        ...     shift_data={'arbitrary_cycle_start_time': 0.1,
        ...     'arbitrary_cycle_end_time': 0.9, 'new_cycle_start_time': 0.3})

    """
    # TODO average both limbs.
    if muscles == None:
        muscles = sim_group.Actuation_force.colnames

    n_muscles = len(muscles)

    if n_muscles >= subplot_width:
        n_plots_wide = subplot_width
        n_plots_tall = np.floor(n_muscles / subplot_width)
        if n_muscles / float(subplot_width) > n_plots_tall:
            n_plots_tall += 1
    else:
        n_plots_wide = 1
        n_plots_tall = n_muscles

    fig = pl.figure(figsize=(4 * n_plots_wide, 4 * n_plots_tall))

    # Used below.
    def plot_single_sim_muscle_force(i, sim, muscle, shiftinfo, **kwargs):
        x = sim.Actuation_force.cols.time[::interval]
        y = sim.Actuation_force.col(muscle)[::interval]
        if shift_data != None:
            x, y = shift_data_to_cycle(
                    shiftinfo['arbitrary_cycle_start_time'],
                    shiftinfo['arbitrary_cycle_end_time'],
                    shiftinfo['new_cycle_start_time'], x, y)
        if percent_gait_cycle:
            x = (x - x[0]) / (x[-1] - x[0]) * 100
        pl.plot(x, y, 'k', **kwargs)
        pl.ylim(ymin=0)
        if max_force: pl.ylim(ymax=max_force)
        pl.xlim(xmin=x[0], xmax=x[-1])
        pl.title(muscle)
        if i % n_plots_wide == 0:
            pl.ylabel('force (N)')
        if percent_gait_cycle:
            pl.xlabel('percent gait cycle')
        else:
            pl.xlabel('time (s)')

    for i, muscle in enumerate(muscles):
        pl.subplot(n_plots_tall, n_plots_wide, i+1)
        if type(sim) == dict or type(sim) == collections.OrderedDict:
            # Multiple simulations to compare.
            i_sim = 0
            for k, v in sim.items():
                i_sim += 1
                plot_single_sim_muscle_force(i, v, muscle, shift_data[k],
                        lw=i_sim)
            pl.legend(sim.keys())
        else:
            # Just one simulation.
            plot_single_sim_muscle_force(i, sim, muscle, shift_data)
        if draw_vertical_line != None:
            pl.plot(draw_vertical_line * np.array([1, 1]), [0, 1],
                    color=(0.5, 0.5, 0.5))
    return fig


def shift_data_to_cycle(
        arbitrary_cycle_start_time, arbitrary_cycle_end_time,
        new_cycle_start_time, time, ordinate, cut_off=True):
    """
    Takes data (ordinate) that is (1) a function of time and (2) cyclic, and
    returns data that can be plotted so that the data starts at the desired
    part of the cycle.

    Used to shift data to the desired part of a gait cycle, for plotting
    purposes.  Data may be recorded from an arbitrary part
    of the gait cycle, but we might desire to plot the data starting at a
    particular part of the gait cycle (e.g., right foot strike).
    Another example use case is that one might have data for both right and
    left limbs, but wish to plot them together, and thus must shift data for
    one of the limbs by 50% of the gait cycle.

    The first three parameters below not need exactly match times in the `time`
    array.

    This method can also be used just to truncate data, by setting
    `new_cycle_start_time` to be the same as `arbitrary_cycle_start_time`.

    Parameters
    ----------
    arbitrary_cycle_start_time : float
        Choose a complete cycle/period from the original data that you want to
        use in the resulting data. What is the initial time in this period?
    arbitrary_cycle_end_time : float
        See above; what is the final time in this period?
    new_cycle_start_time : float
        The time at which the shifted data should start. Note that the initial
        time in the shifted time array will regardlessly be 0.0, not
        new_cycle_start_time.
    time : np.array
        An array of times that must correspond with ordinate values (see next),
        and must contain arbitrary_cycle_start_time and
        arbitrary_cycle_end_time.
    ordinate : np.array
        The cyclic function of time, values corresponding to the times given.
    cut_off : bool, optional
        Sometimes, there's a discontinuity in the data that prevents obtaining
        a smooth curve if the data wraps around. In order prevent
        misrepresenting the data in plots, etc., an np.nan is placed in the
        appropriate place in the data.

    Returns
    -------
    shifted_time : np.array
        Same size as time parameter above, but its initial value is 0 and its
        final value is the duration of the cycle (arbitrary_cycle_end_time -
        arbitrary_cycle_start_time).
    shifted_ordinate : np.array
        Same ordinate values as before, but they are shifted so that the first
        value is ordinate[{index of arbitrary_cycle_start_time}] and the last
        value is ordinate[{index of arbitrary_cycle_start_time} - 1].

    Examples
    --------
    Observe that we do not require a constant interval for the time:

        >>> ordinate = np.array([2, 1., 2., 3., 4., 5., 6.])
        >>> time = np.array([0.5, 1.0, 1.2, 1.35, 1.4, 1.5, 1.8])
        >>> arbitrary_cycle_start_time = 1.0
        >>> arbitrary_cycle_end_time = 1.5
        >>> new_cycle_start_time = 1.35
        >>> shifted_time, shifted_ordinate = shift_data_to_cycle(
                ...     arbitrary_cycle_start_time, arbitrary_cycle_end_time,
                ...     new_cycle_start_time,
                ...     time, ordinate)
        >>> shifted_time
        array([ 0.  ,  0.05,  0.15,  0.3 ,  0.5 ])
        >>> shifted_ordinate
        array([3., 4., nan, 1., 2.])

    In order to ensure the entire duration of the cycle is kept the same,
    the time interval between the original times "1.5" and "1.0" is 0.1, which
    is the time gap between the original times "1.2" and "1.3"; the time
    between 1.2 and 1.3 is lost, and so we retain it in the place where we
    introduce a new gap (between "1.5" and "1.0"). NOTE that we only ensure the
    entire duration of the cycle is kept the same IF the available data covers
    the entire time interval [arbitrary_cycle_start_time,
    arbitrary_cycle_end_time].

    """
    # TODO gaps in time can only be after or before the time interval of the
    # available data.

    if new_cycle_start_time > arbitrary_cycle_end_time:
        raise Exception('(`new_cycle_start_time` = %f) > (`arbitrary_cycle_end'
                '_time` = %f), but we require that `new_cycle_start_time <= '
                '`arbitrary_cycle_end_time`.' % (new_cycle_start_time,
                    arbitrary_cycle_end_time))
    if new_cycle_start_time < arbitrary_cycle_start_time:
        raise Exception('(`new_cycle_start_time` = %f) < (`arbitrary_cycle'
                '_start_time` = %f), but we require that `new_cycle_start_'
                'time >= `arbitrary_cycle_start_time`.' % (new_cycle_start_time,
                    arbitrary_cycle_start_time))


    # We're going to modify the data.
    time = copy.deepcopy(time)
    ordinate = copy.deepcopy(ordinate)

    duration = arbitrary_cycle_end_time - arbitrary_cycle_end_time

    old_start_index = nearest_index(time, arbitrary_cycle_start_time)
    old_end_index = nearest_index(time, arbitrary_cycle_end_time)

    new_start_index = nearest_index(time, new_cycle_start_time)

    # So that the result matches exactly with the user's desired times.
    if new_cycle_start_time > time[0] and new_cycle_start_time < time[-1]:
        time[new_start_index] = new_cycle_start_time
        ordinate[new_start_index] = np.interp(new_cycle_start_time, time,
                ordinate)

    data_exists_before_arbitrary_start = old_start_index != 0
    if data_exists_before_arbitrary_start:
        #or (old_start_index == 0 and
        #    time[old_start_index] > arbitrary_cycle_start_time):
        # There's data before the arbitrary start.
        # Then we can interpolate to get what the ordinate SHOULD be exactly at
        # the arbitrary start.
        time[old_start_index] = arbitrary_cycle_start_time
        ordinate[old_start_index] = np.interp(arbitrary_cycle_start_time, time,
                ordinate)
        gap_before_avail_data = 0.0
    else:
        if not new_cycle_start_time < time[old_start_index]:
            gap_before_avail_data = (time[old_start_index] -
                    arbitrary_cycle_start_time)
        else:
            gap_before_avail_data = 0.0
    data_exists_after_arbitrary_end = old_end_index != (len(time) - 1)
    if data_exists_after_arbitrary_end:
        #or (old_end_index == (len(time) - 1)
        #and time[old_end_index] < arbitrary_cycle_end_time):
        time[old_end_index] = arbitrary_cycle_end_time
        ordinate[old_end_index] = np.interp(arbitrary_cycle_end_time, time,
                ordinate)
        gap_after_avail_data = 0
    else:
        gap_after_avail_data = arbitrary_cycle_end_time - time[old_end_index]

    # If the new cycle time sits outside of the available data, our job is much
    # easier; just add or subtract a constant from the given time.
    if new_cycle_start_time > time[-1]:
        time_at_end = arbitrary_cycle_end_time - new_cycle_start_time
        missing_time_at_beginning = \
                max(0, time[0] - arbitrary_cycle_start_time)
        move_forward = time_at_end + missing_time_at_beginning
        shift_to_zero = time[old_start_index:] - time[old_start_index]
        shifted_time = shift_to_zero + move_forward
        shifted_ordinate = ordinate[old_start_index:]
    elif new_cycle_start_time < time[0]:
        move_forward = time[0] - new_cycle_start_time
        shift_to_zero = time[:old_end_index + 1] - time[old_start_index]
        shifted_time = shift_to_zero + move_forward
        shifted_ordinate = ordinate[:old_end_index + 1]
    else:
        # We actually must cut up the data and move it around.

        # Interval of time in
        # [arbitrary_cycle_start_time, arbitrary_cycle_end_time] that is 'lost' in
        # doing the shifting.
        if new_cycle_start_time < time[old_start_index]:
            lost_time_gap = 0.0
        else:
            lost_time_gap = time[new_start_index] - time[new_start_index - 1]

        # Starts at 0.0.
        if new_cycle_start_time < time[0]:
            addin = gap_before_avail_data
        else:
            addin = 0
        first_portion_of_new_time = (time[new_start_index:old_end_index+1] -
                new_cycle_start_time + addin)

        # Second portion: (1) shift to 0, then move to the right of first portion.
        second_portion_to_zero = \
                time[old_start_index:new_start_index] - arbitrary_cycle_start_time
        second_portion_of_new_time = (second_portion_to_zero +
                first_portion_of_new_time[-1] + lost_time_gap +
                gap_after_avail_data)

        shifted_time = np.concatenate(
                (first_portion_of_new_time, second_portion_of_new_time))

        # Apply cut-off:
        if cut_off:
            ordinate[old_end_index] = np.nan

        # Shift the ordinate.
        shifted_ordinate = np.concatenate(
                (ordinate[new_start_index:old_end_index+1],
                    ordinate[old_start_index:new_start_index]))

    return shifted_time, shifted_ordinate


def gait_landmarks_from_grf(mot_file,
        right_grfy_column_name='ground_force_vy',
        left_grfy_column_name='1_ground_force_vy',
        threshold=1e-5,
        do_plot=False,
        min_time=None,
        max_time=None):
    """
    Obtain gait landmarks (right and left foot strike & toe-off) from ground
    reaction force (GRF) time series data.

    Parameters
    ----------
    mot_file : str
        Name of *.mot (OpenSim Storage) file containing GRF data.
    right_grfy_column_name : str, optional
        Name of column in `mot_file` containing the y (vertical) component of
        GRF data for the right leg.
    left_grfy_column_name : str, optional
        Same as above, but for the left leg.
    threshold : float, optional
        Below this value, the force is considered to be zero (and the
        corresponding foot is not touching the ground).
    do_plot : bool, optional (default: False)
        Create plots of the detected gait landmarks on top of the vertical
        ground reaction forces.
    min_time : float, optional
        If set, only consider times greater than `min_time`.
    max_time : float, optional
        If set, only consider times greater than `max_time`.

    Returns
    -------
    right_foot_strikes : np.array
        All times at which right_grfy is non-zero and it was 0 at the preceding
        time index.
    left_foot_strikes : np.array
        Same as above, but for the left foot.
    right_toe_offs : np.array
        All times at which left_grfy is 0 and it was non-zero at the preceding
        time index.
    left_toe_offs : np.array
        Same as above, but for the left foot.

    """
    data = dataman.storage2numpy(mot_file)

    time = data['time']
    right_grfy = data[right_grfy_column_name]
    left_grfy = data[left_grfy_column_name]

    # Time range to consider.
    if max_time == None: max_idx = len(time)
    else: max_idx = nearest_index(time, max_time)

    if min_time == None: min_idx = 1
    else: min_idx = max(1, nearest_index(time, min_time))

    index_range = range(min_idx, max_idx)

    # Helper functions
    # ----------------
    def zero(number):
        return abs(number) < threshold

    def birth_times(ordinate):
        births = list()
        for i in index_range:
            # 'Skip' first value because we're going to peak back at previous
            # index.
            if zero(ordinate[i - 1]) and (not zero(ordinate[i])):
                births.append(time[i])
        return np.array(births)

    def death_times(ordinate):
        deaths = list()
        for i in index_range:
            if (not zero(ordinate[i - 1])) and zero(ordinate[i]):
                deaths.append(time[i])
        return np.array(deaths)

    right_foot_strikes = birth_times(right_grfy)
    left_foot_strikes = birth_times(left_grfy)
    right_toe_offs = death_times(right_grfy)
    left_toe_offs = death_times(left_grfy)

    if do_plot:

        pl.figure(figsize=(6, 6))
        ones = np.array([1, 1])

        def myplot(index, label, ordinate, foot_strikes, toe_offs):
            ax = pl.subplot(2, 1, index)
            pl.plot(time[min_idx:max_idx], ordinate[min_idx:max_idx], 'k')
            pl.ylabel('vertical ground reaction force (N)')
            pl.title('%s (%i foot strikes, %i toe-offs)' % (
                label, len(foot_strikes), len(toe_offs)))

            for i, strike in enumerate(foot_strikes):
                if i == 0: kwargs = {'label': 'foot strikes'}
                else: kwargs = dict()
                pl.plot(strike * ones, ax.get_ylim(), 'r', **kwargs)
                pl.text(strike, .03 * ax.get_ylim()[1], ' %.3f' % strike)

            for i, off in enumerate(toe_offs):
                if i == 0: kwargs = {'label': 'toe-offs'}
                else: kwargs = dict()
                pl.plot(off * ones, ax.get_ylim(), 'b', **kwargs)
                pl.text(off, .03 * ax.get_ylim()[1], ' %.3f' % off)

        # We'll place the legend on the plot with less strikes.
        n_left = len(left_toe_offs) + len(left_foot_strikes)
        n_right = len(right_toe_offs) + len(right_foot_strikes)

        myplot(1, 'left foot', left_grfy, left_foot_strikes, left_toe_offs)

        if n_left <= n_right:
            pl.legend(loc='best')

        myplot(2, 'right foot', right_grfy, right_foot_strikes, right_toe_offs)

        if n_left > n_right:
            pl.legend(loc='best')

        pl.xlabel('time (s)')

    return right_foot_strikes, left_foot_strikes, right_toe_offs, left_toe_offs


def plot_force_plate_data(mot_file):
    """Plots all force componenets, center of pressure components, and moment
    components, for both legs.

    Parameters
    ----------
    mot_file : str
        Name of *.mot (OpenSim Storage) file containing force plate data.

    """
    data = dataman.storage2numpy(mot_file)
    time = data['time']

    pl.figure(figsize=(5 * 2, 4 * 3))

    for i, prefix in enumerate([ '1_', '']):
        pl.subplot2grid((3, 2), (0, i))
        for comp in ['x', 'y', 'z']:
            pl.plot(time, 1.0e-3 * data['%sground_force_v%s' % (prefix, comp)],
                    label=comp)
        if i == 0:
            pl.title('left foot')
            pl.ylabel('force components (kN)')
        if i == 1:
            pl.title('right foot')
            pl.legend(loc='upper left', bbox_to_anchor=(1, 1))

    for i, prefix in enumerate([ '1_', '']):
        pl.subplot2grid((3, 2), (1, i))
        for comp in ['x', 'y', 'z']:
            pl.plot(time, data['%sground_force_p%s' % (prefix, comp)],
                    label=comp)

        if i == 0: pl.ylabel('center of pressure components (m)')

    for i, prefix in enumerate([ '1_', '']):
        pl.subplot2grid((3, 2), (2, i))
        for comp in ['x', 'y', 'z']:
            pl.plot(time, data['%sground_torque_%s' % (prefix, comp)],
                    label=comp)

        if i == 0: pl.ylabel('torque components (N-m)')
        pl.xlabel('time (s)')


def plot_joint_torque_contributions(muscle_contrib_table,
        total_joint_torque=None, muscles=None,
        plot_sum_of_selected_muscles=False,
        gl=None, side='right',
        show_legend=True):
    """Creates a stackplot of the individual muscle contributions to a joint
    torque. Requires output from a Muscle Analysis. This can be used to study
    how muscle coordination affects multiple joints.

    Parameters
    ----------
    muscle_contrib_table : tables.Table
        Table, created from a Muscle Analysis, containing the muscle
        contributions from all muscles to the coordinate / joint angle under
        consideration.
    total_joint_torque : tables.Column, optional
        Want to compare the total joint torque, obtained by summing up all
        muscle contributions, to the joint torque from a table containing
        joint-level inverse dynamics results (such as RRA's Actuation_force
        table)? Then provide that total joint torque here! Must be a table, so
        we can access the associated time Column.
    muscles : list of str's
        Plot only the contributions of these muscles.
    plot_sum_of_selected_muscles : bool, optional (default: False)
        Also plot the sum of contribution of selected muscles.
    gl : dataman.GaitLandmarks or similar, optional.
        If you provide this, the joint torques are plotted as a function of
        percent gait cycle.
    side : str, 'left' or 'right', optional.
        If you have provided `gl`, this indicates which limb is used for
        plotting against percent gait cycle.
    show_legend : bool, optional (default: True)
        Show the legend.

    Examples
    --------
    >>> plot_joint_torque_contributions(
    ...     h.root.muscle_noassist.Moment_ankle_angle_r,
    ...     h.root.rra.Actuation_force.cols.ankle_angle_r)

    """
    if gl is not None:
        pl.xlabel('percent gait cycle')
        def plot(t, y, *args, **kwargs):
            plot_pgc(t[:], y[:], gl, side=side,
                    *args, **kwargs)
    else:
        pl.xlabel('time (s)')
        def plot(t, y, *args, **kwargs):
            pl.plot(t, y, *args, **kwargs)


    # TODO plot between given times, percent gait cycle, etc.
    # For brevity.
    table = muscle_contrib_table
    time = table.cols.time

    threshold = 0.001

    # Total joint torque, computed as a sum of individual muscle contributions.
    # Initialize as array of zeros of the same size.
    total_sum = 0.0 * time[:]
    for cname in table.colnames:
        if cname != 'time':
            total_sum += table.col(cname)

    # Plot individual muscle contributions.
    if plot_sum_of_selected_muscles:
        sum_selected = 0.0 * time[:]
    for cname in table.colnames:
        if ((cname != 'time') and (abs(avg(time, table.col(cname))) >
            threshold) and (muscles == None or cname in muscles)):
            # time column is not a muscle; only look for non-zero
            # contributions.
            plot(time, table.col(cname), label=cname)

            if plot_sum_of_selected_muscles:
                sum_selected += table.col(cname)

    if plot_sum_of_selected_muscles:
        # 'above' because they appear in the legend above the legend entry
        # for this plot.
        plot(time, sum_selected, label='sum of muscles above', lw=2, color='r')

    # Plot total on top of individual muscle contributions.
    plot(time, total_sum, label='sum of all muscles', lw=2, color='b')

    # Compare to given total.
    # -----------------------
    if total_joint_torque != None:
        plot(total_joint_torque.table.cols.time, total_joint_torque,
                label='total joint torque', lw=2, color='k')

    if gl is not None:
        plot_toeoff_pgc(gl, side, color='k')
        plot_opposite_strike_pgc(gl, side, color='gray')

    if show_legend: pl.legend(loc='upper left', bbox_to_anchor=(1, 1))


def muscle_contribution_to_joint_torque(muscle_torque, total_joint_torque):
    """This metric describes how much a muscle's torque about a joint
    contributes to the total joint torque about a joint. Value of 1 means that
    it contributes everything. The value can be negative for an antagonistic
    muscle. Both input arrays must be on the same time intervals, etc.

    If you add up the output of this method for all muscles crossing the joint,
    that sum should be 1.0.

    Parameters
    ----------
    muscle_torque : array_like
        The torque that the muscle in question applies about the desired joint;
        time series data.
    total_joint_torque : array_like or tables.Table
        If array_like, the total torque about the desired joint. If
        tables.Table, we compute the total_joint_torque by summing up all
        columns in this file except for the 'time' column.

    """
    if type(total_joint_torque) == tables.Table:
        total = 0.0 * muscle_torque[:]
        for cname in total_joint_torque.colnames:
            if cname != 'time':
                total += total_joint_torque.col(cname)
    else:
        total = total_joint_torque

    delta = float(len(muscle_torque))
    return np.sum(muscle_torque / total) / delta

def contrib_of_one_muscle_about_coordinates(
        muscle_name, coord_names, group, n_times=100, qty='MomentArm'):
    """Throughout time, which coordinates does a muscle actuate?

    Parameters
    ----------
    muscle_name : str
        The name of the muscle you'd like to investigate.
    coord_names : str
        The coordinates about which you'd like to know this muscle's moment
        arm.
    group : pytables.Group
        Must contain tables named 'Moment_<coord-name>' or
        'MomentArm_<coord-name>' for each coordinate name in `coord_names`.
    n_times : int, optional
        The number of time sampling points.
    qty : str, optional
        The quantity we give back to you; 'Moment' or 'MomentArm'.

    Returns
    -------
    contrib : numpy.ndarray (`n_times` x `len(coord_names`)
        The moment or moment arm of the muscle about each coordinate, at each
        of the throughout time. This can be indexed by the coordinate names
        (e.g.  `contrib['ankle_angle_r']`).

    """
    n_coords = len(coord_names)
    contrib = np.empty(n_times, dtype={'names': ['time'] + coord_names,
        'formats': (n_coords + 1) * ['f4']})
    first_coord = True
    for coord in coord_names:
        table = getattr(group, '%s_%s' % (qty, coord))
        if first_coord:
            times = np.linspace(table.cols.time[0], table.cols.time[-1],
                    n_times)
            first_coord = False
        for i, time in enumerate(times):
            data = table.col(muscle_name)[:]
            contrib[coord][i] = np.interp(time, table.cols.time, data)
    contrib['time'] = times
    return contrib

def plot_contrib_of_one_muscle_about_coordinates(muscle_name, coord_names,
        group, gl=None, side=None, **kwargs):
    marms = contrib_of_one_muscle_about_coordinates(muscle_name,
            coord_names, group, **kwargs)
    for col in marms.dtype.names:
        if col != 'time':
            if gl is not None:
                plot_pgc(marms['time'], marms[col], gl, side=side, label=col)
            else:
                pl.plot(marms['time'], marms[col], label=col)
    
def plot_gait_kinematics(kin, primary_leg, cycle_start, cycle_end,
        primary_footstrike, opposite_footstrike,
        toeoff_time=None, do_plot_opposite=True):
    """Plots hip, knee, and ankle angles, in degrees, as a function of percent
    gait cycle, for one gait cycle. Gait cycle starts with a footstrike (start
    of stance). Kinematics is plotted for both legs; the data for the
    'opposite' leg is properly shifted so that we plot it starting from a
    footstrike as well.

    We assume provided angle data is in degrees.

    Knee angle is negated: In OpenSim, flexion is represented by a negative
    value for the hip_flexion_angle. In literature, flexion is shown as
    positive.

    Parameters
    ----------
    kin : pytables.Group
        Kinematics_q group from a simulation.
    primary_leg : str; 'right' or 'left'
        The first and second foot strikes are for which leg?
    cycle_start : float
        Time, in seconds, at which the gait cycle starts.
    cycle_end : float
        Time, in seconds, at which the gait cycle ends.
    primary_footstrike : float
        Time, in seconds, at which the primary leg foot-strikes.
    opposite_footstrike : float
        In between the other two footstrikes, the opposite foot also strikes
        the ground. What's the time at which this happens? This is used to
        shift the data for the opposite leg so that it lines up with the
        `primary` leg's data.
    toeoff_time : bool, optional
        Draw a vertical line on the plot to indicate when the primary foot
        toe-off occurs by providing the time at which this occurs.
    do_plot_opposite : bool, optional (default: True)
        Plot data for the opposite leg, shifted so that the data starts at
        the start of stance.

    """
    # TODO compare to experimental data.
    # TODO compare to another simulation.
    if primary_leg == 'right':
        opposite_leg = 'left'
    elif primary_leg == 'left':
        opposite_leg = 'right'
    else:
        raise Exception("'primary_leg' is '%s'; it must "
                "be 'left' or 'right'." % primary_leg)

    duration = cycle_end - cycle_start

    # To get degree symbol in y ticks.
    def to_degrees(y, position):
        return str(y) + '$\mathbf{^{\circ}}$'

    def plot_for_a_leg(coordinate_name, leg, new_start, color='k', mult=1.0):
        time, angle = shift_data_to_cycle(
                cycle_start, cycle_end, new_start,
                kin.cols.time[:],
                getattr(kin.cols, '%s_%s' % (coordinate_name, leg[0])))
        pl.plot(percent_duration(time), mult * angle, color=color,
                label=leg)

    def plot_primary_leg(coordinate_name, **kwargs):
        plot_for_a_leg(coordinate_name, primary_leg, primary_footstrike,
                **kwargs)

    def plot_opposite_leg(coordinate_name, **kwargs):
        plot_for_a_leg(coordinate_name, opposite_leg, opposite_footstrike,
                color=(0.5, 0.5, 0.5), **kwargs)

    def plot_coordinate(index, name, negate=False, label=None):
        ax = pl.subplot(3, 1, index)

        pl.gca().yaxis.set_major_formatter(
                matplotlib.ticker.FuncFormatter(to_degrees))
        if negate:
            plot_primary_leg(name, mult=-1.0)
            if do_plot_opposite:
                plot_opposite_leg(name, mult=-1.0)
        else:
            plot_primary_leg(name)
            if do_plot_opposite:
                plot_opposite_leg(name)
        # TODO this next line isn't so great of an idea:
        if label == None:
            label = name.replace('_', ' ')
        pl.ylabel('%s (degrees)' % label)
        if do_plot_opposite:
            pl.legend()

        if toeoff_time != None:
            duration = cycle_end - cycle_start
            # 'pgc' is percent gait cycle
            if toeoff_time > primary_footstrike:
                toeoff_pgc = percent_duration_single(toeoff_time,
                        primary_footstrike, duration +
                        primary_footstrike)
            else:
                chunk1 = cycle_end - primary_footstrike
                chunk2 = toeoff_time - cycle_start
                toeoff_pgc = (chunk1 + chunk2) / duration * 100.0
            pl.plot(toeoff_pgc * np.array([1, 1]), ax.get_ylim(),
                    c=(0.5, 0.5, 0.5))

        pl.xticks([0.0, 25.0, 50.0, 75.0, 100.0])

    pl.figure(figsize=(4, 12))
    plot_coordinate(1, 'hip_flexion')
    plot_coordinate(2, 'knee_angle', negate=True, label='knee flexion')
    plot_coordinate(3, 'ankle_angle', label='ankle dorsiflexion')
    pl.xlabel('percent gait cycle')


def plot_gait_torques(actu, primary_leg, cycle_start, cycle_end,
        primary_footstrike, opposite_footstrike,
        toeoff_time=None):
    """Plots hip, knee, and ankle torques, as a function of percent
    gait cycle, for one gait cycle. Gait cycle starts with a footstrike (start
    of stance). Torques are plotted for both legs; the data for the
    'opposite' leg is properly shifted so that we plot it starting from a
    footstrike as well.

    We assume torques are in units of N-m.

    Knee torque is negated: In OpenSim, a flexion torque is represented by a
    negative value for the torque. In literature, flexion is shown
    as positive.

    Parameters
    ----------
    actu : pytables.Group or dict
        Actuation_force group from a simulation containing joint torques
        (probably an RRA output). If dict, must have fields 'time',
        'hip_flexion_r', 'knee_angle_r' (extension), 'ankle_angle_r'
        (dorsiflexion), 'hip_flexion_l', 'knee_angle_l' (extension),
        'ankle_angle_l' (dorsiflexion).
    primary_leg : str; 'right' or 'left'
        The first and second foot strikes are for which leg?
    cycle_start : float
        Time, in seconds, at which the gait cycle starts.
    cycle_end : float
        Time, in seconds, at which the gait cycle ends.
    primary_footstrike : float
        Time, in seconds, at which the primary leg foot-strikes.
    opposite_footstrike : float
        In between the other two footstrikes, the opposite foot also strikes
        the ground. What's the time at which this happens? This is used to
        shift the data for the opposite leg so that it lines up with the
        `primary` leg's data.
    toeoff_time : bool, optional
        Draw a vertical line on the plot to indicate when the primary foot
        toe-off occurs by providing the time at which this occurs.

    """
    # TODO compare to experimental data.
    # TODO compare to another simulation.
    if primary_leg == 'right':
        opposite_leg = 'left'
    elif primary_leg == 'left':
        opposite_leg = 'right'
    else:
        raise Exception("'primary_leg' is '%s'; it must "
                "be 'left' or 'right'." % primary_leg)

    def plot_for_a_leg(coordinate_name, leg, new_start, color='k', mult=1.0):
        if type(actu) == np.ndarray:
            raw_time = actu['time']
            raw_y = actu[coordinate_name + '_%s_moment' % leg[0]] # TODO uhoh
        else:
            raw_time = actu.cols.time[:]
            raw_y = getattr(actu.cols, '%s_%s' % (coordinate_name, leg[0]))
        time, angle = shift_data_to_cycle( cycle_start, cycle_end, new_start,
                raw_time, raw_y)
        pl.plot(percent_duration(time), mult * angle, color=color,
                label=leg)

    def plot_primary_leg(coordinate_name, **kwargs):
        plot_for_a_leg(coordinate_name, primary_leg, primary_footstrike,
                **kwargs)

    def plot_opposite_leg(coordinate_name, **kwargs):
        plot_for_a_leg(coordinate_name, opposite_leg, opposite_footstrike,
                color=(0.5, 0.5, 0.5), **kwargs)

    def plot_coordinate(index, name, negate=False, label=None):
        ax = pl.subplot(3, 1, index)
        if negate:
            plot_primary_leg(name, mult=-1.0)
            plot_opposite_leg(name, mult=-1.0)
        else:
            plot_primary_leg(name)
            plot_opposite_leg(name)
        # TODO this next line isn't so great of an idea:
        if label == None:
            label = name.replace('_', ' ')
        pl.ylabel('%s (N-m)' % label)
        pl.legend()

        if toeoff_time != None:
            duration = cycle_end - cycle_start
            # 'pgc' is percent gait cycle
            if toeoff_time > primary_footstrike:
                toeoff_pgc = percent_duration_single(toeoff_time,
                        primary_footstrike, duration +
                        primary_footstrike)
            else:
                chunk1 = cycle_end - primary_footstrike
                chunk2 = toeoff_time - cycle_start
                toeoff_pgc = (chunk1 + chunk2) / duration * 100.0
            pl.plot(toeoff_pgc * np.array([1, 1]), ax.get_ylim(),
                    c=(0.5, 0.5, 0.5))

        pl.xticks([0.0, 25.0, 50.0, 75.0, 100.0])

    pl.figure(figsize=(4, 12))
    plot_coordinate(1, 'hip_flexion', label='hip flexion moment')
    plot_coordinate(2, 'knee_angle', negate=True, label='knee flexion moment')
    plot_coordinate(3, 'ankle_angle', label='ankle dorsiflexion moment')
    pl.xlabel('percent gait cycle')


class GaitScrutinyReport:
    """A PDF report that exhaustively compares differences between
    two simulations. The following are compared:

    - joint angles
    - muscle activations
    - muscle forces
    - muscle metabolics
    - verification
        - kinematics errors
        - residuals
        - joint torque reserves

    Expects the gait2392 muscle set, and muscle names.

    """
    def __init__(self, title='OpenSim gait2392 simulation comparison',
            sim_name=None, comp_name=None, do_plot_opposite=True,
            do_plot_joint_torques=False, max_muscle_force=None,
            max_muscle_power=None, max_metabolic_rate=None, muscles=None,
            met_suffix=None,
            sim_primary_color='k', sim_opposite_color='k',
            comp_primary_color='r', comp_opposite_color='r'):
        """

        Parameters
        ----------
        title : str
            Title of the report; will appear at the top of the report.
        sim_name : str, optional
            A brief name (e.g., 'faster') that identifies what this sim data is
            of.
        comp_name : str, optional
            A brief name (e.g., 'control') that identifies what this data is of.
        do_plot_opposite : bool, optional (default: True)
            Plot data for the opposite leg, shifted so that the data starts at
            the start of stance. It is expected then that opposite_strike is
            provided for all sims added to the report.
        do_plot_joint_torques : bool, optional (default: False)
            Plot inverse dynamics joint torques; requires that a joint torque
            table is passed in for all simulations.
        max_muscle_force : float, optional
            Set the ymax value for all muscle force plots to this value. This
            is to make it easier to make comparisons across muscles.
        max_muscle_power : float, optional
            Set the ymax value for all muscle work rate (power) to this value.
        max_metabolic_rate : float, optional
            Set the ymax value for all muscle metabolic consumption rate to
            this value.
        muscles : optional
            TODO
        met_suffix : str, optional
            The suffix in the column names of the metabolic probes. If None,
            then metabolic probe reporter column names are like
            'metabolic_power_soleus_r'. With met_suffix as `umb`, column names
            are expected to be `metabolic_power_umb_soleus_r`.
        sim_primary_color : matplotlib color, optional
            Color of plots for `simulation`, primary leg.
        sim_opposite_color : matplotlib color, optional
        comp_primary_color : matplotlib color, optional
        comp_opposite_color : matplotlib color, optional

        """
        self._title = title
        self._sim_name = sim_name
        self._comp_name = comp_name
        self._sims = list()
        self._comps = list()
        self._do_plot_opposite = do_plot_opposite
        self._do_plot_joint_torques = do_plot_joint_torques
        self._max_muscle_force = max_muscle_force
        self._max_muscle_power = max_muscle_power
        self._max_metabolic_rate = max_metabolic_rate
        self._muscles = muscles
        self.met_suffix = met_suffix
        self.sim_primary_color = sim_primary_color
        self.sim_opposite_color = sim_opposite_color
        self.comp_primary_color = comp_primary_color
        self.comp_opposite_color = comp_opposite_color
        from matplotlib.colors import ColorConverter
        self.cconv = ColorConverter()

    def add_simulation(self, name, sim, primary_leg, cycle_start, cycle_end,
            primary_strike, opposite_strike, toeoff=None, description=None,
            joint_torques=None):
        """Add a simulation to the report. These are the simulations we want to
        compare to the 'control' simulations.

        Parameters
        ----------
        name : str
            Short name (e.g., 'subject 1'), possibly used in legends in future
            versions.
        sim : tables.Group
            The CMC output of the simulation under investigation.
        primary_leg : str; 'left' or 'right'
            Leg for which `first_strike` and `second_strike` are specified.
        cycle_start : float
            The time at which a complete gait cycle starts, within the
            simulation.
        cycle_end : float
            The time at which a complete gait cycle ends, within the simulation.
        primary_strike : float
            A time, between cycle_start and cycle_end, at which the
            primary leg enters stance,
        opposite_strike : float
            A time, between cycle_start and cycle_end, at which the
            opposite (not primary) leg enters stance,
        toeoff : float, optional
            Time at which the primary foot enters swing. If given, a vertical
            line is shown on plots to separate stance and swing.
        joint_torques : tables.Table
            Table containing inverse dynamics joint torques; used if
            `do_plot_joint_torques` is True in the constructor.

        """
        simdict = dict()
        simdict['name'] = name
        simdict['sim'] = sim
        simdict['primary_leg'] = primary_leg
        simdict['cycle_start'] = cycle_start
        simdict['cycle_end'] = cycle_end
        simdict['primary_strike'] = primary_strike
        simdict['opposite_strike'] = opposite_strike
        simdict['toeoff'] = toeoff
        simdict['description'] = description
        simdict['joint_torques'] = joint_torques
        self._sims.append(simdict)

    def add_comparison_simulation(self, name, sim, primary_leg, cycle_start,
            cycle_end, primary_strike, opposite_strike, toeoff=None,
            description=None, joint_torques=None):
        """Add a simulation to the report. These are 'control' simulations to
        which you want to compare another simulation.

        Parameters
        ----------
        name : str
            Short name (e.g., 'subject 1'), possibly used in legends in future
            versions.
        sim : tables.Group
            The CMC output of the simulation under investigation.
        primary_leg : str; 'left' or 'right'
            Leg for which `first_strike` and `second_strike` are specified.
        cycle_start : float
            The time at which a complete gait cycle starts, within the
            simulation.
        cycle_end : float
            The time at which a complete gait cycle ends, within the simulation.
        primary_strike : float
            A time, between cycle_start and cycle_end, at which the
            primary leg enters stance,
        opposite_strike : float
            A time, between cycle_start and cycle_end, at which the
            opposite (not primary) leg enters stance,
        toeoff : float, optional
            Time at which the primary foot exits stance and enters swing. If
            given, a vertical line is shown on plots to separate stance and
            swing.
        joint_torques : tables.Table
            Table containing inverse dynamics joint torques; used if
            `do_plot_joint_torques` is True in the constructor.

        """
        simdict = dict()
        simdict['name'] = name
        simdict['sim'] = sim
        simdict['primary_leg'] = primary_leg
        simdict['cycle_start'] = cycle_start
        simdict['cycle_end'] = cycle_end
        simdict['primary_strike'] = primary_strike
        simdict['opposite_strike'] = opposite_strike
        simdict['toeoff'] = toeoff
        simdict['description'] = description
        simdict['joint_torques'] = joint_torques
        self._comps.append(simdict)

    def generate(self, fname):
        """Generates the report, and prints it to a PDF.

        Parameters
        ----------
        fname : str
            Name of the file to save the PDF to.
        """
        # TODO joint torques

        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(fname)

        # Helper methods
        # --------------
        def plot_for_a_leg(table, landmarks, coordinate_name, leg, new_start,
                color='k', mult=None, interval=1, cut_off=False, **kwargs):
            if landmarks['primary_leg'] == 'right': 
                right_strike = landmarks['primary_strike']
                right_toeoff = landmarks['toeoff']
                left_strike = landmarks['opposite_strike']
                left_toeoff = np.nan
            elif landmarks['primary_leg'] == 'left':
                left_strike = landmarks['primary_strike']
                left_toeoff = landmarks['toeoff']
                right_strike = landmarks['opposite_strike']
                right_toeoff = np.nan
            gl = dataman.GaitLandmarks(
                    primary_leg=landmarks['primary_leg'],
                    cycle_start=landmarks['cycle_start'],
                    cycle_end=landmarks['cycle_end'],
                    left_strike=left_strike,
                    left_toeoff=left_toeoff,
                    right_strike=right_strike,
                    right_toeoff=right_toeoff)
            time, ordinate = shift_data_to_cycle(
                    landmarks['cycle_start'],
                    landmarks['cycle_end'],
                    landmarks[new_start],
                    table.cols.time[::interval],
                    getattr(table.cols,
                        coordinate_name.replace('!', leg[0]))[::interval],
                    )#cut_off=True)

            pgc, ordinate = data_by_pgc(
                    table.cols.time[::interval],
                    getattr(table.cols,
                        coordinate_name.replace('!', leg[0]))[::interval],
                    gl, side=leg)

            if mult != None: ordinate *= mult

            pl.plot(pgc, ordinate, color=color, label=leg, **kwargs)

        def plot_primary_leg(table, series, coordinate_name, **kwargs):
            plot_for_a_leg(table, series, coordinate_name,
                    series['primary_leg'], 'primary_strike', **kwargs)

        def plot_opposite_leg(table, series, coordinate_name, **kwargs):
            if series['primary_leg'] == 'right': opposite_leg = 'left'
            elif series['primary_leg'] == 'left': opposite_leg = 'right'
            else: raise Exception(
                    "'primary_leg' is %s; must be 'right' or 'left'." % (
                        series['primary_leg']))

            if 'color' in kwargs:
                rgba = list(self.cconv.to_rgba(kwargs['color']))
                # Edit alpha value.
                rgba[3] = 0.5 * rgba[3]
                kwargs['color'] = rgba
            else:
                kwargs['color'] = (0.5, 0.5, 0.5)
            plot_for_a_leg(table, series, coordinate_name, opposite_leg,
                    'opposite_strike', **kwargs)

        def plot_coordinate(grid, loc, table, name, units='-', negate=False,
                label=None, title=None, ylims=None, **kwargs):

            def plot_a_series(series, label, **kwargs):

                # TODO this next line isn't so great of an idea:
                if label: pl.ylabel('%s (%s)' % (label, units))
                if title: pl.title(title, fontsize=12)

                if table == 'joint_torques':
                    thetable = series['joint_torques']
                else:
                    thetable = getattr(series['sim'], table)

                mult = -1.0 if negate else None

                plot_primary_leg(thetable, series, name, mult=mult, **kwargs)
                if self._do_plot_opposite:
                    plot_opposite_leg(thetable, series, name, mult=mult,
                            **kwargs)

            def plot_landmarks(series, ylims, **kwargs):
                if series['toeoff'] != None:
                    duration = series['cycle_end'] - series['cycle_start']
                    # 'pgc' is percent gait cycle
                    if series['toeoff'] > series['primary_strike']:
                        toeoff_pgc = percent_duration_single(series['toeoff'],
                                series['primary_strike'], duration +
                                series['primary_strike'])
                    else:
                        chunk1 = series['cycle_end'] - series['primary_strike']
                        chunk2 = series['toeoff'] - series['cycle_start']
                        toeoff_pgc = (chunk1 + chunk2) / duration * 100.0

                    if ylims == None: ylims = ax.get_ylim()
                    pl.plot(toeoff_pgc * np.array([1, 1]), ylims,
                            c=(0.8, 0.8, 0.8), zorder=0, **kwargs)

            if type(grid) == tuple:
                ax = pl.subplot2grid(grid, loc)
            else:
                ax = pl.subplot(grid[loc[0], loc[1]])

            # Add a curve to this plot for each sim and each comp.
            for sim in self._sims:
                plot_a_series(sim, label, color=self.sim_primary_color,
                        lw=1.5, **kwargs)
            for comp in self._comps:
                plot_a_series(comp, label, color=self.comp_primary_color,
                        lw=1.5, **kwargs)

            # Must be done after all other plotting, so that we use the correct
            # ylims.
            for sim in self._sims:
                plot_landmarks(sim, ylims, color=self.sim_primary_color,
                        lw=1.5)
            for comp in self._comps:
                plot_landmarks(comp, ylims, color=self.comp_primary_color,
                        lw=1.5)

            pl.xticks([0.0, 25.0, 50.0, 75.0, 100.0])
            pl.xlim(0, 100)

        # Title page
        # ----------
        ftitle = pl.figure()
        ax = pl.subplot(111)
        pl.axis('off')
        ftitle.suptitle(self._title, fontweight='bold')
        desc = str()
        if self._sim_name and self._comp_name:
            desc = ('Comparison between:\n- %s (black lines), '
                    'and\n- %s (red lines). \n' % (
                    self._sim_name, self._comp_name))
        desc += """
        Black lines are for the primary limb; gray lines are for the
        opposite limb, and that data is wrapped around so that it starts
        with stance."""

        # TODO def specific_metabolic_cost(subject_mass,
        # TODO         time, value, cycle_duration, cycle_start=None):

        pl.text(0, 0.7, desc)
        pp.savefig(ftitle)

        # Joint angles
        # ------------
        print 'Processing joint angles.'
        fkin = pl.figure(figsize=(3.5, 7.5))
        gs = pl.GridSpec(3, 1)
        pl.suptitle('JOINT ANGLES', weight='bold')

        def plot_kin(index, name, **kwargs):
            plot_coordinate(gs, (index, 0), 'Kinematics_q',
                    '%s_!' % name, 'degrees', **kwargs)

        plot_kin(0, 'hip_flexion', title='hip', label='hip flexion')
        plot_kin(1, 'knee_angle', negate=True, label='knee flexion',
                title='knee')
        plot_kin(2, 'ankle_angle', label='ankle dorsiflexion', title='ankle')
        pl.xlabel('percent gait cycle')
        pl.xlim(0, 100)

        gs.tight_layout(fkin, rect=[0, 0, 1, 0.95])
        pp.savefig(fkin)
        pl.close(fkin)

        # Joint torques!
        # --------------
        if self._do_plot_joint_torques:
            print 'Processing joint torques.'
            fjt = pl.figure(figsize=(3.5, 7.5))
            gs = pl.GridSpec(3, 1)
            pl.suptitle('JOINT TORQUES', weight='bold')

            def plot_jt(index, name, **kwargs):
                plot_coordinate(gs, (index, 0), 'joint_torques',
                        '%s_!' % name, 'N-m', **kwargs)

            plot_jt(0, 'hip_flexion', title='hip', label='hip flexion moment')
            plot_jt(1, 'knee_angle', negate=True, label='knee flexion moment',
                    title='knee')
            plot_jt(2, 'ankle_angle',
                    label='ankle dorsiflexion moment', title='ankle')
            pl.xlabel('percent gait cycle')
            pl.xlim(0, 100)

            gs.tight_layout(fjt, rect=[0, 0, 1, 0.95])
            pp.savefig(fjt)
            pl.close(fjt)

        # Muscles!
        # ========
        def create_plate(subtitle, table, ylabel, pattern, dims, mset,
                yticks=None, **kwargs):
            print 'Processing muscle %s.' % subtitle
            f = pl.figure(figsize=(3.5 * dims[1], 2.5 * dims[0]))
            gs = pl.GridSpec(dims[0], dims[1])
            pl.suptitle('MUSCLE %s' % subtitle.upper(), weight='bold')

            for loc, name in mset.items():
                try:
                    plot_coordinate(gs, loc, table,
                            pattern % name, title=muscle_names[name], **kwargs)
                    if loc[1] == 0: pl.ylabel(ylabel)
                    if loc[0] == (dims[0] - 1): pl.xlabel('percent gait cycle')

                    if yticks: pl.yticks(yticks)
                except Exception as e:
                    print("[GaitScrutiny] Can't plot %s." % name)

            gs.tight_layout(f, rect=[0, 0, 1, 0.95])
            pp.savefig(f)
            pl.close(f)

        # Define muscles to use for the remaining sets of plots.
        muscle_names = dict()
        muscle_names['rect_fem'] = 'rectus femoris'
        muscle_names['vas_med'] = 'vastus medialis'
        muscle_names['vas_int'] = 'vastus intermedius'
        muscle_names['vas_lat'] = 'vastus lateralis'
        muscle_names['semimem'] = 'semimembranosus'
        muscle_names['semiten'] = 'semitendinosus'
        muscle_names['bifemsh'] = 'biceps femoris short head'
        muscle_names['bifemlh'] = 'biceps femoris long head'
        muscle_names['tib_ant'] = 'tibialis anterior'
        muscle_names['ext_dig'] = 'extensor digitorum'
        muscle_names['ext_hal'] = 'extensor hallucis'
        muscle_names['per_tert'] = 'peroneus tertius'
        muscle_names['med_gas'] = 'medial gastrocnemius'
        muscle_names['lat_gas'] = 'lateral gastrocnemius'
        muscle_names['soleus'] = 'soleus'
        muscle_names['tib_post'] = 'tibialis posterior'

        muscle_names['flex_dig'] = 'flexor digitorum'
        muscle_names['flex_hal'] = 'flexor hallucis'
        muscle_names['per_brev'] = 'peroneus brevis'
        muscle_names['per_long'] = 'peroneus longus'
        muscle_names['ercspn'] = 'erector spinae'
        muscle_names['extobl'] = 'external obliques'
        muscle_names['intobl'] = 'internal obliques'
        muscle_names['pect'] = 'pectineus'
        muscle_names['quad_fem'] = 'quadratus femoris'
        muscle_names['gem'] = 'gemellus'
        muscle_names['peri'] = 'periformis'
        muscle_names['grac'] = 'gracilis'
        muscle_names['sar'] = 'sartorius'
        muscle_names['tfl'] = 'tensor fascia latae'
        muscle_names['iliacus'] = 'iliacus'
        muscle_names['psoas'] = 'psoas major'

        muscle_names['glut_max1'] = 'gluteus maximus 1'
        muscle_names['glut_max2'] = 'gluteus maximus 2'
        muscle_names['glut_max3'] = 'gluteus maximus 3'
        muscle_names['glut_med1'] = 'gluteus medius 1'
        muscle_names['glut_med2'] = 'gluteus medius 2'
        muscle_names['glut_med3'] = 'gluteus medius 3'
        muscle_names['glut_min1'] = 'gluteus minimus 1'
        muscle_names['glut_min2'] = 'gluteus minimus 2'
        muscle_names['glut_min3'] = 'gluteus minimus 3'
        muscle_names['add_mag1'] = 'adductor magnus 1'
        muscle_names['add_mag2'] = 'adductor magnus 2'
        muscle_names['add_mag3'] = 'adductor magnus 3'
        muscle_names['add_mag4'] = 'adductor magnus 4'
        muscle_names['add_long'] = 'adductor longus'
        muscle_names['add_brev'] = 'adductor brevis'

        if self._muscles == None:
            mset1 = dict()

            # Quadriceps.
            mset1[(0, 0)] = 'rect_fem'
            mset1[(0, 1)] = 'vas_med'
            mset1[(0, 2)] = 'vas_int'
            mset1[(0, 3)] = 'vas_lat'

            # Hamstrings.
            mset1[(1, 0)] = 'semimem'
            mset1[(1, 1)] = 'semiten'
            mset1[(1, 2)] = 'bifemsh'
            mset1[(1, 3)] = 'bifemlh'

            # Dorsiflexors.
            mset1[(2, 0)] = 'tib_ant'
            mset1[(2, 1)] = 'ext_dig'
            mset1[(2, 2)] = 'ext_hal'
            mset1[(2, 3)] = 'per_tert'

            # Plantarflexors.
            mset1[(3, 0)] = 'med_gas'
            mset1[(3, 1)] = 'lat_gas'
            mset1[(3, 2)] = 'soleus'
            mset1[(3, 3)] = 'tib_post'

            # Define muscles to use for the remaining sets of plots.
            mset2 = dict()

            # Torso.
            mset2[(0, 0)] = 'ercspn'
            mset2[(0, 1)] = 'extobl'
            mset2[(0, 2)] = 'intobl'
            mset2[(0, 3)] = 'psoas'

            # Butt muscles.
            mset2[(1, 0)] = 'pect'
            mset2[(1, 1)] = 'quad_fem'
            mset2[(1, 2)] = 'gem'
            mset2[(1, 3)] = 'peri'

            # Thigh muscles.
            mset2[(2, 0)] = 'grac'
            mset2[(2, 1)] = 'sar'
            mset2[(2, 2)] = 'tfl'
            mset2[(2, 3)] = 'iliacus'

            # Plantarflexors.
            mset2[(3, 0)] = 'flex_dig'
            mset2[(3, 1)] = 'flex_hal'
            mset2[(3, 2)] = 'per_brev'
            mset2[(3, 3)] = 'per_long'

            # Define muscles to use for the remaining sets of plots.
            mset3 = dict()

            mset3[(0, 0)] = 'glut_max1'
            mset3[(0, 1)] = 'glut_max2'
            mset3[(0, 2)] = 'glut_max3'

            mset3[(1, 0)] = 'glut_med1'
            mset3[(1, 1)] = 'glut_med2'
            mset3[(1, 2)] = 'glut_med3'

            mset3[(2, 0)] = 'glut_min1'
            mset3[(2, 1)] = 'glut_min2'
            mset3[(2, 2)] = 'glut_min3'

            mset3[(0, 3)] = 'add_mag1'
            mset3[(1, 3)] = 'add_mag2'
            mset3[(2, 3)] = 'add_mag3'
            mset3[(3, 3)] = 'add_mag3'

            mset3[(3, 0)] = 'add_long'
            mset3[(3, 1)] = 'add_brev'

            msubnames = ['key locomotion muscles', 'misc muscles',
                    'remaining hip muscles']
            msets = [mset1, mset2, mset3]

            # Activations.
            for i in range(3):
                create_plate('activations (%s)' % msubnames[i],
                        'states', 'activation (-)', '%s_!_activation',
                        (4, 4), msets[i],
                        yticks=[0.0, 0.5, 1.0], interval=10, ylims=(0, 1))

            # Forces.
            for i in range(3):
                if self._max_muscle_force: ylims = (0, self._max_muscle_force)
                else: ylims = None
                create_plate('forces (%s)' % msubnames[i],
                        'Actuation_force', 'force (N)', '%s_!', (4, 4), msets[i],
                        ylims=ylims)

            # Power.
            for i in range(3):
                if self._max_muscle_power: ylims = (0, self._max_muscle_power)
                else: ylims = None
                create_plate('work rate / power (%s)' % msubnames[i],
                        'Actuation_power', 'work rate (W)', '%s_!', (4, 4),
                        msets[i], ylims=ylims)

            # Metabolics.
            if self.met_suffix:
                met_pattern = 'metabolic_power_{}_%s_!'.format(
                        self.met_suffix)
            else:
                met_pattern = 'metabolic_power_%s_!'

            try:
                for i in range(3):
                    if self._max_metabolic_rate:
                        ylims = (0, self._max_metabolic_rate)
                    else: ylims = None
                    create_plate('metabolics (%s)' % msubnames[i],
                            'ProbeReporter_probes',
                            'metabolic consumption rate (W)',
                            met_pattern, (4, 4), msets[i],
                            ylims=ylims)
            except Exception, e:
                print e.message
        else:
            n_rows = 0
            n_cols = 0
            for mk, mv in self._muscles.items():
                if mk[0] > n_rows: n_rows = mk[0]
                if mk[1] > n_cols: n_cols = mk[1]
            grid = (n_rows + 1, n_cols + 1)

            # Activations.
            create_plate('activations',
                    'states', 'activation (-)', '%s_!_activation',
                    grid, self._muscles,
                    yticks=[0.0, 0.5, 1.0], interval=10, ylims=(0, 1))

            # Forces.
            if self._max_muscle_force: ylims = (0, self._max_muscle_force)
            else: ylims = None
            create_plate('forces',
                    'Actuation_force', 'force (N)', '%s_!', grid,
                    self._muscles,
                    ylims=ylims)

            # Power.
            if self._max_muscle_power: ylims = (0, self._max_muscle_power)
            else: ylims = None
            create_plate('work rate / power',
                    'Actuation_power', 'work rate (W)', '%s_!', grid,
                    self._muscles, ylims=ylims)

            # Metabolics.
            try:
                if self._max_metabolic_rate:
                    ylims = (0, self._max_metabolic_rate)
                else: ylims = None
                create_plate('metabolics',
                        'ProbeReporter_probes',
                        'metabolic consumption rate (W)',
                        'metabolic_power_%s_!', grid, self._muscles,
                        ylims=ylims)
            except Exception, e:
                print e.message


        # Verification.
        # -------------
        for all_series in [self._sims, self._comps]:
            for a_sim in all_series:
                fig = plot_simulation_verification(a_sim['sim'])
                pl.suptitle('VERIFICATION OF %s' % (a_sim['name']),
                        fontweight='bold')
                pp.savefig(fig)

        # That's all, folks.
        # ------------------
        pp.close()

def percent_duration_single(time, start, end):
    """Converts a single time value to a percent duration (e.g., percent gait
    cycle) value. The difference between this method and `percent_duration` is
    that this works on a single float, rather than on an array.

    Parameters
    ----------
    time : float
        The time value to convert, with units of time (e.g., seconds).
    start : float
        The start of the duration (0 %), in the same units of time.
    end : float
        The end of the duration (100 %), in the same units of time.

    """
    return (time - start) / (end - start) * 100.0


def percent_duration(time, start=None, end=None):
    """Converts a time array to percent duration (e.g., percent gait cycle).

    Parameters
    ----------
    time : np.array
        The time data to convert, with units of time (e.g., seconds).
    start : float, optional
        Start time of the duration. If not provided, we use time[0].
    end : float, optional
        End time of the duration. If not provided, we use time[-1].

    Returns
    -------
    percent_duration : np.array
        Varies from 0 to 100.

    """
    if start == None: start = time[0]
    if end == None: end = time[-1]
    return (time - start) / (end - start) * 100.0


def plot(column, *args, **kwargs):
    """Plots a column of a pyTables table, against time.

    """
    pl.plot(column.table.cols.time, column, *args, **kwargs)
    pl.xlabel('time (s)')

def data_by_pgc(time, data, gl, side='left'):

    if side == 'left':
        strike = gl.left_strike
    elif side == 'right':
        strike = gl.right_strike
    else:
        raise Exception("side '%s' not recognized." % side)

    cycle_duration = (gl.cycle_end - gl.cycle_start)

    if strike < gl.cycle_start:
        strike += cycle_duration
    if strike > gl.cycle_end:
        strike -= cycle_duration                   

    ts, ys = shift_data_to_cycle(gl.cycle_start,
            gl.cycle_end, strike, time, data)

    pgc = percent_duration(ts, 0, cycle_duration)

    if np.any(pgc > 100.0):
        print('Percent gait cycle greater than 100: %f' % np.max(pgc))
    if np.any(pgc > 100.01) or np.any(pgc < 0.0):
        raise Exception('Percent gait cycle out of range.')

    return pgc, ys

def plot_opposite_strike_pgc(gl, side, axes=None, *args, **kwargs):
    if axes:
        ax = axes
    else:
        ax = pl
    if side == 'left':
        strike = gl['right_strike']
    elif side == 'right':
        strike = gl['left_strike']
    else:
        raise Exception()
    cycle_duration = (gl.cycle_end - gl.cycle_start)
    while strike < gl[side + '_strike']:
        strike += cycle_duration
    while strike > gl[side + '_strike'] + cycle_duration:
        strike -= cycle_duration
    ax.axvline(percent_duration_single(strike,
        gl[side + '_strike'],
        gl[side + '_strike'] + cycle_duration),
        *args, **kwargs)

def toeoff_pgc(gl, side):
    toeoff = gl[side + '_toeoff']
    cycle_duration = (gl.cycle_end - gl.cycle_start)
    while toeoff < gl[side + '_strike']:
        toeoff += cycle_duration
    while toeoff > gl[side + '_strike'] + cycle_duration:
        toeoff -= cycle_duration
    return percent_duration_single(toeoff,
            gl[side + '_strike'],
            gl[side + '_strike'] + cycle_duration)

def plot_toeoff_pgc(gl, side, axes=None, *args, **kwargs):
    if axes:
        ax = axes
    else:
        ax = pl
    ax.axvline(toeoff_pgc(gl, side), *args, **kwargs)

def plot_pgc(time, data, gl, side='left', axes=None, plot_toeoff=False, *args,
        **kwargs):

    pgc, ys = data_by_pgc(time, data, gl, side=side)

    if axes:
        ax = axes
    else:
        ax = pl
    ax.plot(pgc, ys, *args, **kwargs)

    if plot_toeoff:
        if 'color' not in kwargs and 'c' not in kwargs:
            kwargs['color'] = 'gray'
        if 'label' in kwargs: kwargs.pop('label')
        plot_toeoff_pgc(gl, side, ax, *args, zorder=0, **kwargs)


def avg_and_std_time_series_across_gait_trials(list_of_dicts, n_points=400):

    output_pgc = np.linspace(0, 100, n_points)
    data = np.empty((n_points, len(list_of_dicts)))
    for i, item in enumerate(list_of_dicts):
        pgc, ys = data_by_pgc(**item)
        data[:, i] = np.interp(output_pgc, pgc, ys, left=np.nan, right=np.nan)
    return output_pgc, nanmean(data, axis=1), nanstd(data, axis=1)

def avg_and_std_toeoff(list_of_dicts):
    toeoffs = np.empty(len(list_of_dicts))
    for i, item in enumerate(list_of_dicts):
        toeoffs[i] = toeoff_pgc(item['gl'], item['side'])
    return nanmean(toeoffs), nanstd(toeoffs)

def plot_avg_and_std_time_series_across_gait_trials(list_of_dicts,
        n_points=400, lw=1.0, alpha=0.5, label=None, plot_toeoff=False, *args,
        **kwargs):

    pgc, avg, std = avg_and_std_time_series_across_gait_trials(list_of_dicts,
            n_points=n_points)

    if plot_toeoff:
        pl.axvline(avg_and_std_toeoff(list_of_dicts)[0], lw=lw,
                color='lightgray', zorder=0)

    pl.fill_between(pgc, avg + std, avg - std, alpha=alpha, *args, **kwargs)

    pl.plot(pgc, avg, *args, lw=lw, label=label, **kwargs)




