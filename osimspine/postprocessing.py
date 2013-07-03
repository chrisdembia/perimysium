"""Contains methods for postprocessing simulations, or computing derived
quantities, such as sum of squared activations, from simulation (e.g., CMC)
results.

"""

import numpy as np
import pylab as pl
import tables

def nearest_index(array, val):
    return np.abs(array - val).argmin()


def sum_of_squared_activations(states_table):
    """Computes the sum, across all muscles, of the square of muscle
    activations.

    Parameters
    ----------
    states_table : tables.Table
        A pyTables table containing muscle activation time series data.

    Returns
    -------
    SSA : numpy.ndarray
        Time series of the sum of squared muscle activations.
    
    """
    SSA = np.zeros(states_table.col('time').shape)
    for col_name in states_table.colnames:
        if col_name.endswith('activation'):
            SSA += states_table.col(col_name)**2
    return SSA

def avg_sum_of_squared_activations(states_table):
    """Computes the average value of the sum of squared activations.

    Parameters
    ----------
    states_table : tables.Table
        A pyTables table containing muscle activation time series data.

    Returns
    -------
    ASSA : float
        Average value, in time, of the sum of squared muscle activations.

    """
    time = states_table.cols.time
    duration = time[-1] - time[0]
    integral = np.trapz(sum_of_squared_activations(states_table), x=time)
    return integral / duration

def change_in_activations_between_two_simulations(states_table):
    """
    """
    pass

def avg(time, value, init_time=None, final_time=None):
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

    return np.trapz(
            value[init_idx:final_idx],
            x=time[init_idx:final_idx])

def sorted_maxabs(table, init_time=None, final_time=None, exclude=None):
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
    exclude : list of str's, optional
        The names of columns to exclude from the sorting (e.g.,
        ['umb_val_wholebody']).

    Returns
    -------
    sorted_vals : numpy.array
        The max absolute value for all non-excluded columns, sorted in
        ascending order.
    sorted_args : list of str's
        The column names corresponding to the sorted values, in ascending
        order.

    """
    if exclude != None and not type(exclude) == list:
        raise Exception("'exclude' must be None or a list.")

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

    for coln in table.colnames:
        if coln != 'time' and (exclude == None or not coln in exclude):
            maxs[coln] = np.max(np.abs(table.col(coln)[init_idx:final_idx]))

    sorted_vals = np.sort(maxs.values())
    sorted_args = [maxs.keys()[idx] for idx in np.argsort(maxs.values())]

    return sorted_vals, sorted_args


def sorted_avg(table, init_time=None, final_time=None, exclude=None):
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
    exclude : list of str's, optional
        The names of columns to exclude from the sorting (e.g.,
        ['umb_val_wholebody']).

    Returns
    -------
    sorted_vals : numpy.array
        The average value for all non-excluded columns, sorted in ascending
        order.
    sorted_args : list of str's
        The column names corresponding to the sorted values, in ascending
        order.

    """
    if exclude != None and not type(exclude) == list:
        raise Exception("'exclude' must be None or a list.")

    time = table.cols.time

    avgs = dict()

    for coln in table.colnames:
        if coln != 'time' and (exclude == None or not coln in exclude):
            avgs[coln] = avg(time, table.col(coln), init_time, final_time)

    sorted_vals = np.sort(avgs.values())
    sorted_args = [avgs.keys()[idx] for idx in np.argsort(avgs.values())]

    return sorted_vals, sorted_args


def sorted_avg_difference(table1, table2,
        init_time=None, final_time=None, exclude=None):
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
        ['umb_val_wholebody']).

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
        if coln != 'time' and (exclude == None or not coln in exclude):
            diff_avgs[coln] = (
                    avg(table1.cols.time, table1.col(coln),
                        init_time, final_time) - 
                    avg(table2.cols.time, table2.col(coln),
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
    sim_grup : tables.Group
        A pyTables group holding tables docked from a simulation (e.g., CMC).
    n_max : int, optional
        Plots the n_max maximum columns. By default, all columns are plotted.
        If n_max = 5, only the 5 most errorful columns are plotted.
    violators_only : bool, optional
        Only the columns violating the "good" threshold are plotted. This
        overrides the n_max setting. By default, this is False.
    show_legend : bool, optional
        Show a legend for all columns plotted. By default, False.

    """
    plot_kinematics_verification(sim_group.pErr, **kwargs)
    plot_residuals_verification(sim_group.Actuation_force, **kwargs)
    plot_reserves_verification(sim_group.Actuation_force, **kwargs)


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


def plot_kinematics_verification(pErr_table,
        n_max=None, violators_only=False, show_legend=False):
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
    ylim_trans = np.max(trans_okay_thresh, 1.1 * max_trans)

    # --- Translation.
    pl.figure()

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

    pl.xlabel('time (s)')
    pl.ylabel('error (cm)')
    pl.title('Translation error')

    if show_legend: pl.legend(**legend_kwargs)

    # --- Rotation.
    pl.figure()

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


def plot_residuals_verification(actforce_table, n_max=None,
        violators_only=False, show_legend=False):
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
        if coln.startswith('F'):
            max_force = max(max_force, sorted_actf[i])
    ylim_force = max(force_okay_thresh, 1.1 * max_force)

    # --- Translation.
    pl.figure()

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
        if coln.startswith('F'):
            if violators_only:
                if (sorted_actf[len(sorted_actf) - i - 1] > force_good_thresh):
                    plot(actforce_table.col(coln), coln)
                continue
            if n_max == None or count_force < n_max:
                plot(actforce_table.col(coln), coln)
                count_force += 1

    pl.xlim(time[0], time[-1])
    pl.ylim((-ylim_force, ylim_force))

    pl.xlabel('time (s)')
    pl.ylabel('residual force (N)')
    pl.title('Residual forces')

    if show_legend: pl.legend(**legend_kwargs)

    # --- Rotation.
    pl.figure()

    max_moment = -np.inf
    for i, coln in enumerate(sorted_colns):
        if coln.startswith('M'):
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
        if coln.startswith('M'):
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


def plot_reserves_verification(actforce_table, n_max=None,
        violators_only=False, show_legend=False):
    """Plots reserves verification information, with comparisons to "good"
    and "okay", and "bad" thresholds, using matplotlib. Assumes all reserves
    are torque actuators, and their names end with 'reserve'.

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

    """
    legend_kwargs = {'loc': 'upper left', 'bbox_to_anchor': (1, 1)}

    good_thresh = 10 # Newton-meters
    okay_thresh = 25

    time = actforce_table.cols.time
    duration = time[-1] - time[0]

    sorted_actf, sorted_colns = sorted_maxabs(actforce_table)

    pl.figure()

    max_torque = -np.inf
    for i, coln in enumerate(sorted_colns):
        if coln.endswith('reserve'):
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
        if coln.endswith('reserve'):
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


def verify_simulation(sim_group):
    """Prints maximum kinematics errors, residuals, and reserves, and compares
    them to thresholds.

    Parameters
    ----------
    sim_grup : tables.Group
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

    print 'Translation: %s with maximum %f cm' % (trans_eval, max_trans)
    print 'Rotation: %s with maximum %f degrees' % (rot_eval, max_rot)


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
        if coln.startswith('F'):
            max_force = max(max_force, sorted_actf[i])

    max_moment = -np.inf
    for i, coln in enumerate(sorted_colns):
        if coln.startswith('M'):
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
        if coln.endswith('reserve'):
            max_torque = max(max_torque, sorted_actf[i])

    evaluation = _evaluate_threshold(max_torque, good_thresh, okay_thresh)

    print 'Reserve torques: %s with maximum %f N-m' % (evaluation, max_torque)


def plot(column):
    """Plots a column of a pyTables table, against time.

    """
    pl.plot(column.table.cols.time, column, label=column.name)
    pl.xlabel('time (s)')
    pl.title(column.name)

# Detect gait landmarks given GRF data.

