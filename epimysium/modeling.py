# ----------------------------------------------------------------------- #
# The OpenSim API is a toolkit for musculoskeletal modeling and           #
# simulation. See http://opensim.stanford.edu and the NOTICE file         #
# for more information. OpenSim is developed at Stanford University       #
# and supported by the US National Institutes of Health (U54 GM072970,    #
# R24 HD065690) and by DARPA through the Warrior Web program.             #
#                                                                         #
# Copyright (c) 2005-2012 Stanford University and the Authors             #
#                                                                         #
# Licensed under the Apache License, Version 2.0 (the "License");         #
# you may not use this file except in compliance with the License.        #
# You may obtain a copy of the License at                                 #
# http://www.apache.org/licenses/LICENSE-2.0.                             #
#                                                                         #
# Unless required by applicable law or agreed to in writing, software     #
# distributed under the License is distributed on an "AS IS" BASIS,       #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         #
# implied. See the License for the specific language governing            #
# permissions and limitations under the License.                          #
# ----------------------------------------------------------------------- #
#
# Authors: Christopher Dembia, Thomas Uchida
# Stanford University

"""Modifying OpenSim Model's, and any other OpenSim object.

"""

# TODO currently, this module assumes this package is installed via Jython, and
# that the OpenSim Jython wrapping is on the Jython CLASSPATH. When an OpenSim
# Python wrapping is available, we can use CPython instead.

import org.opensim.modeling as osm

def control_set_from_storage_files(sto_list):
    """
    Parameters
    ----------
    sto_list: list of org.opensim.modeling.Storage's
        Each column is written to a single ControlSet, as ControlLinear's.

    Returns
    -------
    cset : org.opensim.modeling.ControlSet

    """
    # TODO documentation.
    cset = osm.ControlSet()
    for sto in sto_list:
        this_cset = osm.ControlSet(sto)
        for i_control in range(this_cset.getSize()):
            this_clin = osm.ControlLinear.safeDownCast(this_cset.get(i_control))
            this_clin.setName(this_clin.getName() + ".excitation")
            this_clin.setExtrapolate(False)
            # Undo the default that makes sense for muscles.
            this_clin.setDefaultParameterMin(-10000.0)
            this_clin.setDefaultParameterMax(10000.0)

            for i_param in range(this_clin.getNumParameters()):
                this_clin.setControlValueMin(
                        this_clin.getParameterTime(i_param),
                        this_clin.getParameterValue(i_param))
                this_clin.setControlValueMax(
                        this_clin.getParameterTime(i_param),
                        this_clin.getParameterValue(i_param))
            cset.cloneAndAppend(this_clin)
    return cset


def storage2piecewise_linear_function(sto, column_name, scale_factor=None):
    """Returns a column from a storage file as an
    org.opensim.modeling.PiecewiseLinearFunction. We advise that, once you get
    the function, you name it.

    Parameters
    ----------
    sto : org.opensim.modeling.Storage
        An OpenSim Storage object.
    column_name : str
        Name of a column in `sto`.
    scale_factor : float, optional
        Scale the column data before placing it in the function.

    Returns
    -------
    plf : org.opensim.modeling.PiecewiseLinearFunction
        Just like you asked for.

    """
    time = osm.ArrayDouble()
    sto.getTimeColumn(time)

    state_index = sto.getStateIndex(column_name)

    if type(scale_factor) == float:
        sto.multiplyColumn(state_index, scale_factor)
    elif scale_factor == None:
        pass
    else:
        raise Exception('scale_factor, if specified, must be a float.')

    ordinate = osm.ArrayDouble()
    sto.getDataColumn(state_index, ordinate)

    return osm.PiecewiseLinearFunction(time.getSize(), time.get(),
            ordinate.get())

def add_metabolics_probes(model, twitch_ratio_set='gait2392'):
    """Adds Umberger2010MuscleMetabolicsProbes to an OpenSim model. Adds a probe
    for each muscle, as well as a whole-body probe that returns cumulative
    energy expended across all muscles in the model. When possible, we use
    published twitch ratios for the probes. For muscles for which we do not
    have data, we use a twitch ratio of 0.5. This method doesn't return
    anything; the model given to the method is simply modified.

    Parameters
    ----------
    model : org.opensim.modeling.Model
        An OpenSim Model that has muscles.
    twitch_ratio_set : str; 'gait2392' or 'gait1018'
        The experimental data set to use for the model, depending on the model
        we are adding probes to.

    """
    # Twitch ratios
    # -------------
    if twitch_ratio_set == 'gait2392':
        twitchRatios = {
            'glut_med1': 0.55, 'glut_med2': 0.55, 'glut_med3': 0.55,
            'glut_min1': 0.55, 'glut_min2': 0.55, 'glut_min3': 0.55,
            'semimem': 0.4925, 'semiten': 0.425,
            'bifemlh': 0.5425, 'bifemsh': 0.529,
            'add_mag1': 0.552, 'add_mag2': 0.552, 'add_mag3': 0.552,
            'glut_max1': 0.55, 'glut_max2': 0.55, 'glut_max3': 0.55,
            'iliacus': 0.5, 'psoas': 0.5, 'rect_fem': 0.3865,
            'vas_med': 0.503, 'vas_int': 0.543, 'vas_lat': 0.455,
            'med_gas': 0.566, 'lat_gas': 0.507, 'soleus': 0.803,
            'tib_post': 0.6, 'flex_dig': 0.6, 'flex_hal': 0.6, 'tib_ant': 0.7,
            'per_brev': 0.6, 'per_long': 0.6, 'per_tert': 0.75,
            'ext_dig': 0.75, 'ext_hal': 0.75,
            'ercspn': 0.6, 'intobl': 0.56, 'extobl': 0.58,
            'sar': -1, 'add_long': -1, 'add_brev': -1,
            'tfl': -1, 'pect': -1, 'grac': -1,
            'quad_fem': -1, 'gem': -1, 'peri': -1}
    elif twitch_ratio_set == 'gait1018':
        twitchRatios = {
            'hamstrings': 0.49, 'bifemsh': 0.53, 'glut_max': 0.55,
            'iliopsoas': 0.50, 'rect_fem': 0.39, 'vasti': 0.50,
            'gastroc': 0.54, 'soleus': 0.80,
            'tib_ant': 0.70}

    # Parameters used for all probes
    # ------------------------------
    # The following booleans are constructor arguments for the Umberger probe.
    # These settings are used for all probes.
    activationMaintenanceRateOn = True
    shorteningRateOn = True
    basalRateOn = False
    mechanicalWorkRateOn = True

    # The mass of each muscle will be calculated using data from the model:
    #   muscleMass = (maxIsometricForce / sigma) * rho * optimalFiberLength
    # where sigma = 0.25e6 is the specific tension of mammalian muscle (in
    # Pascals) and rho = 1059.7 is the density of mammalian muscle (in kg/m^3).

    # The slow-twitch ratio used for muscles that either do not appear in the
    # file, or appear but whose proportion of slow-twitch fibers is unknown.
    defaultTwitchRatio = 0.5

    # Whole-body probe
    # ----------------
    # Define a whole-body probe that will report the total metabolic energy
    # consumption over the simulation.
    wholeBodyProbe = osm.Umberger2010MuscleMetabolicsProbe(
        activationMaintenanceRateOn,
        shorteningRateOn,
        basalRateOn,
        mechanicalWorkRateOn)
    wholeBodyProbe.setOperation("value")
    wholeBodyProbe.set_report_total_metabolics_only(False);

    # Add the probe to the model and provide a name.
    model.addProbe(wholeBodyProbe)
    wholeBodyProbe.setName("metabolic_power")

    # Loop through all muscles, adding parameters for each into the whole-body
    # probe.
    for iMuscle in range(model.getMuscles().getSize()):
        thisMuscle = model.getMuscles().get(iMuscle)

        # Get the slow-twitch ratio from the data we read earlier. Start with
        # the default value.
        slowTwitchRatio = defaultTwitchRatio

        # Set the slow-twitch ratio to the physiological value, if it is known.
        for key, val in twitchRatios.items():
            if thisMuscle.getName().startswith(key) and val != -1:
                slowTwitchRatio = val

        # Add this muscle to the whole-body probe. The arguments are muscle
        # name, slow-twitch ratio, and muscle mass. Note that the muscle mass
        # is ignored unless we set useProvidedMass to True.
        wholeBodyProbe.addMuscle(thisMuscle.getName(),
                                 slowTwitchRatio,
                                 -1)

def enable_probes(model_fpath):
    """Ensures that all probes are enabled (isDisabled is false). Writes over
    the given model file.

    Parameters
    ----------
    model_fpath : str
        Path to a model (.OSIM) file.

    """
    model = osm.Model(model_fpath)
    n_probes = model.getProbeSet().getSize()
    for i_probe in range(n_probes):
        model.updProbeSet().get(i_probe).setDisabled(False)
    model.print(model_fpath)

def strengthen_muscles(model_fpath, new_model_fpath, scale_factor):
    """Scales all muscles' maximum isometric force by `scale_factor`.

    Parameters
    ----------
    model_fpath : str
        Path to model (.OSIM) file for the model to strengthen.
    new_model_fpath : str
        Path to which to save the strengthened model.
    scale_factor : float
        All muscle optimal forces are scaled by this number.

    """
    m = osm.Model(model_fpath)
    for i_m in range(m.getMuscles().getSize()):
        m.updMuscles().get(i_m).setMaxIsometricForce(
                m.getMuscles().get(i_m).getMaxIsometricForce() * scale_factor)
    m.print(new_model_fpath)


def set_model_state_from_storage(model, storage, time, state=None):
    """Set the state of the model from a state described in a states Storage
    (.STO) file, at the specified time in the Storage file. Note that the model
    is not modified in any way; we just use the model to set the state.

    The storage should have beeng generated with a model that has the same
    exact states.

    Parameters
    ----------
    model : str or opensim.Model
        If str, a valid path to an OpenSim model file (.osim).
    storage : str or opensim.Storage
        If str, a valid path to a states Storage file.
    time : float
        A time, within the range of the times in the storage file, at which we
        should extract the state from the storage file.
    state : simtk.State
        If you don't want us to call `initSystem()` on the model, then give us
        a state!

    Returns
    -------
    state : simtk.State
        A state object that represents the state given in `storage` at time
        `time`.

    """
    if type(model) == str:
        model = osm.Model(model)
    if type(storage) == str:
        storage = osm.Storage(storage)

    if state == None:
        state = model.initState()

    n_states = storage.getSize()

    state_names = storage.getColumnLabels()

    # Interpolate the data to obtain the state (placed into sto_state) at the
    # specified time. Grab all the states (n_states). I'm assuming that these
    # state values are in the same order as is given by getStateIndex.
    sto_state = osm.ArrayDouble()
    sto_state.setSize(n_states)
    storage.getDataAtTime(time, n_states, sto_state)

    for i in range(state_names.getSize()):
        # I'm not even assuming that these
        # state values are returned in the same order given by state_names.
        if state_names.getitem(i) != 'time':
            sto_idx = storage.getStateIndex(state_names.getitem(i))
            model.setStateVariable(state, state_names.getitem(i),
                    sto_state.getitem(sto_idx))

    # TODO Maybe CAN rely on this.
    state.setTime(time)

    return state

def compute_state_dependent_quantity_in_time(model, storage, fcn):
    """This basically does the same thing as an OpenSim analysis. Compute the
    result of `fcn` for each time in the states_sto, and return the resulting
    array.

    Parameters
    ----------
    model : str or opensim.Model
        If str, a valid path to an OpenSim model file (.osim).
    storage : str or opensim.Storage
        If str, a valid path to a states Storage file.
    fcn : function
        This function must have a signature like:

            qty = fcn(model, state)

        where model is an opensim.Model, and state is a
        simtk.State. Note that you can grab the time via state.getTime().

    Returns
    -------
    time : list of float's
        All the times in the storage file.
    qty : list of float's
        This is the result of `qty` at all the times in the states Storage. It
        has the same length as a column in `storage`.

    """
    if type(model) == str:
        model = osm.Model(model)
    if type(storage) == str:
        storage = osm.Storage(storage)

    state = model.initSystem()

    sto_times = osm.ArrayDouble()
    storage.getTimeColumn(sto_times)

    time = sto_times.getSize() * [0]
    for i in range(sto_times.getSize()):
        time[i] = sto_times.getitem(i)

    qty = len(time) * [0]
    for i, t in enumerate(time):
        this_state = set_model_state_from_storage(model, storage, t,
                state=state)
        qty[i] = fcn(model, state)

    return time, qty


class Scale:
    """Wrapper of org.opensim.modeling.Scale, that adds a convenience
    constructor.

    """
    def __init__(self, body_name, x, y, z, scale_set=None):
        """
        Parameters
        ----------
        body_name : str
            org.opensim.modeling.Scale.setSegmentName(body_name)
        x, y, z : float
            org.opensim.modeling.Scale.setScaleFactors([x, y, z])
        scaleset : org.opensim.modeling.ScaleSet, optional
            A ScaleSet to adopt-and-append the org.opensim.modeling.Scale to.

        """
        self.scale = osm.Scale()
        self.scale.setSegmentName(body_name)
        self.scale.setScaleFactors([x, y, z])
        if scale_set:
            scale_set.cloneAndAppend(self.scale)

class Measurement:
    """Wrapper of org.opensim.modeling.Measurement with convenience methods.

    """

    def __init__(self, name, measurement_set=None):
        """
        Parameters
        ----------
        name : str
            org.opensim.modeling.Measurement.setName(name)
        measurement_set : org.opensim.modeling.MeasurementSet, optional
            A MeasurementSet to adopt-and-append the
            org.opensim.modeling.Measurement to.

        """
        self.measurement = osm.Measurement()
        self.measurement.setName(name)
        if measurement_set:
            measurement_set.adoptAndAppend(self.measurement)

    def add_bodyscale(self, name, axes='XYZ'):
        """Adds a BodyScale to the Measurement.

        Parameters
        ----------
        name : str
            org.opensim.modeling.BodyScale.setName(name)
        axes : str, optional
            e.g., 'X', 'XY', ...
            org.opensim.modeling.BodyScale.setAxisNames().
            Default is isometric.

        """
        bs = osm.BodyScale()
        bs.setName(name)
        axis_names = osm.ArrayStr()
        for ax in axes:
            axis_names.append(ax)
        bs.setAxisNames(axis_names)
        self.measurement.getBodyScaleSet().cloneAndAppend(bs)

    def add_bodyscale_bilateral(self, name, *args, **kwargs):
        """Adds a BodyScale to both sides of a model. If `name` is 'calf', then
        the same BodyScale is added to the two bodies 'calf_l' and 'calf_r'.

        Parameters
        ----------
        name : str
            Shared prefix of the body.
        axes : list of str's
            See `add_bodyscale`.
        """
        self.add_bodyscale('%s_l' % name, *args, **kwargs)
        self.add_bodyscale('%s_r' % name, *args, **kwargs)

    def add_markerpair(self, marker0, marker1):
        """Adds a MarkerPair to the Measurement's MarkerPairSet.

        Parameters
        ----------
        marker0 : str
            Name of the first marker in the pair.
        marker1 : str
            Name of the second marker in the pair.

        """
        mp = osm.MarkerPair()
        mp.setMarkerName(0, marker0)
        mp.setMarkerName(1, marker1)
        self.measurement.getMarkerPairSet().cloneAndAppend(mp)

    def add_markerpair_bilateral(self, marker0, marker1):
        """Adds two MarkerPair's to the Measurement's MarkerPairSet; assuming
        the name convention: if `marker0` is 'Heel', and `marker1` is 'Toe',
        then we add the following marker pairs: 'RHeel' and 'RToe', and 'LHeel'
        and 'LToe'.

        """
        self.add_markerpair('L%s' % marker0, 'L%s' % marker1)
        self.add_markerpair('R%s' % marker0, 'R%s' % marker1)


class IKTaskSet:
    """Wrapper of org.opensim.modeling.IKTaskSet with convenience methods.

    """
    def __init__(self, iktaskset=None):
        """Creates an org.opensim.modeling.IKTaskSet, or just uses the one
        provided, if provided.

        """
        if iktaskset:
            self.iktaskset = iktaskset
        else:
            self.iktaskset = osm.IKTaskSet()

    def add_ikmarkertask(self, name, do_apply, weight):
        """Creates an IKMarkerTask and appends it to the IKTaskSet.

        Parameters
        ----------
        name : str
            org.opensim.modeling.IKMarkerTask.setName(name)
        do_apply : bool
            org.opensim.modeling.IKMarkerTask.setApply(do_apply)
        weight : float
            org.opensim.modeling.IKMarkerTask.setWeight(weight)

        """
        ikt = osm.IKMarkerTask()
        ikt.setName(name)
        ikt.setApply(do_apply)
        ikt.setWeight(weight)
        self.iktaskset.cloneAndAppend(ikt)

    def add_ikmarkertask_bilateral(self, name, do_apply, weight):
        """Adds two IKMarkerTask's to the IKTaskSet.

        Parameters
        ----------
        name : str
            If 'name' is 'Elbow', then two tasks for markers 'LElbow' and
            'MElbow' will be added.
        do_apply, weight :
            See `add_ikmarkertask`.


        """
        self.add_ikmarkertask('L%s' % name, do_apply, weight)
        self.add_ikmarkertask('R%s' % name, do_apply, weight)









