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

import platform

def running_in_jython():
    return platform.system() == 'Java'

if running_in_jython():
    import org.opensim.modeling as osm
else:
    import opensim as osm

def printobj(obj, fname):
    if running_in_jython():
        exec('obj.print(fname)')
    else:
        obj.printToXML(fname)

def replace_thelen_muscles_with_millardequilibrium_muscles(model):
    """Edits the given model:

    1. Removes all Thelen muscles.
    2. Replaces them with Millard equilibrium muscles.

    """
    # 1) Loop through all the forces in the model. Copy them over to a new
    # ForceSet, unless it's a Thelen muscle. If it's a Thelen muscle, create a
    # similar Millard muscle.
    fset = osm.ForceSet()
    for i_force in range(model.getForceSet().getSize()):
        force = model.getForceSet().get(i_force)
        old = osm.Thelen2003Muscle.safeDownCast(force)

        # If not a Thelen muscle, just clone it.
        if not old:

            fset.cloneAndAppend(force)

        else:

            new = osm.Millard2012EquilibriumMuscle()
    
            new.setName(old.getName())
    
            # GeometryPath.
            # --------------
            # TODO geometry path wrap, visible object, ...
            old_geopath = old.getGeometryPath()
            new_geopath = new.updGeometryPath()
            for i_pt in range(old_geopath.getPathPointSet().getSize()):
                old_pt = old_geopath.getPathPointSet().get(i_pt)
                new_geopath.updPathPointSet().cloneAndAppend(old_pt)
    
            # Parameters.
            # -----------
            def transfer(new, old, name):
                exec('new.set_%s(old.get_%s())' % (name, name))
    
            transfer(new, old, 'max_isometric_force')
            transfer(new, old, 'optimal_fiber_length')
            transfer(new, old, 'tendon_slack_length')
            transfer(new, old, 'pennation_angle_at_optimal')
            transfer(new, old, 'max_contraction_velocity')
            transfer(new, old, 'activation_time_constant')
            transfer(new, old, 'deactivation_time_constant')
    
            transfer(new, old, 'default_activation')
            transfer(new, old, 'default_fiber_length')

            fset.cloneAndAppend(new)

    # 2) clearAndDestroy the model's ForceSet().
    model.updForceSet().clearAndDestroy()

    # 3) Add all forces from the new ForceSet to the model.
    for i_force in range(fset.getSize()):
        model.updForceSet().cloneAndAppend(fset.get(i_force))


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

twitch_ratios_2392 = {
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

twitch_ratios_1018 = {
            'hamstrings': 0.49, 'bifemsh': 0.53, 'glut_max': 0.55,
            'iliopsoas': 0.50, 'rect_fem': 0.39, 'vasti': 0.50,
            'gastroc': 0.54, 'soleus': 0.80,
            'tib_ant': 0.70}

# For the muscles that are divided in OpenSim across multiple paths,
# divide the published mass evenly between them.
# Volumes are in cm^3.
_H2014glut_med_volume = 323.2
_H2014glut_min_volume = 104.5
_H2014add_mag_volume = 559.8
_H2014glut_max_volume = 849.0
_H2014small_ext_rotators = 16.1
_H2014extensors_volume = 102.3
_H2014peroneals_volume = 130.8
_Handsfield2014_muscle_volumes = {
        'glut_med1_': _H2014glut_med_volume / 3.0,
        'glut_med2_': _H2014glut_med_volume / 3.0,
        'glut_med3_': _H2014glut_med_volume / 3.0,
        'glut_min1_': _H2014glut_min_volume / 3.0,
        'glut_min2_': _H2014glut_min_volume / 3.0,
        'glut_min3_': _H2014glut_min_volume / 3.0,
        'semimem_': 245.4,
        'semiten_': 186.0,
        'bifemlh_': 206.5,
        'bifemsh_': 100.1,
        'sar_': 163.7,
        'add_long_': 162.1,
        'add_brev_': 104.0,
        'add_mag1_': _H2014add_mag_volume / 3.0,
        'add_mag2_': _H2014add_mag_volume / 3.0,
        'add_mag3_': _H2014add_mag_volume / 3.0,
        'tfl_': 64.9,
        'pect_': 66.3,
        'grac_': 104.0,
        'glut_max1_': _H2014glut_max_volume / 3.0,
        'glut_max2_': _H2014glut_max_volume / 3.0,
        'glut_max3_': _H2014glut_max_volume / 3.0,
        'iliacus_': 176.8,
        'psoas_': 274.8,
        'quad_fem_': 32.4,
        'gem_': _H2014small_ext_rotators,
        'peri_': 42.8, # piriformis.
        'rect_fem_': 269.0,
        'vas_med_': 423.6,
        'vas_int_': 270.5,
        'vas_lat_': 830.9,
        'med_gas_': 257.4,
        'lat_gas_': 150.0,
        'soleus_': 438.2,
        'tib_post_': 104.8,
        'flex_dig_': 30.0,
        'flex_hal_': 78.8,
        'tib_ant_': 135.2,
        'per_brev_': _H2014peroneals_volume / 2.0,
        'per_long_': _H2014peroneals_volume / 2.0,
        'per_tert_': _H2014extensors_volume / 3.0,
        'ext_dig_': _H2014extensors_volume / 3.0,
        'ext_hal_': _H2014extensors_volume / 3.0,
        #'ercspn_': ,
        #'intobl_': ,
        #'extobl_': ,
        }

# In kg.
Handsfield2014_muscle_masses = dict()
muscle_density = 0.001056 # kg / cm^3
for key, val in _Handsfield2014_muscle_volumes.items():
    Handsfield2014_muscle_masses[key] = muscle_density * val

"""
_W2009glut_max_mass = 547.2
_W2009glut_med_mass = 273.5
_Ward2009_muscle_volumes_grams = {
        'glut_med1_': _W2009glut_med_mass / 3.0,
        'glut_med2_': _W2009glut_med_mass / 3.0,
        'glut_med3_': _W2009glut_med_mass / 3.0,
        'glut_min1_': _W2009glut_,
        'glut_min2_':,
        'glut_min3_':,
        'semimem_':,
        'semiten_':,
        'bifemlh_':,
        'bifemsh_':,
        'sar_':,
        'add_long_':,
        'add_brev_':,
        'add_mag1_':,
        'add_mag2_':,
        'add_mag3_':,
        'tfl_':,
        'pect_':,
        'grac_':,
        'glut_max1_':,
        'glut_max2_':,
        'glut_max3_':,
        'iliacus_': 113.7,
        'psoas_': 97.7,
        'quad_fem_':,
        'gem_':,
        'peri_':, # piriformis.
        'rect_fem_':,
        'vas_med_':,
        'vas_int_':,
        'vas_lat_':,
        'med_gas_':,
        'lat_gas_':,
        'soleus_':,
        'tib_post_':,
        'flex_dig_':,
        'flex_hal_':,
        'tib_ant_':,
        'per_brev_':,
        'per_long_':,
        'per_tert_':,
        'ext_dig_':,
        'ext_hal_':,
        #'ercspn_': ,
        #'intobl_': ,
        #'extobl_': ,
        }

Ward2009_muscle_volumes = dict()
for key, val in _Ward2009_muscle_volumes_grams.items():
    Ward2009_muscle_volumes[key] = 0.001 * val
"""

def add_metabolics_probes(model, twitch_ratio_set='gait2392',
        activationMaintenanceRateOn=True,
        shorteningRateOn=True,
        basalRateOn=False,
        mechanicalWorkRateOn=True,
        muscle_masses=None,
        exclude=[]):
    """Adds Umberger2010MuscleMetabolicsProbes to an OpenSim model. Adds a
    probe for each muscle, as well as a whole-body probe that returns
    cumulative energy expended across all muscles in the model. When possible,
    we use published twitch ratios for the probes. For muscles for which we do
    not have data, we use a twitch ratio of 0.5. This method doesn't return
    anything; the model given to the method is simply modified.

    Parameters
    ----------
    model : org.opensim.modeling.Model
        An OpenSim Model that has muscles.
    twitch_ratio_set : float or str ('gait2392' or 'gait1018')
        The experimental data set to use for the model, depending on the model
        we are adding probes to. If a float, use that constant value for all
        muscles.
    activationMaintenanceRateOn : bool, optional
    shorteningRateOn : bool, optional
    basalRateOn : bool, optional
    mechanicalWorkRateOn : bool, optional
    muscle_masses : str, or dict; optional
        * str: 'Handsfield2014' or 'Ward2009' to use lower body masses from the
          respective paper. NOTE: this set of muscle masses does NOT contain
          ercspn, intobl, or extobl masses. Consider excluding those muscles.
        * dict: For muscles in this dict, use the given value as the muscle's
          mass. If the muscle is not specified in this dict, compute the
          muscle's mass from the model's muscle properties.
    exclude : list of str's, optional
        List of muscle names to exclude.

    """
    # Twitch ratios
    # -------------
    if type(twitch_ratio_set) == float:
        twitchRatios = twitch_ratio_set
    elif twitch_ratio_set == 'gait2392':
        twitchRatios = twitch_ratios_2392
    elif twitch_ratio_set == 'gait1018':
        twitchRatios = twitch_ratios_1018
    else:
        raise Exception("Invalid value for `twitch_ratio_set`.")

    # Muscle masses
    # -------------
    if muscle_masses:
        if muscle_masses == 'Handsfield2014':
            muscle_masses = Handsfield2014_muscle_masses
        else:
            raise Exception("Unexpected muscle_masses {}.".format(
                muscle_masses))
    else:
        muscle_masses = dict()

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

        if not (thisMuscle.getName() in exclude):

            # Get the slow-twitch ratio from the data we read earlier. Start
            # with the default value.
            slowTwitchRatio = defaultTwitchRatio
    
            # Set the slow-twitch ratio to the physiological value, if it is
            # known.
            if type(twitch_ratio_set) == float:
                slowTwitchRatio = twitchRatios
            else:
                for key, val in twitchRatios.items():
                    if thisMuscle.getName().startswith(key) and val != -1:
                        slowTwitchRatio = val
    
            # Add this muscle to the whole-body probe. The arguments are muscle
            # name, slow-twitch ratio, and muscle mass. Note that the muscle
            # mass is ignored unless we set useProvidedMass to True.
            wholeBodyProbe.addMuscle(thisMuscle.getName(),
                                     slowTwitchRatio)
    
            # If we are given a muscle mass, use it in the probe.
            for key, val in muscle_masses.items():
                if thisMuscle.getName().startswith(key):
                    wholeBodyProbe.useProvidedMass(thisMuscle.getName(),
                            val)


def add_bhargava_metabolic_probes(model, twitch_ratio_set='gait2392',
        activationRateOn=True,
        maintenanceRateOn=True,
        shorteningRateOn=True,
        basalRateOn=False,
        workRateOn=True,
        muscle_masses=None,
        exclude=[],
        ):
    """Adds Bhargava2004MuscleMetabolicsProbe's to an OpenSim model. Adds a
    probe for each muscle, as well as a whole-body probe that returns
    cumulative energy expended across all muscles in the model. When possible,
    we use published twitch ratios for the probes. For muscles for which we do
    not have data, we use a twitch ratio of 0.5. This method doesn't return
    anything; the model given to the method is simply modified.

    Parameters
    ----------
    model : org.opensim.modeling.Model
        An OpenSim Model that has muscles.
    twitch_ratio_set : float or str ('gait2392' or 'gait1018')
        The experimental data set to use for the model, depending on the model
        we are adding probes to. If a float, use that constant value for all
        muscles.
    activationRateOn : bool, optional
    maintenanceRateOn : bool, optional
    shorteningRateOn : bool, optional
    basalRateOn : bool, optional
    workRateOn : bool, optional
    muscle_masses : str, or dict; optional
        * str: 'Handsfield2014' or 'Ward2009' to use lower body masses from the
          respective paper. NOTE: this set of muscle masses does NOT contain
          ercspn, intobl, or extobl masses. Consider excluding those muscles.
        * dict: For muscles in this dict, use the given value as the muscle's
          mass. If the muscle is not specified in this dict, compute the
          muscle's mass from the model's muscle properties.
    exclude : list of str's, optional
        List of muscle names to exclude.

    """
    # Twitch ratios
    # -------------
    if type(twitch_ratio_set) == float:
        twitchRatios = twitch_ratio_set
    elif twitch_ratio_set == 'gait2392':
        twitchRatios = twitch_ratios_2392
    elif twitch_ratio_set == 'gait1018':
        twitchRatios = twitch_ratios_1018
    else:
        raise Exception("Invalid value for `twitch_ratio_set`.")

    # Muscle masses
    # -------------
    if muscle_masses:
        if muscle_masses == 'Handsfield2014':
            muscle_masses = Handsfield2014_muscle_masses
        else:
            raise Exception("Unexpected muscle_masses {}.".format(
                muscle_masses))
    else:
        muscle_masses = dict()

    # Parameters used for all probes.
    # -------------------------------
    activation_constant_slow_twitch = 40.0
    activation_constant_fast_twitch = 133.0
    maintenance_constant_slow_twitch = 74.0
    maintenance_constant_fast_twitch = 111.0

    defaultTwitchRatio = 0.5

    # Whole-body probe.
    # -----------------
    wholeBodyProbe = osm.Bhargava2004MuscleMetabolicsProbe(
            activationRateOn,
            maintenanceRateOn,
            shorteningRateOn,
            basalRateOn,
            workRateOn)
    wholeBodyProbe.setOperation("value")
    wholeBodyProbe.set_report_total_metabolics_only(False)

    model.addProbe(wholeBodyProbe)
    wholeBodyProbe.setName("metabolic_power_bhar")

    for iMuscle in range(model.getMuscles().getSize()):
        thisMuscle = model.getMuscles().get(iMuscle)

        if not (thisMuscle.getName() in exclude):

            slowTwitchRatio = defaultTwitchRatio
    
            # Set the slow-twitch ratio to the physiological value, if it is
            # known.
            if type(twitch_ratio_set) == float:
                slowTwitchRatio = twitchRatios
            else:
                for key, val in twitchRatios.items():
                    if thisMuscle.getName().startswith(key) and val != -1:
                        slowTwitchRatio = val
    
            wholeBodyProbe.addMuscle(thisMuscle.getName(), slowTwitchRatio,
                    activation_constant_slow_twitch,
                    activation_constant_fast_twitch,
                    maintenance_constant_slow_twitch,
                    maintenance_constant_fast_twitch)
    
            # If we are given a muscle mass, use it in the probe.
            for key, val in muscle_masses.items():
                if thisMuscle.getName().startswith(key):
                    wholeBodyProbe.useProvidedMass(thisMuscle.getName(),
                            val)

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
    printobj(model, model_fpath)

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
    printobj(m, new_model_fpath)


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

def analysis(model, storage, fcn, times=None):
    """This basically does the same thing as an OpenSim analysis. Compute the
    result of `fcn` for each time in the states_sto, using the model's state,
    and return the resulting array.

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
    times : array_like of float's
        Times at which to evaluate `fcn`.

    Returns
    -------
    times : list of float's
        The times corresponding to the evaluations of `fcn`.
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

    if times == None:
        times = sto_times.getSize() * [0]
        for i in range(sto_times.getSize()):
            times[i] = sto_times.getitem(i)

    qty = len(times) * [0]
    for i, t in enumerate(times):
        this_state = set_model_state_from_storage(model, storage, t,
                state=state)
        qty[i] = fcn(model, state)

    return times, qty


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









