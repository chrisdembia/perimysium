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

"""Modifying OpenSim Model's.

"""

# TODO currently, this module assumes this package is installed via Jython, and
# that the OpenSim Jython wrapping is on the Jython CLASSPATH. When an OpenSim
# Python wrapping is available, we can use CPython instead.

import org.opensim.modeling as osm

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

