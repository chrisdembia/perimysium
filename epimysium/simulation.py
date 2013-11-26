"""Forward-integration of an OpenSim model, using integrators in SciPy.

"""

import numpy as np
from scipy.integrate import ode

class Simulation(ode):
    def __init__(self, opensim_model, controller_func=None):
        """ Initializes all states to zero. The initial states can be changed
        via `set_initial_value`.

        """
        super(Simulation, self).__init__(self.f)
        self.opensim_model = opensim_model
        if self.opensim_model.getWorkingState().getNY() == 0:
            self.opensim_state = self.opensim_model.initSystem()
        else:
            self.opensim_state = self.opensim_model.updWorkingState()
        self._num_states = self.opensim_state.getNY()
        self.set_initial_value(np.zeros(self.num_states))

        self._controller_func = controller_func

    @property
    def num_states(self):
        return self._num_states

    def f(self, t, y):
        """See `scipy.integrate.ode` documentation.
        """
        # Update the OpenSim state to reflect the state given to us by the
        # integrator.
        self.opensim_state.setTime(t)
        for i_state in range(self.num_states):
            self.opensim_state.updY().set(i_state, y[i_state])
    
        # Compute derivatives of the states from the model.
        self.opensim_model.computeStateVariableDerivatives(self.opensim_state)
    
        # Let the user control the model, given its updated state.
        if self._controller_func:
            controls_vector = controller_fcn(opensim_model, self.opensim_state)
            self.opensim_model.setControls(self.opensim_state, controls_vector)
    
        # Update the derivatives again, given the controls.
        self.opensim_model.computeStateVariableDerivatives(self.opensim_state)
    
        # Format the state derivatives for python.
        state_derivatives = np.empty(self.num_states)
        for i_state in range(self.num_states):
            state_derivatives[i_state] = \
                    self.opensim_state.getYDot().get(i_state)
    
        return state_derivatives


