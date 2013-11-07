from scipy.integrate import ode


class Simulation(ode):
    def __init__(self, opensim_model, controller_func=None):
        super(Simulation, self).__init__(self.f)
        self.opensim_model = opensim_model
        if self.opensim_model.getWorkingState().getNY() == 0:
            self.opensim_state = self.opensim_model.initSystem()
        else:
            self.opensim_state = self.opensim_model.updWorkingState()

        self.controller_func = controller_func

    def f(self, t, y):
        """
        Parameters
        ----------
        t : float
            Time in the integration.
        y : array_like of float's
            The state given to us by the integrator. We will set the OpenSim state
            to be this state.
        controller_fcn : function
            TODO
    
        Returns
        -------
        state_derivatives : `numpy.array`
            Time derivatives of the states `y`.
    
        """
        # For looping.
        num_states = self.opensim_state.getNY()
    
        # Update the OpenSim state to reflect the state given to us by the
        # integrator.
        self.opensim_state.setTime(t)
        for i_state in range(num_states):
            self.opensim_state.updY().set(i_state, x[i_state + 1])
    
        # Compute derivatives of the states from the model.
        self.opensim_model.computeStateVariableDerivatives(self.opensim_state)
    
        # Let the user control the model, given its updated state.
        if controller_fcn:
            controls_vector = controller_fcn(opensim_model, self.opensim_state)
            self.opensim_model.setControls(self.opensim_state, controls_vector)
    
        # Update the derivatives again, given the controls.
        self.opensim_model.computeStateVariableDerivatives(self.opensim_state)
    
        # Format the state derivatives for python.
        state_derivatives = np.zeros((num_states, 0))
        for i_state in range(num_states):
            state_derivatives[i_state] = \
                    self.opensim_state.getYDot().get(i_state)
    
        return state_derivatives


