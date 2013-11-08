import os

import numpy as np

from opensim import Model
from epimysium.simulation import Simulation

parentdir = os.path.abspath(os.path.dirname(__file__))

def test_doublependulum_uncontrolled():
    """ TODO
    """
    model = Model(os.path.join(parentdir, 'double_pendulum.osim'))
    model.initSystem()
    sim = Simulation(model)
    sim.set_integrator('vode')
    sim.set_initial_value([np.pi/ 4, 0, 0, 0])
    final_time = 5.0 # seconds.
    time_step = 0.01 # seconds.
    # Excluding the initial time:
    n_time_steps = 500
    time_step = final_time / n_time_steps
    actual_y = np.empty((n_time_steps, sim.num_states))
    for i in range(n_time_steps):
        if sim.successful():
            sim.integrate(sim.t + time_step)
            # sim.y is updated with the new state.
            actual_y[i, :] = sim.y
        else:
            raise Exception("Integration failed.")
    des_y = np.loadtxt(os.path.join(parentdir,
        'doublependulum_uncontrolled_des_states.txt'))
    np.testing.assert_allclose(actual_y, des_y)
