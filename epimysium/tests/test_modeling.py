import os

import org.opensim.modeling as osm

from perimysium import modeling

parentdir = os.path.abspath(os.path.dirname(__file__))

def test_set_model_state_from_storage():

    m = osm.Model(os.path.join(parentdir, 'double_pendulum.osim'))
    state_act = m.initSystem()
    sto = osm.Storage(os.path.join(parentdir, 'double_pendulum_states.sto'))
    state_des = modeling.set_model_state_from_storage(m, sto, 0, state_act)
    q1_des = m.getCoordinateSet().get('q1').getValue(state_des)
    q1_act = osm.ArrayDouble()
    q1_act.setSize(sto.getSize())
    sto.getDataAtTime(0.0, sto.getSize(), q1_act)
    assert q1_act.getitem(0) == q1_des

def test_analysis():

    def fcn(model, time, state):
        return state.getTime()
    
    m = osm.Model(os.path.join(parentdir, 'double_pendulum.osim'))
    sto = osm.Storage(os.path.join(parentdir, 'double_pendulum_states.sto'))
    t, qty = modeling.analysis(m, sto, fcn)

if __name__ == '__main__':
    test_set_model_state_from_storage()
    test_analysis()
