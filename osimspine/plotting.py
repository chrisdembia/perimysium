"""Methods to aid with plotting data from OpenSim.

"""

import pylab as pl
import tables

def bipedal_time_axis(h5fname, h5sim_node, speed, times):
    """Creates a matplotlib plot of a biped at the times specified, walking
    along a time axis at a constant speed. The right leg is colored red.

    Parameters
    ----------

    """
    h5file = tables.openFile(h5fname)
    h5

    h5file.close()

