#!/usr/bin/env python

"""
http://wiki.scipy.org/Cookbook/EmbeddingInTraitsGUI
http://code.enthought.com/projects/traits/docs/html/_static/mpl_figure_editor.py
"""

import sys

from numpy import deg2rad, rad2deg

import wx

import matplotlib
# We want matplotlib to use a wxPython backend
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wx import NavigationToolbar2Wx

from traits.api import HasTraits, Any, Instance
from traitsui.api import View, Item
from traitsui.wx.editor import Editor
from traitsui.wx.basic_editor_factory import BasicEditorFactory

from dataman import storage2numpy

class _MPLFigureEditor(Editor):

    scrollable  = True

    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()
        
    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # The panel lets us add additional controls.
        panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        # matplotlib commands to create a canvas
        mpl_control = FigureCanvas(panel, -1, self.value)
        sizer.Add(mpl_control, 1, wx.LEFT | wx.TOP | wx.GROW)
        toolbar = NavigationToolbar2Wx(mpl_control)
        sizer.Add(toolbar, 0, wx.EXPAND)
        self.value.canvas.SetMinSize((10,10))
        return panel

class MPLFigureEditor(BasicEditorFactory):

    klass = _MPLFigureEditor

class StoragePlotter(HasTraits):

    figure = Instance(Figure, ())

    view = View(Item('figure', editor=MPLFigureEditor(),
                            show_label=False),
                    width=400,
                    height=300,
                    resizable=True)

    def __init__(self, stofpath):
        HasTraits.__init__(self, trait_value=True)
        self.data = storage2numpy(stofpath)
        axes = self.figure.add_subplot(111)
        for name in self.data.dtype.names:
            if name != 'time':
                axes.plot(self.data['time'], self.data[name])

def start_gui(*args, **kwargs):
    '''Start the GUI. The GUI automatically creates a Human, and lets the user
    modify its configuration and observe the resulting change in the human's
    inertia properties.

    Parameters
    ----------
    meas_in : str, optional
        The filename of a measurements file to use for the human.
    '''
    YeadonGUI(*args, **kwargs).configure_traits()

def start_plotter(*args, **kwargs):
    """TODO"""
    plotter = StoragePlotter(*args, **kwargs)
    plotter.configure_traits()

if __name__ == '__main__':
    start_plotter(sys.argv[1])



