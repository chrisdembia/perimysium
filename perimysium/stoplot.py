#!/usr/bin/env python

"""
http://wiki.scipy.org/Cookbook/EmbeddingInTraitsGUI
http://code.enthought.com/projects/traits/docs/html/_static/mpl_figure_editor.py
https://github.com/enthought/traitsui/blob/master/examples/demo/Standard_Editors/CheckListEditor_simple_demo.py
http://stackoverflow.com/questions/16663908/enthought-traits-ui-add-values-dynamically-to-values-trait-of-checklisteditor
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

from traits.api import HasTraits, Any, Instance, List, Str
from traitsui.api import View, Item, HSplit, CheckListEditor
from traitsui.wx.editor import Editor
from traitsui.wx.basic_editor_factory import BasicEditorFactory

from dataman import storage2numpy, TRCFile

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
    columns = List(Str)
    selector = List(editor=CheckListEditor(name='columns'))

    view = View(
            HSplit(
                Item('figure', editor=MPLFigureEditor(),
                    show_label=False
                    ),
                Item('selector', style='custom'),
                ),
            width=400,
            height=300,
            resizable=True,
            )

    def __init__(self, stofpath, *columns):
        HasTraits.__init__(self, trait_value=True)
        if stofpath.endswith('.sto') or stofpath.endswith('.mot'):
            self.data = storage2numpy(stofpath)
        elif stofpath.endswith('.trc'):
            self.data = TRCFile(stofpath).data
        axes = self.figure.add_subplot(111)
        columns = list(self.data.dtype.names)
        columns.remove('time')
        self.columns = columns
        #for name in self.data.dtype.names:
        #    if name != 'time' and (len(columns) == 0 or name in columns):
        #        axes.plot(self.data['time'], self.data[name], label=name)
        #axes.set_xlabel('time (s)')
        #axes.legend(loc='best')

def start_plotter(*args, **kwargs):
    """TODO"""
    plotter = StoragePlotter(*args, **kwargs)
    plotter.configure_traits()

if __name__ == '__main__':
    start_plotter(*sys.argv[1:])



