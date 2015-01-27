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

from traits.api import HasTraits, Any, Instance, List, Str, on_trait_change
from traitsui.api import View, Item, HGroup, CheckListEditor
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
    avail_columns = List(Str)
    columns = List(
            editor=CheckListEditor(name='avail_columns',
                format_func=lambda x: x))

    view = View(
            HGroup(
                Item('figure', editor=MPLFigureEditor(),
                    show_label=False
                    ),
                Item('columns', style='custom'),
    #            scrollable=True,
                ),
            width=700,
            height=400,
            resizable=True,
            )

    def __init__(self, stofpath, *args):
        HasTraits.__init__(self, trait_value=True)
        if stofpath.endswith('.sto') or stofpath.endswith('.mot'):
            if 'StaticOptimization' in stofpath:
                self.data = storage2numpy(stofpath, excess_header_entries=2)
            else:
                self.data = storage2numpy(stofpath)
        elif stofpath.endswith('.trc'):
            self.data = TRCFile(stofpath).data
        avail_columns = list(self.data.dtype.names)
        avail_columns.remove('time')
        self.avail_columns = avail_columns

        self.axes = self.figure.add_subplot(111)

        for arg in args:
            self.columns.append(arg)
            self._columns_changed()

    @on_trait_change('columns')
    def _columns_changed(self):
        self.axes.cla()
        for name in self.columns:
            self.axes.plot(self.data['time'], self.data[name], label=name)
        self.axes.set_xlabel('time (s)')
        self.axes.legend(loc='best')

def start_plotter(*args, **kwargs):
    """TODO"""
    plotter = StoragePlotter(*args, **kwargs)
    plotter.configure_traits()

if __name__ == '__main__':
    start_plotter(*sys.argv[1:])



