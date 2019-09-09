import matplotlib.pyplot as plt
from collections import OrderedDict
from MACE.Visualization.Styles.TrackGroup import *


class SubplotStyle:

    def __init__(self, distance, internal_offset=0, x_multiplier=1.05,
                 y_multiplier=1.10,
                 x_start=0, y_start=0,
                 xaxis_visible=True, yaxis_visible=True,
                 spines_bottom_visible=True, spines_right_visible=True,
                 spines_left_visible=True, spines_top_visible=True,
                 track_group_style=None):

        self.distance = distance
        self.internal_offset = internal_offset
        self.x = x_start
        self.x_multiplier = x_multiplier
        self.y_multiplier = y_multiplier
        self.y = y_start

        self.xaxis_visible = xaxis_visible
        self.yaxis_visible = yaxis_visible

        self.spines_visibility = OrderedDict({'top':    spines_top_visible,
                                              'bottom': spines_bottom_visible,
                                              "right":  spines_right_visible,
                                              "left":   spines_left_visible})

        self.track_group_style = track_group_style

    def apply(self, axes=None):
        if axes is None:
            axes = plt.gca()

        axes.get_yaxis().set_visible(self.yaxis_visible)
        axes.get_xaxis().set_visible(self.xaxis_visible)
        #axes.xaxis.set_major_formatter(x_formatter)

        #subplot.spines['bottom'].set_color('none')

        for edge in self.spines_visibility:
            if self.spines_visibility[edge] is False:
                axes.spines[edge].set_color('none')


default_subplot_style = SubplotStyle(distance=5)
chromosome_subplot_style = SubplotStyle(distance=5, xaxis_visible=True, yaxis_visible=False, spines_bottom_visible=True,
                                        spines_right_visible=False, spines_left_visible=False, spines_top_visible=False)
