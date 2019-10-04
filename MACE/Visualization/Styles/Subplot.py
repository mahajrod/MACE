from collections import OrderedDict

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from MACE.Visualization.Styles.TrackGroup import *


class SubplotStyle:

    def __init__(self, distance, internal_offset=0, x_multiplier=1.05,
                 y_multiplier=1.10,
                 x_start=0, y_start=0,
                 x_logbase=None, y_logbase=None,
                 xaxis_visible=True, yaxis_visible=True,
                 spines_bottom_visible=True, spines_right_visible=True,
                 spines_left_visible=True, spines_top_visible=True,
                 track_group_style=None,
                 x_tick_type=None,
                 y_tick_type=None):

        self.distance = distance
        self.internal_offset = internal_offset
        self.x = x_start
        self.x_multiplier = x_multiplier
        self.y_multiplier = y_multiplier
        self.y = y_start

        self.x_logbase = x_logbase
        self.y_logbase = y_logbase

        self.xaxis_visible = xaxis_visible
        self.yaxis_visible = yaxis_visible

        self.spines_visibility = OrderedDict({'top':    spines_top_visible,
                                              'bottom': spines_bottom_visible,
                                              "right":  spines_right_visible,
                                              "left":   spines_left_visible})

        self.track_group_style = track_group_style

        self.x_tick_type = x_tick_type
        self.y_tick_type = y_tick_type
        self.x_tick_formatter = None
        self.y_tick_formatter = None

    def apply(self, axes=None, x_max=None, y_max=None):
        if axes is None:
            axes = plt.gca()

        axes.get_yaxis().set_visible(self.yaxis_visible)
        axes.get_xaxis().set_visible(self.xaxis_visible)
        #axes.xaxis.set_major_formatter(x_formatter)

        #subplot.spines['bottom'].set_color('none')

        for edge in self.spines_visibility:
            if self.spines_visibility[edge] is False:
                axes.spines[edge].set_color('none')

        if self.x_logbase:
            axes.set_xscale('log', basey=self.x_logbase)
        if self.y_logbase:
            axes.set_yscale('log', basey=self.y_logbase)

        if self.x_tick_type and x_max:
            self.x_tick_formatter = self.create_tick_formatter_function(x_max, tick_type=self.x_tick_type)
            axes.xaxis.set_major_formatter(self.x_tick_formatter)

        if self.y_tick_type and y_max:
            self.y_tick_formatter = self.create_tick_formatter_function(y_max, tick_type=self.y_tick_type)
            axes.yaxis.set_major_formatter(self.y_tick_formatter)

    @staticmethod
    def create_tick_formatter_function(max_value, tick_type="nucleotide"):
        max_val = max_value * 1.1
        if tick_type == "nucleotide":
            if max_val // (10 ** 9) > 2:
                def tick_formater(x, pos):
                    return '%1.1f Gbp' % (x * 1e-9)
            elif max_val // (10 ** 6) > 200:
                def tick_formater(x, pos):
                    return '%.0f Mbp' % (x * 1e-6)
            elif max_val // (10 ** 6) > 2:
                def tick_formater(x, pos):
                    return '%.1f Mbp' % (x * 1e-6)
            elif max_val // (10 ** 3) > 2:
                def tick_formater(x, pos):
                    return '%.1f kbp' % (x * 1e-3)
            else:
                def tick_formater(x, pos):
                    return '%i bp' % (int(x))

            return FuncFormatter(tick_formater)

        else:
            raise ValueError("ERROR!!! Tick formter for %s is not implemented yet!" % tick_type)


default_subplot_style = SubplotStyle(distance=5)
chromosome_subplot_style = SubplotStyle(distance=5, xaxis_visible=True, yaxis_visible=False, spines_bottom_visible=True,
                                        spines_right_visible=False, spines_left_visible=False, spines_top_visible=False,
                                        x_tick_type="nucleotide")

rainfall_subplot_style = SubplotStyle(distance=15, xaxis_visible=True, yaxis_visible=True, spines_bottom_visible=True,
                                      spines_right_visible=False, spines_left_visible=False, spines_top_visible=False,
                                      x_tick_type="nucleotide")
