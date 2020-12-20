
from MACE.Visualization.Styles.Subplot import SubplotStyle, default_subplot_style
from MACE.Visualization.Legends import DensityLegend, CoverageLegend
import math
from collections import Iterable, OrderedDict

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


class Subplot(OrderedDict):

    def __init__(self, track_groups=None, y_start=None, y_end=None, x_start=None, x_end=None,
                 style=default_subplot_style, title=None,
                 legend=None, axes=None, type="track"):
        if track_groups:
            OrderedDict.__init__(self, track_groups)
        else:
            OrderedDict.__init__(self)

        self.type = type
        self.y_start = y_start
        self.y_end = y_end
        self.x_start = x_start
        self.x_end = x_end

        if self.type == "track":
            if self.y_start is None:
                self.y_start = 0
            if self.x_start is None:
                self.x_start = 0
            if self.x_end is None:
                self.x_end = 1

        self.style = style
        self.title = title
        self.legend = legend
        self.axes = axes

    def init_coordinates(self):
        if self.type == "track":
            y = self.y_start + self.style.internal_offset - self.style.distance

            for track_group_name in self:

                self[track_group_name].y_start = y + self.style.distance
                self[track_group_name].init_coordinates()
                y = self[track_group_name].y_end
                self.x_end = max(self.x_end, self[track_group_name].x_end)

            self.x_end = self.x_end * self.style.x_multiplier
            self.y_end = (y + self.style.internal_offset) * self.style.y_multiplier

            if isinstance(CoverageLegend, self.legend) or isinstance(DensityLegend, self.legend):
                legend_height = (len(self.legend.thresholds) + 3) * self.legend.element_size
            else:
                legend_height = None

            self.legend.x_start = self.x_end

            if legend_height:
                self.legend.y_start = (self.y_end - legend_height) / 2
            else:
                self.legend.y_start = self.y_end/2
            self.legend.x_size = self.x_end / self.legend.style.x_size_denominator

        elif self.type == "plot":
            pass

    def draw(self, axes=None):
        axes_to_use = axes if axes else self.axes if self.axes else plt.gca()
        self.axes = axes_to_use
        self.init_coordinates()

        for track_group in self:
            self[track_group].draw(axes=axes_to_use)

        self.style.apply(x_max=self.x_end, y_max=self.y_end, axes=axes_to_use)

        plt.xlim(xmin=self.x_start, xmax=self.x_end*1.1)
        plt.ylim(ymin=self.y_start, ymax=self.y_end*1.1)

        if self.title:
            plt.title(self.title)

        if self.legend:
            self.legend.draw()

    def hide(self, axes=None):
        axes_to_use = axes if axes else self.axes if self.axes else plt.gca()
        axes_to_use.set_axis_off()
        #axes_to_use.get_xaxis().set_visible(False)
        #axes_to_use.get_yaxis().set_visible(False)