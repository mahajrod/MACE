
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
                 legend=None, axes=None, type="track", auto_scale=False, x_scale_factor=1,
                 y_scale_factor=1, figure_x_y_ratio=None,
                 xmax_multiplier=1.05, ymax_multiplier=1.1,
                 xmin_multiplier=0.001, ymin_multiplier=0.001):
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

        # TODO: add autoscale implementation
        #self.auto_scale = auto_scale
        self.x_scale_factor = x_scale_factor
        self.y_scale_factor = y_scale_factor
        self.x_y_ratio = None
        self.figure_x_y_ratio = figure_x_y_ratio

        self.xmax_multiplier = xmax_multiplier
        self.ymax_multiplier = ymax_multiplier

        self.xmin_multiplier = xmin_multiplier
        self.ymin_multiplier = ymin_multiplier

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

            #if self.auto_scale:
            # self.y_scale_factor = 1
            self.x_y_ratio = self.x_end / self.y_end
            for track_group_name in self:
                self[track_group_name].subplot_x_y_ratio = self.x_y_ratio
                self[track_group_name].figure_x_y_ratio = self.figure_x_y_ratio
                for track_name in self[track_group_name]:
                    self[track_group_name][track_name].figure_x_y_ratio = self.figure_x_y_ratio
                    self[track_group_name][track_name].subplot_x_y_ratio = self.x_y_ratio
            #print(self.x_scale_factor)

            # TODO: rewrite coordinate calculation for legends. Create method init_coordinates in Legend class
            if isinstance(self.legend, (CoverageLegend, DensityLegend)):
                legend_height = (len(self.legend.thresholds) + 3) * self.legend.element_size
            else:
                legend_height = None

            if self.legend:
                self.legend.x_start = self.x_end

                if legend_height:
                    self.legend.y_start = (self.y_end - legend_height) / 2
                else:
                    self.legend.y_start = self.y_end/2
                self.legend.x_size = self.x_end / self.legend.style.x_size_denominator

                self.legend.init_coordinates()
                self.x_end = self.legend.x_end
        elif self.type == "plot":
            pass

    def draw(self, axes=None):
        axes_to_use = axes if axes else self.axes if self.axes else plt.gca()
        self.axes = axes_to_use
        self.init_coordinates()

        for track_group in self:
            self[track_group].draw(axes=axes_to_use)

        #for track_group in self:
        #    self[track_group].draw_borders(axes=axes_to_use)

        self.style.apply(x_max=self.x_end, y_max=self.y_end, axes=axes_to_use)

        plt.xlim(xmin=self.x_start - (self.x_end * self.xmin_multiplier), xmax=self.x_end )
        plt.ylim(ymin=self.y_start - (self.y_end * self.ymin_multiplier), ymax=self.y_end * self.ymax_multiplier)

        if self.title:
            plt.title(self.title, fontsize=self.style.title_fontsize, fontweight=self.style.title_fontweight)

        if self.legend:
            self.legend.draw()

    def hide(self, axes=None):
        axes_to_use = axes if axes else self.axes if self.axes else plt.gca()
        axes_to_use.set_axis_off()
        #axes_to_use.get_xaxis().set_visible(False)
        #axes_to_use.get_yaxis().set_visible(False)