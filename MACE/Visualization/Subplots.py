
from MACE.Visualization.Styles.Subplot import SubplotStyle, default_subplot_style
import math
from collections import Iterable, OrderedDict

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


class Subplot(OrderedDict):

    def __init__(self, track_groups=None, y_start=0, x_start=0, x_end=1, style=default_subplot_style, title=None,
                 legend=None):
        if track_groups:
            OrderedDict.__init__(self, track_groups)
        else:
            OrderedDict.__init__(self)

        self.y_start = y_start
        self.y_end = None

        self.x_start = x_start
        self.x_end = x_end
        self.style = style
        self.title = title
        self.legend = legend

    def init_coordinates(self):
        y = self.y_start + self.style.internal_offset - self.style.distance

        for track_group_name in self:

            self[track_group_name].y_start = y + self.style.distance
            self[track_group_name].init_coordinates()
            y = self[track_group_name].y_end
            self.x_end = max(self.x_end, self[track_group_name].x_end)

        self.x_end = self.x_end * self.style.x_multiplier
        self.y_end = (y + self.style.internal_offset) * self.style.y_multiplier

        self.legend.x_start = self.x_end
        self.legend.y_start = self.y_end/2
        self.legend.x_size = self.x_end / self.legend.style.x_size_denominator

    def draw(self, axes=None):

        self.init_coordinates()
        for track_group in self:
            self[track_group].draw(axes=axes)

        self.style.apply()

        plt.xlim(xmin=self.x_start, xmax=self.x_end*1.1)
        plt.ylim(ymin=self.y_start, ymax=self.y_end)

        if self.title:
            plt.title(self.title)

        if self.legend:
            self.legend.draw()
