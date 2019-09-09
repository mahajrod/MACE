
from MACE.Visualization.Styles.TrackGroup import TrackGroupStyle, default_track_group_style
import math
from collections import Iterable, OrderedDict

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


class TrackGroup(OrderedDict):

    def __init__(self, tracks=None, y_start=None, x_start=0, x_end=1, style=default_track_group_style):
        if tracks:
            OrderedDict.__init__(self, tracks)
        else:
            OrderedDict.__init__(self)

        self.y_start = y_start
        self.y_end = None

        self.x_start = x_start
        self.x_end = x_end

        self.style = style

    def init_coordinates(self):
        y = self.y_start + self.style.internal_offset - self.style.distance

        for track_name in self:
            self[track_name].y_start = y + self.style.distance
            self.x_end = max(self.x_end, self[track_name].x_end)
            y += self[track_name].style.height

        self.x_end = self.x_end * self.style.x_multiplier
        self.y_end = y + self.style.internal_offset

    def draw(self, axes=None):
        self.init_coordinates()
        for track_name in self:
            self[track_name].draw(axes=axes)

