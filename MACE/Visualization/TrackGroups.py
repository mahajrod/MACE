
from MACE.Visualization.Styles.TrackGroup import TrackGroupStyle, default_track_group_style
import math
from collections import Iterable, OrderedDict

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


class TrackGroup(OrderedDict):

    def __init__(self, tracks=None, y_start=None, x_start=0, x_end=1, style=default_track_group_style,
                 label=None):
        if tracks:
            OrderedDict.__init__(self, tracks)
        else:
            OrderedDict.__init__(self)

        self.y_start = y_start
        self.y_end = None

        self.x_start = x_start
        self.x_end = x_end

        self.style = style
        self.label = label
        self.track_label_param_list = None

    def init_coordinates(self):
        y = self.y_start + self.style.internal_offset - self.style.distance
        self.track_label_param_list = [[len(self[track_name].label) if (self[track_name].label and self[track_name].style.show_label) else 0,
                                       self[track_name].style.label_fontsize] for track_name in self]
        for track_name in self:
            self[track_name].y_start = y + self.style.distance
            self.x_end = max(self.x_end, self[track_name].x_end)
            y += self[track_name].style.height

        self.x_end = self.x_end * self.style.x_multiplier
        self.y_end = y + self.style.internal_offset

    def draw(self, axes=None, style=None, label_shift=0):
        self.init_coordinates()

        used_style = style if style else self.style
        current_subplot = axes if axes else plt.gca()

        for track_name in self:
            self[track_name].draw(axes=axes)

        #x_text_offset = max(list(map(lambda s: s[0]*s[1], self.track_label_param_list)))
        # TODO: add automatic calculation for x axis

        label_shift_full = (-100 if max(list(map(lambda s: s[0]*s[1], self.track_label_param_list))) else 0) + label_shift
        if self.label and used_style.show_label:
            current_subplot.annotate(self.label, xy=(0, (self.y_start + self.y_end)/2 + self.style.label_y_shift), xycoords='data',
                                     fontsize=self.style.label_fontsize,
                                     fontweight=self.style.label_fontweight,
                                     xytext=(self.style.label_x_shift + label_shift_full,
                                             (self.y_end - self.y_start) / 2), textcoords='offset points',
                                     ha=self.style.label_hor_aln, va=self.style.label_vert_aln)

