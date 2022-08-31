
from MACE.Visualization.Styles.TrackGroup import TrackGroupStyle, default_track_group_style
import math
from collections import Iterable, OrderedDict

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Polygon


class TrackGroup(OrderedDict):

    def __init__(self, tracks=None, y_start=None, x_start=0, x_end=1, style=default_track_group_style,
                 label=None, x_scale_factor=1, y_scale_factor=1, auto_scale=False,
                 subplot_x_y_ratio=None, figure_x_y_ratio=None, highlight=False):
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

        self.x_scale_factor = x_scale_factor
        self.y_scale_factor = y_scale_factor

        # TODO: finish autoscale implementation
        self.auto_scale = auto_scale

        self.x_y_ratio = None
        self.subplot_x_y_ratio = subplot_x_y_ratio
        self.figure_x_y_ratio = figure_x_y_ratio

        self.highlight = highlight

    def init_coordinates(self, style=None):
        used_style = style if style else self.style

        y = self.y_start + self.style.internal_offset - self.style.distance
        self.track_label_param_list = [[len(self[track_name].label) if (self[track_name].label and self[track_name].style.show_label) else 0,
                                       self[track_name].style.label_fontsize] for track_name in self]
        for track_name in self:
            self[track_name].y_start = y + self.style.distance
            self.x_end = max(self.x_end, self[track_name].x_end)
            y += self[track_name].style.height + self.style.distance
        #y -= self.style.distance
        self.x_end = self.x_end #* self.style.x_multiplier
        self.y_end = y + self.style.internal_offset

        #if self.auto_scale:
        #    self.y_scale_factor = 1
        self.x_y_ratio = self.x_end / self.y_end

        for track_name in self:
            self[track_name].track_group_x_y_ratio = self.x_y_ratio

        if self.highlight:
            for track_name in self:
                self[track_name].track_group_highlight = True
                self[track_name].track_group_highlight_color = used_style.highlight_color

    def draw(self, axes=None, style=None, label_shift=0):
        self.init_coordinates()

        used_style = style if style else self.style
        current_subplot = axes if axes else plt.gca()

        if self.highlight and (used_style.highlight_color is not None):
            highlight_patch = Polygon(np.array([[self.x_start, self.y_start],
                                                [self.x_start, self.y_end],
                                                [self.x_end, self.y_end],
                                                [self.x_end, self.y_start]]),
                                      #color=used_style.empty_color if (
                                      #        used_style.fill_empty and self.records is None) else used_style.face_color,
                                      fill=True,
                                      edgecolor=used_style.highlight_color,
                                      facecolor=used_style.highlight_color,
                                      linewidth=used_style.highlight_edge_width, zorder=used_style.zorder['highlight'])
            current_subplot.add_patch(highlight_patch)

        for track_name in self:
            self[track_name].draw(axes=axes)

        #x_text_offset = max(list(map(lambda s: s[0]*s[1], self.track_label_param_list)))
        # TODO: add automatic calculation for x axis

        label_shift_full = (-120 if max(list(map(lambda s: s[0]*s[1], self.track_label_param_list))) else 0) + label_shift
        if self.label and used_style.show_label:
            #print(self.label)
            #print(self.style.label_x_shift + label_shift_full, (self.y_end + self.y_start) / 2)
            #print(self.y_start, self.y_end)
            #print()
            current_subplot.annotate(self.label, xy=(0, (self.y_start + self.y_end)/2 + self.style.label_y_shift + 2),
                                     xycoords='data',
                                     fontsize=self.style.label_fontsize,
                                     fontstyle=self.style.label_fontstyle,
                                     fontweight=self.style.label_fontweight,
                                     xytext=(self.style.label_x_shift + label_shift_full,
                                             0), textcoords='offset points',
                                     ha=self.style.label_hor_aln, va=self.style.label_vert_aln)

    #def draw_borders(self, axes=None, style=None, label_shift=0):
        #self.init_coordinates()

    #    for track_name in self:
    #        self[track_name].draw_borders(axes=axes)