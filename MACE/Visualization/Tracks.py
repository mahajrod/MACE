
from MACE.Visualization.Styles.Track import default_track_style
from MACE.Visualization.Styles.Feature import default_feature_style
import math
from collections import Iterable, OrderedDict

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


class Track:

    def __init__(self, style):
        self.style = style


class WindowTrack(Track):

    def __init__(self, windows_df, window_size, window_step, y_start=None, x_start=0, x_end=None,
                 style=default_track_style, label=None,
                 window_type="stacking", feature_style=default_feature_style, color_expression=None,
                 multiplier=1000, colormap=None, thresholds=np.array((0.0, 0.1, 0.25, 0.5, 1.0)),
                 colors=("#333a97", "green", "yellow", "orange", "red"), background="white",
                 masked="grey"):
        Track.__init__(self, style)

        self.track_type = "window"
        self.windows = windows_df
        self.windows["density"] = self.windows / window_size * multiplier
        self.window_type = window_type
        self.window_size = window_size
        self.window_step = window_step
        self.feature_style = feature_style

        self.y_start = y_start
        self.y_end = None

        self.x_start = x_start
        self.x_end = x_end

        self.label = label

        self.color_expression = color_expression
        self.colormap = colormap
        self.thresholds = thresholds
        self.colors = colors

        if colormap:
            self.cmap = plt.get_cmap(self.colormap, len(self.thresholds))
            self.colors = [self.cmap(i) for i in range(0, len(thresholds))]

        self.background = background
        self.masked = masked

    def set_color(self):
        self.windows["color"] = self.windows["value"].apply()

    def color_threshold_expression(self, value):
                                   # colors=("white", "#333a97", "#3d3795", "#5d3393","#813193", "#9d2d7f", "#b82861",
                                   #         "#d33845", "#ea2e2e", "#f5ae27"))
                                   #thresholds=np.array((0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5)),
                                   #colors=("white", "#333a97", "#3d3795", "#5d3393","#813193", "#9d2d7f", "#b82861",
                                  #         "#d33845", "#ea2e2e", "#f5ae27")):
        if value <= self.thresholds[0]:
            return self.background
        if value > self.thresholds[-1]:
            return self.colors[-1]
        for i in range(0, len(self.thresholds) - 1):
            if self.thresholds[i] < value <= self.thresholds[i+1]:
                return self.colors[i]

    def add_color(self, expression=None):

        self.windows["color"] = map(expression if expression else self.color_threshold_expression,
                                    self.windows.iloc[:, -1].to_list())

    def create_patch_collection(self, y_start=None, style=None, feature_style=None,
                                track_xmin=None, track_xmax=None):
        # TODO: add min and max coordinate filtering
        if self.style is None and style is None:
            raise ValueError("ERROR!!! No style was set!")

        used_style = style if style else self.style
        used_feature_style = feature_style if feature_style else self.feature_style

        y_track = y_start if y_start else self.y_start
        #color_exp = color_expression if color_expression else self.color_expression

        patch_list = []
        if used_feature_style.patch_type == "rectangle":
            if "color" in self.windows.columns:
                for window_tuple in self.windows.itertuples(index=True, name=None):
                    window_start = window_tuple[0] * self.window_size
                    #print window_tuple[-1]
                    #print window_start, y_track, self.window_size, used_style.height

                    patch_list.append(Rectangle((window_start, y_track), self.window_size,
                                                used_style.height,
                                                fill=True,
                                                edgecolor=used_feature_style.edge_color,
                                                facecolor=window_tuple[-1],
                                                linewidth=used_feature_style.edge_width))

            else:
                for window_tuple in self.windows.itertuples(index=True, name=None):
                    window_start = window_tuple[0] * self.window_size

                    patch_list.append(Rectangle((window_start, y_track), self.window_size,
                                                used_style.height,
                                                fill=used_feature_style.fill,
                                                edgecolor=used_feature_style.edge_color,
                                                facecolor=used_feature_style.face_color,
                                                linewidth=used_feature_style.edge_width))

        #return patch_list
        return PatchCollection(patch_list, match_original=True)

    def draw(self, axes=None, style=None):

        used_style = style if style else self.style

        current_subplot = axes if axes else plt.gca()

        if used_style.edge:
            current_subplot.add_patch(Rectangle((self.x_start, self.y_start), self.x_end,
                                      used_style.height,
                                      fill=used_style.fill,
                                      edgecolor=used_style.edge_color,
                                      facecolor=used_style.face_color,
                                      linewidth=used_style.edge_width))

        current_subplot.add_collection(self.create_patch_collection())
        #for patch in self.create_patch_collection():
        #    current_subplot.add_patch(patch)

        if self.label and self.style.show_label:
            current_subplot.annotate(self.label, xy=(0, self.y_start + self.style.label_y_shift), xycoords='data',
                                     fontsize=self.style.label_fontsize,
                                     xytext=(-15, 1.5 * self.style.label_y_shift), textcoords='offset points',
                                     ha=self.style.label_hor_aln, va=self.style.label_vert_aln)

        #




