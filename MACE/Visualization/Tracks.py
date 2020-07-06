
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

    def __init__(self, records, style, y_start=None, x_start=0, x_end=None, label=None, feature_style=default_feature_style,
                 color_expression=None, colormap=None, thresholds=None, colors=None, background=None, masked=None,
                 patch_function=None):

        self.records = records

        self.style = style
        self.feature_style = feature_style

        self.y_start = y_start
        self.y_end = None

        self.x_start = x_start
        self.x_end = x_end

        self.label = label

        if color_expression:
            self.style.color_expression = color_expression
        if thresholds:
            self.style.thresholds = thresholds
        if colors:
            self.style.colors = colors
        if background:
            self.style.background = background
        if masked:
            self.style.masked = masked
        if colormap:
            self.style.colormap = colormap
            if thresholds:
                self.style.cmap = plt.get_cmap(self.style.colormap, len(self.style.thresholds))
                self.style.colors = [self.style.cmap(i) for i in range(0, len(self.style.thresholds))]

        self.patch_function = patch_function

    def color_threshold_expression(self, value):
                                   # colors=("white", "#333a97", "#3d3795", "#5d3393","#813193", "#9d2d7f", "#b82861",
                                   #         "#d33845", "#ea2e2e", "#f5ae27"))
                                   #thresholds=np.array((0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5)),
                                   #colors=("white", "#333a97", "#3d3795", "#5d3393","#813193", "#9d2d7f", "#b82861",
                                  #         "#d33845", "#ea2e2e", "#f5ae27")):
        if value <= self.style.thresholds[0]:
            return self.style.background
        if value > self.style.thresholds[-1]:
            return self.style.colors[-1]
        for i in range(0, len(self.style.thresholds) - 1):
            if self.style.thresholds[i] < value <= self.style.thresholds[i+1]:
                #print(i)
                #print(self.style.colors)
                #print(self.style.thresholds)
                return self.style.colors[i]

    #def set_color(self):
    #    self.records["color"] = self.records["value"].apply()
    def create_patch_function(self, style=None, feature_style=None, y_start=None, *args, **kwargs):
        """
        Track type dependent function
        :param style:
        :param feature_style:
        :param y_start:
        :param args:
        :param kwargs:
        :return:
        """
        raise ValueError("ERROR!!! Call for ancestor abstract class method, i.e this method "
                         "was not implemented in called descendant class")
        pass

    def create_patch_collection(self, y_start=None, style=None, feature_style=None,
                                track_xmin=None, track_xmax=None, patch_function=None, *args, **kwargs):
        # TODO: add min and max coordinate filtering
        if self.style is None and style is None:
            raise ValueError("ERROR!!! No style was set!")
        used_style = style if style else self.style
        used_feature_style = feature_style if feature_style else self.feature_style

        y_track = y_start if y_start else self.y_start

        if patch_function:
            self.patch_function = patch_function
        else:
            if not self.patch_function:
                self.patch_function = self.create_patch_function(style=used_style,
                                                                 feature_style=used_feature_style,
                                                                 y_start=y_track,
                                                                 *args, **kwargs)

        return PatchCollection(self.records.apply(self.patch_function, axis=1), match_original=True,)

    def add_color(self, expression=None, value_column_index=-1, value_column_name=None, masking=True):

        self.records["color"] = map(expression if expression else self.color_threshold_expression,
                                    self.records[value_column_name].to_list() if value_column_name else self.records.iloc[:, value_column_index].to_list())
        self.records["color"].astype('category', copy=False)
        if masking and ("masked" in self.records.columns):
            self.records.loc[self.records["masked"] == True, "color"] = self.style.masked
        print(self.records)

    def add_color_by_dict(self, value_column_name=None, value_column_index=None, default_color='black'):
        if (value_column_name is None) and (value_column_index is None):
            raise ValueError("ERROR!!! Both column name and column index were not set!")
        if value_column_name:
            self.records["color"] = self.records[value_column_name].replace(self.color_dict)
        else:
            self.records["color"] = self.records.iloc[:, value_column_index].replace(self.color_dict)
        self.records["color"][~self.records["color"].isin(list(self.color_dict.keys()))] = default_color

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

        if self.label and self.style.show_label:
            current_subplot.annotate(self.label, xy=(0, self.y_start + self.style.label_y_shift), xycoords='data',
                                     fontsize=self.style.label_fontsize,
                                     xytext=(-15, 1.5 * self.style.label_y_shift), textcoords='offset points',
                                     ha=self.style.label_hor_aln, va=self.style.label_vert_aln)


class WindowTrack(Track):

    def __init__(self, windows_df, window_size, window_step, y_start=None, x_start=0, x_end=None,
                 style=default_track_style, label=None,
                 window_type="stacking", multiplier=1000, feature_style=default_feature_style, color_expression=None,
                 colormap=None, thresholds=None,
                 colors=None, background=None, masked=None):
        """

        :param windows_df:      Example
                                                All
                                CHROM	WINDOW
                                chr4	    0	58
                                            1	61
                                            2	54
        :param window_size:
        :param window_step:
        :param y_start:
        :param x_start:
        :param x_end:
        :param style:
        :param label:
        :param window_type:
        :param multiplier:
        :param feature_style:
        :param color_expression:
        :param colormap:
        :param thresholds:
        :param colors:
        :param background:
        :param masked:
        """

        Track.__init__(self, windows_df, style, y_start=y_start, x_start=x_start, x_end=x_end, label=label,
                       feature_style=feature_style, color_expression=color_expression,
                       colormap=colormap, thresholds=thresholds, colors=colors, background=background, masked=masked)

        self.track_type = "window"
        self.window_type = window_type
        self.records["density"] = self.records[list(self.records.columns).remove("masked")] / window_size * multiplier

        self.window_size = window_size
        self.window_step = window_step
        self.multiplier = multiplier
        self.preprocess_data()

    def preprocess_data(self):
        if self.window_type == "stacking":
            self.records.reset_index(level="WINDOW", inplace=True)
            self.records["start"] = self.records["WINDOW"] * self.window_step
            self.records = self.records[["start"] + list(self.records.columns[:-1])]

        elif self.window_type == "sliding":
            # TODO: implement sliding window drawing
            pass
        else:
            raise ValueError("ERROR!!! Unknown window type: %s" % self.window_type)

    def create_patch_function(self, style=None, feature_style=None, y_start=None, *args, **kwargs):

        if feature_style.patch_type == "rectangle":
            if "color" in self.records.columns:

                def create_patch(row, style=style if style else self.style,
                                 feature_style=feature_style if feature_style else self.feature_style, y_start=y_start,
                                 window_size=self.window_size):

                    return Rectangle((row['start'], y_start), window_size,
                                     style.height,
                                     fill=True,
                                     edgecolor=feature_style.edge_color,
                                     facecolor=row["color"],
                                     linewidth=feature_style.edge_width)

                return create_patch

            else:
                def create_patch(row, style=style, feature_style=feature_style, y_start=y_start,
                                 window_size=self.window_size):

                    return Rectangle((row['start'], y_start), window_size,
                                     style.height,
                                     fill=True,
                                     edgecolor=feature_style.edge_color,
                                     facecolor=feature_style.face_color,
                                     linewidth=feature_style.edge_width)

                return create_patch

    """
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
            if "color" in self.records.columns:
                for window_tuple in self.records.itertuples(index=True, name=None):
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
                for window_tuple in self.records.itertuples(index=True, name=None):
                    window_start = window_tuple[0] * self.window_size

                    patch_list.append(Rectangle((window_start, y_track), self.window_size,
                                                used_style.height,
                                                fill=used_feature_style.fill,
                                                edgecolor=used_feature_style.edge_color,
                                                facecolor=used_feature_style.face_color,
                                                linewidth=used_feature_style.edge_width))

        #return patch_list
        return PatchCollection(patch_list, match_original=True)
    """


class FeatureTrack(Track):

    def __init__(self, feature_df, y_start=None, x_start=0, x_end=None,
                 style=default_track_style, label=None,
                 feature_style=default_feature_style, color_expression=None,
                 colormap=None, thresholds=None,
                 colors=None, background=None, masked=None):
        Track.__init__(self, feature_df, style, y_start=y_start, x_start=x_start, x_end=x_end, label=label,
                       feature_style=feature_style, color_expression=color_expression,
                       colormap=colormap, thresholds=thresholds, colors=colors, background=background, masked=masked)

        self.track_type = "feature"
        self.preprocess_data()

    def preprocess_data(self):
        self.records["length"] = self.records["end"] - self.records["start"]

    def create_patch_function(self, style=None, feature_style=None, y_start=None, *args, **kwargs):

        if feature_style.patch_type == "rectangle":
            if "color" in self.records.columns:

                def create_patch(row, style=style if style else self.style,
                                 feature_style=feature_style if feature_style else self.feature_style, y_start=y_start):

                    return Rectangle((row['start'], y_start), row['length'],
                                     style.height,
                                     fill=True,
                                     edgecolor=feature_style.edge_color,
                                     facecolor=row["color"],
                                     linewidth=feature_style.edge_width)

                return create_patch

            else:
                def create_patch(row, style=style, feature_style=feature_style, y_start=y_start):

                    return Rectangle((row['start'], y_start), row['length'],
                                     style.height,
                                     fill=True,
                                     edgecolor=feature_style.edge_color,
                                     facecolor=feature_style.face_color,
                                     linewidth=feature_style.edge_width)

                return create_patch