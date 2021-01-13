
from MACE.Visualization.Styles.Track import default_track_style
from MACE.Visualization.Styles.Feature import default_feature_style
import math
from collections import Iterable, OrderedDict

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle, Ellipse


class Track:

    def __init__(self, records, style, y_start=None, x_start=0, x_end=None, label=None,
                 feature_style=default_feature_style, color_expression=None, colormap=None, thresholds=None,
                 colors=None, background=None, masked=None, patch_function=None, subplot_scale=False,
                 track_group_scale=False, figure_x_y_ration=None, subplot_x_y_ratio=None, track_group_x_y_ratio=None):

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

        self.figure_x_y_ratio = figure_x_y_ration
        self.subplot_x_y_ratio = subplot_x_y_ratio
        self.track_group_x_y_ratio = track_group_x_y_ratio

        self.subplot_scale = subplot_scale
        self.track_group_scale = track_group_scale

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
    def create_patch_function(self, style=None, feature_style=None, y_start=None,  *args, **kwargs):
        """
        Track type dependent function
        :param track_group_scale:
        :param subplot_scale:
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
        if self.records is None:
            return 0
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

        self.records["color"] = list(map(expression if expression else self.color_threshold_expression,
                                    self.records[value_column_name].to_list() if value_column_name else self.records.iloc[:, value_column_index].to_list()))
        self.records["color"].astype('category', copy=False)
        if masking and ("masked" in self.records.columns):
            self.records.loc[self.records["masked"] == True, "color"] = self.style.masked

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
            #print (self.x_start, self.y_start, self.x_end)
            current_subplot.add_patch(Rectangle((self.x_start, self.y_start), self.x_end,
                                      used_style.height,
                                      color=used_style.empty_color if (used_style.fill_empty and self.records is None) else used_style.face_color,
                                      fill=True if (used_style.fill_empty and self.records is None) else used_style.fill,
                                      edgecolor=used_style.edge_color,
                                      facecolor=used_style.face_color,
                                      linewidth=used_style.edge_width))
        if self.records is not None:
            current_subplot.add_collection(self.create_patch_collection())

        if self.label and self.style.show_label:
            current_subplot.annotate(self.label, xy=(0, self.y_start + self.style.height/2.5), xycoords='data',
                                     fontsize=self.style.label_fontsize,
                                     xytext=(-15, 0), textcoords='offset points',
                                     ha=self.style.label_hor_aln, va=self.style.label_vert_aln)


class WindowTrack(Track):

    def __init__(self, windows_df, window_size, window_step, y_start=None, x_start=0, x_end=None,
                 style=default_track_style, label=None, norm=False,
                 window_type="stacking", multiplier=None, feature_style=default_feature_style, color_expression=None,
                 colormap=None, thresholds=None,
                 colors=None, background=None, masked=None, subplot_scale=False,
                 track_group_scale=False, figure_x_y_ration=None, subplot_x_y_ratio=None, track_group_x_y_ratio=None):
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
                       colormap=colormap, thresholds=thresholds, colors=colors, background=background, masked=masked,
                       subplot_scale=subplot_scale, track_group_scale=track_group_scale,
                       figure_x_y_ration=figure_x_y_ration, subplot_x_y_ratio=subplot_x_y_ratio,
                       track_group_x_y_ratio=track_group_x_y_ratio
                       )

        self.track_type = "window"
        self.window_type = window_type
        #print(self.records)
        #print(self.records.columns.to_list())
        #print(self.records.columns.to_list().remove("masked"))

        count_columns = self.records.columns.to_list()
        if "masked" in count_columns:
            count_columns.remove("masked")

        if norm and multiplier:
            self.records["density"] = self.records.loc[:, count_columns] * (float(multiplier) / float(window_size))
        elif multiplier:
            self.records["density"] = self.records.loc[:, count_columns] * float(multiplier)
        elif norm:
            self.records["density"] = self.records.loc[:, count_columns] / float(window_size)
        else:
            self.records["density"] = self.records.loc[:, count_columns]

        self.window_size = window_size
        self.window_step = window_step
        self.multiplier = multiplier
        self.preprocess_data()

    def preprocess_data(self):
        window_index_name = self.records.index.names[0]
        if self.window_type == "stacking":
            self.records.reset_index(level=window_index_name, inplace=True)
            self.records["start"] = self.records[window_index_name] * self.window_step
            self.records = self.records[["start"] + list(self.records.columns[:-1])]

        elif self.window_type == "sliding":
            # TODO: implement sliding window
            # TODO: verify - sliding window drawing is already implemented  somewhere else:)))
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
                 colors=None, background=None, masked=None,
                 feature_start_column_id="start", 
                 feature_end_column_id="end", 
                 feature_color_column_id="color",
                 feature_length_column_id="length",
                 x_scale_factor=1,
                 y_scale_factor=1,
                 auto_scale=False,
                 subplot_scale=False,
                 track_group_scale=False,
                 figure_x_y_ration=None, subplot_x_y_ratio=None, track_group_x_y_ratio=None
                 ):

        Track.__init__(self, feature_df, style, y_start=y_start, x_start=x_start, x_end=x_end, label=label,
                       feature_style=feature_style, color_expression=color_expression,
                       colormap=colormap, thresholds=thresholds, colors=colors, background=background, masked=masked,
                       subplot_scale=subplot_scale,
                       track_group_scale=track_group_scale,
                       figure_x_y_ration=figure_x_y_ration, subplot_x_y_ratio=subplot_x_y_ratio,
                       track_group_x_y_ratio=track_group_x_y_ratio
                       )

        self.track_type = "feature"

        self.feature_start_column_id = feature_start_column_id
        self.feature_end_column_id = feature_end_column_id
        self.feature_color_column_id = feature_color_column_id
        self.feature_length_column_id = feature_length_column_id

        self.x_scale_factor = x_scale_factor
        self.y_scale_factor = y_scale_factor
        self.auto_scale = auto_scale
        self.preprocess_data()

    def preprocess_data(self):
        if self.records is not None:
            self.records[self.feature_length_column_id] = self.records[self.feature_end_column_id] - self.records[self.feature_start_column_id]

    def create_patch_function(self, style=None, feature_style=None, y_start=None, *args, **kwargs):
        if self.records is None:
            return 0

        if feature_style.patch_type == "rectangle":
            if self.feature_color_column_id in self.records.columns:

                def create_patch(row, style=style if style else self.style,
                                 feature_style=feature_style if feature_style else self.feature_style, y_start=y_start):

                    return Rectangle((row[self.feature_start_column_id], y_start), row[self.feature_length_column_id],
                                     style.height,
                                     fill=True,
                                     edgecolor=feature_style.edge_color,
                                     facecolor=row[self.feature_color_column_id],
                                     linewidth=feature_style.edge_width)

                return create_patch

            else:
                def create_patch(row, style=style, feature_style=feature_style, y_start=y_start):

                    return Rectangle((row[self.feature_start_column_id], y_start), row[self.feature_length_column_id],
                                     style.height,
                                     fill=True,
                                     edgecolor=feature_style.edge_color,
                                     facecolor=feature_style.face_color,
                                     linewidth=feature_style.edge_width)

                return create_patch

        elif feature_style.patch_type == "circle":

            if self.feature_color_column_id in self.records.columns:

                def create_patch(row, style=style if style else self.style,
                                 feature_style=feature_style if feature_style else self.feature_style, y_start=y_start):

                    return Circle(((row[self.feature_end_column_id] + row[self.feature_start_column_id]) / 2 , y_start + style.height / 2), radius=feature_style.height / 2,
                                  fill=True,
                                  edgecolor=feature_style.edge_color,
                                  facecolor=row[self.feature_color_column_id],
                                  linewidth=feature_style.edge_width)

                return create_patch

            else:
                def create_patch(row, style=style, feature_style=feature_style, y_start=y_start):
                    return Circle(((row[self.feature_end_column_id] + row[self.feature_start_column_id]) / 2, y_start + style.height / 2),
                                  radius=feature_style.height / 2,
                                  fill=True,
                                  edgecolor=feature_style.edge_color,
                                  facecolor=feature_style.face_color,
                                  linewidth=feature_style.edge_width)
                
                return create_patch

        elif feature_style.patch_type == "ellipse":

            if self.subplot_scale:
                x_scale_factor = self.subplot_x_y_ratio / self.figure_x_y_ratio
            elif self.track_group_scale:
                x_scale_factor = self.track_group_x_y_ratio / self.figure_x_y_ratio
            else:
                x_scale_factor = 1
            y_scale_factor = 1

            if self.feature_color_column_id in self.records.columns:

                def create_patch(row, style=style if style else self.style,
                                 feature_style=feature_style if feature_style else self.feature_style, y_start=y_start):

                    return Ellipse(((row[self.feature_end_column_id] + row[self.feature_start_column_id]) / 2,
                                   y_start + style.height / 2),
                                   height=feature_style.height * y_scale_factor,
                                   width=feature_style.height * x_scale_factor if self.subplot_scale else (row[self.feature_end_column_id] - row[self.feature_start_column_id] + 1) * x_scale_factor,
                                   fill=True,
                                   edgecolor=feature_style.edge_color,
                                   facecolor=row[self.feature_color_column_id],
                                   linewidth=feature_style.edge_width)

                return create_patch

            else:
                def create_patch(row, style=style, feature_style=feature_style, y_start=y_start):
                    return Ellipse(((row[self.feature_end_column_id] + row[self.feature_start_column_id]) / 2,
                                   y_start + style.height / 2),
                                   height=feature_style.height * y_scale_factor,
                                   width=feature_style.height * x_scale_factor if self.subplot_scale else (row[self.feature_end_column_id] - row[self.feature_start_column_id] + 1) * x_scale_factor,
                                   fill=True,
                                   edgecolor=feature_style.edge_color,
                                   facecolor=feature_style.face_color,
                                   linewidth=feature_style.edge_width)

                return create_patch