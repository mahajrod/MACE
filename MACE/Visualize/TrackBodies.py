
import sys
import pandas as pd
import logging
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Iterable

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle, Ellipse, Polygon
from matplotlib.lines import Line2D

from MACE.Visualize.Polygons import LinearChromosome
from MACE.Visualize.Styles.RecordStyles import default_record_style
from MACE.Visualize.Styles.TrackBodyStyles import default_track_body_style


class TrackBody:  # TODO: add breaks
    def __init__(self, record_dict, record_style_dict, style=default_track_body_style,
                 y_start=None, x_start=0, x_end=None, y_end=None, subplot_scale=False,
                 track_group_scale=False, figure_x_y_ratio=None, subplot_x_y_ratio=None, track_group_x_y_ratio=None,
                 stranded=None, rounded=None, centromere_start=None, centromere_end=None,
                 record_x_start_column_id="start", record_x_end_column_id="end", record_strand_column_id="strand",
                 default_record_style=default_record_style):

        self.record_dict = record_dict
        self.style = style
        self.record_style_dict = record_style_dict
        self.default_record_style = default_record_style
        self.y_start = y_start
        self.y_end = y_end

        self.x_start = x_start
        self.x_end = x_end

        self.figure_x_y_ratio = figure_x_y_ratio
        self.subplot_x_y_ratio = subplot_x_y_ratio
        self.track_group_x_y_ratio = track_group_x_y_ratio

        #----------------- Track Body options-------------------
        self.subplot_scale = subplot_scale
        self.track_group_scale = track_group_scale
        self.stranded = stranded
        self.rounded = rounded
        #self.middle_break = middle_break

        self.centromere_start = centromere_start
        self.centromere_end = centromere_end

        self.x_scale_factor = None
        self.recognized_track_types = ["window", "hist", "plot", "marker"]
        #----------------------------------------------

        #----------------- Record options -----------------

        self.record_x_start_column_id = record_x_start_column_id
        self.record_x_end_column_id = record_x_end_column_id
        self.record_strand_column_id = record_strand_column_id
        #--------------------------------------------------

        self.polygon_dict = {}

        if self.record_dict is None:
            self.empty = True
        elif not self.record_dict:
            self.empty = True
        elif sum(map(lambda r: not r.empty, self.record_dict.values())) == 0:
            self.empty = True
        else:
            self.empty = False

    @staticmethod
    def check_na(value):
        na_check = pd.isnull(value)
        return sum(na_check) if isinstance(na_check, Iterable) else na_check

    def create_patch_function_dict(self, y_start, y_end, style=None, record_style_dict=None, stranded=None):
        if self.empty:
            return None
        used_style = style if style is None else self.style
        patch_function_dict = OrderedDict()

        for track_name in self.record_dict:
            track_type = self.record_dict[track_name].track_type_list[0]
            if track_type not in self.recognized_track_types:
                logging.error(f"ERROR!!! Unrecognized track type - '{track_type}'. Allowed track types: {','.join(self.recognized_track_types)}")
                raise ValueError(f"ERROR!!! Unrecognized track type - '{track_type}'. Allowed track types: {','.join(self.recognized_track_types)}")

            print(used_style)
            #----------- detect style of the track ------------
            if track_name in record_style_dict:
                used_record_style = record_style_dict[track_name]
            elif track_name in self.record_style_dict:
                used_record_style = self.record_style_dict[track_name]
            elif self.default_record_style is not None:
                used_record_style = self.default_record_style
            else:
                logging.error(f"ERROR!!! Record style was not set for track {track_name}")
                ValueError(f"ERROR!!! Record style was not set for track {track_name}")
            #-------------------------------------------------

            #----------- detect strandness of the track -----------
            if stranded is not None:
                strandedness = stranded
            elif self.stranded is not None:
                strandedness = self.stranded
            else:
                strandedness = None
            #------------------------------------------------------

            zorder_shift = 0
            if (self.record_dict[track_name].track_feature_dict is not None) and self.record_dict[track_name].track_feature_dict:
                if track_name in self.record_dict[track_name].track_feature_dict:
                    if "zorder" in self.record_dict[track_name].track_feature_dict[track_name]:
                        zorder_shift = self.record_dict[track_name].track_feature_dict[track_name]["zorder"]

            if track_type == "window":
                patch_function_dict[track_name] = self.create_patch_function_window(used_style,
                                                                                    used_record_style[track_type],
                                                                                    y_start, y_end,
                                                                                    stranded=strandedness,
                                                                                    zorder_shift=zorder_shift,
                                                                                    record_color_column_id=self.record_dict[track_name].track_metadata_df[self.record_dict[track_name].track_metadata_df == "type"].loc[track_name, "column_name"])
            elif track_type == "hist":
                patch_function_dict[track_name] = self.create_patch_function_hist(used_style,
                                                                                  used_record_style[track_type],
                                                                                  y_start, y_end,
                                                                                  stranded=strandedness,
                                                                                  zorder_shift=zorder_shift)
            elif track_type == "plot":
                patch_function_dict[track_name] = self.create_patch_function_plot(used_style,
                                                                                  used_record_style[track_type],
                                                                                  y_start, y_end,
                                                                                  stranded=strandedness,
                                                                                  zorder_shift=zorder_shift)
            elif track_type == "marker":
                patch_function_dict[track_name] = self.create_patch_function_marker(used_style,
                                                                                    used_record_style[track_type],
                                                                                    y_start, y_end,
                                                                                    stranded=strandedness,
                                                                                    zorder_shift=zorder_shift)

        return patch_function_dict

    def create_patch_function_window(self, style, patch_style, y_start, y_end, stranded=None,
                                     record_color_column_id=None, record_edge_color_column_id=None,
                                     record_edge_width_column_id=None, zorder_shift=0):
        print(type(style))
        if stranded:
            def forward_patch_function(row):
                print(row)

                if self.check_na(row[record_color_column_id]):
                    return pd.NA

                return Rectangle((row[self.record_x_start_column_id], (y_start + y_end)/2),
                                 row[self.record_x_start_column_id] + row[self.record_x_end_column_id],
                                 (y_end - y_start)/2,
                                 fill=True if record_color_column_id and (row[record_color_column_id]) != "default" else True if patch_style.default_color else False,
                                 edgecolor=row[record_edge_color_column_id] if record_edge_color_column_id and (row[record_edge_color_column_id]) != "default" else patch_style.default_edge_color,
                                 facecolor=row[record_color_column_id] if record_color_column_id and (row[record_color_column_id]) != "default" else patch_style.default_color,
                                 linewidth=row[record_edge_width_column_id] if record_edge_width_column_id and (row[record_edge_width_column_id]) != "default" else patch_style.default_edge_width,
                                 zorder=style.zorder_dict["records"] + patch_style.zorder_shift + zorder_shift)

            def reverse_patch_function(row):
                print(row)
                if self.check_na(row[record_color_column_id]):
                    return pd.NA

                return Rectangle((row[self.record_x_start_column_id], y_start),
                                 row[self.record_x_start_column_id] + row[self.record_x_end_column_id],
                                 (y_end - y_start)/2,
                                 fill=True if record_color_column_id and (row[record_color_column_id]) != "default" else True if patch_style.default_color else False,
                                 edgecolor=row[record_edge_color_column_id] if record_edge_color_column_id and (row[record_edge_color_column_id]) != "default" else patch_style.default_edge_color,
                                 facecolor=row[record_color_column_id] if record_color_column_id and (row[record_color_column_id]) != "default" else patch_style.default_color,
                                 linewidth=row[record_edge_width_column_id] if record_edge_width_column_id and (row[record_edge_width_column_id]) != "default" else patch_style.default_edge_width,
                                 zorder=style.zorder_dict["records"] + patch_style.zorder_shift + zorder_shift)

            return forward_patch_function, reverse_patch_function

        else:
            def patch_function(row):
                print(row)
                #print(Rectangle((row[self.record_x_start_column_id], y_start),
                #                 row[self.record_x_start_column_id] + row[self.record_x_end_column_id],
                #                 y_end - y_start,
                #                 fill=True if record_color_column_id and (row[record_color_column_id]) != "default" else True if patch_style.default_color else False,
                #                 edgecolor=row[record_edge_color_column_id] if record_edge_color_column_id and (row[record_edge_color_column_id]) != "default" else patch_style.default_edge_color,
                #                 facecolor=row[record_color_column_id] if record_color_column_id and (row[record_color_column_id]) != "default" else patch_style.default_color,
                #                 linewidth=row[record_edge_width_column_id] if record_edge_width_column_id and (row[record_edge_width_column_id]) != "default" else patch_style.default_edge_width,
                #                 zorder=style.zorder_dict["records"] + patch_style.zorder_shift))
                if self.check_na(row[record_color_column_id]):
                    return pd.NA
                return Rectangle((row[self.record_x_start_column_id], y_start),
                                 row[self.record_x_start_column_id] + row[self.record_x_end_column_id],
                                 y_end - y_start,
                                 fill=True if record_color_column_id and (row[record_color_column_id]) != "default" else True if patch_style.default_color else False,
                                 edgecolor=row[record_edge_color_column_id] if record_edge_color_column_id and (row[record_edge_color_column_id]) != "default" else patch_style.default_edge_color,
                                 facecolor=row[record_color_column_id] if record_color_column_id and (row[record_color_column_id]) != "default" else patch_style.default_color,
                                 linewidth=row[record_edge_width_column_id] if record_edge_width_column_id and (row[record_edge_width_column_id]) != "default" else patch_style.default_edge_width,
                                 zorder=style.zorder_dict["records"] + patch_style.zorder_shift)
            return patch_function, None

    def create_patch_function_hist(self, style, record_style, y_start, y_end, stranded=None, zorder_shift=0):
        pass

    def create_patch_function_plot(self, style, record_style, y_start, y_end, stranded=None, zorder_shift=0):
        pass

    def create_patch_function_marker(self, style, record_style, y_start, y_end, stranded=None, zorder_shift=0):
        pass

    def create_patch_collection(self, y_start=None, y_end=None, style=None, record_style_dict=None, stranded=None):
        if self.empty:
            return None
        #if self.style is None and style is None:
        #    raise ValueError("ERROR!!! No style was set!")
        used_style = style if style is not None else self.style
        used_record_style_dict = record_style_dict if record_style_dict else self.record_style_dict

        start = y_start if y_start else self.y_start
        end = y_end if y_end else self.y_end

        patch_function_dict = self.create_patch_function_dict(start, end, style=used_style,
                                                              record_style_dict=used_record_style_dict,
                                                              stranded=stranded)
        if patch_function_dict is None:
            return []

        patch_collection_list = []

        for track_name in patch_function_dict:
            if patch_function_dict[track_name][1] is None:     # consider records as unstranded:
                tmp = self.record_dict[track_name].apply(patch_function_dict[track_name][0], axis=1)
                patch_collection_list.append(PatchCollection(tmp[~pd.isnull(tmp)],
                                                             match_original=True,
                                                             antialiased=False))
            else:
                if self.record_strand_column_id in self.record_dict[track_name].columns:
                    if not self.record_dict[track_name][self.record_dict[track_name][self.record_strand_column_id] == "+"].empty:
                        tmp = self.record_dict[track_name][self.record_dict[track_name][self.record_strand_column_id] == "+"].apply(patch_function_dict[track_name][0], axis=1)
                        patch_collection_list.append(PatchCollection(tmp[~pd.isnull(tmp)],
                                                                     match_original=True,
                                                                     antialiased=False))

                        tmp = self.record_dict[track_name][self.record_dict[track_name][self.record_strand_column_id] == "-"].apply(patch_function_dict[track_name][0], axis=1)

                        patch_collection_list.append(PatchCollection(tmp[~pd.isnull(tmp)],
                                                                     match_original=True,
                                                                     antialiased=False))
                    else:
                        tmp = self.record_dict[track_name].apply(patch_function_dict[track_name][0], axis=1)
                        patch_collection_list.append(PatchCollection(tmp[~pd.isnull(tmp)],
                                                                     match_original=True,
                                                                     antialiased=False))
                else:
                    tmp = self.record_dict[track_name].apply(patch_function_dict[track_name][0], axis=1)
                    patch_collection_list.append(PatchCollection(tmp[~pd.isnull(tmp)],
                                                                 match_original=True,
                                                                 antialiased=False))
        return patch_collection_list

    def draw(self, axes=None, style=None, y_start=None, y_end=None, record_style_dict=None, stranded=None):

        current_subplot = axes if axes else plt.gca()

        if not self.polygon_dict:
            self.calculate_service_polygons(style=style)

        # add service polygons
        for patch in self.polygon_dict.values():
            current_subplot.add_patch(patch)

        # add features\windows\etc first
        end_y = None
        if y_end is not None:
            end_y = y_end
        elif style is not None:
            end_y = self.y_start + style.height
        elif self.style is not None:
            end_y = self.y_start + self.style.height
        else:
            logging.error(f"ERROR!!! Track body style was not set!")
            ValueError(f"ERROR!!! Track body style was not set!")

        if (self.record_dict is not None) and (not self.empty):
            print("BBBBB")
            for collection in self.create_patch_collection(y_start=y_start if y_start is not None else self.y_start,
                                                           y_end=end_y,
                                                           record_style_dict=record_style_dict,
                                                           stranded=stranded):
                if collection:
                    current_subplot.add_collection(collection)

    def calculate_service_polygons(self, style=None):

        used_style = style if style else self.style
        if used_style is None:
            raise ValueError("ERROR!!! No track body style was set")

        if not self.empty:
            background_color = used_style.background_color
            background_alpha = used_style.background_alpha
        else:
            if used_style.highlight_empty:
                background_color = used_style.empty_color if used_style.empty_color else used_style.background_color
                background_alpha = used_style.empty_alpha if used_style.empty_alpha else used_style.background_alpha
            else:
                background_color = used_style.background_color
                background_alpha = used_style.background_alpha

        self.x_scale_factor = self.subplot_x_y_ratio / (self.figure_x_y_ratio if self.figure_x_y_ratio is not None else 1)  # same scalling for all tracks

        left_middle_point = np.array([self.x_start, self.y_start + used_style.height / 2])
        right_middle_point = np.array([self.x_end, self.y_start + used_style.height / 2])

        self.polygon_dict["outer"] = Rectangle((self.x_start, self.y_start), self.x_end - self.x_start, used_style.height,
                                               edgecolor=used_style.outer_edge_color if used_style.outer_edge_color is not None else used_style.outer_color,
                                               fill=True if used_style.outer_color else False,
                                               facecolor=used_style.outer_color,
                                               alpha=used_style.outer_alpha,
                                               linewidth=used_style.outer_edge_width,
                                               zorder=used_style.zorder_dict["outer"])

        self.polygon_dict["background"] = LinearChromosome(self.x_start, self.y_start, self.x_end - self.x_start, used_style.height,
                                                           stranded=used_style.stranded, rounded=used_style.rounded,
                                                           stranded_end=used_style.stranded_end,
                                                           centromere_start=self.centromere_start, centromere_end=self.centromere_end,
                                                           show_centromere=used_style.show_centromere,
                                                           arc_point_number=used_style.arc_point_number,
                                                           x_scale_factor=self.x_scale_factor,
                                                           edgecolor=background_color,
                                                           facecolor=background_color,
                                                           alpha=background_alpha,
                                                           fill=True if used_style.background_color else False,
                                                           zorder=used_style.zorder_dict["background"],
                                                           outer_zorder=used_style.zorder_dict["masking_patches"],
                                                           outercolor=used_style.outer_color,
                                                           linewidth=0)
        self.polygon_dict = {
                             **self.polygon_dict,
                             **self.polygon_dict["background"].masking_point_array_dict
                             }

        if used_style.show_middle_line:
            self.polygon_dict["middle_line"] = Line2D((left_middle_point[0], right_middle_point[0]),
                                                      (left_middle_point[1], right_middle_point[1]),
                                                      color=used_style.middle_line_color,
                                                      linewidth=used_style.middle_line_width,
                                                      zorder=used_style.zorder["middle_line"])

        if (not self.empty) or used_style.empty_edge:
            self.polygon_dict["border"] = Polygon(self.polygon_dict["background"].point_array,
                                                  edgecolor=used_style.edge_color, linewidth=used_style.edge_width, facecolor=None,
                                                  fill=False, alpha=used_style.edge_alpha, zorder=used_style.zorder_dict["border"])

