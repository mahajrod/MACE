#!/usr/bin/env python
import math
import datetime

from copy import deepcopy
from collections.abc import Iterable
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from RouToolPa.Collections.General import TwoLvlDict
from RouToolPa.Routines.Drawing import DrawingRoutines

from MACE.Data.DataTypes import *

from MACE.Visualize.Tracks import *
from MACE.Visualize.VerticalTrackGroups import *
from MACE.Visualize.HorizontalTrackGroups import *
from MACE.Visualize.Subplots import *
from MACE.Visualize.Figures import *

from MACE.Visualize.Legends import *
from MACE.Functions.Generators import recursive_generator
from MACE.Visualize.Styles.Subplot import *
from MACE.Visualize.Styles.Figure import *
from MACE.Visualize.Styles.Feature import *
from MACE.Visualize.Styles.Track import *


class Visualization(DrawingRoutines):

    def __init__(self):
        DrawingRoutines.__init__(self)
        self.colormap_list = plt.colormaps()
        """
        self.colormap_dict = OrderedDict({'sequential_uniform': ['viridis', 'plasma', 'inferno', 'magma', 'cividis'],

                                          'sequential': ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                                                         'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                                                         'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'],
                                          'sequential_2': ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
                                                          'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                                                          'hot', 'afmhot', 'gist_heat', 'copper'],
                                          'diverging': ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                                                        'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'],
                                          'cyclic': ['twilight', 'twilight_shifted', 'hsv'],
                                          'qualitative': ['Pastel1', 'Pastel2', 'Paired', 'Accent',
                                                          'Dark2', 'Set1', 'Set2', 'Set3',
                                                          'tab10', 'tab20', 'tab20b', 'tab20c'],
                                          'miscellaneous': ['flag', 'prism', 'ocean', 'gist_earth', 'terrain',
                                                            'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix',
                                                            'brg', 'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral',
                                                            'gist_ncar']})
        """

    @staticmethod
    def color_threshold_expression(value, thresholds, colors, background):
        # TODO: Needs at least partial implementation as ColorStyle
        """
        :param value:
        :param thresholds:
        :param background:
        :return:
        """
        # colors=("white", "#333a97", "#3d3795", "#5d3393","#813193", "#9d2d7f", "#b82861",
        #         "#d33845", "#ea2e2e", "#f5ae27"))
        # thresholds=np.array((0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5)),
        # colors=("white", "#333a97", "#3d3795", "#5d3393","#813193", "#9d2d7f", "#b82861",
        #         "#d33845", "#ea2e2e", "#f5ae27")):
        if value <= thresholds[0]:
            return background
        if value > thresholds[-1]:
            return colors[-1]
        for i in range(0, len(thresholds) - 1):
            if thresholds[i] < value <= thresholds[i + 1]:
                # print(i)
                # print(self.style.colors)
                # print(self.style.thresholds)
                return colors[i]

    @staticmethod
    def add_color_to_track_df(track_df, expression, value_column_index=-1, value_column_name=None,
                                masking=False, masked_color="grey"):
        output_df = track_df.copy()

        output_df["color"] = list(map(expression, output_df[value_column_name].to_list() if value_column_name else output_df.iloc[:, value_column_index].to_list()))
        output_df["color"].astype('category', copy=False)
        if masking and ("masked" in output_df):
            output_df.loc[output_df["masked"] == True, "color"] = masked_color

        return output_df

    @staticmethod
    def apply_masking_to_track_df(track_df, masking_df, masking_color=None):
        #print(track_df)
        output_df_columns = list(track_df.columns)
        output_df = track_df.copy().reset_index(drop=False).set_index(["scaffold", "start", "end"])
        tmp_masking_df = masking_df.reset_index(drop=False).set_index(["scaffold", "start", "end"])
        if "color" not in output_df.columns:
            raise ValueError("ERROR!!! 'color' column is absent in track dataframe. Mask must be added after coloring dataframe!")

        #check index overlap
        overlap_sr = tmp_masking_df.index.isin(output_df.index)
        if sum(overlap_sr) != len(overlap_sr):
            print(overlap_sr)
            print(overlap_sr[overlap_sr])
            raise ValueError("ERROR!!! Not all windows from masking track are present in track df!")

        # apply masking
        output_df.loc[tmp_masking_df.index, "color"] = tmp_masking_df["color"] if masking_color is None else masking_color
        output_df["color"].astype('category', copy=False)
        output_df = output_df.reset_index(drop=False).set_index("scaffold")
        output_df = output_df[output_df_columns]
        #print(output_df)

        return output_df

    @staticmethod
    def density_legend(colors, thresholds, colormap=None, feature_name="SNPs"):
        return DensityLegend(colors=colors, colormap=colormap, thresholds=thresholds, feature_name=feature_name)

    @staticmethod
    def coverage_legend(colormap, thresholds):
        return CoverageLegend(colormap=colormap, thresholds=thresholds)

    @staticmethod
    def chromosome_legend(species_color_df_dict, reference_scaffold_order_list):
        return ChromosomeLegend(chromosome_df_dict=species_color_df_dict,
                                scaffold_order_list=reference_scaffold_order_list)

    @staticmethod
    def feature_legend(legend_df, colormap):
        return FeatureLegend(legend_df, colormap=colormap, ) if legend_df is not None else None

    def draw_tracks(self, track_dict, track_length_df_dict, track_order_dict,
                    output_prefix=None,



                    centromere_df=None,
                    highlight_df_dict=None,
                    legend=None,


                    #colormap=None, thresholds=None, colors=None, background=None,
                    default_color="red", # TODO: check if it is possible to remove it
                    title=None,
                    extensions=("png",),
                    feature_start_column_id="start",
                    feature_end_column_id="end",
                    feature_color_column_id="color",
                    feature_length_column_id="length",
                    feature_strand_column_id="strand",
                    feature_shape="rectangle",
                    feature_height_fraction=0.7,
                    stranded_tracks=False,
                    rounded_tracks=False,
                    stranded_end_tracks=False,
                    middle_break=False,
                    fill_empty_tracks=True,
                    empty_color="lightgrey",
                    subplot_scale=False,
                    track_group_scale=False,
                    track_group_label_fontstyle='normal',
                    track_group_distance=2,
                    x_tick_type="nucleotide",
                    xmax_multiplier=1.1, ymax_multiplier=1.1,
                    xtick_fontsize=None,
                    subplot_title_fontsize=None,
                    subplot_title_fontweight='bold',


                    show_vertical_track_group_label=True,
                    show_horizontal_track_group_label=True,
                    show_track_label=True,
                    close_figure=False,
                    save_figure=False,
                    figure=None,
                    axes=None,
                    figure_width=15,
                    figure_height_per_scaffold=0.5,
                    figure_header_height=0,
                    dpi=300,
                    subplots_adjust_left=None,
                    subplots_adjust_bottom=None,
                    subplots_adjust_right=None,
                    subplots_adjust_top=None,
                    autoscale_figure=True
                    ):
        vertical_track_group_number = len(track_dict)
        vertical_track_group_number_dict = {vertical_track_group: len(track_dict[vertical_track_group]) for vertical_track_group in track_dict}
        horizontal_track_group_number = sum(vertical_track_group_number_dict.values())


        track_group_dict = OrderedDict()
        #print(scaffold_order_list)
        scaffolds = scaffold_order_list.to_list() if isinstance(scaffold_order_list, (pd.Series, pd.Index)) else scaffold_order_list  # scaffold_order_list[::-1] if scaffold_order_list else collection_gff.records.index.get_level_values(level=0).unique().to_list()
        scaffold_number = len(scaffolds)

        synteny_feature_track_style = TrackStyle(height=10, colormap=None, background="white",
                                                 masked="grey", fill_empty=fill_empty_tracks, empty_color=empty_color,
                                                 stranded=stranded_tracks,
                                                 rounded=rounded_tracks,
                                                 stranded_end=stranded_end_tracks,
                                                 middle_break=middle_break,
                                                 centromere=True if centromere_df is not None else False)

        feature_height = 5 if stranded_tracks else 10

        if feature_shape == "rectangle":
            feature_style = FeatureStyle(patch_type="rectangle", height=feature_height, label_fontsize=10,
                                         face_color=default_color)
        elif feature_shape == "circle":
            feature_style = FeatureStyle(patch_type="circle", height=feature_height * feature_height_fraction,
                                         label_fontsize=10, face_color=default_color)
        elif feature_shape == "ellipse":
            feature_style = FeatureStyle(patch_type="ellipse", height=feature_height * feature_height_fraction,
                                         label_fontsize=10, face_color=default_color)
        elif feature_shape == "hist":
            feature_style = FeatureStyle(patch_type="hist", height=feature_height, label_fontsize=10, face_color=default_color)
        else:
            raise ValueError("ERROR!!! Unknown feature style")

        #feature_style = FeatureStyle(patch_type="rectangle", height=feature_height, label_fontsize=10)
        track_number = 0
        for chr in scaffolds:  # count_df.index.get_level_values(level=0).unique():

            highlight = False
            highlight_color = None
            if (highlight_df is not None) and (not highlight_df.empty):
                if chr in highlight_df.index:
                    highlight = True
                    highlight_color = highlight_df.loc[chr, "color"]

            track_group_style = TrackGroupStyle(label_fontstyle=track_group_label_fontstyle,
                                                distance=track_group_distance,
                                                highlight_color=highlight_color)

            track_group_dict[chr] = TrackGroup(label=chr if show_trackgroup_label else None,
                                               style=track_group_style,
                                               highlight=highlight
                                               )
            centromere_start = None
            centromere_end = None
            if (centromere_df is not None) and (not centromere_df.empty):
                if chr in centromere_df.index:
                    #print(chr)
                    centromere_start = centromere_df.loc[chr, "start"]
                    centromere_end = centromere_df.loc[chr, "end"]

            for species in bed_collection_dict:
                records = bed_collection_dict[species].records if hasattr(bed_collection_dict[species], "records") else bed_collection_dict[species]

                # print(species)
                #print(scaffold_length_df)
                #print(scaffold_length_df)
                #print(scaffold_length_df.loc[chr])
                #print(scaffold_length_df.loc[chr][0])
                track_group_dict[chr][species] = FeatureTrack(
                    records.loc[[chr]] if chr in records.index else None, x_end=scaffold_length_df.loc[chr].iloc[0],
                    label=species if show_track_label else None, #colormap=colormap, thresholds=thresholds,
                    style=synteny_feature_track_style,
                    #colors=colors, background=background,
                    feature_style=feature_style,
                    feature_start_column_id=feature_start_column_id,
                    feature_end_column_id=feature_end_column_id,
                    feature_color_column_id=feature_color_column_id,
                    feature_length_column_id=feature_length_column_id,
                    feature_strand_column_id=feature_strand_column_id,
                    subplot_scale=subplot_scale,
                    track_group_scale=track_group_scale,
                    stranded=stranded_tracks,
                    middle_break=middle_break,
                    centromere_start=centromere_start,
                    centromere_end=centromere_end)
                track_number += 1
                # print(track_group_dict[chr][species].records)
                #if feature_color_column_id not in records.columns:
                #    track_group_dict[chr][species].add_color_by_dict(default_color=default_color) if default_color else \
                #    track_group_dict[chr][species].add_color_by_dict()
        #print(track_number)
        subplot_style = SubplotStyle(distance=5, xaxis_visible=True, yaxis_visible=False, spines_bottom_visible=True,
                                     spines_right_visible=False, spines_left_visible=False, spines_top_visible=False,
                                     x_tick_type=x_tick_type,
                                     title_fontsize=subplot_title_fontsize,
                                     title_fontweight=subplot_title_fontweight,
                                     x_tick_major_fontsize=xtick_fontsize,
                                     x_tick_minor_fontsize=(xtick_fontsize - 1) if xtick_fontsize is not None else None)
        chromosome_subplot = Subplot(track_group_dict, title=title, style=subplot_style,
                                     legend=legend, #ChromosomeLegend(chromosome_df_dict=species_color_df_dict, scaffold_order_list=scaffold_order_list),
                                     auto_scale=True,
                                     figure_x_y_ratio=figure_width / max(1, int(track_number * figure_height_per_scaffold + figure_header_height)),
                                     xmax_multiplier=xmax_multiplier, ymax_multiplier=ymax_multiplier)
        #print((figure_width,
        #                       max(1, int(scaffold_number * figure_height_per_scaffold + figure_header_height))))
        #print((scaffold_number, figure_height_per_scaffold, figure_header_height))

        if figure is not None:
            current_figure = figure
        else:
            if axes is None:
                current_figure = plt.figure(1, figsize=(figure_width,
                                            max(1, int(horizontal_track_group_number * figure_height_per_scaffold + figure_header_height))), dpi=dpi)
            else:
                current_figure = None

        chromosome_subplot.draw()

        if (subplots_adjust_left is not None) or (subplots_adjust_right is not None) or (subplots_adjust_top is not None) or (subplots_adjust_bottom is not None):
            plt.subplots_adjust(left=subplots_adjust_left, right=subplots_adjust_right,
                                top=subplots_adjust_top, bottom=subplots_adjust_bottom)
        if save_figure:
            if output_prefix is None:
                raise ValueError("ERROR!!! Can't save figure - output prefix is not set")
            for ext in extensions:
                if autoscale_figure:
                    plt.savefig("%s.%s" % (output_prefix, ext), bbox_inches="tight", dpi=dpi)
                else:
                    plt.savefig(f"{output_prefix}.{ext}")

        if (current_figure is not None) and close_figure:
            plt.close(current_figure.number)

        else:
            return current_figure

