#!/usr/bin/env python
import math
import datetime

from copy import deepcopy
from collections import Iterable, OrderedDict

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

#from RouToolPa.Routines import FileRoutines
from RouToolPa.Collections.General import TwoLvlDict
from RouToolPa.Routines.Drawing import DrawingRoutines

from MACE.Visualization.Tracks import WindowTrack, FeatureTrack
from MACE.Visualization.TrackGroups import TrackGroup
from MACE.Visualization.Subplots import Subplot
from MACE.Visualization.Figures import Figure

from MACE.Visualization.Legends import DensityLegend, FeatureLegend, CoverageLegend
from MACE.Functions.Generators import recursive_generator
from MACE.Visualization.Styles.Subplot import chromosome_subplot_style
from MACE.Visualization.Styles.Figure import plot_figure_style, rainfall_figure_style, chromosome_figure_style, one_plot_figure_style


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
    def zygoty_bar_plot(zygoty_counts, output_prefix, extension_list=("png",),
                        figsize=(5, 5), dpi=200, title=None, color_dict=None):

        default_color_dict = OrderedDict({
                                          "homo": "orange",
                                          "hetero": "blue",
                                          "ref": "green",
                                          "absent": "red",
                                          })

        colors = color_dict if color_dict else default_color_dict

        #zygoty_counts = self.count_zygoty(outfile="%s.counts" % output_prefix)
        df_shape = np.shape(zygoty_counts)
        fig = plt.figure(1, figsize=figsize, dpi=dpi)

        bar_width = 1.0 / (df_shape[0] + 1)
        bin_coord = np.arange(df_shape[1])

        for i in range(0, df_shape[0]):
            plt.bar(bin_coord + i * bar_width,
                    zygoty_counts.loc[zygoty_counts.index[i]],
                    width=bar_width, edgecolor='white',
                    color=default_color_dict[zygoty_counts.index[i]],
                    label=zygoty_counts.index[i])

        plt.ylabel('Variants', fontweight='bold')
        plt.xlabel('Sample', fontweight='bold')
        plt.xticks([coord + bar_width for coord in range(len(bin_coord))], zygoty_counts.columns,
                   rotation=45)
        if title:
            plt.title(title, fontweight='bold')
        plt.legend()
        for extension in extension_list:
            plt.savefig("%s.%s" % (output_prefix, extension), bbox_inches='tight')
        plt.close()

        return zygoty_counts

    @staticmethod
    def singleton_bar_plot(singleton_counts, output_prefix, extension_list=("png",),
                           figsize=(5, 5), dpi=200, title=None, color="blue"):

        fig = plt.figure(1, figsize=figsize, dpi=dpi)

        bar_width = 0.5
        bin_coord = np.arange(len(singleton_counts))

        for i in range(0, df_shape[0]):
            plt.bar(bin_coord + i * bar_width,
                    zygoty_counts.loc[zygoty_counts.index[i]],
                    width=bar_width, edgecolor='white',
                    color=default_color_dict[zygoty_counts.index[i]],
                    label=zygoty_counts.index[i])

        plt.ylabel('Variants', fontweight='bold')
        plt.xlabel('Sample', fontweight='bold')
        plt.xticks([coord + bar_width for coord in range(len(bin_coord))], zygoty_counts.columns,
                   rotation=45)
        if title:
            plt.title(title, fontweight='bold')
        plt.legend()
        for extension in extension_list:
            plt.savefig("%s.%s" % (output_prefix, extension), bbox_inches='tight')
        plt.close()

        return zygoty_counts

    def draw_variant_window_densities(self, count_df, window_size, window_step, scaffold_length_df,
                                      output_prefix,
                                      figure_width=15, figure_height_per_scaffold=0.5, dpi=300,
                                      colormap=None, thresholds=None, colors=None, background=None, masked=None,
                                      title=None,
                                      extensions=("png", ),
                                      scaffold_order_list=None,
                                      test_colormaps=False,
                                      masking=True,
                                      multiplier=1000,
                                      subplots_adjust_left=None,
                                      subplots_adjust_bottom=None,
                                      subplots_adjust_right=None,
                                      subplots_adjust_top=None,
                                      ):

        self.draw_windows(count_df, window_size, window_step, scaffold_length_df,
                          output_prefix,
                          figure_width=figure_width,
                          figure_height_per_scaffold=figure_height_per_scaffold, dpi=dpi,
                          colormap=colormap, thresholds=thresholds,
                          colors=colors, background=background, masked=masked,
                          title=title,
                          extensions=extensions,
                          scaffold_order_list=scaffold_order_list,
                          test_colormaps=test_colormaps,
                          masking=masking,
                          multiplier=multiplier,
                          norm=True,
                          subplots_adjust_left=subplots_adjust_left,
                          subplots_adjust_bottom=subplots_adjust_bottom,
                          subplots_adjust_right=subplots_adjust_right,
                          subplots_adjust_top=subplots_adjust_top,
                          )

    def draw_coverage_windows(self, count_df, window_size, window_step, scaffold_length_df,
                              mean_coverage,
                              output_prefix,
                              figure_width=15, figure_height_per_scaffold=0.5, dpi=300,
                              colormap=None, thresholds=(0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5),
                              colors=None, background=None,
                              title=None,
                              extensions=("png", ),
                              scaffold_order_list=None,
                              test_colormaps=False,
                              multiplier=None,
                              absolute_coverage_values=False,
                              subplots_adjust_left=None,
                              subplots_adjust_bottom=None,
                              subplots_adjust_right=None,
                              subplots_adjust_top=None,
                              ):
        if absolute_coverage_values:
            cov_df = count_df
            final_thresholds = list(np.array(thresholds) * float(mean_coverage))

        else:
            cov_df = count_df.copy(deep=True)
            cov_df = cov_df / mean_coverage
            final_thresholds = list(np.array(thresholds))

        self.draw_windows(cov_df, window_size, window_step, scaffold_length_df,
                          output_prefix,
                          figure_width=figure_width,
                          plot_type="coverage",
                          figure_height_per_scaffold=figure_height_per_scaffold, dpi=dpi,
                          colormap=colormap, thresholds=final_thresholds,
                          colors=colors, background=background, masked=False,
                          title=title,
                          extensions=extensions,
                          scaffold_order_list=scaffold_order_list,
                          test_colormaps=test_colormaps,
                          masking=False,
                          multiplier=multiplier,
                          norm=False,
                          subplots_adjust_left=subplots_adjust_left,
                          subplots_adjust_bottom=subplots_adjust_bottom,
                          subplots_adjust_right=subplots_adjust_right,
                          subplots_adjust_top=subplots_adjust_top,
                          )

    def draw_windows(self, count_df, window_size, window_step, scaffold_length_df,
                     output_prefix,
                     plot_type="densities",
                     figure_width=15, figure_height_per_scaffold=0.5, dpi=300,
                     colormap=None, thresholds=None, colors=None, background=None, masked=None,
                     title=None,
                     extensions=("png", ),
                     scaffold_order_list=None,
                     test_colormaps=False,
                     masking=True,
                     multiplier=None,
                     norm=False,
                     subplots_adjust_left=None,
                     subplots_adjust_bottom=None,
                     subplots_adjust_right=None,
                     subplots_adjust_top=None,):

        track_group_dict = OrderedDict()
        window_step_final = window_step if window_step else window_size
        scaffolds = scaffold_order_list[::-1] if scaffold_order_list else count_df.index.get_level_values(level=0).unique().to_list()
        scaffold_number = len(scaffolds)

        if test_colormaps:
            for colormap_entry in self.colormap_list:
                print("%s\tDrawing using %s colormap..." % (str(datetime.datetime.now()), colormap_entry))
                if plot_type == "densities":
                    legend = DensityLegend(colormap=colormap_entry,
                                           thresholds=thresholds)
                elif plot_type == "coverage":
                    legend = CoverageLegend(colormap=colormap_entry,
                                            thresholds=thresholds)
                # TODO: replace color recalculation for whole dataframe by replacenments in category
                for chr in scaffolds: # count_df.index.get_level_values(level=0).unique():
                    track_group_dict[chr] = TrackGroup(
                        {chr: WindowTrack(count_df.xs(chr), window_size, window_step_final, x_end=scaffold_length_df.loc[chr][0],
                                          multiplier=multiplier, label=chr, colormap=colormap_entry, thresholds=thresholds,
                                          colors=colors, background=background, masked=masked, norm=norm)})
                    track_group_dict[chr][chr].add_color(masking=masking)
                # track_group_dict
                # track_group_dict["chr13"]
                chromosome_subplot = Subplot(track_group_dict,
                                             title=(title + " (colormap %s)" % colormap_entry) if title else "Colormap %s" % colormap_entry,
                                             style=chromosome_subplot_style,
                                             legend=legend)

                plt.figure(1, figsize=(figure_width, int(scaffold_number*figure_height_per_scaffold)), dpi=dpi)

                chromosome_subplot.draw()
                plt.subplots_adjust(left=subplots_adjust_left, right=subplots_adjust_right,
                                    top=subplots_adjust_top, bottom=subplots_adjust_bottom)

                for ext in extensions:
                    plt.savefig("%s.%s.%s" % (output_prefix, colormap_entry, ext))
                plt.close(1)
        else:
            if plot_type == "densities":
                legend = DensityLegend(colormap=colormap,
                                       thresholds=thresholds)
            elif plot_type == "coverage":
                legend = CoverageLegend(colormap=colormap,
                                        thresholds=thresholds)

            for chr in scaffolds:  # count_df.index.get_level_values(level=0).unique():
                print(scaffold_length_df)
                print(scaffold_length_df[chr])
                track_group_dict[chr] = TrackGroup(
                    {chr: WindowTrack(count_df.xs(chr), window_size, window_step_final, x_end=scaffold_length_df.loc[chr][0],
                                      multiplier=multiplier, label=chr, colormap=colormap, thresholds=thresholds,
                                      colors=colors, background=background, masked=masked, norm=norm)})
                track_group_dict[chr][chr].add_color(masking=masking)
            # track_group_dict
            # track_group_dict["chr13"]
            chromosome_subplot = Subplot(track_group_dict, title=title, style=chromosome_subplot_style,
                                         legend=legend)

            plt.figure(1, figsize=(figure_width, int(scaffold_number * figure_height_per_scaffold)), dpi=dpi)

            chromosome_subplot.draw()
            plt.subplots_adjust(left=subplots_adjust_left, right=subplots_adjust_right,
                                top=subplots_adjust_top, bottom=subplots_adjust_bottom)
            for ext in extensions:
                plt.savefig("%s.%s" % (output_prefix, ext))

    # ----------------------- In progress ------------------------------
    @staticmethod
    def plot_clustering_threshold_tests(cluster_df, output_prefix, scaffold_order_list=None,
                                        extensions=("png", ), suptitle="Test of clustering thresholds",
                                        figure_width=16, figure_height=16, y_logbase=10):
        threshold_number = len(cluster_df.columns)

        figure_style = deepcopy(plot_figure_style)
        figure_style.width = figure_width
        figure_style.height = figure_height

        cluster_number_df = cluster_df.groupby(level=0).nunique()

        cluster_count_list = [cluster_df.groupby(level=0).agg(lambda s: sum(s.value_counts() == i)) for i in (1, 2, 3, 5)]

        triple_plus_count_df = cluster_df.groupby(level=0).agg(lambda s: sum(s.value_counts() > 2))
        five_plus_count_df = cluster_df.groupby(level=0).agg(lambda s: sum(s.value_counts() > 4))

        scaffolds = scaffold_order_list if scaffold_order_list else cluster_df.index.get_level_values(
            level=0).unique().to_list()
        scaffold_number = len(scaffolds)

        subplot_dict = OrderedDict([(chr, None) for chr in scaffolds])

        figure = Figure(subplots=subplot_dict, style=figure_style, suptitle=suptitle)
        figure.draw()
        df_list = [cluster_number_df] + cluster_count_list + [triple_plus_count_df, five_plus_count_df]
        label_list = ["All", "1", "2", "3", "5", "3+", "5+"]
        color_list = ["blue", "red", "orange", "magenta", "green", "black", "brown"]
        for (scaffold, subplot_index) in zip(subplot_dict.keys(), range(0, scaffold_number)):
            subplot_hor_index = subplot_index % figure.horizontal_subplot_number
            subplot_vert_index = subplot_index // figure.horizontal_subplot_number
            axes = figure.axes_array[subplot_vert_index][subplot_hor_index]
            for data, label, color in zip(df_list, label_list, color_list):
                #print data
                if scaffold in data.index:
                    #print scaffold
                    #print data.loc[scaffold]
                    axes.plot(data.columns, data.loc[scaffold], label=label,
                              color=color)
                    axes.grid()
            if (subplot_hor_index == (figure.horizontal_subplot_number - 1)) and (subplot_vert_index == 0):
                axes.legend()

            axes.set_title(scaffold)

        if output_prefix:
            for i, j in zip((0, 1, 2, 3), (1, 2, 3, 5)):
                cluster_count_list[i].to_csv("%s.cluster.%i.counts" % (output_prefix, j), sep="\t", index_label="scaffold")
            cluster_df.to_csv("%s.cluster" % output_prefix, sep="\t", index_label="scaffold")
            cluster_number_df.to_csv("%s.cluster.all.counts" % output_prefix, sep="\t", index_label="scaffold")
            triple_plus_count_df.to_csv("%s.cluster.3plus.counts" % output_prefix, sep="\t", index_label="scaffold")
            five_plus_count_df.to_csv("%s.cluster.5plus.counts" % output_prefix, sep="\t", index_label="scaffold")

        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

        plt.yscale('log', basey=y_logbase)
        plt.ylim(ymin=0.1)
        for ext in extensions:
            plt.savefig("%s.log%i.%s" % (output_prefix, y_logbase, ext))

        plt.close()

        all_chr = map(lambda s: s.sum(), df_list)
        figure = Figure(subplots={"all": []}, style=one_plot_figure_style, suptitle=suptitle)
        figure.draw()

        axes = plt.gca()

        for data, label, color in zip(all_chr, label_list, color_list):
            # print data
            axes.plot(data.index.tolist(), data, label=label,
                      color=color)
        axes.grid()
        axes.legend()

        axes.set_title("All")

        for ext in extensions:
            plt.savefig("%s.all.%s" % (output_prefix, ext))

        plt.yscale('log', basey=y_logbase)
        plt.ylim(ymin=0.1)

        for ext in extensions:
            plt.savefig("%s.all.log.%i.%s" % (output_prefix, y_logbase, ext))

    def rainfall_plot(self, distance_df, scaffold_length_df,
                      output_prefix,
                      figure_width=15,
                      figure_height_per_scaffold=0.5,
                      dpi=300,
                      colormap='jet', title=None,
                      extensions=("png",),
                      scaffold_order_list=None
                      ):
        subplot_group_dict = OrderedDict()

        scaffolds = scaffold_order_list[::-1] if scaffold_order_list else distance_df.index.get_level_values(
            level=0).unique().to_list()
        scaffold_number = len(scaffolds)

        pass

    def draw_gff(self, collection_gff, scaffold_length_df,
                 output_prefix,
                 figure_width=15, figure_height_per_scaffold=0.5, dpi=300,
                 colormap=None, thresholds=None, colors=None, background=None,
                 title=None,
                 extensions=("png", ),
                 scaffold_order_list=None,
                 ):

        track_group_dict = OrderedDict()

        scaffolds = scaffold_order_list[::-1] if scaffold_order_list else collection_gff.records.index.get_level_values(level=0).unique().to_list()
        scaffold_number = len(scaffolds)

        for chr in scaffolds:  # count_df.index.get_level_values(level=0).unique():
            track_group_dict[chr] = TrackGroup(
                {chr: FeatureTrack(collection_gff.records.xs(chr), x_end=scaffold_length_df.loc[chr][0],
                                   label=chr, colormap=colormap, thresholds=thresholds,
                                   colors=colors, background=background)})
            track_group_dict[chr][chr].add_color_by_dict()
        # track_group_dict
        # track_group_dict["chr13"]
        chromosome_subplot = Subplot(track_group_dict, title=title, style=chromosome_subplot_style,
                                     legend=FeatureLegend(colormap=colormap,
                                                          thresholds=thresholds))

        plt.figure(1, figsize=(figure_width, int(scaffold_number * figure_height_per_scaffold)), dpi=dpi)

        chromosome_subplot.draw()

        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))



    # ------------------ Not rewritten yet --------------------------------

    """
    self.get_filtered_scaffold_list(last_collection.target_scaffold_list,
                                                               scaffold_black_list=target_black_list,
                                                               sort_scaffolds=False,
                                                               scaffold_ordered_list=target_ordered_list,
                                                               scaffold_white_list=target_white_list)
    """

    def draw_window_density_distribution(self, count_dict, window_size, output_prefix=None, suptitle="SNP density distribution",
                                         density_multiplicator=1000,
                                         number_of_bins=None, width_of_bins=None,
                                         max_threshold=None, min_threshold=None,
                                         scaffold_black_list=[], scaffold_white_list=[],
                                         sort_scaffolds=False, scaffold_ordered_list=None, subplot_size=4,
                                         per_scaffold_histo_dir="per_scaffold_histo_dir/",
                                         subplot_tuple=None, share_x_axis=True, share_y_axis=True,
                                         extensions=("png",), show_mean_and_median=True):
        """
        scaffold_threshold: if number of scaffolds is higher draw only separated_histograms
        """
        samples_list = count_dict.keys()
        final_scaffold_list = self.get_filtered_scaffold_list(count_dict,
                                                              scaffold_black_list=scaffold_black_list,
                                                              sort_scaffolds=sort_scaffolds,
                                                              scaffold_ordered_list=scaffold_ordered_list,
                                                              scaffold_white_list=scaffold_white_list)

        scaffold_number = len(final_scaffold_list)

        self.safe_mkdir(per_scaffold_histo_dir)

        xlabel = "Number of SNPs"
        ylabel = "Number of windows"

        scaled_count_dict = OrderedDict()

        empty_windows_scaffold_dict = OrderedDict()
        for scaffold in final_scaffold_list:
            for sample in count_dict:
                if scaffold not in count_dict[sample]:
                    continue
                empty_windows_scaffold_dict[scaffold] = np.zeros(len(count_dict[sample][scaffold]))
                break

        for sample in samples_list:
            scaled_count_dict[sample] = OrderedDict()
            for scaffold in final_scaffold_list:
                if scaffold not in count_dict[sample]:
                    scaled_count_dict[sample][scaffold] = empty_windows_scaffold_dict[scaffold]
                scaled_count_dict[sample][scaffold] = np.array(map(float, count_dict[sample][scaffold])) * density_multiplicator / window_size

        print("Drawing separated histograms for each scaffold...")
        extended_label_dict = OrderedDict()
        for scaffold in final_scaffold_list:
            print("Drawing histogram for scaffold %s" % scaffold)
            #scaffold_data = [scaled_count_dict[sample][scaffold] if scaffold in scaled_count_dict[sample] else empty_windows_scaffold_dict[scaffold] for sample in samples_list]
            scaffold_data = [scaled_count_dict[sample][scaffold] for sample in samples_list]

            out_prefix = "%s/%s.%s" % (per_scaffold_histo_dir, output_prefix, scaffold) if output_prefix else "%s/%s" % (per_scaffold_histo_dir, scaffold)
            for sample in samples_list:
                median = np.median(scaled_count_dict[sample][scaffold])
                mean = np.mean(scaled_count_dict[sample][scaffold])
                extended_label = "%s: Med. %.2f, Avg: %.2f" % (sample, float(median), float(mean))
                print(extended_label)
                if scaffold in extended_label_dict:
                    extended_label_dict[scaffold].append(extended_label)
                else:
                    extended_label_dict[scaffold] = [extended_label]
            #print scaffold_data
            self.draw_histogram(scaffold_data, output_prefix=out_prefix, number_of_bins=number_of_bins,
                                width_of_bins=width_of_bins,
                                max_threshold=max_threshold, min_threshold=min_threshold,
                                xlabel=xlabel, ylabel=ylabel,
                                title=scaffold, extensions=extensions, ylogbase=None, subplot=None, suptitle=None,
                                close_figure=True,
                                data_label_list=extended_label_dict[scaffold] if show_mean_and_median else samples_list)
        #print scaled_count_dict
        print("Drawing histograms for all scaffolds on same figure...")
        data = list(recursive_generator(scaled_count_dict))
        min_value = min(data) if data else 0
        max_value = max(data) if data else 0
        #print len(scaled_count_dict)
        #print data
        bin_array = self.generate_bin_array(data, y_list=None, bin_number=number_of_bins,
                                            bin_width=width_of_bins, bin_array=None,
                                            min_x_value=min_threshold, max_x_value=max_threshold, min_y_value=None,
                                            max_y_value=None, add_max_value=True)

        plt.suptitle(suptitle)
        if subplot_tuple is None:
            side = math.sqrt(scaffold_number)
            rounded_side = int(side)
            side = rounded_side + 1 if side % rounded_side else rounded_side
            subplot_tupleee = (side, side)
            #print subplot_tupleee
        else:
            subplot_tupleee = subplot_tuple

            if len(subplot_tupleee) != 2:
                raise ValueError("Subplot tuple should contain exactly two values, not %i!" % len(subplot_tuple))
            if not (isinstance(subplot_tuple[0], int) and isinstance(subplot_tuple[1], int)):
                raise ValueError("Subplot tuple should contain two values, not (%s, %s)!" % (str(type(subplot_tuple[0])),
                                                                                             str(type(subplot_tuple[1]))))

        figure = plt.figure(256, figsize=(subplot_size * subplot_tupleee[0],
                                          subplot_size * subplot_tupleee[1]),
                            dpi=200)
        print (subplot_size * subplot_tupleee[0], subplot_size * subplot_tupleee[1])
        number_of_subplots = subplot_tupleee[0] * subplot_tupleee[1]
        subplot_list = []
        for dataset_index in range(0, len(final_scaffold_list)):
            scaffold = final_scaffold_list[dataset_index]
            if dataset_index > 0:
                if share_x_axis and share_y_axis:
                    subplot_list.append(figure.add_subplot(subplot_tupleee[0],
                                                           subplot_tupleee[1],
                                                           dataset_index + 1,
                                                           sharex=subplot_list[0],
                                                           sharey=subplot_list[0]))
                elif share_x_axis:
                    subplot_list.append(figure.add_subplot(subplot_tupleee[0],
                                                           subplot_tupleee[1],
                                                           dataset_index + 1,
                                                           sharex=subplot_list[0]))
                elif share_y_axis:
                    subplot_list.append(figure.add_subplot(subplot_tupleee[0],
                                                           subplot_tupleee[1],
                                                           dataset_index + 1,
                                                           sharex=subplot_list[0],
                                                           sharey=subplot_list[0]))
                else:
                    subplot_list.append(figure.add_subplot(subplot_tupleee[0],
                                                           subplot_tupleee[1],
                                                           dataset_index + 1))
            else:
                subplot_list.append(figure.add_subplot(subplot_tupleee[0],
                                                       subplot_tupleee[1],
                                                       dataset_index + 1))
            #print dataset_index + 1
            #print subplot_tupleee[0] * (subplot_tupleee[1] - 1)
            #print ((dataset_index + 1) > (subplot_tupleee[0] * (subplot_tupleee[1] - 1)))

            histo = self.draw_histogram([scaled_count_dict[sample][scaffold] for sample in samples_list],
                                        number_of_bins=None,
                                        width_of_bins=None, max_threshold=None,
                                        min_threshold=None,
                                        xlabel=xlabel if ((dataset_index + 1) > (subplot_tupleee[0] * (subplot_tupleee[1] - 1))) else None,
                                        ylabel=ylabel if ((dataset_index + 1) % subplot_tupleee[0]) == 1 else None,
                                        title=scaffold, extensions=("png",), ylogbase=None,
                                        subplot=subplot_list[dataset_index],
                                        suptitle=None,
                                        data_label_list=extended_label_dict[scaffold] if show_mean_and_median else samples_list,
                                        bin_array=bin_array)
            plt.xlim(xmin=min_threshold if min_threshold and (min_threshold >= min_value) else min_value,
                     xmax=max_threshold if max_threshold and (max_threshold <= max_value) else max_value)
            #print histo
            """
            if output_prefix:
                output_histo_file = "%s.%s.%shisto" % (output_prefix,
                                                       dataset_index if parameters[8] is None else parameters[10],
                                                       ("log%i." % parameters[7]) if parameters[7] else "")
                np.savetxt(output_histo_file, histo, fmt="%i\t%i")
            """
        plt.tight_layout()
        if output_prefix:
            for ext in extensions:
                plt.savefig("%s.%s" % (output_prefix, ext))
        plt.close(figure)

        print("Drawing combined histogram for all scaffolds...")

        combined_count_dict = OrderedDict()
        extended_combined_label_list = []
        for sample in samples_list:
            combined_count_dict[sample] = []
            for scaffold in count_dict[sample]:
                combined_count_dict[sample] = combined_count_dict[sample] + count_dict[sample][scaffold]

            combined_count_dict[sample] = np.array(map(float, combined_count_dict[sample]))* density_multiplicator / window_size

            median = np.median(combined_count_dict[sample])
            mean = np.mean(combined_count_dict[sample])
            extended_label = "%s: Med. %.2f, Avg: %.2f" % (sample, float(median), float(mean))
            print(extended_label)
            extended_combined_label_list.append(extended_label)

        #print combined_count_dict
        figure = plt.figure(384, figsize=(8,8))
        self.draw_histogram([combined_count_dict[sample] for sample in combined_count_dict],
                            output_prefix="%s.combined" % output_prefix if output_prefix else "combined",
                            number_of_bins=number_of_bins,
                            width_of_bins=width_of_bins,
                            max_threshold=max_threshold, min_threshold=min_threshold,
                            xlabel=xlabel, ylabel=ylabel,
                            title="SNP density distribution(all scaffolds)",
                            extensions=extensions, ylogbase=None, subplot=None, suptitle=None,
                            close_figure=True, data_label_list=extended_combined_label_list if show_mean_and_median else samples_list)

    def generate_bin_array(self, x_list, y_list=None, bin_number=20, bin_width=None, bin_array=None,
                           min_x_value=None, max_x_value=None, min_y_value=None,
                           max_y_value=None, add_max_value=True):
        if (bin_width is not None) and (bin_array is not None):
            raise ValueError("Both bin width and bin array were set")
        #print x_list
        min_x, max_x = min(map(min, x_list) if isinstance(x_list[0], Iterable) else x_list), \
                       max(map(max, x_list) if isinstance(x_list[0], Iterable) else x_list)
        if y_list:
            min_y, max_y = min(map(min, y_list) if isinstance(y_list[0], Iterable) else y_list), \
                           max(map(max, y_list) if isinstance(y_list[0], Iterable) else y_list)

        if bin_width:
            xbins = self.generate_bin_array_by_width(min_x_value if min_x_value is not None else min_x,
                                                     max_x_value if max_x_value is not None else max_x,
                                                     bin_width if isinstance(bin_width, int) else bin_width[0],
                                                     add_max_value=add_max_value)
            if y_list:
                ybins = self.generate_bin_array_by_width(min_y_value if min_y_value is not None else min_y,
                                                         max_y_value if max_y_value is not None else max_y,
                                                         bin_width if isinstance(bin_width, int) else bin_width[1],
                                                         add_max_value=add_max_value)
            bins = (xbins, ybins) if y_list else xbins

        elif bin_array:
            bins = bin_array
        else:
            if bin_number is None:
                print("WARNINNG!!! No bin_number or bin_width or bin_array were set. "
                      "Instead default value(20) for bin number is used.")
            xbins = np.linspace(min_x_value if min_x_value is not None else min_x,
                                max_x_value if max_x_value is not None else max_x,
                                20 if bin_number is None else bin_number if isinstance(bin_number, int) else bin_number[0])
            if y_list:
                ybins = np.linspace(min_y_value if min_y_value is not None else min_y,
                                    max_y_value if max_y_value is not None else max_y,
                                    20 if bin_number is None else bin_number if isinstance(bin_number, int) else bin_number[1])
            bins = (xbins, ybins) if y_list else xbins

        return bins

    @staticmethod
    def draw_histogram(data_array, output_prefix=None, number_of_bins=None, width_of_bins=None, bin_array=None,
                       max_threshold=None, min_threshold=None, xlabel=None, ylabel=None,
                       title=None, extensions=("png",), ylogbase=None, subplot=None,
                       suptitle=None, show_legend=True, close_figure=False, data_label_list=None):

        if (not(number_of_bins is None)) and (not(width_of_bins is None)) and (not (bin_array is None)):
            raise AttributeError("Options -w/--width_of_bins and -b/--number_of_bins mustn't be set simultaneously")

        if max_threshold and min_threshold:
            if max_threshold < min_threshold:
                raise ValueError("Maximum threshold (%s) is lower than minimum threshold(%s)" % (str(max_threshold),
                                                                                                 str(min_threshold)))

        if not (data_label_list is None):
            if len(data_array) != len(data_label_list):
                raise ValueError("Length of sample_list is different from number of sample arrays")
        #print "UUUUUUUUUU"
        #print data_array
        if isinstance(data_array[0], Iterable):
            max_lenn = max([max(sample) for sample in data_array])
            min_lenn = min([min(sample) for sample in data_array])
        else:
            max_lenn = max(data_array)
            min_lenn = min(data_array)

        data_names_list = data_label_list if data_label_list else ["S%i" % i for i in range(1, len(data_array) + 1)]

        max_len = max_threshold if (not(max_threshold is None)) and (max_threshold < max_lenn) else max_lenn
        min_len = min_threshold if (not(min_threshold is None)) and (min_lenn < min_threshold) else min_lenn
        filtered = []

        if (max_len < max_lenn) and (min_len > min_lenn):
            for entry in data_array:
                if min_len <= entry <= max_len:
                    filtered.append(entry)
        elif max_len < max_lenn:
            for entry in data_array:
                if entry <= max_len:
                    filtered.append(entry)
        elif min_len > min_lenn:
            for entry in data_array:
                if min_len <= entry:
                    filtered.append(entry)
        else:
            filtered = data_array
        if subplot is None:
            #print "aaaaaaaaaa"
            figure = plt.figure(1, figsize=(6, 6),)
            subplott = figure.add_subplot(1, 1, 1)
        else:
            plt.axes(subplot)
        if number_of_bins:
            bins = number_of_bins
        elif width_of_bins:
            bins = np.arange(min_len, max_len, width_of_bins)
            #print bins
            #bins[0] += 1
            bins = np.append(bins, [max_len])
        elif bin_array is not None:
            #print bin_array
            bins = bin_array
        else:
            bins = 30

        n, bins, patches = plt.hist(filtered, bins=bins, label=data_names_list) # , log=False if ylogbase is None else True)
        #print n, bins, patches
        bin_centers = (bins + ((bins[1] - bins[0])/2))[:-1]
        #print bin_centers
        #print len(n)
        #print len(bin_centers)
        #print min_len, max_len
        plt.xlim(xmin=min_len, xmax=max_len)
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)
        if title:
            plt.title(title)
        if suptitle:
            plt.suptitle(suptitle)

        if ylogbase:
            subplot.set_yscale('log', basey=ylogbase)
            plt.ylim(ymin=1)

        if show_legend:
            plt.legend(loc="best")
        if output_prefix:
            for ext in extensions:
                plt.savefig("%s.%s" % (output_prefix, ext))

                # save histo values

                #np.savetxt("%s.histo" % output_prefix, zip(bin_centers, n), fmt="%i\t%i")

        if subplot is None:
            if close_figure:
                plt.close(figure)

        return zip(bin_centers, n)

    def draw_multi_histogram_picture(self, list_of_data_arrays, subplot_tuple, output_prefix=None,
                                     figsize=(10, 10), number_of_bins_list=None, width_of_bins_list=None,
                                     bin_array_list=None,
                                     max_threshold_list=None, min_threshold_list=None, xlabel_list=None, ylabel_list=None,
                                     title_list=None, ylogbase_list=None, label_list=None,
                                     extensions=("png",), suptitle=None, share_y_axis=False,
                                     share_x_axis=False):
        figure = plt.figure(1, figsize=figsize)
        if suptitle:
            plt.suptitle(suptitle)
        if len(subplot_tuple) != 2:
            raise ValueError("Subplot tuple should contain exactly two values, not %i!" % len(subplot_tuple))
        if not (isinstance(subplot_tuple[0], int) and isinstance(subplot_tuple[1], int)):
            raise ValueError("Subplot tuple should contain two values, not (%s, %s)!" % (str(type(subplot_tuple[0])),
                                                                                         str(type(subplot_tuple[1]))))

        number_of_subplots = subplot_tuple[0] * subplot_tuple[1]
        number_of_datasets = len(list_of_data_arrays)

        parameters_list = [number_of_bins_list, width_of_bins_list, max_threshold_list, min_threshold_list,
                           xlabel_list, ylabel_list, title_list, ylogbase_list, label_list, bin_array_list]
        """
        parameter index:
        0   number_of_bins_list
        1   width_of_bins_list
        2   max_threshold_list
        3   min_threshold_list
        4   xlabel_list
        5   ylabel_list
        6   title_list
        7   ylogbase_list
        8   label_list
        9   bin_array_list
        """

        subplot_list = []
        for dataset_index in range(0, number_of_datasets):
            parameters = [None, None, None, None, None, None, None, None, None, None, None]
            for parameter_index in range(0, 10):
                if parameters_list[parameter_index]:
                    if dataset_index < len(parameters_list[parameter_index]):
                        parameters[parameter_index] = parameters_list[parameter_index][dataset_index]
            if dataset_index > 0:
                if share_x_axis and share_y_axis:
                    subplot_list.append(figure.add_subplot(subplot_tuple[0],
                                                           subplot_tuple[1],
                                                           dataset_index + 1,
                                                           sharex=subplot_list[0],
                                                           sharey=subplot_list[0]))
                elif share_x_axis:
                    subplot_list.append(figure.add_subplot(subplot_tuple[0],
                                                           subplot_tuple[1],
                                                           dataset_index + 1,
                                                           sharex=subplot_list[0]))
                elif share_y_axis:
                    subplot_list.append(figure.add_subplot(subplot_tuple[0],
                                                           subplot_tuple[1],
                                                           dataset_index + 1,
                                                           sharex=subplot_list[0],
                                                           sharey=subplot_list[0]))
                else:
                    subplot_list.append(figure.add_subplot(subplot_tuple[0],
                                                           subplot_tuple[1],
                                                           dataset_index + 1))
            else:
                subplot_list.append(figure.add_subplot(subplot_tuple[0],
                                                       subplot_tuple[1],
                                                       dataset_index + 1))

            histo = self.draw_histogram(list_of_data_arrays[dataset_index], number_of_bins=parameters[0],
                                        width_of_bins=parameters[1], max_threshold=parameters[2],
                                        min_threshold=parameters[3], xlabel=parameters[4], ylabel=parameters[5],
                                        title=parameters[6], extensions=("png",), ylogbase=parameters[7],
                                        subplot=subplot_list[dataset_index],
                                        suptitle=None, data_label_list=parameters[8], bin_array=bin_array_list[9])
            #print histo
            if output_prefix:
                output_histo_file = "%s.%s.%shisto" % (output_prefix,
                                                       dataset_index if parameters[8] is None else parameters[10],
                                                       ("log%i." % parameters[7]) if parameters[7] else "")
                np.savetxt(output_histo_file, histo, fmt="%i\t%i")

        if output_prefix:
            for ext in extensions:
                plt.savefig("%s.%s" % (output_prefix, ext))

        return figure
