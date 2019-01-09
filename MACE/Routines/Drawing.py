#!/usr/bin/env python

import os
import math
from collections import Iterable, OrderedDict

import numpy as np

import matplotlib
matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
from matplotlib import colors
from matplotlib import text

from MACE.General import FileRoutines
from MACE.Functions.Generators import recursive_generator

from MACE.General.GeneralCollections import TwoLvlDict

class DrawingRoutines:

    @staticmethod
    def get_filtered_scaffold_list(count_dict,
                                   scaffold_black_list=[],
                                   sort_scaffolds=False,
                                   scaffold_ordered_list=None,
                                   scaffold_white_list=[]):
        white_set = set(scaffold_white_list)
        black_set = set(scaffold_black_list)

        scaffold_set = set()

        for sample in count_dict:
            scaffold_set |= set(count_dict[sample])

        if white_set:
            #print "XXXXXXXXXXXXXX"
            #print scaffold_set
            #print "YYYYYYYYYYYYYYYY"
            #print white_set
            scaffold_set = scaffold_set & white_set
            #print "ZZZZZZZZZZZZZZZZ"
            #print scaffold_set
        if black_set:
            #print "WWWWWWWWWWWW"
            #print black_set
            scaffold_set = scaffold_set - black_set
            #print "QQQQQQQQQQQQQQQQ"
            #print black_set

        #print "RRRRRRRRRRRRRRRRRR"
        #print scaffold_set
        scaffold_list = list(scaffold_set)
        #print "PPPPPPPPPPP"
        #print scaffold_list
        if sort_scaffolds:
            scaffold_list.sort()

        final_scaffold_list = []

        if scaffold_ordered_list:
            for entry in scaffold_ordered_list:
                if entry in scaffold_list:
                    final_scaffold_list.append(entry)
                    scaffold_list.remove(entry)
                else:
                    print("WARNING!!!Entry(%s) from order list is absent in list of scaffolds!" % entry)
            final_scaffold_list = final_scaffold_list + scaffold_list
        else:
            final_scaffold_list = scaffold_list

        return final_scaffold_list

    def draw_variant_window_densities(self, count_dict, scaffold_length_dict, window_size, window_step, output_prefix,
                                      masking_dict=None,
                                      gap_fraction_threshold=0.4,
                                      record_style=None, ext_list=("svg", "png"),
                                      label_fontsize=13, left_offset=0.2, figure_width=12,
                                      figure_height_scale_factor=0.5, scaffold_synonym_dict=None,
                                      id_replacement_mode="partial", suptitle=None, density_multiplicator=1000,
                                      scaffold_black_list=[], sort_scaffolds=False, scaffold_ordered_list=None,
                                      scaffold_white_list=[], add_sample_name_to_labels=False,
                                      dist_between_scaffolds_scaling_factor=1,
                                      gap_color="grey",
                                      masked_color="grey", no_snp_color="white",
                                      colormap=None,
                                      colors=("#333a97", "#3d3795","#5d3393", "#813193", "#9d2d7f", "#b82861",
                                                        "#d33845", "#ea2e2e", "#f5ae27"),
                                      thresholds=(0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5),
                                      colormap_tuple_list=((0.0, "#333a97"), (0.1, "#3d3795"), (0.5, "#5d3393"),
                                                           (0.75, "#813193"), (1.0, "#9d2d7f"), (1.25, "#b82861"),
                                                           (1.5, "#d33845"), (2.0, "#ea2e2e"), (2.5, "#f5ae27"))):
        """ cont_dict = {sample: {scaffold: }}"""

        if dist_between_scaffolds_scaling_factor < 1:
            raise ValueError("Scaling factor for distance between scaffolds have to be >=1.0")

        """
        white_set = set(scaffold_white_list)
        black_set = set(scaffold_black_list)

        scaffold_set = set()
        for sample in count_dict:
            scaffold_set |= set(count_dict[sample])

        if white_set:
            scaffold_set = scaffold_set & white_set

        if black_set:
            scaffold_set = scaffold_set - black_set

        scaffold_list = list(scaffold_set)

        if sort_scaffolds:
            scaffold_list.sort()

        final_scaffold_list = []
        if scaffold_ordered_list:
            for entry in scaffold_ordered_list:
                final_scaffold_list.append(entry)
                scaffold_list.remove(entry)
            final_scaffold_list = final_scaffold_list + scaffold_list
        else:
            final_scaffold_list = scaffold_list
        """

        final_scaffold_list = self.get_filtered_scaffold_list(count_dict,
                                                              scaffold_black_list=scaffold_black_list,
                                                              sort_scaffolds=sort_scaffolds,
                                                              scaffold_ordered_list=scaffold_ordered_list,
                                                              scaffold_white_list=scaffold_white_list)
        scaffold_number = len(final_scaffold_list)
        max_scaffold_length = max([scaffold_length_dict[scaf] for scaf in final_scaffold_list])
        #max_scaffold_length = max(scaffold_length_dict.values())

        figure = plt.figure(figsize=(figure_width,
                                     int(figure_height_scale_factor * scaffold_number * len(count_dict))))
        subplot = plt.subplot(1, 1, 1)

        subplot.get_yaxis().set_visible(False)
        #subplot.get_xaxis().set_visible(False)
        #axes.xaxis.set_major_formatter(x_formatter)

        #subplot.spines['bottom'].set_color('none')
        subplot.spines['right'].set_color('none')
        subplot.spines['left'].set_color('none')
        subplot.spines['top'].set_color('none')

        scaffold_height = 10

        dist_between_scaffolds = 5
        start_x = 0
        start_y = - dist_between_scaffolds

        label_line_y_shift = int(scaffold_height/2)
        label_line_y_jump = int(scaffold_height/2)

        #normalize_color_func = LinearSegmentedColormap.from_list("Densities_custom", colormap_tuple_list)
        #plt.register_cmap(cmap=colormap)
        #colormap = cm.get_cmap(name="plasma", lut=None)
        #normalize_colors = colors.BoundaryNorm(boundaries_for_colormap, len(boundaries_for_colormap) - 1) * int(256/(len(boundaries_for_colormap) - 1))
        #normalize_colors = colors.Normalize(vmin=boundaries_for_colormap[0], vmax=boundaries_for_colormap[-1])

        masked_windows_count_dict = TwoLvlDict()
        no_snps_windows_count_dict = TwoLvlDict()

        for sample in count_dict:
            masked_windows_count_dict[sample] = OrderedDict()

        if colormap:
            cmap = plt.get_cmap(colormap, len(thresholds))

        for scaffold in final_scaffold_list:

            sample_index = 0
            for sample in count_dict:
                masked_windows_count_dict[sample][scaffold] = 0
                no_snps_windows_count_dict[sample][scaffold] = 0
                #if scaffold in scaffold_black_list:
                #    continue
                #print gap_coords_list, gap_len_list

                start_y += scaffold_height + dist_between_scaffolds * (dist_between_scaffolds_scaling_factor if sample_index == 0 else 1)
                label_y_start = label_line_y_shift + start_y
                gap_y_jump = label_y_start + label_line_y_jump
                prev_x = 0

                #figure.text(0, start_y, scaffold, rotation=0, fontweight="bold", transform=subplot.transAxes, fontsize=9,
                #             horizontalalignment='center',
                #             verticalalignment='center')

                if scaffold_synonym_dict:
                    if id_replacement_mode == "exact":
                        if scaffold in scaffold_synonym_dict:
                            scaffold_label = scaffold_synonym_dict[scaffold]
                        else:
                            scaffold_label = scaffold
                            print("WARNING!!! Synonym for %s was not found" % scaffold)
                    elif id_replacement_mode == "partial":

                        partial_syn_list = []
                        for partial_syn in scaffold_synonym_dict:
                            if partial_syn in scaffold:
                                partial_syn_list.append(partial_syn)

                        if len(partial_syn_list) > 1:
                            print("WARNING!!! More than one possible replacement for %s was found: %s. No replacement then." % (scaffold, ",".join(partial_syn_list)))
                            scaffold_label = scaffold
                        elif not partial_syn_list:
                            scaffold_label = scaffold
                            print("WARNING!!! Synonym for %s was not found" % scaffold)
                        else:
                            scaffold_label = scaffold_synonym_dict[partial_syn_list[0]]
                    else:
                        raise ValueError("Unknown id replacement mode")

                else:
                    scaffold_label = scaffold

                subplot.annotate(("%s (%s)" % (scaffold, sample))if add_sample_name_to_labels else scaffold_label,
                                 xy=(0, label_y_start), xycoords='data', fontsize=16,
                                 xytext=(-15, 1.5 * label_line_y_shift), textcoords='offset points',
                                 ha='right', va='top')
                if scaffold in count_dict[sample]:
                    for i in range(0, len(count_dict[sample][scaffold])):

                        window_start = i * window_step
                        window_end = window_start + window_size - 1  # TODO: check end coordinate
                        if masking_dict:
                            if scaffold in masking_dict:
                                unmasked_length = window_size - masking_dict[scaffold][i]
                                if unmasked_length > 0:
                                    variant_density = float(count_dict[sample][scaffold][i] * density_multiplicator) / float(unmasked_length)
                                else:
                                    variant_density = None
                        else:
                            variant_density = float(count_dict[sample][scaffold][i] * density_multiplicator) / float(window_size)

                        if variant_density:
                            if colormap:
                                if variant_density <= thresholds[0]:
                                    window_color = no_snp_color
                                else:
                                    for i in range(0, len(thresholds) - 1):
                                        if thresholds[i] < variant_density <= thresholds[i+1]:
                                            window_color = cmap[i]
                                            break
                                    else:
                                        window_color = cmap[i+1]

                            else:
                                if variant_density <= colormap_tuple_list[0][0]:
                                    window_color = no_snp_color
                                else:
                                    for lower_boundary, color in colormap_tuple_list:
                                        if variant_density <= lower_boundary:
                                            break
                                        if variant_density > lower_boundary:
                                            prev_color = color
                                    else:
                                        prev_color = color
                                    window_color = prev_color
                        else:
                            window_color = masked_color

                        if masking_dict:
                            if scaffold in masking_dict:
                                if float(masking_dict[scaffold][i]) / float(window_size) > gap_fraction_threshold:
                                    window_color = masked_color
                        #print scaffold
                        #print i, variant_density, window_color

                        if window_color == masked_color:
                            masked_windows_count_dict[sample][scaffold] += 1
                        elif window_color == no_snp_color:
                            no_snps_windows_count_dict[sample][scaffold] += 1

                        window = Rectangle((window_start, start_y), window_size, scaffold_height, fill=True,
                                           edgecolor=None, facecolor=window_color, linewidth=0.0000000000001)
                        #print prev_x
                        #print gap_coords[0] - prev_x

                        subplot.add_patch(window)

                # draw_chromosome

                fragment = Rectangle((0, start_y), scaffold_length_dict[scaffold], scaffold_height, fill=False,
                                     edgecolor="black", facecolor=None, linewidth=0.5)
                subplot.add_patch(fragment)
                sample_index += 1

        legend_y_position = int(start_y/2)
        legend_x_position = int(max_scaffold_length * 1.05)
        legend_element_side = scaffold_height

        square_y_pos = legend_y_position - legend_element_side

        for color, legend_label in zip((masked_color, no_snp_color), ("masked", "no SNPs")):
            square_y_pos += legend_element_side
            fragment = Rectangle((legend_x_position, square_y_pos), max_scaffold_length/64, legend_element_side, fill=True,
                                 edgecolor="black", facecolor=color, linewidth=0.5)

            subplot.add_patch(fragment)
            subplot.annotate(legend_label,
                             xy=(legend_x_position + 2 * max_scaffold_length/64, square_y_pos), xycoords='data', fontsize=13,
                             xytext=(legend_x_position + 2 * max_scaffold_length/64, square_y_pos),)

        for i in range(0, len(colormap_tuple_list)):
            square_y_pos += legend_element_side
            #print (colormap_tuple_list[i][1])
            fragment = Rectangle((legend_x_position, square_y_pos), max_scaffold_length/64, legend_element_side, fill=True,
                                 edgecolor="black", facecolor=colormap_tuple_list[i][1], linewidth=0.5)

            subplot.add_patch(fragment)
            if i == (len(colormap_tuple_list) - 1):
                legend_element_label = "> %.2f" % colormap_tuple_list[i][0]
            else:
                legend_element_label = "%.2f - %.2f" % (colormap_tuple_list[i][0], colormap_tuple_list[i + 1][0])

            subplot.annotate(legend_element_label,
                             xy=(legend_x_position + 2 * max_scaffold_length/64, square_y_pos), xycoords='data', fontsize=13,
                             xytext=(legend_x_position + 2 * max_scaffold_length/64, square_y_pos),)

        plt.xlim(xmin=0, xmax=int(1.2 * max_scaffold_length))
        plt.ylim(ymin=0, ymax=start_y + 2 * scaffold_height)
        #plt.colorbar(subplot)
        #plt.tight_layout()

        plt.subplots_adjust(left=left_offset, right=0.95)#bottom=0.1, right=0.8, top=0.9)
        if suptitle:
            plt.suptitle(suptitle)
        for extension in ext_list:
            plt.savefig("%s.%s" % (output_prefix, extension))
        plt.close()

        no_snps_windows_count_dict.write("%s.no_snps.windows.count" % output_prefix)
        masked_windows_count_dict.write("%s.masked.windows.count" % output_prefix)

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

        FileRoutines.safe_mkdir(per_scaffold_histo_dir)

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





    """
    @staticmethod
    def draw_variant_window_densities(count_dict, scaffold_length_dict, window_size, window_step, output_prefix,
                                      record_style=None, ext_list=("svg", "png"),
                                      label_fontsize=13, left_offset=0.2, figure_width=8, scaffold_synonym_dict=None,
                                      id_replacement_mode="partial", suptitle=None, density_multiplicator=1000,
                                      scaffold_black_list=[], sort_scaffolds=False, scaffold_ordered_list=None,
                                      scaffold_white_list=[],
                                      gap_color="grey",
                                      colormap_tuple_list=((0.0, "#333a97"), (0.1, "#3d3795"), (0.5, "#5d3393"),
                                                           (0.75, "#813193"), (1.0, "#9d2d7f"), (1.25, "#b82861"),
                                                           (1.5, "#d33845"), (2.0, "#ea2e2e"), (2.5, "#f5ae27"))):

        white_set = set(scaffold_white_list)
        black_set = set(scaffold_black_list)
        scaffold_set = set(count_dict)

        if white_set:
            scaffold_set = scaffold_set & white_set

        if black_set:
            scaffold_set = scaffold_set - black_set

        scaffold_list = list(scaffold_set)

        if sort_scaffolds:
            scaffold_list.sort()

        final_scaffold_list = []
        if scaffold_ordered_list:
            for entry in scaffold_ordered_list:
                final_scaffold_list.append(entry)
                scaffold_list.remove(entry)
            final_scaffold_list = final_scaffold_list + scaffold_list
        else:
            final_scaffold_list = scaffold_list

        scaffold_number = len(final_scaffold_list)
        max_scaffold_length = max([scaffold_length_dict[scaf] for scaf in final_scaffold_list])
        #max_scaffold_length = max(scaffold_length_dict.values())

        figure = plt.figure(figsize=(scaffold_number, figure_width))
        subplot = plt.subplot(1, 1, 1)

        subplot.get_yaxis().set_visible(False)
        #subplot.get_xaxis().set_visible(False)
        #axes.xaxis.set_major_formatter(x_formatter)

        #subplot.spines['bottom'].set_color('none')
        subplot.spines['right'].set_color('none')
        subplot.spines['left'].set_color('none')
        subplot.spines['top'].set_color('none')

        scaffold_height = 10

        dist_between_scaffolds = 5
        start_x = 0
        start_y = - dist_between_scaffolds

        label_line_y_shift = int(scaffold_height/2)
        label_line_y_jump = int(scaffold_height/2)

        #normalize_color_func = LinearSegmentedColormap.from_list("Densities_custom", colormap_tuple_list)
        #plt.register_cmap(cmap=colormap)
        #colormap = cm.get_cmap(name="plasma", lut=None)
        #normalize_colors = colors.BoundaryNorm(boundaries_for_colormap, len(boundaries_for_colormap) - 1) * int(256/(len(boundaries_for_colormap) - 1))
        #normalize_colors = colors.Normalize(vmin=boundaries_for_colormap[0], vmax=boundaries_for_colormap[-1])

        for scaffold in final_scaffold_list:
            #if scaffold in scaffold_black_list:
            #    continue
            #print gap_coords_list, gap_len_list

            start_y += scaffold_height + dist_between_scaffolds
            label_y_start = label_line_y_shift + start_y
            gap_y_jump = label_y_start + label_line_y_jump
            prev_x = 0

            #figure.text(0, start_y, scaffold, rotation=0, fontweight="bold", transform=subplot.transAxes, fontsize=9,
            #             horizontalalignment='center',
            #             verticalalignment='center')

            if scaffold_synonym_dict:
                if id_replacement_mode == "exact":
                    if scaffold in scaffold_synonym_dict:
                        scaffold_label = scaffold_synonym_dict[scaffold]
                    else:
                        scaffold_label = scaffold
                        print("WARNING!!! Synonym for %s was not found" % scaffold)
                elif id_replacement_mode == "partial":

                    partial_syn_list = []
                    for partial_syn in scaffold_synonym_dict:
                        if partial_syn in scaffold:
                            partial_syn_list.append(partial_syn)

                    if len(partial_syn_list) > 1:
                        print("WARNING!!! More than one possible replacement for %s was found: %s. No replacement then." % (scaffold, ",".join(partial_syn_list)))
                        scaffold_label = scaffold
                    elif not partial_syn_list:
                        scaffold_label = scaffold
                        print("WARNING!!! Synonym for %s was not found" % scaffold)
                    else:
                        scaffold_label = scaffold_synonym_dict[partial_syn_list[0]]
                else:
                    raise ValueError("Unknown id replacement mode")

            else:
                scaffold_label = scaffold

            subplot.annotate(scaffold_label, xy=(0, label_y_start), xycoords='data', fontsize=16,
                             xytext=(-15, 1.5 * label_line_y_shift), textcoords='offset points',
                             ha='right', va='top')

            for i in range(0, len(count_dict[scaffold])):

                window_start = i * window_step
                window_end = window_start + window_size - 1  # TODO: check end coordinate

                variant_density = float(count_dict[scaffold][i] * density_multiplicator) / float(window_size)

                if variant_density <= colormap_tuple_list[0][0]:
                    window_color = "white"
                else:
                    for lower_boundary, color in colormap_tuple_list:
                        if variant_density <= lower_boundary:
                            break
                        if variant_density > lower_boundary:
                            prev_color = color
                    else:
                        prev_color = color
                    window_color = prev_color

                #print scaffold
                #print i, variant_density, window_color

                window = Rectangle((window_start, start_y), window_size, scaffold_height, fill=True,
                                   edgecolor=None, facecolor=window_color, linewidth=0.0000000000001)
                #print prev_x
                #print gap_coords[0] - prev_x

                subplot.add_patch(window)

            # draw_chromosome

            fragment = Rectangle((0, start_y), scaffold_length_dict[scaffold], scaffold_height, fill=False,
                                 edgecolor="black", facecolor=None, linewidth=0.5)
            subplot.add_patch(fragment)

        legend_y_position = int(start_y/2)
        legend_x_position = int(max_scaffold_length * 1.05)
        legend_element_side = scaffold_height

        square_y_pos = legend_y_position - legend_element_side

        for i in range(0, len(colormap_tuple_list)):
            square_y_pos = square_y_pos + legend_element_side
            #print (colormap_tuple_list[i][1])
            fragment = Rectangle((legend_x_position, square_y_pos), max_scaffold_length/64, legend_element_side, fill=True,
                                 edgecolor="black", facecolor=colormap_tuple_list[i][1], linewidth=0.5)

            subplot.add_patch(fragment)
            if i == (len(colormap_tuple_list) - 1):
                legend_element_label = "> %.2f" % colormap_tuple_list[i][0]
            else:
                legend_element_label = "%.2f - %.2f" % (colormap_tuple_list[i][0], colormap_tuple_list[i + 1][0])

            subplot.annotate(legend_element_label,
                             xy=(legend_x_position + 2 * max_scaffold_length/64, square_y_pos), xycoords='data', fontsize=13,
                             xytext=(legend_x_position + 2 * max_scaffold_length/64, square_y_pos),)

        plt.xlim(xmin=0, xmax=int(1.1 * max_scaffold_length))
        plt.ylim(ymin=0, ymax=start_y + 2 * scaffold_height)
        #plt.colorbar(subplot)
        #plt.tight_layout()

        plt.subplots_adjust(left=left_offset, right=0.95)#bottom=0.1, right=0.8, top=0.9)
        if suptitle:
            plt.suptitle(suptitle)
        for extension in ext_list:
            plt.savefig("%s.%s" % (output_prefix, extension))
        """