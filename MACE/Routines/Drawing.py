#!/usr/bin/env python

import os

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


class DrawingRoutines():

    @staticmethod
    def draw_variant_window_densities(count_dict, scaffold_length_dict, window_size, window_step, output_prefix,
                                      record_style=None, ext_list=("svg", "png"),
                                      label_fontsize=13, left_offset=0.2, figure_width=8, scaffold_synonym_dict=None,
                                      id_replacement_mode="partial", suptitle=None, density_multiplicator=1000,
                                      colormap_tuple_list=((0.0, "#333a97"), (0.1, "#3d3795"), (0.5, "#5d3393"),
                                                           (0.75, "#813193"), (1.0, "#9d2d7f"), (1.25, "#b82861"),
                                                           (1.5, "#d33845"), (2.0, "#ea2e2e"), (2.5, "#f5ae27"))):

        scaffold_number = len(count_dict)
        max_scaffold_length = max(scaffold_length_dict.values())

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
        for scaffold in count_dict:
            #print gap_coords_list, gap_len_list

            start_y += scaffold_height + dist_between_scaffolds
            label_y_start = label_line_y_shift + start_y
            gap_y_jump = label_y_start + label_line_y_jump
            prev_x = 0
            """
            figure.text(0, start_y, scaffold, rotation=0, fontweight="bold", transform=subplot.transAxes, fontsize=9,
                         horizontalalignment='center',
                         verticalalignment='center')
            """
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


        plt.xlim(xmin=0, xmax=int(1.05 * max_scaffold_length))
        plt.ylim(ymin=0, ymax=start_y + 2 * scaffold_height)
        #plt.colorbar(subplot)
        #plt.tight_layout()
        plt.subplots_adjust(left=left_offset, right=0.95)#bottom=0.1, right=0.8, top=0.9)
        if suptitle:
            plt.suptitle(suptitle)
        for extension in ext_list:
            plt.savefig("%s.%s" % (output_prefix, extension))
