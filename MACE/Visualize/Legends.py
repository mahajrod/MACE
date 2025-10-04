
from MACE.Visualization.Styles.Legend import default_legend_style


from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


class Legend:
    def __init__(self, y_start=0, y_end=None, x_start=0, x_end=None, x_size=10, element_size=10, style=default_legend_style,
                 colormap=None, fontsize=13):
        self.style = style
        self.y_start = y_start
        self.y_end = y_end
        self.x_start = x_start
        self.x_end = x_end

        self.element_size = element_size
        self.x_size = x_size

        self.colormap = colormap
        self.fontsize = fontsize

    def init_coordinates(self):
        """
        Abstract method to be defined in each legend class
        :return:
        """
        pass


class DensityLegend(Legend):

    def __init__(self, y_start=0, y_end=None, x_start=0, x_end=None, x_size=10, element_size=10, style=default_legend_style,
                 colormap=None, thresholds=np.array((0.0, 0.1, 0.25, 0.5, 1.0)),
                 colors=("#333a97", "green", "yellow", "orange", "red"), background="white", feature_name="SNPs",
                 masked="grey", fontsize=13):

        Legend.__init__(self, y_start=y_start, y_end=y_end, x_start=x_start, x_end=x_end, x_size=x_size, element_size=element_size,
                        style=style, colormap=colormap, fontsize=fontsize)

        self.colors = colors
        self.thresholds = thresholds

        if colormap:
            self.cmap = plt.get_cmap(self.colormap, len(self.thresholds))
            self.colors = [self.cmap(i) for i in range(0, len(thresholds))]

        self.background = background
        self.masked = masked
        self.feature_name = feature_name

    def init_coordinates(self, style=None):
        self.x_end = self.x_start + (2 + 5) * self.x_size
        self.y_end = self.y_start + (len(self.thresholds) + 3) * self.element_size

    def draw(self, axes=None, style=None):

        used_style = style if style else self.style

        current_subplot = axes if axes else plt.gca()

        square_y_pos = self.y_start - self.element_size

        for color, legend_label in zip((self.masked, self.background), ("masked", "no {0}".format(self.feature_name))):
            square_y_pos += self.element_size
            #print (self.x_start, square_y_pos), self.x_size, self.element_size, color
            fragment = Rectangle((self.x_start , square_y_pos), self.x_size, self.element_size,
                                 fill=True,
                                 edgecolor="black", facecolor=color, linewidth=0.5)

            current_subplot.add_patch(fragment)

            current_subplot.annotate(legend_label,
                                     xy=(self.x_start + 2 * self.x_size, square_y_pos + self.element_size*0.25), xycoords='data',
                                     fontsize=self.fontsize,
                                     xytext=(self.x_start + 2 * self.x_size, square_y_pos+ self.element_size*0.25), )

        for i in range(0, len(self.thresholds)):
            square_y_pos += self.element_size
            # print (colormap_tuple_list[i][1])
            fragment = Rectangle((self.x_start, square_y_pos), self.x_size, self.element_size,
                                 fill=True,
                                 edgecolor="black", facecolor=self.colors[i], linewidth=0.5)

            current_subplot.add_patch(fragment)
            if i == (len(self.thresholds) - 1):
                legend_element_label = "> %.2f" % self.thresholds[i]
            else:
                legend_element_label = "%.2f - %.2f" % (self.thresholds[i], self.thresholds[i + 1])

            current_subplot.annotate(legend_element_label,
                                     xy=(self.x_start + 2 * self.x_size, square_y_pos+ self.element_size*0.25), xycoords='data',
                                     fontsize=self.fontsize,
                                     xytext=(self.x_start + 2 * self.x_size, square_y_pos + self.element_size*0.25), )


class CoverageLegend(Legend):

    def __init__(self, y_start=0, y_end=None, x_start=0, x_end=None, x_size=10, element_size=10, style=default_legend_style,
                 colormap=None, thresholds=np.array((0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5)),
                 colors=("#333a97", "green", "yellow", "orange", "red"), background="white",
                 masked="grey", fontsize=13):

        Legend.__init__(self, y_start=y_start, y_end=y_end, x_start=x_start, x_end=x_end, x_size=x_size, element_size=element_size,
                        style=style, colormap=colormap, fontsize=fontsize)

        self.colors = colors
        self.thresholds = thresholds

        if colormap:
            self.cmap = plt.get_cmap(self.colormap, len(self.thresholds))
            self.colors = [self.cmap(i) for i in range(0, len(thresholds))]

        self.background = background
        self.masked = masked

    def init_coordinates(self):
        self.x_end = self.x_start + (2 + 5) * self.x_size

    def draw(self, axes=None, style=None):

        used_style = style if style else self.style

        current_subplot = axes if axes else plt.gca()

        square_y_pos = self.y_start - self.element_size

        for color, legend_label in zip((self.masked, self.background), ("masked", "no coverage")):
            square_y_pos += self.element_size
            #print (self.x_start, square_y_pos), self.x_size, self.element_size, color
            fragment = Rectangle((self.x_start, square_y_pos), self.x_size, self.element_size,
                                 fill=True,
                                 edgecolor="black", facecolor=color, linewidth=0.5)

            current_subplot.add_patch(fragment)

            current_subplot.annotate(legend_label,
                                     xy=(self.x_start + 2 * self.x_size, square_y_pos + self.element_size*0.25), xycoords='data',
                                     fontsize=self.fontsize,
                                     xytext=(self.x_start + 2 * self.x_size, square_y_pos+ self.element_size*0.25), )

        for i in range(0, len(self.thresholds)):
            square_y_pos += self.element_size
            # print (colormap_tuple_list[i][1])
            fragment = Rectangle((self.x_start, square_y_pos), self.x_size, self.element_size,
                                 fill=True,
                                 edgecolor="black", facecolor=self.colors[i], linewidth=0.5)

            current_subplot.add_patch(fragment)
            if i == (len(self.thresholds) - 1):
                legend_element_label = "> %.2fx" % self.thresholds[i]
            else:
                legend_element_label = "%.2fx - %.2fx" % (self.thresholds[i], self.thresholds[i + 1])

            current_subplot.annotate(legend_element_label,
                                     xy=(self.x_start + 2 * self.x_size, square_y_pos + self.element_size*0.25), xycoords='data',
                                     fontsize=self.fontsize,
                                     xytext=(self.x_start + 2 * self.x_size, square_y_pos + self.element_size*0.25), )


class FeatureLegend(Legend):

    def __init__(self, legend_df=None, featuretype_list=None, y_start=0, y_end=None, x_start=0, x_end=None,
                 x_size=10, element_size=10, style=default_legend_style,
                 colormap=None, fontsize=13):

        Legend.__init__(self, y_start=y_start, y_end=y_end, x_start=x_start, x_end=x_end, x_size=x_size, element_size=element_size,
                        style=style, colormap=colormap, fontsize=fontsize)

        self.legend_df = legend_df
        self.featuretype_list = featuretype_list
        if legend_df is None:
            if colormap and featuretype_list:
                self.cmap = plt.get_cmap(self.colormap, len(self.featuretype_list))

                self.legend_df = pd.DataFrame.from_dict(OrderedDict([(featuretype_list[i], self.cmap(i)) for i in range(0, len(featuretype_list))]))

    def init_coordinates(self, style=None):
        self.x_end = self.x_start + (2 + 5) * self.x_size
        self.y_end = self.y_start + len(style.legend_df if (style and (style.legend_df is not None)) else self.legend_df) * self.element_size

    def draw(self, axes=None, style=None):
        l_df = style.legend_df if (style and (style.legend_df is not None)) else self.legend_df

        current_subplot = axes if axes else plt.gca()

        square_y_pos = self.y_start - self.element_size
        #print(l_df)
        for color in l_df.index:
            square_y_pos += self.element_size
            # print (colormap_tuple_list[i][1])
            fragment = Rectangle((self.x_start, square_y_pos), self.x_size, self.element_size,
                                 fill=True,
                                 edgecolor="black", facecolor=color, linewidth=0.5)

            current_subplot.add_patch(fragment)
            #print(l_df.loc[color].iloc[0])
            current_subplot.annotate(str(l_df.loc[color].iloc[0]),
                                     xy=(self.x_start + 2 * self.x_size, square_y_pos + self.element_size * 0.25), xycoords='data',
                                     fontsize=self.fontsize,
                                     xytext=(self.x_start + 2 * self.x_size, square_y_pos + self.element_size * 0.25), )


class ChromosomeLegend(Legend):

    def __init__(self, chromosome_df_dict, scaffold_order_list, y_start=0, y_end=None, x_start=0, x_end=None,
                 x_size=10, element_size=10, style=default_legend_style,
                 colormap=None, background="white",
                 fontsize=13):

        Legend.__init__(self, y_start=y_start, y_end=y_end, x_start=x_start, x_end=x_end, x_size=x_size, element_size=element_size,
                        style=style, colormap=colormap, fontsize=fontsize)

        self.chromosome_df_dict = chromosome_df_dict
        self.scaffold_order_list = scaffold_order_list
        self.max_chromosomes = max([len(self.chromosome_df_dict[species]) for species in self.chromosome_df_dict])
        self.reference_chromosome_number = len(self.scaffold_order_list)
        if colormap:
            self.cmap = plt.get_cmap(self.colormap, self.max_chromosomes)
            self.colors = [self.cmap(i) for i in range(0, self.max_chromosomes)]

        self.background = background

    def init_coordinates(self):
        self.x_end = self.x_start + (2 + len(self.chromosome_df_dict) * 7) * self.x_size

    def draw(self, axes=None, style=None):

        used_style = style if style else self.style

        current_subplot = axes if axes else plt.gca()

        last_x_start = self.x_start + 2 * self.x_size
        square_y_pos_reference = (self.max_chromosomes + 3) * self.element_size # TODO: make coordinate calculation more clear
        #square_y_pos_reference = self.y_start + self.reference_chromosome_number * self.element_size
        for species in self.chromosome_df_dict:
            chromosome_number = len(self.chromosome_df_dict[species])
            square_y_pos = square_y_pos_reference
            current_subplot.annotate(species, xy=(last_x_start, square_y_pos + self.element_size * 0.25),
                                     xycoords='data',
                                     fontsize=self.fontsize,
                                     fontweight='bold',
                                     xytext=(last_x_start, square_y_pos + self.element_size * 0.25), )
            square_y_pos -= 2 * self.element_size
            for i in range(0, chromosome_number)[::-1]:
                square_y_pos -= self.element_size
                # print (colormap_tuple_list[i][1])

                fragment = Rectangle((last_x_start - 2 * self.x_size, square_y_pos), self.x_size, self.element_size,
                                     fill=True,
                                     edgecolor="black", facecolor=self.chromosome_df_dict[species]["color"].iloc[i],
                                     linewidth=0.5)

                current_subplot.add_patch(fragment)

                current_subplot.annotate(self.chromosome_df_dict[species].index[i],
                                         xy=(last_x_start, square_y_pos + self.element_size * 0.25), xycoords='data',
                                         fontsize=self.fontsize,
                                         xytext=(last_x_start, square_y_pos + self.element_size * 0.25), )

            last_x_start = last_x_start + 7 * self.x_size
