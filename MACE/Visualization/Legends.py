
from MACE.Visualization.Styles.Legend import default_legend_style


from collections import Iterable, OrderedDict

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


class DensityLegend:

    def __init__(self, y_start=0, x_start=0, x_size=10, element_size=10, style=default_legend_style,
                 colormap=None, thresholds=np.array((0.0, 0.1, 0.25, 0.5, 1.0)),
                 colors=("#333a97", "green", "yellow", "orange", "red"), background="white",
                 masked="grey", fontsize=13):
        self.style = style
        self.y_start = y_start
        self.x_start = x_start

        self.element_size = element_size
        self.x_size = x_size

        self.colormap = colormap
        self.colors = colors
        self.thresholds = thresholds

        if colormap:
            self.cmap = plt.get_cmap(self.colormap, len(self.thresholds))
            self.colors = [self.cmap(i) for i in range(0, len(thresholds))]

        self.background = background
        self.masked = masked

        self.fontsize = fontsize

    def draw(self, axes=None, style=None):

        used_style = style if style else self.style

        current_subplot = axes if axes else plt.gca()

        square_y_pos = self.y_start - self.element_size

        for color, legend_label in zip((self.masked, self.background), ("masked", "no SNPs")):
            square_y_pos += self.element_size
            print (self.x_start, square_y_pos), self.x_size, self.element_size, color
            fragment = Rectangle((self.x_start - 2 * self.x_size, square_y_pos), self.x_size, self.element_size,
                                 fill=True,
                                 edgecolor="black", facecolor=color, linewidth=0.5)

            current_subplot.add_patch(fragment)

            current_subplot.annotate(legend_label,
                                     xy=(self.x_start, square_y_pos + self.element_size*0.25), xycoords='data',
                                     fontsize=self.fontsize,
                                     xytext=(self.x_start, square_y_pos+ self.element_size*0.25), )

        for i in range(0, len(self.thresholds)):
            square_y_pos += self.element_size
            # print (colormap_tuple_list[i][1])
            fragment = Rectangle((self.x_start - 2 * self.x_size, square_y_pos), self.x_size, self.element_size,
                                 fill=True,
                                 edgecolor="black", facecolor=self.colors[i], linewidth=0.5)

            current_subplot.add_patch(fragment)
            if i == (len(self.thresholds) - 1):
                legend_element_label = "> %.2f" % self.thresholds[i]
            else:
                legend_element_label = "%.2f - %.2f" % (self.thresholds[i], self.thresholds[i + 1])

            current_subplot.annotate(legend_element_label,
                                     xy=(self.x_start, square_y_pos+ self.element_size*0.25), xycoords='data',
                                     fontsize=self.fontsize,
                                     xytext=(self.x_start, square_y_pos+ self.element_size*0.25), )


