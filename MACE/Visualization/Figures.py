
from MACE.Visualization.Styles.Figure import FigureStyle, default_figure_style
import math
from collections import Iterable, OrderedDict

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


class Figure(OrderedDict):

    def __init__(self, subplots=None, figure_id=1, width=None, height=None,
                 subplot_number=None,
                 horizontal_subplot_number=None,
                 vertical_subplot_number=None,
                 suptitle=None,
                 style=default_figure_style,
                 dpi=None):
        if subplots:
            OrderedDict.__init__(self, subplots)
        else:
            OrderedDict.__init__(self)

        self.id = figure_id
        self.style = style

        self.suptitle = suptitle
        self.figure = None
        self.axes_array = None

        self.subplot_number = subplot_number if subplot_number else None

        if subplots:
            self.subplot_number = len(subplots)

        self.width = width if width else self.style.width
        self.height = height if height else self.style.height

        self.horizontal_subplot_number = horizontal_subplot_number
        self.vertical_subplot_number = vertical_subplot_number

        if subplots:
            if self.horizontal_subplot_number and self.vertical_subplot_number:
                pass
            elif self.horizontal_subplot_number:
                ratio = self.subplot_number // self.horizontal_subplot_number
                redundant = self.subplot_number % self.horizontal_subplot_number
                if ratio < 1:
                    self.horizontal_subplot_number = self.subplot_number
                    self.vertical_subplot_number = 1
                elif redundant == 0:
                    self.vertical_subplot_number = ratio
                else:
                    self.vertical_subplot_number = ratio + 1
            elif self.vertical_subplot_number:
                ratio = self.subplot_number // self.vertical_subplot_number
                redundant = self.subplot_number % self.vertical_subplot_number
                if ratio < 1:
                    self.vertical_subplot_number = self.subplot_number
                    self.horizontal_subplot_number = 1
                elif redundant == 0:
                    self.horizontal_subplot_number = ratio
                else:
                    self.horizontal_subplot_number = ratio + 1
            else:
                sqr = int(math.sqrt(self.subplot_number))
                if self.subplot_number == sqr:
                    self.horizontal_subplot_number = sqr
                    self.vertical_subplot_number = sqr
                elif self.subplot_number <= sqr * (sqr + 1):
                    self.horizontal_subplot_number = sqr + 1
                    self.vertical_subplot_number = sqr
                else:
                    self.horizontal_subplot_number = sqr + 1
                    self.vertical_subplot_number = sqr + 1

        if dpi:
            self.style.dpi = dpi

        if horizontal_subplot_number and vertical_subplot_number:
            self.axes_number = self.horizontal_subplot_number * self.vertical_subplot_number

        if self.style.width_per_subplot and horizontal_subplot_number:

            self.width = self.style.width_per_subplot * horizontal_subplot_number

        if self.style.height_per_subplot and vertical_subplot_number:
            self.height = self.style.height_per_subplot * vertical_subplot_number

    def draw(self):

        if self.horizontal_subplot_number and self.vertical_subplot_number:
            self.figure, self.axes_array = plt.subplots(self.vertical_subplot_number, self.horizontal_subplot_number,
                                                        sharex=self.style.share_x_axis, sharey=self.style.share_y_axis,
                                                        num=self.id,
                                                        dpi=self.style.dpi,
                                                        figsize=(self.horizontal_subplot_number * self.style.width_per_subplot,
                                                                 self.vertical_subplot_number * self.style.height_per_subplot),
                                                        squeeze=False)
        else:
            self.figure = plt.figure(self.id, dpi=self.style.dpi,
                                     figsize=(self.width, self.height))

        if self:
            for (subplot, subplot_index) in zip(self.keys(), range(0, self.subplot_number)):
                if self[subplot]:
                    subplot_hor_index = subplot_index % self.horizontal_subplot_number
                    subplot_vert_index = subplot_index // self.horizontal_subplot_number
                    self[subplot].draw(axes=self.axes_array[subplot_hor_index][subplot_vert_index])

        if self.suptitle:
            plt.suptitle(self.suptitle)

