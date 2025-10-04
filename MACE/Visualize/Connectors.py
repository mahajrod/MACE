
from matplotlib.patches import Polygon
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import matplotlib.pyplot as plt
import numpy as np


class Connector(PathPatch):
    def __init__(self, top_start: tuple = (0, 1),
                 top_end: tuple = (1, 2),
                 bottom_start: tuple = (3, 0),
                 bottom_end: tuple = (4, +0.5),
                 x_fraction_parameter: float = 4,
                 y_fraction_parameter: float = 4,
                 y_shift: [None, float] = None,
                 edgecolor: str = "grey",
                 facecolor: str = "grey",
                 alpha: float = 0.3,
                 fill: bool = True,
                 zorder=None,
                 linewidth=0
                 ):
        self.top_start = top_start
        self.top_end = top_end
        self.bottom_start = bottom_start
        self.bottom_end = bottom_end
        self.x_fraction_parameter = x_fraction_parameter
        self.y_fraction_parameter = y_fraction_parameter
        self. y_shift = y_shift

        self.path_coordinates, self.path_commands = self.calculate_path()
        PathPatch.__init__(self,
                           Path(self.path_coordinates, self.path_commands),
                           facecolor=facecolor,
                           edgecolor=edgecolor,
                           linewidth=linewidth,
                           fill=fill,
                           alpha=alpha,
                           zorder=zorder)

    def calculate_path(self) -> (list, list):
        """
        Abstract method to be defined in each connector class
        :return:
        """
        
        return None, None


class CubicBezierConnector(Connector):

    def calculate_path(self):
        if self.y_shift is None:
            return (
                    [self.top_start,

                     (self.top_start[0] + (self.bottom_start[0] - self.top_start[0]) / self.x_fraction_parameter,
                      self.top_start[1] + (self.top_end[1] - self.top_start[1]) / self.y_fraction_parameter),
                     (self.bottom_start[0] - (self.bottom_start[0] - self.top_start[0]) / self.x_fraction_parameter,
                      self.bottom_start[1] + (self.bottom_end[1] - self.bottom_start[1]) / self.y_fraction_parameter),
                     self.bottom_start,

                     self.bottom_end,

                     (self.bottom_end[0] - (self.bottom_end[0] - self.top_end[0]) / self.x_fraction_parameter,
                      self.bottom_end[1] - (self.bottom_end[1] - self.bottom_start[1]) / self.y_fraction_parameter),
                     (self.top_end[0] + (self.bottom_end[0] - self.top_end[0]) / self.x_fraction_parameter,
                      self.top_end[1] - (self.top_end[1] - self.top_start[1]) / self.y_fraction_parameter),
                     self.top_end,

                     (0, 0), ],

                    [Path.MOVETO,
                     Path.CURVE4, Path.CURVE4, Path.CURVE4,
                     Path.LINETO,
                     Path.CURVE4, Path.CURVE4, Path.CURVE4,
                     Path.CLOSEPOLY, ]
                     )
        else:
            return (
                    [self.top_start,

                     (self.top_start[0],
                      self.top_start[1] - self.y_shift / self.y_fraction_parameter),
                     (self.bottom_start[0],
                      self.bottom_start[1] + self.y_shift / self.y_fraction_parameter),
                     self.bottom_start,

                     self.bottom_end,

                     (self.bottom_end[0],
                      self.bottom_end[1] + self.y_shift / self.y_fraction_parameter),
                     (self.top_end[0],
                      self.top_end[1] - self.y_shift / self.y_fraction_parameter),
                     self.top_end,

                     (0, 0), ],

                    [Path.MOVETO,
                     Path.CURVE4, Path.CURVE4, Path.CURVE4,
                     Path.LINETO,
                     Path.CURVE4, Path.CURVE4, Path.CURVE4,
                     Path.CLOSEPOLY, ]
                     )