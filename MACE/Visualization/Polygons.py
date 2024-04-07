
from matplotlib.patches import Polygon

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()


class ChromosomePolygon(Polygon):
    def __init__(self, x_start: float, y_start: float, length: float, height: float,
                 stranded: bool = False, rounded: bool = False,
                 centromere_start: (None, float) = None, centromere_end: (None, float) = None,
                 show_centromere: (None, bool) = None,
                 arc_point_number: int = 100, x_scale_factor: float = 1,
                 edgecolor: str = "grey",
                 facecolor: str = "grey",
                 alpha: float = 0.3,
                 fill: bool = True,
                 zorder=None,
                 linewidth=None
                 ):
        self.x_start = x_start
        self.y_start = y_start
        self.length = length
        self.height = height
        self.centromere_start = centromere_start
        self.centromere_end = centromere_end

        #self.edgecolor = edgecolor
        #self.facecolor = facecolor
        #self.alpha = alpha
        #self.fill = fill
        #self.zorder = zorder

        self.y_end = self.y_start + height
        self.x_end = self.x_start + length
        self.centromere_x_start = (self.x_start + self.centromere_start) if self.centromere_start is not None else None
        self.centromere_x_end = (self.x_start + self.centromere_end) if self.centromere_start is not None else None
        
        self.stranded = stranded
        self.rounded = rounded
        self.show_centromere = show_centromere

        self.x_scale_factor = x_scale_factor

        self.general_x_smooth_element_len = None
        self.left_x_smooth_element_len = None
        self.right_x_smooth_element_len = None
        self.centromere_x_smooth_element_len = None

        self.y_radius = None
        self.left_x_radius = None
        self.right_x_radius = None

        self.left_center_point = None
        self.right_center_point = None

        self.left_bottom_outer_point = None
        self.left_top_outer_point = None
        self.right_top_outer_point = None
        self.right_bottom_outer_point = None

        self.left_middle_point = None
        self.left_top_point = None
        self.left_bottom_point = None

        self.right_middle_point = None
        self.right_top_point = None
        self.right_bottom_point = None

        self.arc_point_number = arc_point_number
        self.arc_angles_dict = {}
        self.arc_center_dict = {}
        self.x_radius_dict = {}
        self.y_radius_dict = {}
        self.arc_point_dict = {}

        self.centromere_middle_point = None
        self.centromere_left_top_point = None
        self.centromere_right_top_point = None

        self.centromere_right_bottom_point = None
        self.centromere_left_bottom_point = None

        self.left_right_overlap = None
        self.left_centromere_middle_overlap = None
        self.right_centromere_middle_overlap = None
        self.left_centromere_overlap = None
        self.right_centromere_overlap = None
        self.show_centromere = show_centromere

        self.point_array = self.init_point_array()
        
        Polygon.__init__(self, self.point_array, edgecolor=edgecolor, linewidth=linewidth, facecolor=facecolor,
                         fill=fill, alpha=alpha, zorder=zorder)
        
    def init_coordinates(self):
        """
        Abstract method to be defined in each custom polygon class
        :return:
        """
        pass


class LinearChromosome(ChromosomePolygon):

    def init_point_array(self):
        #print(self.show_centromere, self.centromere_start, self.centromere_end)
        #height = self.y_end - self.y_start
        # coordinates of outer track rectangle
        self.left_bottom_outer_point = np.array([self.x_start, self.y_start])
        self.left_top_outer_point = np.array([self.x_start, self.y_start + self.height])
        self.right_top_outer_point = np.array([self.x_end, self.y_start + self.height])
        self.right_bottom_outer_point = np.array([self.x_end, self.y_start])

        self.general_x_smooth_element_len = self.height / 2 * self.x_scale_factor
        self.left_x_smooth_element_len = self.general_x_smooth_element_len
        self.right_x_smooth_element_len = self.general_x_smooth_element_len
        self.centromere_x_smooth_element_len = self.general_x_smooth_element_len

        self.y_radius = float(self.height) / 2
        self.left_x_radius = self.left_x_smooth_element_len
        self.right_x_radius = self.right_x_smooth_element_len

        self.left_center_point = np.array([self.x_start + self.left_x_smooth_element_len, self.y_start + self.height / 2])  # (x, y)
        self.right_center_point = np.array([self.x_end - self.right_x_smooth_element_len, self.y_start + self.height / 2])

        self.left_middle_point = np.array([self.x_start, self.y_start + self.height / 2])
        self.left_top_point = np.array([self.x_start + self.left_x_smooth_element_len, self.y_start + self.height])
        self.left_bottom_point = np.array([self.x_start + self.left_x_smooth_element_len, self.y_start])

        self.right_middle_point = np.array([self.x_end, self.y_start + self.height / 2])
        self.right_top_point = np.array([self.x_end - self.right_x_smooth_element_len, self.y_start + self.height])
        self.right_bottom_point = np.array([self.x_end - self.right_x_smooth_element_len, self.y_start])

        # verify of left/right overlaps and adjust x coordinates
        if self.left_top_point[0] > self.right_top_point[0]:
            self.left_right_overlap = True
            self.left_top_point[0] = (self.left_top_point[0] + self.right_top_point[0]) / 2

            self.right_top_point[0] = self.left_top_point[0]
            self.right_bottom_point[0] = self.left_top_point[0]
            self.left_bottom_point[0] = self.left_top_point[0]
        else:
            self.left_right_overlap = False

        # print(self.right_top_point, self.right_middle_point, self.right_bottom_point)
        #print(self.show_centromere, self.centromere_x_start, self.centromere_x_start is not None, self.centromere_x_end, self.centromere_x_end is not None)
        if self.show_centromere and (self.centromere_x_start is not None) and (self.centromere_x_end is not None):
            #print("AAAA")
            centromere_middle = float(self.centromere_x_start + self.centromere_x_end) / 2

            self.centromere_middle_point = np.array([centromere_middle, self.y_start + self.height / 2])

            self.centromere_left_top_point = np.array([centromere_middle - self.centromere_x_smooth_element_len,
                                                       self.y_start + self.height])
            self.centromere_right_top_point = np.array([centromere_middle + self.centromere_x_smooth_element_len,
                                                        self.y_start + self.height])

            self.centromere_right_bottom_point = np.array([centromere_middle + self.centromere_x_smooth_element_len,
                                                           self.y_start])
            self.centromere_left_bottom_point = np.array([centromere_middle - self.centromere_x_smooth_element_len,
                                                          self.y_start])

            # verify and adjust centromere coordinates
            #if self.left_right_overlap:
            #   pass # self.show_centromere = False
            #else:

            # check overlaps with centromere
            self.left_centromere_middle_overlap = True if centromere_middle < self.left_top_point[0] else False
            self.right_centromere_middle_overlap = True if centromere_middle > self.right_top_point[0] else False
            self.left_centromere_overlap = True if self.centromere_left_top_point[0] < self.left_top_point[
                0] else False
            self.right_centromere_overlap = True if self.centromere_right_top_point[0] > self.right_top_point[
                0] else False

            if self.left_centromere_middle_overlap:
                self.left_x_radius = (centromere_middle - self.left_middle_point[0]) / 2
                self.left_top_point[0] = self.left_middle_point[0] + self.left_x_radius
                self.left_bottom_point[0] = self.left_top_point[0]
                self.left_center_point[0] = self.left_top_point[0]

                self.centromere_left_top_point[0] = self.left_top_point[0]
                self.centromere_left_bottom_point[0] = self.left_top_point[0]
            elif self.left_centromere_overlap:
                self.left_x_radius = (self.left_top_point[0] + self.centromere_left_top_point[0]) / 2 - self.x_start
                self.left_top_point[0] = self.left_middle_point[0] + self.left_x_radius
                self.left_bottom_point[0] = self.left_top_point[0]
                self.left_center_point[0] = self.left_top_point[0]

                self.centromere_left_top_point[0] = self.left_top_point[0]
                self.centromere_left_bottom_point[0] = self.left_top_point[0]

            if self.right_centromere_middle_overlap:
                self.right_x_radius = (self.right_middle_point[0] - centromere_middle) / 2
                self.right_top_point[0] = self.right_middle_point[0] - self.right_x_radius
                self.right_bottom_point[0] = self.right_top_point[0]
                self.right_center_point[0] = self.right_top_point[0]

                self.centromere_right_top_point[0] = self.right_top_point[0]
                self.centromere_right_bottom_point[0] = self.right_top_point[0]
            elif self.right_centromere_overlap:
                self.right_x_radius = self.x_end - (
                            self.centromere_right_top_point[0] + self.right_top_point[0]) / 2
                self.right_top_point[0] = self.right_middle_point[0] - self.right_x_radius
                self.right_bottom_point[0] = self.right_top_point[0]
                self.right_center_point[0] = self.right_top_point[0]

                self.centromere_right_top_point[0] = self.right_top_point[0]
                self.centromere_right_bottom_point[0] = self.right_top_point[0]

        # print(self.right_top_point, self.right_middle_point, self.right_bottom_point, self.x_end, self.right_x_radius)
        self.arc_angles_dict = {"left_bottom": np.linspace(1.5 * np.pi, np.pi, self.arc_point_number),
                                "left_top": np.linspace(np.pi, np.pi / 2, self.arc_point_number),
                                "right_top": np.linspace(np.pi / 2, 0, self.arc_point_number),
                                "right_bottom": np.linspace(2 * np.pi, 1.5 * np.pi, self.arc_point_number),
                                }
        self.arc_center_dict = {"left_bottom": self.left_center_point,
                                "left_top": self.left_center_point,
                                "right_top": self.right_center_point,
                                "right_bottom": self.right_center_point,
                                }
        self.x_radius_dict = {"left_bottom": self.left_x_radius,
                              "left_top": self.left_x_radius,
                              "right_top": self.right_x_radius,
                              "right_bottom": self.right_x_radius,
                              }
        self.y_radius_dict = {"left_bottom": self.y_radius,
                              "left_top": self.y_radius,
                              "right_top": self.y_radius,
                              "right_bottom": self.y_radius,
                              }

        self.arc_point_dict = {}

        for arc in self.arc_angles_dict:
            self.arc_point_dict[arc] = np.column_stack(
                [self.x_radius_dict[arc] * np.cos(self.arc_angles_dict[arc]) + self.arc_center_dict[arc][0],
                 self.y_radius_dict[arc] * np.sin(self.arc_angles_dict[arc]) + self.arc_center_dict[arc][1]])

        self.masking_point_array_dict = {}
        # print (self.x_start, self.y_start, self.x_end)
        if self.stranded and self.rounded:
            if self.stranded_end:
                left_point_list = [[self.left_bottom_point],
                                   [self.left_middle_point],
                                   self.arc_point_dict["left_top"],
                                   [self.left_top_point]
                                   ]
                right_point_list = [[self.right_top_point],
                                    [self.right_middle_point],
                                    self.arc_point_dict["right_bottom"],
                                    [self.right_bottom_point]
                                    ]

                self.masking_point_array_dict = {"left": np.concatenate([[self.left_bottom_point],
                                                                         [self.left_middle_point],
                                                                         self.arc_point_dict["left_top"],
                                                                         [self.left_top_point],
                                                                         [self.left_top_outer_point],
                                                                         [self.left_bottom_outer_point]
                                                                         ]),
                                                 "right": np.concatenate([[self.right_top_point],
                                                                          [self.right_middle_point],
                                                                          self.arc_point_dict["right_bottom"],
                                                                          [self.right_bottom_point],
                                                                          [self.right_bottom_outer_point],
                                                                          [self.right_top_outer_point]
                                                                          ])
                                                 }
            else:
                left_point_list = [[self.left_bottom_point],
                                   self.arc_point_dict["left_bottom"],
                                   [self.left_middle_point],
                                   self.arc_point_dict["left_top"],
                                   [self.left_top_point],
                                   ]
                right_point_list = [[self.right_top_point],
                                    self.arc_point_dict["right_top"],
                                    [self.right_middle_point],
                                    self.arc_point_dict["right_bottom"],
                                    [self.right_bottom_point]
                                    ]
                self.masking_point_array_dict = {"left": np.concatenate([[self.left_bottom_point],
                                                                         self.arc_point_dict["left_bottom"],
                                                                         [self.left_middle_point],
                                                                         self.arc_point_dict["left_top"],
                                                                         [self.left_top_point],
                                                                         [self.left_top_outer_point],
                                                                         [self.left_bottom_outer_point]
                                                                         ]),
                                                 "right": np.concatenate([[self.right_top_point],
                                                                          self.arc_point_dict["right_top"],
                                                                          [self.right_middle_point],
                                                                          self.arc_point_dict["right_bottom"],
                                                                          [self.right_bottom_point],
                                                                          [self.right_bottom_outer_point],
                                                                          [self.right_top_outer_point]
                                                                          ])
                                                 }
        elif self.rounded:
            left_point_list = [[self.left_bottom_point],
                               self.arc_point_dict["left_bottom"],
                               [self.left_middle_point],
                               self.arc_point_dict["left_top"],
                               [self.left_top_point]
                               ]
            right_point_list = [[self.right_top_point],
                                self.arc_point_dict["right_top"],
                                [self.right_middle_point],
                                self.arc_point_dict["right_bottom"],
                                [self.right_bottom_point]
                                ]
            self.masking_point_array_dict = {"left": np.concatenate([[self.left_bottom_point],
                                                                     self.arc_point_dict["left_bottom"],
                                                                     [self.left_middle_point],
                                                                     self.arc_point_dict["left_top"],
                                                                     [self.left_top_point],
                                                                     [self.left_top_outer_point],
                                                                     [self.left_bottom_outer_point]
                                                                     ]),
                                             "right": np.concatenate([[self.right_top_point],
                                                                      self.arc_point_dict["right_top"],
                                                                      [self.right_middle_point],
                                                                      self.arc_point_dict["right_bottom"],
                                                                      [self.right_bottom_point],
                                                                      [self.right_bottom_outer_point],
                                                                      [self.right_top_outer_point]
                                                                      ])
                                             }
        elif self.stranded:
            if self.stranded_end:
                left_point_list = [[self.left_bottom_point],
                                   [self.left_middle_point],
                                   [self.left_top_outer_point]
                                   ]
                right_point_list = [[self.right_top_point],
                                    [self.right_middle_point],
                                    [self.right_bottom_outer_point],
                                    ]
                self.masking_point_array_dict = {"left": np.concatenate([[self.left_bottom_point],
                                                                         [self.left_middle_point],
                                                                         [self.left_bottom_outer_point]
                                                                         ]),
                                                 "right": np.concatenate([[self.right_top_point],
                                                                          [self.right_middle_point],
                                                                          [self.right_top_outer_point]
                                                                          ])
                                                 }
            else:
                left_point_list = [[self.left_bottom_outer_point],
                                   [self.left_top_outer_point]
                                   ]
                right_point_list = [[self.right_top_outer_point],
                                    [self.right_bottom_outer_point]
                                    ]
        else:
            left_point_list = [[self.left_bottom_outer_point],
                               [self.left_top_outer_point]
                               ]
            right_point_list = [[self.right_top_outer_point],
                                [self.right_bottom_outer_point]
                                ]

            self.masking_point_array_dict = {}

        top_middle_point_list = []
        bottom_middle_point_list = []

        if self.show_centromere:
            # do not draw centromere if rounding points from left and right overlap
            if self.show_centromere and (self.centromere_x_start is not None) and (self.centromere_x_end is not None):
                top_middle_point_list = [[self.centromere_left_top_point],
                                         [self.centromere_middle_point],
                                         [self.centromere_right_top_point]]
                bottom_middle_point_list = [[self.centromere_right_bottom_point],
                                            [self.centromere_middle_point],
                                            [self.centromere_left_bottom_point]]
                self.masking_point_array_dict["top_centromere"] = np.concatenate(top_middle_point_list)
                self.masking_point_array_dict["bottom_centromere"] = np.concatenate(bottom_middle_point_list)

        point_array = np.concatenate(
            left_point_list + top_middle_point_list + right_point_list + bottom_middle_point_list)
        
        return point_array
