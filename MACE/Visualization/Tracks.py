
from MACE.Visualization.Styles.Track import default_track_style
from MACE.Visualization.Styles.Feature import default_feature_style
import math
from collections import Iterable, OrderedDict

import numpy as np

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle, Ellipse, Polygon
from matplotlib.lines import Line2D


class Track:

    def __init__(self, records, style, y_start=None, x_start=0, x_end=None, label=None,
                 feature_style=default_feature_style, color_expression=None, colormap=None, thresholds=None,
                 colors=None, background=None, masked=None, patch_function=None,
                 forward_patch_function=None, reverse_patch_function=None, subplot_scale=False,
                 track_group_scale=False, figure_x_y_ratio=None, subplot_x_y_ratio=None, track_group_x_y_ratio=None,
                 stranded=False, rounded=False, middle_break=False, centromere_start=None, centromere_end=None,
                 track_group_highlight=False, track_group_highlight_color=None, highlight=False):

        self.records = records

        self.style = style
        self.feature_style = feature_style

        self.y_start = y_start
        self.y_end = None

        self.x_start = x_start
        self.x_end = x_end

        self.label = label

        if color_expression:
            self.style.color_expression = color_expression
        if thresholds:
            self.style.thresholds = thresholds
        if colors:
            self.style.colors = colors
        if background:
            self.style.background = background
        if masked:
            self.style.masked = masked
        if colormap:
            self.style.colormap = colormap
            if thresholds:
                self.style.cmap = plt.get_cmap(self.style.colormap, len(self.style.thresholds))
                self.style.colors = [self.style.cmap(i) for i in range(0, len(self.style.thresholds))]

        self.patch_function = patch_function

        self.forward_patch_function = forward_patch_function
        self.reverse_patch_function = reverse_patch_function

        self.figure_x_y_ratio = figure_x_y_ratio
        self.subplot_x_y_ratio = subplot_x_y_ratio
        self.track_group_x_y_ratio = track_group_x_y_ratio

        self.subplot_scale = subplot_scale
        self.track_group_scale = track_group_scale
        self.stranded = stranded
        self.rounded = rounded
        self.middle_break = middle_break

        self.x_scale_factor = None

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

        self.arc_angles_dict = {}
        self.arc_center_dict = {}
        self.x_radius_dict = {}
        self.y_radius_dict = {}
        self.arc_point_dict = {}

        self.point_array = None
        self.masking_point_array_dict = {}
        self.masking_patch_dict = {}
        self.track_background_patch = None
        self.track_border_patch = None
        self.centromere_start = centromere_start
        self.centromere_end = centromere_end

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
        self.show_centromere = None

        self.middle_break_left_top_point = None
        self.middle_break_right_top_point = None
        self.middle_break_right_bottom_point = None
        self.middle_break_left_bottom_point = None

        self.middle_break_array = None
        self.middle_break_patch = None
        self.middle_break_left_line = None
        self.middle_break_right_line = None

        self.track_group_highlight = track_group_highlight
        self.track_group_highlight_color = track_group_highlight_color
        self.highlight = highlight

    #def set_color(self):
    #    self.records["color"] = self.records["value"].apply()
    def create_patch_function(self, style=None, feature_style=None, y_start=None,  *args, **kwargs):
        """
        Track type dependent function
        :param track_group_scale:
        :param subplot_scale:
        :param style:
        :param feature_style:
        :param y_start:
        :param args:
        :param kwargs:
        :return: (forward_patch_function, reverse_patch_function) or (patch_function, None)
        """
        raise ValueError("ERROR!!! Call for ancestor abstract class method, i.e this method "
                         "was not implemented in called descendant class")
        pass

    def create_patch_collection(self, y_start=None, style=None, feature_style=None,
                                track_xmin=None, track_xmax=None, patch_function=None,
                                forward_patch_function=None, reverse_patch_function=None,
                                strand_column="strand", *args, **kwargs):
        if self.records is None:
            return 0
        # TODO: add min and max coordinate filtering
        if self.style is None and style is None:
            raise ValueError("ERROR!!! No style was set!")
        used_style = style if style else self.style
        used_feature_style = feature_style if feature_style else self.feature_style

        y_track = y_start if y_start else self.y_start

        # TODO: simplify following code, self.path function is already tuple of two functions or Nones.
        if not self.stranded:
            if patch_function:
                self.patch_function = patch_function
            else:
                if not self.patch_function:
                    self.patch_function = self.create_patch_function(style=used_style,
                                                                     feature_style=used_feature_style,
                                                                     y_start=y_track,
                                                                     *args, **kwargs)[0]
            return [PatchCollection(self.records.apply(self.patch_function, axis=1),
                                    match_original=True,
                                    antialiased=False,
                                    zorder=used_style.zorder["element"]), None]
        else:
            if (forward_patch_function is not None) and (reverse_patch_function is not None):
                self.forward_patch_function = forward_patch_function
                self.reverse_patch_function = reverse_patch_function
            else:
                self.forward_patch_function, self.reverse_patch_function = self.create_patch_function(style=used_style,
                                                                                                      feature_style=used_feature_style,
                                                                                                      y_start=y_track,
                                                                                                      *args, **kwargs)

            if strand_column in self.records.columns:
                #print(self.records)
                forward_patch_collection = PatchCollection(self.records[self.records[strand_column] == "+"].apply(self.forward_patch_function, axis=1),
                                                           antialiased=False,
                                                           match_original=True,
                                                           zorder=used_style.zorder["element"]) if not self.records[self.records[strand_column] == "+"].empty else None
                reverse_patch_collection = PatchCollection(self.records[self.records[strand_column] == "-"].apply(self.reverse_patch_function, axis=1),
                                                           antialiased=False,
                                                           match_original=True,
                                                           zorder=used_style.zorder["element"]) if not self.records[self.records[strand_column] == "-"].empty else None
            else:
                forward_patch_collection = PatchCollection(self.records.apply(self.forward_patch_function, axis=1),
                                                           antialiased=False,
                                                           match_original=True,
                                                           zorder=used_style.zorder["element"]) if not self.records.empty else None
                reverse_patch_collection = None
            return forward_patch_collection, reverse_patch_collection

    def draw(self, axes=None, style=None):
        """Element drawing was separated from border drawing because of overlap of masking pathes and borders
        resulting in visual bugs. Now all elements are drawn first, and all borders are added only after that"""

        used_style = style if style else self.style

        current_subplot = axes if axes else plt.gca()

        if self.highlight and (used_style.highlight_color is not None):
            highlight_patch = Polygon(np.array([[self.x_start, self.y_start],
                                                [self.x_start, self.y_end],
                                                [self.x_end, self.y_end],
                                                [self.x_end, self.y_start]]),
                                      # color=used_style.empty_color if (
                                      #        used_style.fill_empty and self.records is None) else used_style.face_color,
                                      fill=True,
                                      edgecolor=used_style.highlight_edge_color,
                                      facecolor=used_style.highlight_color,
                                      linewidth=used_style.highlight_edge_width,
                                      zorder=used_style.zorder['highlight'])
            masking_color = used_style.highlight_color
            axes.add_patch(highlight_patch)
        elif self.track_group_highlight and (self.track_group_highlight_color is not None):
            highlight_patch = Polygon(np.array([[self.x_start, self.y_start],
                                                [self.x_start, self.y_end],
                                                [self.x_end, self.y_end],
                                                [self.x_end, self.y_start]]),
                                      # color=used_style.empty_color if (
                                      #        used_style.fill_empty and self.records is None) else used_style.face_color,
                                      fill=True,
                                      edgecolor=self.track_group_highlight_color,
                                      facecolor=self.track_group_highlight_color,
                                      linewidth=used_style.highlight_edge_width,
                                      zorder=used_style.zorder['highlight'])
            masking_color = self.track_group_highlight_color
            axes.add_patch(highlight_patch)
        else:
            masking_color = None


        # calculate coordinates for masking patches and border patch
        #print("AAAAAA")
        #print(self.subplot_x_y_ratio)
        #print(self.figure_x_y_ratio)
        self.x_scale_factor = self.subplot_x_y_ratio / (self.figure_x_y_ratio if self.figure_x_y_ratio is not None else 1)  # same scalling for all tracks

        # coordinates of outer track rectangle
        self.left_bottom_outer_point = np.array([self.x_start, self.y_start])
        self.left_top_outer_point = np.array([self.x_start, self.y_start + used_style.height])
        self.right_top_outer_point = np.array([self.x_end, self.y_start + used_style.height])
        self.right_bottom_outer_point = np.array([self.x_end, self.y_start])

        self.general_x_smooth_element_len = used_style.height / 2 * self.x_scale_factor
        self.left_x_smooth_element_len = self.general_x_smooth_element_len
        self.right_x_smooth_element_len = self.general_x_smooth_element_len
        self.centromere_x_smooth_element_len = self.general_x_smooth_element_len

        self.y_radius = float(used_style.height) / 2
        self.left_x_radius = self.left_x_smooth_element_len
        self.right_x_radius = self.right_x_smooth_element_len

        self.left_center_point = np.array([self.x_start + self.left_x_smooth_element_len, self.y_start + used_style.height / 2])  # (x, y)
        self.right_center_point = np.array([self.x_end - self.right_x_smooth_element_len, self.y_start + used_style.height / 2])

        self.left_middle_point = np.array([self.x_start, self.y_start + used_style.height / 2])
        self.left_top_point = np.array([self.x_start + self.left_x_smooth_element_len, self.y_start + used_style.height])
        self.left_bottom_point = np.array([self.x_start + self.left_x_smooth_element_len, self.y_start])

        self.right_middle_point = np.array([self.x_end, self.y_start + used_style.height / 2])
        self.right_top_point = np.array([self.x_end - self.right_x_smooth_element_len, self.y_start + used_style.height])
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

        if used_style.centromere and (self.centromere_start is not None) and (self.centromere_end is not None):
            centromere_middle = float(self.centromere_start + self.centromere_end) / 2

            self.centromere_middle_point = np.array([centromere_middle, self.y_start + used_style.height / 2])

            self.centromere_left_top_point = np.array([centromere_middle - self.centromere_x_smooth_element_len,
                                                       self.y_start + used_style.height])
            self.centromere_right_top_point = np.array([centromere_middle + self.centromere_x_smooth_element_len,
                                                        self.y_start + used_style.height])

            self.centromere_right_bottom_point = np.array([centromere_middle + self.centromere_x_smooth_element_len,
                                                           self.y_start])
            self.centromere_left_bottom_point = np.array([centromere_middle - self.centromere_x_smooth_element_len,
                                                          self.y_start])

            # verify and adjust centromere coordinates
            if self.left_right_overlap:
                self.show_centromere = False
            else:
                self.show_centromere = True
                # check overlaps with centromere
                self.left_centromere_middle_overlap = True if centromere_middle < self.left_top_point[0] else False
                self.right_centromere_middle_overlap = True if centromere_middle > self.right_top_point[0] else False
                self.left_centromere_overlap = True if self.centromere_left_top_point[0] < self.left_top_point[0] else False
                self.right_centromere_overlap = True if self.centromere_right_top_point[0] > self.right_top_point[0] else False

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
                    self.right_x_radius = self.x_end - (self.centromere_right_top_point[0] + self.right_top_point[0]) / 2
                    self.right_top_point[0] = self.right_middle_point[0] - self.right_x_radius
                    self.right_bottom_point[0] = self.right_top_point[0]
                    self.right_center_point[0] = self.right_top_point[0]

                    self.centromere_right_top_point[0] = self.right_top_point[0]
                    self.centromere_right_bottom_point[0] = self.right_top_point[0]

        self.arc_angles_dict = {"left_bottom": np.linspace(1.5 * np.pi, np.pi, used_style.arc_point_number),
                                "left_top": np.linspace(np.pi, np.pi / 2, used_style.arc_point_number),
                                "right_top": np.linspace(np.pi / 2, 0, used_style.arc_point_number),
                                "right_bottom": np.linspace(2 * np.pi, 1.5 * np.pi, used_style.arc_point_number),
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
            self.arc_point_dict[arc] = np.column_stack([self.x_radius_dict[arc] * np.cos(self.arc_angles_dict[arc]) + self.arc_center_dict[arc][0],
                                                        self.y_radius_dict[arc] * np.sin(self.arc_angles_dict[arc]) + self.arc_center_dict[arc][1]])

        self.masking_point_array_dict = {}
        # print (self.x_start, self.y_start, self.x_end)
        if used_style.stranded and used_style.rounded:
            if used_style.stranded_end:
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
        elif used_style.rounded:
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
        elif used_style.stranded:
            if used_style.stranded_end:
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
            """
            self.track_patch = Rectangle((self.x_start, self.y_start), self.x_end - self.x_start,
                                         used_style.height,
                                         color=used_style.empty_color if (
                                                 used_style.fill_empty and self.records is None) else used_style.face_color,
                                         fill=True if (
                                                     used_style.fill_empty and self.records is None) else used_style.fill,
                                         edgecolor=used_style.edge_color,
                                         facecolor=used_style.face_color,
                                         linewidth=used_style.edge_width)
            """
        top_middle_point_list = []
        bottom_middle_point_list = []

        if self.show_centromere:
            # do not draw centromere if rounding points from left and right overlap
            if used_style.centromere and (self.centromere_start is not None) and (self.centromere_end is not None):
                top_middle_point_list = [[self.centromere_left_top_point],
                                         [self.centromere_middle_point],
                                         [self.centromere_right_top_point]]
                bottom_middle_point_list = [[self.centromere_right_bottom_point],
                                            [self.centromere_middle_point],
                                            [self.centromere_left_bottom_point]]
                self.masking_point_array_dict["top_centromere"] = np.concatenate(top_middle_point_list)
                self.masking_point_array_dict["bottom_centromere"] = np.concatenate(bottom_middle_point_list)

        if self.middle_break:
            x_middle = (self.x_end + self.x_start) / 2
            self.middle_break_left_top_point = np.array([x_middle, self.y_start + used_style.height * (1 + used_style.middle_break_y_overhang)])
            self.middle_break_right_top_point = np.array([x_middle + 2 * self.centromere_x_smooth_element_len,
                                                          self.y_start + used_style.height * (1 + used_style.middle_break_y_overhang)])

            self.middle_break_right_bottom_point = np.array([x_middle, self.y_start - used_style.height * used_style.middle_break_y_overhang])
            self.middle_break_left_bottom_point = np.array([x_middle - 2 * self.centromere_x_smooth_element_len,
                                                            self.y_start - used_style.height * used_style.middle_break_y_overhang])

            self.middle_break_array = np.concatenate([[self.middle_break_left_top_point],
                                                      [self.middle_break_right_top_point],
                                                      [self.middle_break_right_bottom_point],
                                                      [self.middle_break_left_bottom_point]])
            #print(self.middle_break_array )
            self.middle_break_patch = Polygon(self.middle_break_array,
                                              #color=used_style.empty_color if (
                                              #        used_style.fill_empty and self.records is None) else used_style.background,
                                              fill=True,
                                              edgecolor=used_style.background,
                                              facecolor=used_style.background,
                                              linewidth=used_style.edge_width,
                                              zorder=used_style.zorder["middle_break"])

            self.middle_break_left_line = Line2D((self.middle_break_left_bottom_point[0], self.middle_break_left_top_point[0]),
                                                 (self.middle_break_left_bottom_point[1], self.middle_break_left_top_point[1]),
                                                 color=used_style.empty_color if self.records is None else used_style.edge_color,
                                                 linewidth=used_style.edge_width,
                                                 zorder=used_style.zorder["middle_break"])

            self.middle_break_right_line = Line2D((self.middle_break_right_bottom_point[0], self.middle_break_right_top_point[0]),
                                                 (self.middle_break_right_bottom_point[1], self.middle_break_right_top_point[1]),
                                                 color=used_style.empty_color if self.records is None else used_style.edge_color,
                                                 linewidth=used_style.edge_width,
                                                 zorder=used_style.zorder["middle_break"])
            current_subplot.add_patch(self.middle_break_patch)
            current_subplot.add_line(self.middle_break_left_line)
            current_subplot.add_line(self.middle_break_right_line)

        self.point_array = np.concatenate(left_point_list + top_middle_point_list + right_point_list + bottom_middle_point_list)
        self.track_background_patch = Polygon(self.point_array,
                                              #color=used_style.empty_color if (
                                              #        used_style.fill_empty and self.records is None) else used_style.background,
                                              fill=True,
                                              edgecolor=used_style.empty_color if (
                                                      used_style.fill_empty and self.records is None) else used_style.background,
                                              facecolor=used_style.empty_color if (
                                                      used_style.fill_empty and self.records is None) else used_style.background,
                                              linewidth=used_style.edge_width,
                                              zorder=used_style.zorder["background"])
        self.track_border_patch = Polygon(self.point_array,
                                          #color=used_style.empty_color if (
                                          #        used_style.fill_empty and self.records is None) else used_style.face_color,
                                          fill=True if (
                                                  used_style.fill_empty and self.records is None) else used_style.fill,
                                          edgecolor=used_style.empty_color if (
                                                  used_style.fill_empty and self.records is None) else used_style.edge_color,
                                          facecolor=used_style.empty_color if (
                                                  used_style.fill_empty and self.records is None) else used_style.face_color,
                                          linewidth=used_style.edge_width,
                                          zorder=used_style.zorder["border"])

        self.masking_patch_dict = {masking_path: Polygon(self.masking_point_array_dict[masking_path],
                                                         # color=used_style.background,
                                                         fill=True,
                                                         edgecolor=masking_color if masking_color else used_style.background,
                                                         facecolor=masking_color if masking_color else used_style.background,
                                                         linewidth=used_style.edge_width,
                                                         zorder=used_style.zorder["masking_patches"]) for masking_path in self.masking_point_array_dict}

        # add middle (strand) line if necessary
        if used_style.stranded:
            current_subplot.add_line(Line2D((self.left_middle_point[0], self.right_middle_point[0]),
                                            (self.left_middle_point[1], self.right_middle_point[1]),
                                            color=used_style.middle_line_color,
                                            linewidth=used_style.middle_line_width,
                                            zorder=used_style.zorder["strand_line"]))

        current_subplot.add_patch(self.track_background_patch)
        # add features\windows\etc first
        if self.records is not None:
            for collection in self.create_patch_collection():
                if collection is not None:
                    current_subplot.add_collection(collection)
        # add masking patches:
        for patch in self.masking_patch_dict:
            current_subplot.add_patch(self.masking_patch_dict[patch])

        if self.label and self.style.show_label:
            current_subplot.annotate(self.label, xy=(0, self.y_start + self.style.height/2.5), xycoords='data',
                                     fontsize=self.style.label_fontsize,
                                     xytext=(-15, 0), textcoords='offset points',
                                     ha=self.style.label_hor_aln, va=self.style.label_vert_aln)

        if used_style.edge:
            # add track
            current_subplot.add_patch(self.track_border_patch)

    #def draw_borders(self, axes=None, style=None):
    #    """Border drawing was separated from element drawing because of overlap of masking pathes and borders
    #    resulting in visual bugs. Now all elements are drawn first, and all borders are added only after that"""
    #    used_style = style if style else self.style
    #
    #    current_subplot = axes if axes else plt.gca()
    #
    #    if used_style.edge:
    #        # add track
    #        current_subplot.add_patch(self.track_patch)


class WindowTrack(Track):

    def __init__(self, windows_df, window_size, window_step, y_start=None, x_start=0, x_end=None,
                 style=default_track_style, label=None, norm=False,
                 window_type="stacking", multiplier=None, feature_style=default_feature_style, color_expression=None,
                 colormap=None, thresholds=None,
                 colors=None, background=None, masked=None, subplot_scale=False,
                 track_group_scale=False, figure_x_y_ratio=1, subplot_x_y_ratio=None, track_group_x_y_ratio=None,
                 stranded=False, rounded=False, centromere_start=None, centromere_end=None):
        """

        :param windows_df:      Example
                                                All
                                CHROM	WINDOW
                                chr4	    0	58
                                            1	61
                                            2	54
        :param window_size:
        :param window_step:
        :param y_start:
        :param x_start:
        :param x_end:
        :param style:
        :param label:
        :param window_type:
        :param multiplier:
        :param feature_style:
        :param color_expression:
        :param colormap:
        :param thresholds:
        :param colors:
        :param background:
        :param masked:
        """

        Track.__init__(self, windows_df, style, y_start=y_start, x_start=x_start, x_end=x_end, label=label,
                       feature_style=feature_style, color_expression=color_expression,
                       colormap=colormap, thresholds=thresholds, colors=colors, background=background, masked=masked,
                       subplot_scale=subplot_scale, track_group_scale=track_group_scale,
                       figure_x_y_ratio=figure_x_y_ratio, subplot_x_y_ratio=subplot_x_y_ratio,
                       track_group_x_y_ratio=track_group_x_y_ratio,
                       stranded=stranded, rounded=rounded,
                       centromere_start=centromere_start, centromere_end=centromere_end
                       )

        self.track_type = "window"
        self.window_type = window_type

        count_columns = self.records.columns.to_list()
        if "masked" in count_columns:
            count_columns.remove("masked")

        if norm and multiplier:
            self.records["density"] = self.records.loc[:, count_columns] * (float(multiplier) / float(window_size))
        elif multiplier:
            self.records["density"] = self.records.loc[:, count_columns] * float(multiplier)
        elif norm:
            self.records["density"] = self.records.loc[:, count_columns] / float(window_size)
        else:
            self.records["density"] = self.records.loc[:, count_columns]

        self.window_size = window_size
        self.window_step = window_step
        self.multiplier = multiplier
        self.preprocess_data()

    def preprocess_data(self):
        window_index_name = self.records.index.names[0]
        if self.window_type == "stacking":
            self.records.reset_index(level=window_index_name, inplace=True)
            self.records["start"] = self.records[window_index_name] * self.window_step
            self.records = self.records[["start"] + list(self.records.columns[:-1])]

        elif self.window_type == "sliding":
            # TODO: implement sliding window
            # TODO: verify - sliding window drawing is already implemented  somewhere else:)))
            pass
        else:
            raise ValueError("ERROR!!! Unknown window type: %s" % self.window_type)

    def create_patch_function(self, style=None, feature_style=None, y_start=None, *args, **kwargs):

        if feature_style.patch_type == "rectangle":
            if "color" in self.records.columns:

                def create_patch(row, style=style if style else self.style,
                                 feature_style=feature_style if feature_style else self.feature_style, y_start=y_start,
                                 window_size=self.window_size):

                    return Rectangle((row['start'], y_start), window_size,
                                     style.height,
                                     fill=True,
                                     edgecolor=feature_style.edge_color,
                                     facecolor=row["color"],
                                     linewidth=feature_style.edge_width)

                return create_patch, None

            else:
                def create_patch(row, style=style, feature_style=feature_style, y_start=y_start,
                                 window_size=self.window_size):

                    return Rectangle((row['start'], y_start), window_size,
                                     style.height,
                                     fill=True,
                                     edgecolor=feature_style.edge_color,
                                     facecolor=feature_style.face_color,
                                     linewidth=feature_style.edge_width)

                return create_patch, None


class FeatureTrack(Track):

    def __init__(self, feature_df, y_start=None, x_start=0, x_end=None,
                 style=default_track_style, label=None,
                 feature_style=default_feature_style, color_expression=None,
                 colormap=None, thresholds=None,
                 colors=None, background=None, masked=None,
                 patch_function=None,
                 forward_patch_function=None, reverse_patch_function=None,
                 feature_start_column_id="start", 
                 feature_end_column_id="end", 
                 feature_color_column_id="color",
                 feature_length_column_id="length",
                 feature_strand_column_id="strand",
                 feature_value_column_id=None,
                 x_scale_factor=1,
                 y_scale_factor=1,
                 auto_scale=False,
                 subplot_scale=False,
                 track_group_scale=False,
                 stranded=False,
                 rounded=False,
                 middle_break=False,
                 figure_x_y_ratio=None, subplot_x_y_ratio=None, track_group_x_y_ratio=None,
                 centromere_start=None, centromere_end=None):

        Track.__init__(self, feature_df, style, y_start=y_start, x_start=x_start, x_end=x_end, label=label,
                       feature_style=feature_style, color_expression=color_expression,
                       colormap=colormap, thresholds=thresholds, colors=colors, background=background, masked=masked,
                       subplot_scale=subplot_scale,
                       patch_function=patch_function,
                       forward_patch_function=forward_patch_function, reverse_patch_function=reverse_patch_function,
                       track_group_scale=track_group_scale,
                       figure_x_y_ratio=figure_x_y_ratio, subplot_x_y_ratio=subplot_x_y_ratio,
                       track_group_x_y_ratio=track_group_x_y_ratio,
                       stranded=stranded, rounded=rounded, middle_break=middle_break,
                       centromere_start=centromere_start, centromere_end=centromere_end
                       )

        self.track_type = "feature"

        self.feature_start_column_id = feature_start_column_id
        self.feature_end_column_id = feature_end_column_id
        self.feature_color_column_id = feature_color_column_id
        self.feature_length_column_id = feature_length_column_id
        self.feature_strand_column_id = feature_strand_column_id
        self.feature_value_column_id = feature_value_column_id

        self.x_scale_factor = x_scale_factor
        self.y_scale_factor = y_scale_factor
        self.auto_scale = auto_scale
        self.preprocess_data()

    def preprocess_data(self):
        if self.records is not None:
            self.records[self.feature_length_column_id] = self.records[self.feature_end_column_id] - self.records[self.feature_start_column_id]

    def create_patch_function(self, style=None, feature_style=None, y_start=None, stranded=None, *args, **kwargs):
        if self.records is None:
            return 0

        stranded_feature = stranded if stranded is not None else self.stranded

        if feature_style.patch_type == "rectangle":
            if self.feature_color_column_id in self.records.columns:
                if stranded_feature:
                    def create_forward_patch(row, style=style if style else self.style,
                                     feature_style=feature_style if feature_style else self.feature_style, y_start=y_start):

                        return Rectangle((row[self.feature_start_column_id], y_start + (style.height / 2)),
                                         row[self.feature_length_column_id],
                                         style.height / 2,
                                         fill=True,
                                         edgecolor=None,
                                         facecolor=row[self.feature_color_column_id],
                                         linewidth=feature_style.edge_width)

                    def create_reverse_patch(row, style=style if style else self.style,
                                     feature_style=feature_style if feature_style else self.feature_style,
                                     y_start=y_start):

                        return Rectangle((row[self.feature_start_column_id], y_start),
                                         row[self.feature_length_column_id],
                                         style.height / 2,
                                         fill=True,
                                         edgecolor=feature_style.edge_color,
                                         facecolor=row[self.feature_color_column_id],
                                         linewidth=feature_style.edge_width)

                    return create_forward_patch, create_reverse_patch

                else:
                    def create_patch(row, style=style if style else self.style,
                                     feature_style=feature_style if feature_style else self.feature_style, y_start=y_start):

                        return Rectangle((row[self.feature_start_column_id], y_start),
                                         row[self.feature_length_column_id],
                                         style.height,
                                         fill=True,
                                         edgecolor=feature_style.edge_color,
                                         facecolor=row[self.feature_color_column_id],
                                         linewidth=feature_style.edge_width)

                    return create_patch, None

            else:
                if stranded_feature:
                    def create_forward_patch(row, style=style if style else self.style,
                                     feature_style=feature_style if feature_style else self.feature_style, y_start=y_start):

                        return Rectangle((row[self.feature_start_column_id], y_start + (style.height / 2)),
                                         row[self.feature_length_column_id],
                                         style.height / 2,
                                         fill=True,
                                         edgecolor=feature_style.edge_color,
                                         facecolor=feature_style.face_color,
                                         linewidth=feature_style.edge_width)

                    def create_reverse_patch(row, style=style if style else self.style,
                                     feature_style=feature_style if feature_style else self.feature_style,
                                     y_start=y_start):

                        return Rectangle((row[self.feature_start_column_id], y_start),
                                         row[self.feature_length_column_id],
                                         style.height / 2,
                                         fill=True,
                                         edgecolor=feature_style.edge_color,
                                         facecolor=feature_style.face_color,
                                         linewidth=feature_style.edge_width)

                    return create_forward_patch, create_reverse_patch

                else:
                    def create_patch(row, style=style, feature_style=feature_style, y_start=y_start):

                        return Rectangle((row[self.feature_start_column_id], y_start),
                                         row[self.feature_length_column_id],
                                         style.height,
                                         fill=True,
                                         edgecolor=feature_style.edge_color,
                                         facecolor=feature_style.face_color,
                                         linewidth=feature_style.edge_width)

                    return create_patch, None

        elif feature_style.patch_type == "circle":

            if self.feature_color_column_id in self.records.columns:

                def create_patch(row, style=style if style else self.style,
                                 feature_style=feature_style if feature_style else self.feature_style, y_start=y_start):

                    return Circle(((row[self.feature_end_column_id] + row[self.feature_start_column_id]) / 2 , y_start + style.height / 2), radius=feature_style.height / 2,
                                  fill=True,
                                  edgecolor=feature_style.edge_color,
                                  facecolor=row[self.feature_color_column_id],
                                  linewidth=feature_style.edge_width)

                return create_patch, None

            else:
                def create_patch(row, style=style, feature_style=feature_style, y_start=y_start):
                    return Circle(((row[self.feature_end_column_id] + row[self.feature_start_column_id]) / 2, y_start + style.height / 2),
                                  radius=feature_style.height / 2,
                                  fill=True,
                                  edgecolor=feature_style.edge_color,
                                  facecolor=feature_style.face_color,
                                  linewidth=feature_style.edge_width)
                
                return create_patch, None

        elif feature_style.patch_type == "ellipse":

            if self.subplot_scale:
                x_scale_factor = self.subplot_x_y_ratio / self.figure_x_y_ratio
            elif self.track_group_scale:
                x_scale_factor = self.track_group_x_y_ratio / self.figure_x_y_ratio
            else:
                x_scale_factor = 1
            y_scale_factor = 1

            if self.feature_color_column_id in self.records.columns:

                def create_patch(row, style=style if style else self.style,
                                 feature_style=feature_style if feature_style else self.feature_style, y_start=y_start):

                    return Ellipse(((row[self.feature_end_column_id] + row[self.feature_start_column_id]) / 2,
                                   y_start + style.height / 2),
                                   height=feature_style.height * y_scale_factor,
                                   width=feature_style.height * x_scale_factor if self.subplot_scale else (row[self.feature_end_column_id] - row[self.feature_start_column_id] + 1) * x_scale_factor,
                                   fill=True,
                                   edgecolor=feature_style.edge_color,
                                   facecolor=row[self.feature_color_column_id],
                                   linewidth=feature_style.edge_width)

                return create_patch, None

            else:
                def create_patch(row, style=style, feature_style=feature_style, y_start=y_start):
                    return Ellipse(((row[self.feature_end_column_id] + row[self.feature_start_column_id]) / 2,
                                   y_start + style.height / 2),
                                   height=feature_style.height * y_scale_factor,
                                   width=feature_style.height * x_scale_factor if self.subplot_scale else (row[self.feature_end_column_id] - row[self.feature_start_column_id] + 1) * x_scale_factor,
                                   fill=True,
                                   edgecolor=feature_style.edge_color,
                                   facecolor=feature_style.face_color,
                                   linewidth=feature_style.edge_width)

                return create_patch, None

        elif feature_style.patch_type == "hist":

            def create_patch(row, style=style, feature_style=feature_style, y_start=y_start):

                bin_values = np.array(map(float, row[self.feature_value_column_id].split(",")))
                bin_heights = bin_values / sum(bin_values) * style.height
                bin_cumulative_bottoms = np.zeros(len(bin_values))
                bin_cumulative_bottoms[1:] = np.cumsum(bin_heights)[0:-1]
                bin_colors = [style.hist_colors[i] for i in range(0, len(bin_values))]

                return [Rectangle((row[self.feature_start_column_id], y_start + y_bottom),
                                  row[self.feature_length_column_id],
                                  style.height,
                                  fill=True,
                                  edgecolor=feature_style.edge_color,
                                  facecolor=feature_style.face_color,
                                  linewidth=0) for y_bottom, y_height, color in zip(bin_cumulative_bottoms,
                                                                                    bin_heights,
                                                                                    bin_colors)]

            return create_patch, None