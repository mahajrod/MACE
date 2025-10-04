from MACE.Visualize.Styles.Label import default_x_axis_label_style, default_y_axis_label_style


class AxisStyle:

    def __init__(self, axis_width=0.01, major_tick_height=0.02, minor_tick_height=0.01,
                 major_tick_width=0.02, minor_tick_width=0.01, axis_label_position='bottom',
                 tick_label_style=None, tick_label_distance=0.02,
                 axis_label_style=None, axis_label_distance=0.02, orientation="horizontal",
                 show_axis=False, show_ticks=False, show_tick_labels=False, show_label=False,
                 zorder=100, zorder_shift=0):
        """

        :param axis_width:
        :param major_tick_height:
        :param minor_tick_height:
        :param major_tick_width:
        :param minor_tick_width:
        :param axis_label_position:
        :param label_style:
        :param label_distance:
        :param orientation:
        :param zorder:
        """

        self.axis_width = axis_width
        self.major_tick_height = major_tick_height
        self.minor_tick_height = minor_tick_height
        self.major_tick_width = major_tick_width
        self.minor_tick_width = minor_tick_width
        self.tick_label_distance = tick_label_distance
        self.tick_label_style = tick_label_style
        self.axis_label_distance = axis_label_distance
        self.axis_label_style = axis_label_style
        self.axis_label_distance = axis_label_distance
        self.axis_label_position = axis_label_position
        self.orientation = orientation
        self.show_axis = show_axis
        self.show_ticks = show_ticks
        self.show_tick_labels = show_tick_labels
        self.show_label = show_label
        self.zorder = zorder
        self.zorder_shift = zorder_shift


default_x_axis_style = AxisStyle(orientation="horizontal", axis_label_style=default_x_axis_label_style,
                                 show_axis=False, show_ticks=False, show_label=False)

default_y_axis_style = AxisStyle(orientation="vertical", axis_label_style=default_y_axis_label_style,
                                 show_axis=False, show_ticks=False, show_label=False,)
