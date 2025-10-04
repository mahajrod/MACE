from MACE.Visualize.Styles.AxisStyles import default_x_axis_style, default_y_axis_style
from MACE.Visualize.Styles.LabelStyles import default_track_label_style
from MACE.Visualize.Styles.TrackBodyStyles import default_track_body_style


class TrackStyle:

    def __init__(self, edge=False, background_color=None, edge_color="black", edge_width=None,
                 dimh_space_left=0,
                 dimh_y_axis_space=0.02,
                 dimh_space_right=0,
                 dimv_space_top=0.2,
                 dimv_label_space_bottom=0.2,
                 dimv_xaxis_space_top=0.1,
                 dimv_space_bottom=0.2,
                 body_style=None,
                 label_style=None,
                 x_axis_style=None,
                 y_axis_style=None,
                 zorder_shift=0,
                 show_label=True,
                 show_x_axis=False,
                 show_y_axis=False):
        # track
        self.dimh_space_left = dimh_space_left                        # Tsl
        self.dimh_y_axis_space = dimh_y_axis_space                    # Tys
        self.dimh_space_right = dimh_space_right                      # Tsr

        self.dimv_space_top = dimv_space_top                          # Tst
        self.dimv_label_space_bottom = dimv_label_space_bottom        # Tlsb
        self.dimv_xaxis_space_top = dimv_xaxis_space_top              # Txst
        self.dimv_space_bottom = dimv_space_bottom                    # Tsb

        self.edge = edge
        self.edge_color = edge_color
        self.edge_width = edge_width
        self.background_color = background_color

        # styles
        self.body_style = body_style
        self.label_style = label_style
        self.x_axis_style = x_axis_style
        self.y_axis_style = y_axis_style

        #self.middle_break = middle_break
        #self.middle_break_y_overhang = middle_break_y_overhang
        self.zorder_dict = {
                            "polygon": 10,
                           }
        self.zorder_shift = zorder_shift

        self.show_label = show_label
        self.show_x_axis = show_x_axis
        self.show_y_axis = show_y_axis


default_track_style = TrackStyle(body_style=default_track_body_style,
                                 label_style=default_track_label_style,
                                 x_axis_style=default_x_axis_style,
                                 y_axis_style=default_y_axis_style,
                                 show_label=True,
                                 show_x_axis=False,
                                 show_y_axis=False)
