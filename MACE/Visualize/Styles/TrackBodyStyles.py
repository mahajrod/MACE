
class TrackBodyStyle:

    def __init__(self, height=10,
                 edge=True, outer_color="white", outer_alpha=None, outer_edge_color=None, outer_edge_width=None,
                 background_color=None, background_alpha=None, edge_color="black", edge_alpha=None, edge_width=None,
                 highlight_empty=True, empty_color="lightgrey", empty_edge=False, empty_alpha=None,
                 middle_line_color="black", middle_line_width=0.5, show_middle_line=False, stranded=False, rounded=True,
                 stranded_end=False, show_centromere=True, arc_point_number=100,
                 zorder_shift=0,
                 ):
        # track polygon
        self.height = height
        self.highlight_empty = highlight_empty
        self.empty_color = empty_color
        self.empty_edge = empty_edge
        self.empty_alpha = empty_alpha
        self.outer_color = outer_color
        self.outer_alpha = outer_alpha
        self.outer_edge_color = outer_edge_color
        self.outer_edge_width = outer_edge_width
        self.background_color = background_color
        self.background_alpha = background_alpha
        self.edge = edge
        self.edge_color = edge_color
        self.edge_alpha = edge_alpha
        self.edge_width = edge_width

        # track body parameters
        self.stranded = stranded
        self.rounded = rounded
        self.stranded_end = stranded_end
        self.arc_point_number = arc_point_number
        self.show_centromere = show_centromere
        self.middle_line_color = middle_line_color
        self.middle_line_width = middle_line_width
        self.show_middle_line = show_middle_line
        if self.stranded:
            self.show_middle_line = True

        #self.middle_break = middle_break
        #self.middle_break_y_overhang = middle_break_y_overhang
        self.zorder_dict = {
                            'outer': 20,
                            'background': 30,
                            'records': 50,
                            'masking_patches': 70,
                            'middle_line': 90,
                            'border': 100,
                            'breaks': 110
                           }
        self.zorder_shift = zorder_shift
        for element in self.zorder_dict:
            self.zorder_dict[element] += self.zorder_shift


default_track_body_style = TrackBodyStyle(height=10, background_color="white", highlight_empty=True,
                                          empty_color="lightgrey", empty_edge=False)

