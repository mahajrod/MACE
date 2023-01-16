

class TrackGroupStyle:

    def __init__(self, distance=2, internal_offset=2, x_multiplier=1.05, show_label=True,
                 label_fontsize=16, label_fontweight="normal",
                 label_fontstyle='normal',
                 label_hor_aln='right', label_vert_aln='center',
                 label_y_shift=0,
                 label_x_shift=-15,
                 highlight_color=None,
                 highlight_edge_color=None,
                 highlight_edge_width=None
                 ):
        self.distance = distance
        self.internal_offset = internal_offset
        self.x_multiplier = x_multiplier

        self.show_label = show_label
        self.label_fontweight = label_fontweight
        self.label_fontstyle = label_fontstyle
        self.label_fontsize = label_fontsize
        self.label_hor_aln = label_hor_aln
        self.label_vert_aln = label_vert_aln

        self.label_y_shift = label_y_shift
        self.label_x_shift = label_x_shift

        self.highlight_color = highlight_color
        self.highlight_edge_color = highlight_edge_color
        self.highlight_edge_width = highlight_edge_width

        self.zorder = {
                       'highlight': 10,
                       'background': 20,
                       }


default_track_group_style = TrackGroupStyle()
italic_track_group_style = TrackGroupStyle(label_fontstyle='italic')
