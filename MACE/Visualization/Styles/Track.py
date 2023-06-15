from numpy import array

from matplotlib.pyplot import get_cmap


class TrackStyle:

    def __init__(self, height, edge=True, fill=False, face_color=None, edge_color="black", edge_width=None,
                 show_label=True, label_fontsize=16, label_fontstyle='normal', label_fontweight='normal',
                 label_hor_aln='right', label_vert_aln='center',
                 label_y_shift=0, colormap=None, hist_colormap="Set1",
                 hist_colormap_level_number=10,
                 thresholds=array((0.0, 0.1, 0.25, 0.5, 1.0)),
                 colors=("#333a97", "green", "yellow", "orange", "red"), background="white",
                 masked="grey", color_expression=None, fill_empty=False, empty_color="lightgrey",
                 middle_line_color="black", middle_line_width=0.5, stranded=False, rounded=False,
                 stranded_end=False,
                 centromere=False, arc_point_number=100,
                 highlight_color=None,
                 highlight_edge_color=None,
                 highlight_edge_width=None,
                 middle_break=False,
                 middle_break_y_overhang=0.2
                 ):
        self.height = height
        self.fill = fill
        self.fill_empty = fill_empty
        self.empty_color = empty_color
        self.face_color = face_color
        self.edge = edge
        self.edge_color = edge_color
        self.edge_width = edge_width

        self.show_label = show_label
        self.label_fontsize = label_fontsize
        self.label_fontstyle = label_fontstyle
        self.label_fontweight = label_fontweight
        self.label_hor_aln = label_hor_aln
        self.label_vert_aln = label_vert_aln

        self.label_y_shift = label_y_shift

        self.colormap = colormap

        self.color_expression = color_expression
        self.thresholds = thresholds
        self.colors = colors
        self.background = background
        self.masked = masked

        self.hist_colormap = hist_colormap
        self.hist_colors = get_cmap(hist_colormap, hist_colormap_level_number)

        self.stranded = stranded
        self.rounded = rounded
        self.stranded_end = stranded_end
        self.arc_point_number = arc_point_number
        self.centromere = centromere

        self.middle_line_color = middle_line_color
        self.middle_line_width = middle_line_width

        self.highlight_color = highlight_color
        self.highlight_edge_color = highlight_edge_color
        self.highlight_edge_width = highlight_edge_width

        self.middle_break = middle_break
        self.middle_break_y_overhang = middle_break_y_overhang
        self.zorder = {
                       'highlight': 20,
                       'background': 30,
                       'element': 50,
                       'masking_patches': 70,
                       'strand_line': 90,
                       'border': 100,
                       'middle_break': 110
                       }
        if colormap:
            self.cmap = get_cmap(self.colormap, len(self.thresholds))
            self.colors = [self.cmap(i) for i in range(0, len(thresholds))]


default_track_style = TrackStyle(height=10, colormap=None, thresholds=array((0.0, 0.1, 0.25, 0.5, 1.0)),
                                 colors=("#333a97", "green", "yellow", "orange", "red"), background="white",
                                 masked="grey", fill_empty=False)

feature_track_style = TrackStyle(height=10, colormap=None, background="white",
                                 masked="grey", fill_empty=True, empty_color="lightgrey")
