from numpy import array

from matplotlib.pyplot import get_cmap


class TrackStyle:

    def __init__(self, height, edge=True, fill=False, face_color=None, edge_color="black", edge_width=None,
                 show_label=True, label_fontsize=16, label_hor_aln='right', label_vert_aln='center',
                 label_y_shift=None, colormap=None, thresholds=array((0.0, 0.1, 0.25, 0.5, 1.0)),
                 colors=("#333a97", "green", "yellow", "orange", "red"), background="white",
                 masked="grey", color_expression=None):
        self.height = height
        self.fill = fill
        self.face_color = face_color
        self.edge = edge
        self.edge_color = edge_color
        self.edge_width = edge_width

        self.show_label = show_label
        self.label_fontsize = label_fontsize
        self.label_hor_aln = label_hor_aln
        self.label_vert_aln = label_vert_aln

        self.label_y_shift = label_y_shift if label_y_shift else height / 2

        self.colormap = colormap
        self.color_expression = color_expression
        self.thresholds = thresholds
        self.colors = colors
        self.background = background
        self.masked = masked

        if colormap:
            self.cmap = get_cmap(self.colormap, len(self.thresholds))
            self.colors = [self.cmap(i) for i in range(0, len(thresholds))]


default_track_style = TrackStyle(height=10, colormap=None, thresholds=array((0.0, 0.1, 0.25, 0.5, 1.0)),
                                 colors=("#333a97", "green", "yellow", "orange", "red"), background="white",
                                 masked="grey")
