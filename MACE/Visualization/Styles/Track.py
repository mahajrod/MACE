

class TrackStyle():

    def __init__(self, height, edge=True, fill=False, face_color=None, edge_color="black", edge_width=None,
                 show_label=True, label_fontsize=16, label_hor_aln='right', label_vert_aln='top',
                 label_y_shift=None):
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


default_track_style = TrackStyle(height=10)
