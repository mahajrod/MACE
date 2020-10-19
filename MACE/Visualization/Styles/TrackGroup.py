

class TrackGroupStyle:

    def __init__(self, distance=5, internal_offset=2, x_multiplier=1.05, show_label=True,
                 label_fontsize=16, label_hor_aln='right', label_vert_aln='top',
                 label_y_shift=0,
                 label_x_shift=-15):
        self.distance = distance
        self.internal_offset = internal_offset
        self.x_multiplier = x_multiplier

        self.show_label = show_label
        self.label_fontsize = label_fontsize
        self.label_hor_aln = label_hor_aln
        self.label_vert_aln = label_vert_aln

        self.label_y_shift = label_y_shift
        self.label_x_shift = label_x_shift

default_track_group_style = TrackGroupStyle()

