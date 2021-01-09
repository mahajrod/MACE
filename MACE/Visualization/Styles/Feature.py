

class FeatureStyle:

    def __init__(self, patch_type, height, fill=True, edge=True, edge_color=None, face_color=None, edge_width=0.0000001,
                 show_label=True, label_fontsize=16, label_hor_aln='right', label_vert_aln='center',
                 label_y_shift=None, x_scale_factor=1,
                 y_scale_factor=1):
        self.patch_type = patch_type
        self.height = height

        self.fill = fill
        self.edge = edge
        self.edge_color = edge_color
        self.face_color = face_color
        self.edge_width = edge_width

        self.show_label = show_label
        self.label_fontsize = label_fontsize
        self.label_hor_aln = label_hor_aln
        self.label_vert_aln = label_vert_aln

        self.label_y_shift = label_y_shift if label_y_shift else height / 2

        self.x_scale_factor = x_scale_factor
        self.y_scale_factor = y_scale_factor


default_feature_style = FeatureStyle(patch_type="rectangle", height=9,)
circle_feature_style = FeatureStyle(patch_type="circle", height=5,)
ellipse_feature_style = FeatureStyle(patch_type="ellipse", height=5,)
