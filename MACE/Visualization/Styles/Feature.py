

class FeatureStyle():

    def __init__(self, patch_type, fill=True, edge_color=None, face_color=None, edge_width=0.0000001):
        self.patch_type = patch_type
        self.fill = fill
        self.edge_color = edge_color
        self.face_color = face_color
        self.edge_width = edge_width


default_feature_style = FeatureStyle(patch_type="rectangle")
