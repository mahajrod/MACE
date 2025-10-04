

class LegendStyle():

    def __init__(self, patch_type="rectangle", x_size_denominator=32, legend_df=None):

        self.patch_type = patch_type
        self.x_size_denominator = x_size_denominator
        self.legend_df = legend_df


default_legend_style = LegendStyle()
