from collections import OrderedDict


class RecordStyle(OrderedDict):

    def __init__(self, *args, **kwargs):
        OrderedDict.__init__(self, *args, **kwargs)


class AbstractStyle:
    def __init__(self, y_bottom_offset: float = 0.0, y_top_offset: float = 0.0,
                 default_color: str or None = None, default_edge_color: str or None = None,
                 stranded: bool or None = None,
                 default_edge_width: float or None = None,
                 zorder_shift=0):
        self.y_bottom_offset = y_bottom_offset
        self.y_top_offset = y_top_offset
        self.default_color = default_color
        self.default_edge_color = default_edge_color
        self.stranded = stranded
        self.default_edge_width = default_edge_width
        self.zorder_shift = zorder_shift


class WindowStyle(AbstractStyle):
    def __init__(self, y_bottom_offset: float = 0.0, y_top_offset: float = 0.0,
                 default_color: str or None = None, default_edge_color: str or None = None,
                 stranded: bool or None = None, default_edge_width=0.00001,
                 zorder_shift=0):
        AbstractStyle.__init__(self, y_bottom_offset=y_bottom_offset, y_top_offset=y_top_offset,
                               default_color=default_color, default_edge_color=default_edge_color,
                               stranded=stranded, default_edge_width=default_edge_width,
                               zorder_shift=zorder_shift)


class HistStyle:
    pass


class PlotStyle:
    pass


class MarkerStyle:
    pass


default_record_style = RecordStyle({
                                    "window": WindowStyle(),
                                    "hist": HistStyle(),
                                    "plot": PlotStyle(),
                                    "marker": MarkerStyle()
                                    })

#default_feature_style = FeatureStyle(patch_type="rectangle", height=9,)
#circle_feature_style = FeatureStyle(patch_type="circle", height=5,)
#ellipse_feature_style = FeatureStyle(patch_type="ellipse", height=5,)
