import matplotlib.pyplot as plt


class FigureStyle:

    def __init__(self, dpi=300, width=None,
                 height=None, width_per_subplot=None,
                 height_per_subplot=None,
                 share_x_axis=False, share_y_axis=False):
        """

        :param figure_id:
        :param dpi:
        :param figure_width:
        :param figure_height:
        :param figure_width_per_subplot:
        :param figure_height_per_subplot:
        :param horizontal_subplot_number:
        :param vertical_subplot_number:
        :param share_x_axis:
        :param share_y_axis:
        Controls sharing of properties among x (sharex) or y (sharey) axes:
            True or 'all': x- or y-axis will be shared among all subplots.
            False or 'none': each subplot x- or y-axis will be independent.
            'row': each subplot row will share an x- or y-axis.
            'col': each subplot column will share an x- or y-axis.
        """

        self.width = width
        self.height = height

        self.width_per_subplot = width_per_subplot
        self.height_per_subplot = height_per_subplot

        self.share_x_axis = share_x_axis
        self.share_y_axis = share_y_axis

        self.dpi = dpi

    def apply(self):
        pass


default_figure_style = FigureStyle(dpi=300, width_per_subplot=None,
                                   height_per_subplot=None,
                                   share_x_axis=False, share_y_axis=False)

chromosome_figure_style = FigureStyle(dpi=300, width=15, height_per_subplot=None)

rainfall_figure_style = FigureStyle(dpi=300, width_per_subplot=None,
                                    height_per_subplot=2,
                                    share_x_axis=True, share_y_axis=True)
plot_figure_style = FigureStyle(dpi=300, width_per_subplot=None,
                                height_per_subplot=2,
                                share_x_axis=True, share_y_axis=True,
                                width=8, height=8)
