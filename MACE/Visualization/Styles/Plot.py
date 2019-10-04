from numpy import array

from matplotlib.pyplot import get_cmap


class PlotStyle:

    def __init__(self, color=None, marker=None, linestyle='solid', linewidth=2, markersize=12, **kwargs):
        """

        :param color:
        :param marker:
            .   point
            ,   pixel
            o   circle
            v   triangle_down
            ^   triangle up
            +   plus
            x   x
            s   square
            D   diamond
            d   thin diamonds
            |   vertica line
            _   horizontal line
            *   star
            p   pentagon

            for other markers check
                https://matplotlib.org/3.1.1/api/markers_api.html#module-matplotlib.markers
        :param linestyle:
            -   or  'solid'
            --  or  'dased'
            -.  or  'dashdot'
            :   or  dotted
        :param linewidth:
        :param markersize:
        :param kwargs:
        """
        self.color = color
        self.marker = marker
        self.linestyle = linestyle
        self.linewidth = linewidth
        self.markersize = markersize
        self.kwargs = kwargs


default_plot_style = PlotStyle(color=None, marker=None, linestyle='solid', linewidth=2, markersize=12)
