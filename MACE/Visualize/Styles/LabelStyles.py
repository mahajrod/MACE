
class LabelStyle:

    def __init__(self, fontsize=7, rotation="horizontal", fontstyle='normal', fontweight='normal',
                 horizontal_alignment='center', vertical_alignment='center',
                 fontname=None,  fontfamily=None, zorder=120, zorder_shift=0):
        """
        :param fontname:                    one of system fonts
        :param fontfamily:                  'serif' | 'sans-serif' | 'cursive' | 'fantasy' | 'monospace'
        :param fontsize:                    int
        :param rotation:                    angle in degrees | 'vertical' | 'horizontal'
        :param fontstyle:                   'normal' | 'italic' | 'oblique'
        :param fontweight:                  'normal' | 'bold' | 'heavy' | 'light' | 'ultrabold' | 'ultralight'
        :param horizontal_alignment:        'center' | 'right' | 'left'
        :param vertical_alignment:          'center' | 'top' | 'bottom' | 'baseline'
        """
        # track
        self.fontsize = fontsize
        self.fontstyle = fontstyle
        self.fontweight = fontweight
        self.horizontal_alignment = horizontal_alignment
        self.vertical_alignment = vertical_alignment
        self.fontname = fontname
        self.fontfamily = fontfamily
        self.rotation = rotation
        self.zorder = zorder
        self.zorder_shift = zorder_shift


default_track_label_style = LabelStyle(fontsize=7, fontstyle='normal', fontweight='normal',
                                       horizontal_alignment='center', vertical_alignment='bottom',
                                       rotation="horizontal",
                                       fontname=None,  fontfamily=None, zorder=120)

default_htg_label_style = LabelStyle(fontsize=9, fontstyle='normal', fontweight='normal',
                                     horizontal_alignment='right', vertical_alignment='center',
                                     fontname=None,  fontfamily=None, zorder=120,
                                     rotation="horizontal")

default_vtg_label_style = LabelStyle(fontsize=11, fontstyle='normal', fontweight='normal',
                                     horizontal_alignment='right', vertical_alignment='center',
                                     fontname=None,  fontfamily=None, zorder=120,
                                     rotation="horizontal")

default_x_axis_label_style = LabelStyle(fontsize=7, fontstyle='normal', fontweight='normal',
                                        horizontal_alignment='center', vertical_alignment='top',
                                        rotation="horizontal",
                                        fontname=None,  fontfamily=None, zorder=120)

default_y_axis_label_style = LabelStyle(fontsize=7, fontstyle='normal', fontweight='normal',
                                        horizontal_alignment='right', vertical_alignment='center',
                                        rotation="vertical",
                                        fontname=None,  fontfamily=None, zorder=120)

default_x_tick_label_style = LabelStyle(fontsize=5, fontstyle='normal', fontweight='normal',
                                        horizontal_alignment='center', vertical_alignment='top',
                                        rotation="horizontal",
                                        fontname=None,  fontfamily=None, zorder=120)

default_y_tick_label_style = LabelStyle(fontsize=5, fontstyle='normal', fontweight='normal',
                                        horizontal_alignment='right', vertical_alignment='center',
                                        rotation="vertical",
                                        fontname=None,  fontfamily=None, zorder=120)
