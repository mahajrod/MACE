__author__ = 'mahajrod'

import re
import os
import sys
import bz2
import gzip
import shutil
from collections import Iterable, OrderedDict

if sys.version_info[0] == 3:
    from io import TextIOWrapper as file


def metaopen(filename, flags, buffering=None, compresslevel=5):
    if not isinstance(filename, str): # or isinstance(filename, gzip.GzipFile) or isinstance(filename, bz2.BZ2File):
        if isinstance(filename, file):
            return filename
        else:
            raise ValueError("ERROR!!! Not file object or str: {}".format(str(filename)))
    elif filename[-3:] == ".gz":
        return gzip.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
    elif filename[-4:] == ".bz2":
        return bz2.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
    else:
        if buffering is not None:
            return open(filename, flags, buffering=buffering)
        else:
            return open(filename, flags)
