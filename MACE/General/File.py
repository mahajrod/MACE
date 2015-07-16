#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os


def split_filename(filepath):
    directory, basename = os.path.split(filepath)
    prefix, extension = os.path.splitext(basename)
    return directory, prefix, extension


def make_list_of_path_to_files(list_of_dirs_and_files, expression):

    pathes_list = []
    for entry in list_of_dirs_and_files:
        if os.path.isdir(entry):
            files_in_dir = sorted(filter(expression, os.listdir(entry)))
            for filename in files_in_dir:
                pathes_list.append("%s%s" % (entry, filename))
        elif os.path.exists(entry):
            pathes_list.append(entry)
        else:
            print("%s does not exist" % entry)

    return pathes_list
