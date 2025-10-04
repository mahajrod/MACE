__author__ = 'mahajrod'

import os
from pathlib import Path
from setuptools import setup, find_packages

from os.path import join, dirname

"""
scripts_dir_path = Path("scripts")

scripts = []
for script in scripts_dir_path.glob("*"):
    if script.is_dir():  # ignore all dirs
        continue
    else:
        scripts.append(script)

scripts = list(map(str, scripts)) # pip and setuptools doesnt work with Path objects
"""


def get_list_of_files(list_of_dirs_and_files, expression=None, recursive=True):
    file_list = []

    for entry in [list_of_dirs_and_files] if isinstance(list_of_dirs_and_files, str) else list_of_dirs_and_files:

        if os.path.isdir(entry):

            files_in_dir = ["%s%s" % (entry if entry[-1] == "/" else (entry+"/"), filename)
                            for filename in sorted(filter(expression, os.listdir(entry))
                                                   if expression else os.listdir(entry))]
            if recursive:
                for filename in files_in_dir:
                    if os.path.isdir(filename):
                        file_list += get_list_of_files([filename], expression=expression, recursive=True)
                    else:
                        file_list.append(filename)
            else:
                file_list += files_in_dir
        elif os.path.exists(entry):
            file_list.append(os.path.abspath(entry))
        else:
            print("%s does not exist" % entry)
    return file_list


setup(name='MACE',
      version='1.1.35',
      packages=find_packages(),
      author='Sergei F. Kliver',
      author_email='mahajrod@gmail.com',
      install_requires=['scipy', 'pandas', 'numpy', 'matplotlib', 'biopython', 'bcbio-gff'],
      long_description=open(join(dirname(__file__), 'README.md')).read(),
      scripts=get_list_of_files("scripts/"))
