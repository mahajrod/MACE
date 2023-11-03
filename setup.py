__author__ = 'mahajrod'

import os
from pathlib import Path
from setuptools import setup, find_packages

from os.path import join, dirname

scripts_dir_path = Path("scripts")

scripts = []
for script in scripts_dir_path.glob("*"):
    if script.is_dir():  # ignore all dirs
        continue
    else:
        scripts.append(script)
scripts = list(map(str, scripts)) # pip and setuptools doesnt work with Path objects
setup(name='MACE',
      version='1.1.27',
      packages=find_packages(),
      author='Sergei F. Kliver',
      author_email='mahajrod@gmail.com',
      install_requires=['scipy', 'pandas', 'numpy', 'matplotlib', 'biopython', 'bcbio-gff'],
      long_description=open(join(dirname(__file__), 'README.md')).read(),
      scripts=scripts)
