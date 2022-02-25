__author__ = 'mahajrod'

import os
from pathlib import Path
from setuptools import setup, find_packages


setup(name='MACE',
      version='1.1.4',
      packages=find_packages(),
      author='Sergei F. Kliver',
      author_email='mahajrod@gmail.com',
      install_requires=['scipy', 'numpy', 'matplotlib', 'biopython', 'bcbio-gff'],
      long_description=open(join(dirname(__file__), 'README.md')).read(),
      scripts=list(map(str, sorted(Path('scripts/').rglob("*.py")))))
