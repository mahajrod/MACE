__author__ = 'mahajrod'

import os

from setuptools import setup, find_packages

from os.path import join, dirname

scripts = os.listdir("scripts/")

scripts = ["scripts/%s" % script for script in scripts]


setup(name='MACE',
      version='1.1.1',
      packages=find_packages(),
      author='Sergei F. Kliver',
      author_email='mahajrod@gmail.com',
      install_requires=['scipy', 'numpy', 'matplotlib', 'biopython', 'bcbio-gff'],
      long_description=open(join(dirname(__file__), 'README.md')).read(),
      scripts=scripts)
