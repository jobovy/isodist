from setuptools import setup #, Extension
import os, os.path
import re

longDescription= ""


setup(name='isodist',
      version='1.',
      description='spectro-photometric distances to stars',
      author='Jo Bovy',
      author_email='bovy@ias.edu',
      license='New BSD',
      long_description=longDescription,
      url='http://code.google.com/p/isodist/',
      package_dir = {'isodist/': ''},
      packages=['isodist'],
#      dependency_links = ['https://github.com/dfm/MarkovPy/tarball/master#egg=MarkovPy'],
      install_requires=['numpy','scipy']
      )
