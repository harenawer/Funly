#!/usr/bin/env python

from setuptools import setup

setup(name='funly',
      version='0.1.0',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      license='GPLv3',
      package_dir={'': 'lib'},
      packages=['funly'],
      entry_points={'console_scripts': ['funly = funly.user:main2'],},
      )
