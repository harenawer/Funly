#!/usr/bin/env python

from setuptools import setup

setup(name='funly',
      version='0.1.0',
      author='Sergio Pascual',
      author_email='sergiopr@astrax.fis.ucm.es',
      license='GPLv3',
      package_dir={'': 'src'},
      py_modules=['funly'],
      entry_points={'console_scripts': ['funly = funly:main2'],},
      )
