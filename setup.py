# -*- coding: utf-8 -*-
""" Setup configuration file
"""

from setuptools import setup

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

setup(
    name='pathoscopeutils',
    version='0.1',
    description='Utility scripts for analyzing PathoScope output',
    url='https://github.com/gwcbi/PathoScopeUtils',
    author='Matthew L. Bendall',
    license='MIT',
    packages=[
        'pathoscopeutils',
    ],
    zip_safe=False,
    # entry_points = {
    #     'console_scripts': [
    #         'somescript = ',
    #     ],
    # }
)
