import sys

try:
    from setuptools import setup
except:
    from distutils.core import setup


#!/usr/bin/env python


import sys
if sys.version_info < (3, 4):
    sys.exit('Python 3.4 or greater is required.')

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup



with open('RELEASE') as f:
    lines = f.readlines()

version = lines[0]
version = version.strip().split()[-1]

VERSION = version

LICENSE = "MIT"


setup(
    name='LiBis',
    version=VERSION,
    description=(
        'Low input Bisulfite sequencing alignment'
    ),
    long_description='',
    author='Yue Yin',
    author_email='dangertrip@tamu.edu',
    maintainer='Deqiang Sun',
    maintainer_email='dsun@tamu.edu',
    license=LICENSE,
    packages=find_packages(),
    platforms=["all"],
    url='https://github.com/Dangertrip/LiBis',
    install_requires=[
        "matplotlib",
        "numpy",
        "pandas",
        "scikit-learn",
        "scipy",
        "seaborn",
        "pysam"
    ],
    scripts=[
        "bin/LiBis",
    ],
    include_package_data=True,
)
