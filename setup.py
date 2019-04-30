#!/usr/bin/env python

from setuptools import find_packages
from numpy.distutils.core import setup, Extension

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = ['scitools-iris', 'numpy']

setup_requirements = []

test_requirements = []

fortran = Extension('fortran', libraries=['inter'],
                    sources=['source/lagranto/caltra.f90',
                             'source/lagranto/inter.f90',
                             'source/lagranto/trace.f90'])

setup(
    author="Leo Saffin",
    author_email='leo.saffin@physics.ox.ac.uk',
    classifiers=[
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 2 - Pre-Alpha',
        'Natural Language :: English',
        "Programming Language :: Python",
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Atmospheric Science'],
    description="Functionality built on iris I have used for work.",
    install_requires=requirements,
    license="GNU General Public License v3 or later (GPLv3+)",
    long_description=readme,
    include_package_data=True,
    keywords='pylagranto',
    name='pylagranto',
    packages=['pylagranto'],
    package_dir={'': 'source'},
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/LSaffin/pylagranto',
    version='0.1',

    libraries=[('inter', dict(sources=['source/lagranto/inter.f90']))],
    ext_modules=[fortran]
)
