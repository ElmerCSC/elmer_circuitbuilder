#!/usr/bin/env python3

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = ['pytest>=3', ]

setup(
    author="Jonathan Velasco",
    author_email='jvelasco@csc.fi',
    maintainer="CSC - IT Center for Science :  ElmerFEM Developers",
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="Python module for creating the circuit simulation definitions for Elmer FEM. The circuit definitions enable easy setup of coils (e.g. massive, stranded, and foil) in 2D and 3D for magnetodynamics applications.",
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='elmer_circuitbuilder',
    name='elmer_circuitbuilder',
    #packages=find_packages(include=['elmer_circuitbuilder', 'elmer_circuitbuilder.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/ElmerCSC/elmer_circuitbuilder.git',
    version='v0.0.5',
    zip_safe=False,
    extras_require={"dev":["pytest>=3.7"], },
    py_modules=["elmer_circuitbuilder"],
    package_dir={'': 'elmer_circuitbuilder'},

)
