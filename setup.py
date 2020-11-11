# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

setup(
    name                 = 'pyMotif',
    version              = '0.0.1',
    #TODO: написать description
    description          = '''An utility for recognition of consensus motifs and other amino acid enrichments in large scale proteomics''',
    long_description     = (''.join(open('README.MD').readlines())),
    author               = 'Victoria Lineva & Julia Bubis',
    author_email         = 'lineva.vi@phystech.edu',
    install_requires     = ['pyteomics>=4.2', 'pandas','argparse','logging', 'os', 'collections', 're', 'scipy', 'numpy'],
    classifiers          = ['Intended Audience :: Science/Research',
                            'Programming Language :: Python :: 3',
                            'Topic :: Scientific/Engineering :: Bio-Informatics',
                            'Topic :: Scientific/Engineering :: Chemistry',
                            'Topic :: Scientific/Engineering :: Physics'],
    license              = 'License :: OSI Approved :: Apache Software License',
    packages             = find_packages(),
    entry_points         = {'console_scripts': ['pyMotif=pyMotif.main:main']}
    )
