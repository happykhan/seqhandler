#!/usr/bin/env python


import os
import sys

import seqhandler.__init__ as meta

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

packages = [
    'to-gff',
]

requires = []
with open('requirements.txt') as fin:
    lines = fin.readlines()
for l in lines:
    requires.append(l.strip())

setup(
    name='to-gff',
    version=meta.__version__,
    description='to-gff - Generate gff file from EMBL/Genbank for QUAST',
    long_description=open('README.rst').read(),
    author='Mitchell Stanton-Cook',
    author_email='m.stantoncook@gmail.com',
    url='https://github.com/mscook/to-gff',
    packages=packages,
    scripts = [],
    package_data={'': ['LICENSE']},
    package_dir={'to-gff': 'to-gff'},
    include_package_data=True,
    install_requires=requires,
    license=open('LICENSE').read(),
    zip_safe=False,
    classifiers=(
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 2 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ),
)
