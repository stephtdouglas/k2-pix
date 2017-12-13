#!/usr/bin/env python
import sys
from setuptools import setup
# To use a consistent encoding
from codecs import open
import os
from os import path

# Prepare and send a new release to PyPI
if "release" in sys.argv[-1]:
    os.system("python setup.py sdist")
    os.system("twine upload dist/*")
    os.system("rm -rf dist/k2-pix*")
    sys.exit()

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README'), encoding='utf-8') as f:
    long_description = f.read()

entry_points = {'console_scripts': ['k2pix = main:k2pix']}

setup(
    name='k2-pix',
    packages=['k2-pix'],
    version='0.1',
    description="Overlaying a Sky View survey image's contours onto a K2 pixel stamp",
    long_description=long_description,
    url='https://github.com/gully/k2-pix',
    author='Stephanie Douglas and GitHub contributors',
    author_email='',
    license='MIT',
    entry_points=entry_points,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest', 'pytest-cov'],
    keywords='astronomy astrophysics',
    install_requires=['pyketools', 'astroquery', 'numpy', 'matplotlib', 'tqdm', 'astropy']
)
