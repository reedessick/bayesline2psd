#!/usr/bin/env python
__usage__ = "setpy.py command [--options]"
__description__ = "standard install script"
__author__ = "reed.essick@ligo.org"

#-------------------------------------------------

from distutils.core import setup

setup(
    name = 'bayesline2psd',
    version = '0.0',
    url = 'https://github.com/reedessick/bayesline2psd.git',
    author = __author__,
    author_email = 'reed.essick@ligo.org',
    description = __description__,
    liscence = 'MIT License',
    scripts = [
        'bin/chains2psd',
        'bin/psd2process',
    ],
    packages = [
        'bayesline2psd',
    ],
    data_files = [
    ],
    requires = [
    ],
)
