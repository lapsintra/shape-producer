# -*- coding: utf-8 -*-

from cutstring import *
import logging
logger = logging.getLogger(__name__)
"""
"""


# Class to store a variable and its binning
class Variable(object):
    def __init__(self, name, binning):
        self._name = name
        self._binning = binning

    @property
    def name(self):
        return self._name

    @property
    def binning(self):
        return self._binning
