# -*- coding: utf-8 -*-

import os

from histogram import *
from cutstring import *
from systematics import *
from systematic_variations import *

import logging
logger = logging.getLogger(__name__)


class EstimationMethod(object):
    def __init__(self, name, folder, era, directory, channel, mc_campaign):
        self._directory = directory
        self._folder = folder
        self._name = name
        self._mc_campaign = mc_campaign
        self._channel = channel
        self._era = era

    def get_path(self, systematic, folder):
        return systematic.category.channel.name + "_" + folder + "/ntuple"

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def era(self):
        return self._era

    @property
    def channel(self):
        return self._channel

    def get_weights(self):
        return Weights(Weight("1.0", "constant"))

    def get_cuts(self):
        return Cuts()

    # function parsing the datasets helper to return the files
    # overwrite this function
    # TODO: Error handling for missing files!
    def get_files(self):
        raise NotImplementedError

    def artus_file_names(self, files):
        return [os.path.join(self._directory, f, "%s.root" % f) for f in files]

    # wrapper function for the Histogram creation performing the systematic shifts
    def apply_systematic_variations(self, systematic, settings):
        return systematic.variation.shifted_root_objects(settings)

    def define_root_objects(self, systematic):
        histogram_settings = []
        histogram_settings.append({
            "name":
            systematic.name,
            "inputfiles":
            self.get_files,
            "folder": [self.get_path, systematic, self._folder],
            "cuts":
            systematic.category.cuts +
            self.get_cuts(),  # TODO: get_cuts() is correct here?
            "weights":
            self.get_weights,
            "variable":
            systematic.category.variable
        })
        return histogram_settings

    # TODO: Make this less magic
    def create_root_objects(self, systematic):
        root_object_settings = self.define_root_objects(systematic)

        # execute the underlying functions
        root_object_settings = self.apply_systematic_variations(
            systematic, root_object_settings)
        for setting in root_object_settings:
            for key, value in setting.iteritems():
                if isinstance(value, list):
                    setting[key] = value[0](*value[1:])
                elif callable(value):
                    setting[key] = value()

        root_objects = []
        for setting in root_object_settings:
            root_objects.append(create_root_object(**setting))
        return root_objects

    # doing nothing, shape is exactly the histogram as default
    def do_estimation(self, systematic):
        if len(systematic.root_objects) != 1:
            logger.fatal(
                "There are %d histograms associated to the systematic with name %s, but not exactly 1.",
                len(systematic.root_objects), systematic.name)
            raise Exception
        return systematic.root_objects[0]
