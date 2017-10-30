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


class SStoOSEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process):
        super(QCDEstimation, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign=None)
        self._bg_processes = [copy.deepcopy(p) for p in bg_processes]
        self._data_process = copy.deepcopy(data_process)

    def create_root_objects(self, systematic):
        ss_category = copy.deepcopy(systematic.category)
        ss_category.cuts.get("os").name = "ss"
        ss_category.cuts.get("ss").invert()
        ss_category.name = ss_category._name + "_ss"

        root_objects = []
        systematic._qcd_systematics = []
        for process in [self._data_process] + self._bg_processes:
            s = Systematic(
                category=ss_category,
                process=process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125)
            systematic._qcd_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects
        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_qcd_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _qcd_systematics needed for QCD estimation.",
                systematic.name)
            raise Exception

        # Create shapes
        for s in systematic._qcd_systematics:
            s.do_estimation()

        # Data shape
        shape = systematic._qcd_systematics[0].shape

        # Subtract MC shapes from data shape
        for s in systematic._qcd_systematics[1:]:
            shape.result.Add(s.shape.result, -1.0)

        # Test that not a single bin in TH1F shape.result is negative
        if shape.has_negative_entries():
            logger.fatal(
                "Subtraction of Monte Carlo from data results in negative number of events."
            )
            raise Exception

        # Rename root object accordingly
        shape.name = systematic.name
        return shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError