# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)

from cutstring import Constant
from datasets_helper import DatasetsHelperLight as DatasetsHelper

import os
"""
"""


def log_query(name, query, files):
    logger.debug("Query for %s: %s", name, query)
    if len(files) < 1:
        logger.critical("Query for %s did not return any files.", name)
        raise Exception
    for i_file, file_ in enumerate(files):
        logger.debug("File %d: %s", i_file + 1, file_)


class Era(object):
    def __init__(self, name, luminosity, database_path):
        self._name = name
        self._luminosity = luminosity
        self._database_path = database_path
        self._datasets_helper = DatasetsHelper(self._database_path)

    @property
    def datasets_helper(self):
        return self._datasets_helper

    @property
    def name(self):
        return self._name

    @property
    def lumi_weight(self):
        return Constant(str(self._luminosity), "lumi")


class Run2016(Era):
    def __init__(self, database_path):
        super(Run2016, self).__init__("Run2016", 35.87 * 1000.0, database_path)

    def data_files(self, channel):
        query = {
            "data": True,
            "campaign": "Run2016(B|C|D|E|F|G|H)",
            "scenario": "03Feb2017.*"
        }
        if channel.name == "mt":
            query["process"] = "SingleMuon"
        elif channel.name == "et":
            query["process"] = "SingleElectron"
        elif channel.name == "tt":
            query["process"] = "Tau"
        else:
            logger.critical("Channel %s is not implemented.", channel.name)
        files = self.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return files


class Run2017(Era):
    def __init__(self, database_path):
        super(Run2017, self).__init__("Run2017", 41.96 * 1000.0, database_path)
        #super(Run2017, self).__init__("Run2017", 18.90 * 1000.0, database_path) # for B, C, D only
        #super(Run2017, self).__init__("Run2017", 10.22 * 1000.0, database_path) # for Rereco B, C equivalent

    def data_files(self, channel):
        query = {
            "data": True,
            "campaign": "Run2017(B|C|D|E|F)",
            #"campaign": "Run2017(B|C|D)",
            #"campaign": "Run2017(B|C)",
            "scenario": "PromptRecov(1|2|3)"
        }
        if channel.name == "mt":
            query["process"] = "SingleMuon"
        elif channel.name == "et":
            query["process"] = "SingleElectron"
        elif channel.name == "tt":
            query["process"] = "Tau"
        elif channel.name == "em":
            query["process"] = "MuonEG"
        else:
            logger.critical("Channel %s is not implemented.", channel.name)
        files = self.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return files


class Run201712SepRereco(Era):
    def __init__(self, database_path):
        super(Run201712SepRereco, self).__init__("Run201712SepRereco", 10.22 * 1000.0, database_path)

    def data_files(self, channel):
        query = {
            "data": True,
            "campaign": "Run2017(B|C)",
            "scenario": "12Sep2017v1"
        }
        if channel.name == "mt":
            query["process"] = "SingleMuon"
        elif channel.name == "et":
            query["process"] = "SingleElectron"
        elif channel.name == "tt":
            query["process"] = "Tau"
        elif channel.name == "em":
            query["process"] = "MuonEG"
        else:
            logger.critical("Channel %s is not implemented.", channel.name)
        files = self.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return files
