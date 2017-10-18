# -*- coding: utf-8 -*-

# TODO: Factorize this step into new python module? This messes up the independence
# of the systematics generation.
import CombineHarvester.CombineTools.ch as ch

import os
import re

import logging
logger = logging.getLogger(__name__)


class CombineHarvester(object):
    def __init__(self, systematics):
        if not systematics.is_produced:
            logger.fatal("Systematics have not been produced.")
            raise Exception
        self._systematics = systematics

        self._ch = ch.CombineHarvester()
        if logger.isEnabledFor(logging.DEBUG):
            self._ch.SetVerbosity(1)

        self._is_produced = False

    def produce(self):
        self.add_observations()
        self.add_processes()
        self.add_shape_systematics()
        self._is_produced = True
        logger.debug("Successfully produced datacard.")

    def write(self, output_file_prefix):
        if not self._is_produced:
            logger.fatal("Datacard has not been produced.")
            raise Exception
        logger.info("Create datacard files %s.txt and %s.root.",
                    output_file_prefix, output_file_prefix)
        writer = ch.CardWriter("{}.txt".format(output_file_prefix),
                               "{}.root".format(output_file_prefix))
        if logger.isEnabledFor(logging.DEBUG):
            writer.SetVerbosity(1)
        writer.CreateDirectories(
            False)  # TODO: FIXME: Does not work without this?
        # writer.SetWildcardMasses([]) # TODO: What is this doing?
        writer.WriteCards("", self._ch)

    def print_datacard(self):
        self._ch.PrintAll()

    def add_observations(self):
        map_obs = self.create_systematics_map("observed")
        if map_obs == {}:
            logger.warning("No observation process found.")
        for name in map_obs:
            mass = [map_obs[name]["mass"]]
            analysis = [map_obs[name]["analysis"]]
            era = [map_obs[name]["era"]]
            channel = [map_obs[name]["channel"]]
            processes = [k for k in map_obs[name]["processes"]]
            categories = [
                (i, c)
                for i, c in enumerate(map_obs[name]["processes"][processes[0]])
            ]

            logger.debug(
                "Add for %s observation with name %s and categories %s", name,
                processes, categories)
            self._ch.AddObservations(mass, analysis, era, channel, categories)

    def add_processes(self):
        flag = {"signal": True, "background": False}
        for type in flag:
            map = self.create_systematics_map(type)
            if map == {}:
                logger.warning("No %s processes found.", type)
            for name in map:
                mass = [map[name]["mass"]]
                analysis = [map[name]["analysis"]]
                era = [map[name]["era"]]
                channel = [map[name]["channel"]]
                processes = [k for k in map[name]["processes"]]
                categories = [
                    (i, c)
                    for i, c in enumerate(map[name]["processes"][processes[0]])
                ]

                logger.debug("Add for %s %s processes %s and categories %s",
                             name, type, processes, categories)
                self._ch.AddProcesses(mass, analysis, era, channel, processes,
                                      categories, flag[type])

    def add_shape_systematics(self):
        # Add shape systematics to datacard
        # NOTE: It seems that this has to be done after extracting all shapes.
        for s in self._systematics.systematics:
            # TODO: We assume silently that the Down shift exists, make an error handling!
            if s.variation.direction == "Up":
                name = re.search(
                    "{PROCESS}_{BIN}_{ANALYSIS}_{ERA}_{VARIABLE}_(.*)Up".
                    format(
                        PROCESS=s.process.name,
                        BIN=s.category.name,
                        ANALYSIS=s.analysis,
                        ERA=s.era.name,
                        VARIABLE=s.category.variable.name), s.name)
                if len(name.groups()) != 1:
                    logger.fatal(
                        "Regex did not return basename of variation %s.",
                        name.groups())
                    raise Exception
                basename = name.groups(
                )[0]  # TODO: Rename systematics so that the trailing underscore goes away.
                process = s.process.name
                channel = s.channel.name
                strength = s.variation.strength
                logger.debug(
                    "Add systematic %s for process %s and channel %s with strength %.4f.",
                    basename, process, channel, strength)
                self._ch.cp().process([process]).channel([channel]).AddSyst(
                    self._ch, basename, "shape", ch.SystMap()(strength))

        # Extract shape systematics (includes nominal shapes)
        for s in self._systematics.systematics:
            process = s.process.name
            era = s.era.name
            channel = s.channel.name
            logger.debug("Extract shape %s.", s.name)
            self._ch.cp().process([
                process
            ]).era([era]).channel([channel]).ExtractShapes(
                self._systematics.output_file, s.name,
                "$PROCESS_$BIN_{ANALYSIS}_{ERA}_{VARIABLE}_$SYSTEMATIC".format(
                    ANALYSIS=s.analysis,
                    ERA=s.era.name,
                    VARIABLE=s.category.variable.name))

    def add_bin_by_bin(self, threshold, merge_threshold, fix_norm):
        if not self._is_produced:
            logger.fatal("Datacard has not been produced.")
            raise Exception
        logger.debug("Add bin-by-bin uncertainties.")
        bbb = ch.BinByBinFactory()
        if logger.isEnabledFor(logging.DEBUG):
            bbb.SetVerbosity(1)
        bbb.SetAddThreshold(threshold).SetMergeThreshold(
            merge_threshold).SetFixNorm(fix_norm)
        # TODO: FIXME: Only backgrounds?
        bbb.MergeBinErrors(self._ch.cp().backgrounds())
        bbb.AddBinByBin(self._ch.cp().backgrounds(), self._ch)

    # TODO: Make this easier?
    def create_systematics_map(self, type):
        map = {}
        for s in self._systematics.systematics:
            if s.variation.is_nominal():
                if s.process.type == type:
                    mass = str(s.mass)
                    analysis = s.analysis
                    era = s.era.name
                    channel = s.channel.name
                    name = "_".join([mass, analysis, era, channel])

                    category = s.category.name
                    process = s.process.name
                    if not name in map:
                        map[name] = {
                            "mass": mass,
                            "analysis": analysis,
                            "era": era,
                            "channel": channel,
                            "processes": {}
                        }
                    if not process in map[name]["processes"]:
                        map[name]["processes"][process] = []
                    map[name]["processes"][process].append(category)

        for name in map:
            channels = None
            for process in map[name]["processes"]:
                if channels == None:
                    channels = map[name]["processes"][process]
                else:
                    if not (channels == map[name]["processes"][process]):
                        logger.fatal(
                            "Found processes in the same analysis with different categories: %s",
                            map)
                        raise Exception

        if type == "observed":
            for name in map:
                if len(map[name]["processes"]) != 1:
                    logger.fatal(
                        "Found for the observation more than exactly one process is found: %s",
                        map[name]["processes"])
                    raise Exception

        return map
