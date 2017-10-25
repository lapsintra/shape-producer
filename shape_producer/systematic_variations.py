# -*- coding: utf-8 -*-
"""
"""
import logging
logger = logging.getLogger(__name__)


# this helper function can be used in case the systematic variation's name ends with "Down" and "Up"
def create_systematic_variations(name, property_name, systematic_variation):
    results = []
    results.append(systematic_variation(name, property_name, "Down"))
    results.append(systematic_variation(name, property_name, "Up"))
    return results


# class performing the systematic variation
class SystematicVariation(object):
    def __init__(self, name, direction):
        self._name = name
        self._direction = direction

    # TODO: init name vs getter name is different, confusing!
    @property
    def name(self):
        return self._name + self._direction

    def change_histogram_name(self, h_settings, direction):
        if isinstance(h_settings["name"], list):
            h_settings["name"].append(direction)
        else:
            h_settings["name"] = [h_settings["name"], direction]
        return h_settings

    def shifted_root_objects(self, h_settings):
        return h_settings

    def is_nominal(self):
        return False


class Nominal(SystematicVariation):
    # TODO: Do this with super?
    def __init__(self, direction=None):
        self._name = "Nominal"
        self._direction = direction

    @property
    def name(self):
        name = self._name
        if self._direction:
            name += "_" + self._direction
        return name

    def change_histogram_name(self, h_settings, direction):
        return h_settings

    def shifted_root_objects(self, h_settings):
        return h_settings

    def is_nominal(self):
        return True


class DifferentPipeline(SystematicVariation):
    def __init__(self, name, pipeline, direction):
        super(DifferentPipeline, self).__init__(name, direction)
        self._pipeline = pipeline

    def shifted_root_objects(self, h_settings):
        for index in range(len(h_settings)):
            h_settings[index]["folder"][2] = self._pipeline + self._direction
        return h_settings


class SquareAndRemoveWeight(SystematicVariation):
    def __init__(self, name, weight_name, direction):
        super(SquareAndRemoveWeight, self).__init__(name, direction)
        self._weight_name = weight_name

    def shifted_root_objects(self, h_settings):
        for index in range(len(h_settings)):
            if self._direction == "Up":
                h_settings[index]["weights"] = h_settings[index][
                    "weights"]().square(self._weight_name)
            elif self._direction == "Down":
                h_settings[index]["weights"] = h_settings[index][
                    "weights"]().remove(self._weight_name)
        return h_settings


class ReplaceWeight(SystematicVariation):
    def __init__(self, name, weight_name, new_weight, direction):
        super(ReplaceWeight, self).__init__(name, direction)
        self._weight_name = weight_name
        self._new_weight = new_weight

    def shifted_root_objects(self, h_settings):
        for index in range(len(h_settings)):
            h_settings[index]["weights"] = h_settings[index][
                "weights"]().remove(self._weight_name)  #.add(self._new_weight)
            h_settings[index]["weights"].add(self._new_weight)
        return h_settings
