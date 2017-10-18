# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)
"""
"""


# non-functional base class, just to express that the latter ones belong together
class Binning(object):
    pass


# TODO: Rewrite me
"""
class Variable_Binning(Binning):
    def __init__(self, *args):
        if not sorted(list(set(args))) == list(args):
            logger.fatal(
                "An invalid variable binning has been requested with a wrong bin border ordering or a repetition of bins."
            )
            logger.fatal("The binning was: %s.", args)
            raise Exception
        self.bin_borders = args

    def get_nbinsx(self):
        return len(self.bin_borders) - 1

    def extract(self):
        return str(self.bin_borders)
"""


class ConstantBinning(Binning):
    def __init__(self, nbinsx, xlow, xhigh):
        self._nbinsx = int(nbinsx)
        self._xlow = float(xlow)
        self._xhigh = float(xhigh)

    def extract(self):
        return "".join([
            "(",
            str(self._nbinsx), ",",
            str(self._xlow), ",",
            str(self._xhigh), ")"
        ])

    @property
    def nbinsx(self):
        return self._nbinsx
