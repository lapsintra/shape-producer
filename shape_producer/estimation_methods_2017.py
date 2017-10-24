# -*- coding: utf-8 -*-

from cutstring import *
from estimation_methods import EstimationMethod
from estimation_methods_2016 import DataEstimation as DataEstimation2016
from estimation_methods_2016 import QCDEstimation as QCDEstimation2016
from estimation_methods_2016 import VVEstimation as VVEstimation2016
from era import log_query


class DataEstimation(DataEstimation2016):
    pass


class QCDEstimation(QCDEstimation2016):
    pass


class VVEstimation(VVEstimation2016):
    pass


class ZttEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(ZttEstimation, self).__init__(
            name="Ztt",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer17MiniAOD")

    def get_weights(self):
        return Weights(

            # MC related weights
            Weight("generatorWeight","generatorWeight"),
            #Weight("numberGeneratedEventsWeight","numberGeneratedEventsWeight"), # to be used only for one inclusive sample
            #Weight("crossSectionPerEventWeight","crossSectionPerEventWeight"), # to be used only for one inclusive sample
            Weight("(((genbosonmass >= 50.0 && (npartons == 0 || npartons >= 4))*5.75970078e-5) + ((genbosonmass >= 50.0 && npartons == 1)*1.36277241e-5) + ((genbosonmass >= 50.0 && npartons == 2)*7.42888435e-6) + ((genbosonmass >= 50.0 && npartons == 3)*1.62808443e-5) + ((genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight))", "z_stitching_weight"),

            # Weights for corrections
            #Weight("zPtReweightWeight", "zPtReweightWeight"),
            #Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))", "hadronic_tau_sf"),

            # Data related scale-factors
            self.era.lumi_weight
            )

    def get_cuts(self):

        ztt_genmatch_cut = Cut("1 == 1","ztt_genmatch")
        if self.channel.name in ["mt","et"]:
            ztt_genmatch_cut = Cut("gen_match_2==5", "ztt_genmatch")
        elif self.channel.name == "tt":
            ztt_genmatch_cut = Cut("(gen_match_1==5) && (gen_match_2==5)", "ztt_genmatch")
        elif self.channel.name == "tt":
            ztt_genmatch_cut = Cut("(gen_match_1==3) && (gen_match_2==4)", "ztt_genmatch")
        return Cuts(ztt_genmatch_cut)

    def get_files(self):
        query = {
            "process": "(DYJetsToLL_M10to50|DY.?JetsToLL_M50)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            #"version": "v1"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ZllEstimation(ZttEstimation):
    def __init__(self, era, directory, channel):
        super(ZttEstimation, self).__init__(
            name="Zll",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer17MiniAOD")

    def get_cuts(self):
        zll_genmatch_cut = Cut("1 == 1","zll_genmatch_mt")
        if self.channel.name in ["mt","et"]:
            zll_genmatch_cut = Cut("gen_match_2!=5", "zll_genmatch")
        elif self.channel.name == "tt":
            zll_genmatch_cut = Cut("(gen_match_1!=5) || (gen_match_2!=5)", "zll_genmatch")
        elif self.channel.name == "em":
            zll_genmatch_cut = Cut("(gen_match_1!=3) || (gen_match_2!=4)", "zll_genmatch")
        return Cuts(zll_genmatch_cut)

class WJetsEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(WJetsEstimation, self).__init__(
            name="WJets",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer17MiniAOD")

    def get_weights(self):
        return Weights(

            # MC related weights
            Weight("generatorWeight","generatorWeight"),
            #Weight("numberGeneratedEventsWeight","numberGeneratedEventsWeight"), # to be used only for one inclusive sample
            #Weight("crossSectionPerEventWeight","crossSectionPerEventWeight"), # to be used only for one inclusive sample
            Weight("(((npartons == 0 || npartons >= 5)*2.36006270e-3) + ((npartons == 1)*2.34817764e-4) + ((npartons == 2)*1.31144867e-4) + ((npartons == 3)*1.39177532e-4) + ((npartons == 4)*6.46064804e-5))", "wj_stitching_weight"),

            # Weights for corrections
            #Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))", "hadronic_tau_sf"),

            # Data related scale-factors
            self.era.lumi_weight
            )

    def get_files(self):
        query = {
            "process": "W.?JetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class TTEstimation(EstimationMethod):
    def __init__(self, era, directory, channel):
        super(TTEstimation, self).__init__(
            name="TT",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign="RunIISummer17MiniAOD")

    def get_weights(self):
        return Weights(

            # MC related weights
            Weight("generatorWeight","generatorWeight"),
            Weight("numberGeneratedEventsWeight","numberGeneratedEventsWeight"), # to be used only for one inclusive sample
            Weight("crossSectionPerEventWeight","crossSectionPerEventWeight"), # to be used only for one inclusive sample

            # Weights for corrections
            Weight("topPtReweightWeight", "topPtReweightWeight"),
            #Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))", "hadronic_tau_sf"),

            # Data related scale-factors
            self.era.lumi_weight
            )

    def get_files(self):
        query = {
            "process": "^TT$",
            "data": False,
            "campaign": self._mc_campaign
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)
