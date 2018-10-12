# -*- coding: utf-8 -*-

import copy
import os

from estimation_methods import EstimationMethod, SStoOSEstimationMethod, ABCDEstimationMethod, SumUpEstimationMethod
from histogram import *
from cutstring import *
from process import *
from systematics import *
from systematic_variations import *
from era import log_query

import logging
logger = logging.getLogger(__name__)


class DataEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(DataEstimation, self).__init__(
            name="data_obs",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign=None)
        self._channel = channel

    def get_files(self):
        return self.artus_file_names(self.era.data_files(self._channel))

    def get_cuts(self):
        return Cuts()


class FakeEstimationLT(DataEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(DataEstimation, self).__init__(
            name="fakes",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign=None)
        self._channel = channel

    def get_weights(self):
        return Weights(Weight("ff2_nom", "fake_factor"))

    def create_root_objects(self, systematic):
        aiso_systematic = copy.deepcopy(systematic)
        aiso_systematic.category.cuts.remove("tau_iso")
        aiso_systematic.category.cuts.add(
            Cut(
                "byTightIsolationMVArun2v1DBoldDMwLT_2<0.5&&byVLooseIsolationMVArun2v1DBoldDMwLT_2>0.5",
                "tau_aiso"))
        return super(FakeEstimationLT,
                     self).create_root_objects(aiso_systematic)


class FakeEstimationTT(DataEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(DataEstimation, self).__init__(
            name="fakes",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign=None)
        self._channel = channel

    def get_weights(self):
        return Weights(
            Weight(
                "(0.5*ff1_nom*(byTightIsolationMVArun2v1DBoldDMwLT_1<0.5)+0.5*ff2_nom*(byTightIsolationMVArun2v1DBoldDMwLT_2<0.5))",
                "fake_factor"))

    def create_root_objects(self, systematic):
        aiso_systematic = copy.deepcopy(systematic)
        aiso_systematic.category.cuts.remove("tau_1_iso")
        aiso_systematic.category.cuts.remove("tau_2_iso")
        aiso_systematic.category.cuts.add(
            Cut(
                "(byTightIsolationMVArun2v1DBoldDMwLT_2>0.5&&byTightIsolationMVArun2v1DBoldDMwLT_1<0.5&&byVLooseIsolationMVArun2v1DBoldDMwLT_1>0.5)||(byTightIsolationMVArun2v1DBoldDMwLT_1>0.5&&byTightIsolationMVArun2v1DBoldDMwLT_2<0.5&&byVLooseIsolationMVArun2v1DBoldDMwLT_2>0.5)",
                "tau_aiso"))
        return super(FakeEstimationTT,
                     self).create_root_objects(aiso_systematic)


class HTTEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="HTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight(
                "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
                "hadronic_tau_sf"), Weight("eventWeight", "eventWeight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "(^GluGluHToTauTau.*125.*|^VBFHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ggHEstimation(HTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight(
                "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
                "hadronic_tau_sf"), Weight("ggh_NNLO_weight", "gghNNLO"),
            Weight("eventWeight", "eventWeight"), self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^GluGluHToTauTau.*125.*",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ggHEstimation_VBFTOPO_JET3VETO(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH_VBFTOPO_JET3VETO",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==101", "htxs_match"))


class ggHEstimation_VBFTOPO_JET3(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH_VBFTOPO_JET3",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==102", "htxs_match"))


class ggHEstimation_0J(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH_0J",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==103", "htxs_match"))


class ggHEstimation_1J_PTH_0_60(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH_1J_PTH_0_60",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==104", "htxs_match"))


class ggHEstimation_1J_PTH_60_120(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH_1J_PTH_60_120",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==105", "htxs_match"))


class ggHEstimation_1J_PTH_120_200(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH_1J_PTH_120_200",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==106", "htxs_match"))


class ggHEstimation_1J_PTH_GT200(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH_1J_PTH_GT200",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==107", "htxs_match"))


class ggHEstimation_GE2J_PTH_0_60(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH_GE2J_PTH_0_60",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==108", "htxs_match"))


class ggHEstimation_GE2J_PTH_60_120(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH_GE2J_PTH_60_120",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==109", "htxs_match"))


class ggHEstimation_GE2J_PTH_120_200(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH_GE2J_PTH_120_200",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==110", "htxs_match"))


class ggHEstimation_GE2J_PTH_GT200(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH_GE2J_PTH_GT200",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==111", "htxs_match"))


class qqHEstimation(HTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="qqH",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_files(self):
        query = {
            "process": "^VBFHToTauTau.*125.*",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class qqHEstimation_VBFTOPO_JET3VETO(qqHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="qqH_VBFTOPO_JET3VETO",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==201", "htxs_match"))


class qqHEstimation_VBFTOPO_JET3(qqHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="qqH_VBFTOPO_JET3",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==202", "htxs_match"))


class qqHEstimation_VH2JET(qqHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="qqH_VH2JET",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==203", "htxs_match"))


class qqHEstimation_REST(qqHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="qqH_REST",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==204", "htxs_match"))


class qqHEstimation_PTJET1_GT200(qqHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="qqH_PTJET1_GT200",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1cat==205", "htxs_match"))


class VHEstimation(HTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="VH",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_files(self):
        query = {
            "process": "(^W(minus|plus)HToTauTau.*125.*|^ZHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ZTTEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(ZTTEstimation, self).__init__(
            name="ZTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight("eventWeight", "eventWeight"),
            Weight("zPtReweightWeight", "zPtReweightWeight"),
            Weight(
                "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
                "hadronic_tau_sf"),
            Weight(
                "((((genbosonmass >= 150.0 && (npartons == 0 || npartons >= 5))*3.95423374e-5) + ((genbosonmass >= 150.0 && npartons == 1)*1.27486147e-5) + ((genbosonmass >= 150.0 && npartons == 2)*1.3012785e-5) + ((genbosonmass >= 150.0 && npartons == 3)*1.33802133e-5) + ((genbosonmass >= 150.0 && npartons == 4)*1.09698723e-5)+((genbosonmass >= 50.0 && genbosonmass < 150.0 && (npartons == 0 || npartons >= 5))*3.95423374e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 1)*1.27486147e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 2)*1.3012785e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 3)*1.33802133e-5) + ((genbosonmass >= 50.0 && genbosonmass < 150.0 && npartons == 4)*1.09698723e-5)+((genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight))/(numberGeneratedEventsWeight*crossSectionPerEventWeight*sampleStitchingWeight))",
                "z_stitching_weight"), self.era.lumi_weight)

    def get_cuts(self):
        return Cuts(Cut("gen_match_2==5", "ztt_genmatch_mt"))

    def get_files(self):
        query = {
            "process":
            "(DYJetsToLL_M10to50|DYJetsToLL_M50|DY1JetsToLL_M50|DY2JetsToLL_M50|DY3JetsToLL_M50|DY4JetsToLL_M50)",
            "data":
            False,
            "campaign":
            self._mc_campaign,
            "generator":
            "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ZTTEstimationTT(ZTTEstimation):
    def get_cuts(self):
        return Cuts(Cut("(gen_match_1==5&&gen_match_2==5)", "ztt_genmatch_tt"))


class ZTTEstimationLL(ZTTEstimation):
    def get_cuts(self):
        return Cuts(
            Cut(
                "(gen_match_1==3||gen_match_1==4)&&(gen_match_2==3||gen_match_2==4)",
                "ztt_genmatch_ll"))


class ZTTEmbeddedEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(ZTTEmbeddedEstimation, self).__init__(
            name="EMB",
            folder="nominal",
            era=era,
            friend_directory=friend_directory,
            directory=directory,
            channel=channel,
            mc_campaign=None)

    def embedding_stitchingweight(self):
        if self.channel.name == 'mt':
            comp_eff_B = "(1.0/0.899)"
            comp_eff_C = "(1.0/0.881)"
            comp_eff_D = "(1.0/0.877)"
            comp_eff_E = "(1.0/0.939)"
            comp_eff_F = "(1.0/0.936)"
            comp_eff_G = "(1.0/0.908)"
            comp_eff_H = "(1.0/0.962)"
            runB = "((run >= 272007) && (run < 275657))*" + comp_eff_B
            runC = "+((run >= 275657) && (run < 276315))*" + comp_eff_C
            runD = "+((run >= 276315) && (run < 276831))*" + comp_eff_D
            runE = "+((run >= 276831) && (run < 277772))*" + comp_eff_E
            runF = "+((run >= 277772) && (run < 278820))*" + comp_eff_F
            runG = "+((run >= 278820) && (run < 280919))*" + comp_eff_G
            runH = "+((run >= 280919) && (run < 284045))*" + comp_eff_H
            return "(" + runB + runC + runD + runE + runF + runG + runH + ")"
        elif self.channel.name == 'et':
            comp_eff_B = "(1.0/0.902)"
            comp_eff_C = "(1.0/0.910)"
            comp_eff_D = "(1.0/0.945)"
            comp_eff_E = "(1.0/0.945)"
            comp_eff_F = "(1.0/0.915)"
            comp_eff_G = "(1.0/0.903)"
            comp_eff_H = "(1.0/0.933)"
            runB = "((run >= 272007) && (run < 275657))*" + comp_eff_B
            runC = "+((run >= 275657) && (run < 276315))*" + comp_eff_C
            runD = "+((run >= 276315) && (run < 276831))*" + comp_eff_D
            runE = "+((run >= 276831) && (run < 277772))*" + comp_eff_E
            runF = "+((run >= 277772) && (run < 278820))*" + comp_eff_F
            runG = "+((run >= 278820) && (run < 280919))*" + comp_eff_G
            runH = "+((run >= 280919) && (run < 284045))*" + comp_eff_H
            return "(" + runB + runC + runD + runE + runF + runG + runH + ")"
        elif self.channel.name == 'tt':
            comp_eff_B = "(1.0/0.897)"
            comp_eff_C = "(1.0/0.908)"
            comp_eff_D = "(1.0/0.950)"
            comp_eff_E = "(1.0/0.861)"
            comp_eff_F = "(1.0/0.941)"
            comp_eff_G = "(1.0/0.908)"
            comp_eff_H = "(1.0/0.949)"
            runB = "((run >= 272007) && (run < 275657))*" + comp_eff_B
            runC = "+((run >= 275657) && (run < 276315))*" + comp_eff_C
            runD = "+((run >= 276315) && (run < 276831))*" + comp_eff_D
            runE = "+((run >= 276831) && (run < 277772))*" + comp_eff_E
            runF = "+((run >= 277772) && (run < 278820))*" + comp_eff_F
            runG = "+((run >= 278820) && (run < 280919))*" + comp_eff_G
            runH = "+((run >= 280919) && (run < 284045))*" + comp_eff_H
            return "(" + runB + runC + runD + runE + runF + runG + runH + ")"
        elif self.channel.name == 'em':
            comp_eff_B = "(1.0/0.891)"
            comp_eff_C = "(1.0/0.910)"
            comp_eff_D = "(1.0/0.953)"
            comp_eff_E = "(1.0/0.947)"
            comp_eff_F = "(1.0/0.942)"
            comp_eff_G = "(1.0/0.906)"
            comp_eff_H = "(1.0/0.950)"
            runB = "((run >= 272007) && (run < 275657))*" + comp_eff_B
            runC = "+((run >= 275657) && (run < 276315))*" + comp_eff_C
            runD = "+((run >= 276315) && (run < 276831))*" + comp_eff_D
            runE = "+((run >= 276831) && (run < 277772))*" + comp_eff_E
            runF = "+((run >= 277772) && (run < 278820))*" + comp_eff_F
            runG = "+((run >= 278820) && (run < 280919))*" + comp_eff_G
            runH = "+((run >= 280919) && (run < 284045))*" + comp_eff_H
            return "(" + runB + runC + runD + runE + runF + runG + runH + ")"
        else:
            log.error("Embedded currently not implemented for channel \"%s\"!"
                      % self.channel.name)

    def get_weights(self):
        if self.channel.name == "mt":
            return Weights(
                Weight("generatorWeight*(generatorWeight<=1.0)",
                       "simulation_sf"),
                Weight("muonEffTrgWeight", "scale_factor"),
                Weight(self.embedding_stitchingweight(),
                       "2016 stitching weight"),
                Weight("idWeight_1*triggerWeight_1*isoWeight_1", "lepton_sf"),
                Weight("1.0", "mutau_crosstriggerweight"),
                Weight("(gen_match_2==5)*1.02+(gen_match_2!=5)", "emb_tau_id"),
                Weight("embeddedDecayModeWeight", "decayMode_SF"))
        if self.channel.name == "et":
            return Weights(
                Weight("generatorWeight*(generatorWeight<=1.0)",
                       "simulation_sf"),
                Weight("muonEffTrgWeight", "scale_factor"),
                Weight(self.embedding_stitchingweight(),
                       "2016 stitching weight"),
                Weight("idWeight_1*triggerWeight_1*isoWeight_1", "lepton_sf"),
                Weight("(gen_match_2==5)*1.02+(gen_match_2!=5)", "emb_tau_id"),
                Weight("embeddedDecayModeWeight", "decayMode_SF"))
        elif self.channel.name == "tt":
            return Weights(
                Weight("generatorWeight*(generatorWeight<=1.0)",
                       "simulation_sf"),
                Weight("muonEffTrgWeight", "scale_factor"),
                Weight(self.embedding_stitchingweight(),
                       "2016 stitching weight"),
                Weight(
                    "TriggerDataEfficiencyWeight_1*TriggerDataEfficiencyWeight_2*doubleTauTrgWeight",
                    "trg_sf"),
                Weight("((gen_match_1==5)*1.02+(gen_match_1!=5))*((gen_match_2==5)*1.02+(gen_match_2!=5))", "emb_tau_id"),
                Weight("embeddedDecayModeWeight",
                    "decayMode_SF"))
        elif self.channel.name == "em":
            return Weights(
                Weight("generatorWeight", "simulation_sf"),
                Weight("muonEffTrgWeight", "scale_factor"),
                # no trigger sf yet
                Weight("idWeight_1*isoWeight_1*idWeight_2*isoWeight_2",
                       "leptopn_sf"))

    def get_files(self):
        query = {"process": "Embedding2016(B|C|D|E|F|G|H)", "embedded": True}
        if self.channel.name == "mt":
            query["campaign"] = "MuTauFinalState"
            query["scenario"] = ".*v2"
        elif self.channel.name == "et":
            query["campaign"] = "ElTauFinalState"
            query["scenario"] = ".*v2"
        elif self.channel.name == "tt":
            query["campaign"] = "TauTauFinalState"
            query["scenario"] = ".*(v2|v3)"
        elif self.channel.name == "em":
            query["campaign"] = "ElMuFinalState"
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

    def get_cuts(self):
        ztt_genmatch_cut = Cut("1 == 1", "ztt_genmatch")
        if self.channel.name in ["mt", "et"]:
            ztt_genmatch_cut = Cut("gen_match_2==5", "ztt_genmatch")
        elif self.channel.name == "tt":
            ztt_genmatch_cut = Cut("(gen_match_1==5) && (gen_match_2==5)",
                                   "ztt_genmatch")
        elif self.channel.name == "em":
            ztt_genmatch_cut = Cut("(gen_match_1>2) && (gen_match_2>3)",
                                   "ztt_genmatch")
        return Cuts(ztt_genmatch_cut)


class ZLLEstimation(ZTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(ZTTEstimation, self).__init__(
            name="ZLL",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("(gen_match_2<5||gen_match_2==6)", "zll_genmatch_mt"))

    def get_weights(self):
        ztt_weights = super(ZLLEstimation, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*1.0) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.0) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


class ZLLEstimationMTSM(ZLLEstimation):
    def get_weights(self):
        ztt_weights = super(ZLLEstimation, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*0.75) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.0) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


class ZLLEstimationETSM(ZLLEstimation):
    def get_weights(self):
        ztt_weights = super(ZLLEstimation, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*0.98) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.2) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


class ZLEstimationMT(ZTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(ZTTEstimation, self).__init__(
            name="ZL",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("gen_match_2<5", "zl_genmatch_mt"))


class ZLEstimationMTSM(ZLEstimationMT):
    def get_weights(self):
        ztt_weights = super(ZLEstimationMT, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*0.75) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.0) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


class ZLEstimationETSM(ZLEstimationMT):
    def get_weights(self):
        ztt_weights = super(ZLEstimationMT, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*0.98) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.2) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


class ZLEstimationLL(ZTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(ZTTEstimation, self).__init__(
            name="ZL",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(
            Cut("(gen_match_1==1||gen_match_1==2)&&(gen_match_1==gen_match_2)",
                "zl_genmatch_ll"))


class ZJEstimationMT(ZTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(ZTTEstimation, self).__init__(
            name="ZJ",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("gen_match_2==6", "zj_genmatch_mt"))


# et is equivalent to mt
class ZJEstimationET(ZJEstimationMT):
    pass


class ZJEstimationLL(ZTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(ZTTEstimation, self).__init__(
            name="ZJ",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(
            Cut(
                "!(((gen_match_1==3||gen_match_1==4)&&(gen_match_2==3||gen_match_2==4))||((gen_match_1==1||gen_match_1==2)&&(gen_match_1==gen_match_2)))",
                "zj_genmatch_ll"))


class ZLEstimationET(ZLEstimationMT):
    pass


class ZJEstimationTT(ZJEstimationMT):
    def get_cuts(self):
        return Cuts(
            Cut("(gen_match_2 == 6 || gen_match_1 == 6)", "zj_genmatch_tt"))


class ZLEstimationTT(ZLEstimationMT):
    def get_cuts(self):
        return Cuts(
            Cut(
                "(gen_match_1<6&&gen_match_2<6&&!(gen_match_1==5&&gen_match_2==5))",
                "zl_genmatch_tt"))


class EWKWpEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(EWKWpEstimation, self).__init__(
            name="EWKWp",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight(
                "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
                "hadronic_tau_sf"), Weight("eventWeight", "eventWeight"),
            Weight(
                "(5.190747826298e-6)/(numberGeneratedEventsWeight*crossSectionPerEventWeight)",
                "EWKWp_stitching_weight"), self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^EWKWPlus",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, query, files)
        return self.artus_file_names(files)


class EWKWmEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(EWKWmEstimation, self).__init__(
            name="EWKWm",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight(
                "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
                "hadronic_tau_sf"), Weight("eventWeight", "eventWeight"),
            Weight(
                "(4.200367267668e-6)/(numberGeneratedEventsWeight*crossSectionPerEventWeight)",
                "EWKW_stitching_weight"), self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^EWKWMinus",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, query, files)
        return self.artus_file_names(files)


class WEstimationRaw(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(WEstimationRaw, self).__init__(
            name="W",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight(
                "(((npartons == 0 || npartons >= 5)*7.09390278348407e-4) + ((npartons == 1)*1.90063898596475e-4) + ((npartons == 2)*5.8529964471165e-5) + ((npartons == 3)*1.9206444928444e-5) + ((npartons == 4)*1.923548021385e-5))/(numberGeneratedEventsWeight*crossSectionPerEventWeight*sampleStitchingWeight)",
                "wj_stitching_weight"),
            Weight(
                "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
                "hadronic_tau_sf"), Weight("eventWeight", "eventWeight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "W.*JetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class WEstimation(SumUpEstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(WEstimation, self).__init__(
            name="W",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            processes=[
                Process(
                    "W",
                    WEstimationRaw(
                        era,
                        directory,
                        channel,
                        friend_directory=friend_directory)),
                Process(
                    "EWKWp",
                    EWKWpEstimation(
                        era,
                        directory,
                        channel,
                        friend_directory=friend_directory)),
                Process(
                    "EWKWm",
                    EWKWmEstimation(
                        era,
                        directory,
                        channel,
                        friend_directory=friend_directory))
            ])


class WTEstimation(WEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(WEstimation, self).__init__(
            name="W",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("gen_match_1==3||gen_match_1==4", "wt_genmatch"))


class WLEstimation(WEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(WEstimation, self).__init__(
            name="W",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("!(gen_match_1==3||gen_match_1==4)", "wl_genmatch"))


class TTEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight("topPtReweightWeight", "topPtReweightWeight"),
            Weight("eventWeight", "eventWeight"),
            Weight(
                "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
                "hadronic_tau_sf"), self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^TT$",
            "data": False,
            "campaign": self._mc_campaign
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class TTTEstimationMT(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("gen_match_2==5", "ttt_genmatch_mt"))


class TTLEstimationMT(TTEstimation):
    # L refering to a prompt t-quark to lepton decay as opposed to t->tau->lepton (important for embedded events)
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTL",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(
            Cut("gen_match_2==5", "genmatch"),
            Cut("!((gen_match_1==4)&&(gen_match_2==5))",
                "ttbar->tau tau veto for embedded events"))


class TTTTEstimationMT(TTEstimation):
    # true tt->tautau
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(
            Cut("((gen_match_1 == 4) && (gen_match_2 == 5))",
                "select ttbar->tau tau events"))


class TTJEstimationMT(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTJ",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("gen_match_2!=5", "ttj_genmatch_mt"))


class TTTEstimationET(TTTEstimationMT):
    pass


class TTLEstimationET(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTL",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(
            Cut("gen_match_2==5", "genmatch"),
            Cut("!((gen_match_1==3)&&(gen_match_2==5))",
                "ttbar->tau tau veto for embedded events"))


class TTTTEstimationET(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(
            Cut("((gen_match_1 == 3) && (gen_match_2 == 5))",
                "select ttbar->tau tau events"))


class TTJEstimationET(TTJEstimationMT):
    pass


class TTTEstimationTT(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(
            Cut("((gen_match_1 == 5) && (gen_match_2 == 5))",
                "select ttbar->tau tau events"))


class TTLEstimationTT(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTL",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(
            Cut(
                "((gen_match_1 == 5) && (gen_match_2 == 5))*!((gen_match_1 == 5) && (gen_match_2 == 5))",
                "Empty Process")
        )  # All ttbar->real tau events are vetoed for embedded events


class TTJEstimationTT(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(
            Cut("!((gen_match_1==5)&&(gen_match_2==5))",
                "ttbar->tau tau veto"))


class EWKZllEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(EWKZllEstimation, self).__init__(
            name="EWKZll",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight(
                "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
                "hadronic_tau_sf"),
            Weight(
                "(3.989190065346e-6)/(numberGeneratedEventsWeight*crossSectionPerEventWeight)",
                "EWKZll_stitching_weight"),
            Weight("eventWeight", "eventWeight"), self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^EWKZ2Jets.",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, query, files)
        return self.artus_file_names(files)


class EWKZnnEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(EWKZnnEstimation, self).__init__(
            name="EWKZnn",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight(
                "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
                "hadronic_tau_sf"),
            Weight(
                "(3.35561920393e-6)/(numberGeneratedEventsWeight*crossSectionPerEventWeight)",
                "EWKZnn_stitching_weight"),
            Weight("eventWeight", "eventWeight"), self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^EWKZ2Jets$",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, query, files)
        return self.artus_file_names(files)


class EWKZEstimation(SumUpEstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(EWKZEstimation, self).__init__(
            name="EWKZ",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            processes=[
                Process(
                    "EWKZll",
                    EWKZllEstimation(
                        era,
                        directory,
                        channel,
                        friend_directory=friend_directory)),
                Process(
                    "EWKZnn",
                    EWKZnnEstimation(
                        era,
                        directory,
                        channel,
                        friend_directory=friend_directory))
            ])


class VVEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(VVEstimation, self).__init__(
            name="VV",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_weights(self):
        return Weights(
            Weight(
                "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
                "hadronic_tau_sf"), Weight("eventWeight", "eventWeight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process":
            "(WWTo1L1Nu2Q|" + "WZTo1L1Nu2Q|" + "WZTo1L3Nu|" + "WZTo2L2Q|" +
            "ZZTo2L2Q" + ")",
            "data":
            False,
            "campaign":
            self._mc_campaign,
            "generator":
            "amcatnlo-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process": "(VVTo2L2Nu|ZZTo4L)",
            "extension": "ext1",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo-pythia8"
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process": "WZJToLLLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "pythia8"
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process":
            "(STt-channelantitop4finclusiveDecays|STt-channeltop4finclusiveDecays|STtWantitop5finclusiveDecays|STtWtop5finclusiveDecays)",
            "data":
            False,
            "campaign":
            self._mc_campaign
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, "<optimzed out>", files)
        return self.artus_file_names(files)


class VVTEstimationLT(VVEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(VVEstimation, self).__init__(
            name="VVT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("gen_match_2==5", "vvt_genmatch_lt"))


class VVJEstimationLT(VVEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(VVEstimation, self).__init__(
            name="VVJ",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("gen_match_2!=5", "vvj_genmatch_lt"))


class VVTEstimationTT(VVEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(VVEstimation, self).__init__(
            name="VVT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(
            Cut("((gen_match_1 == 5) && (gen_match_2 == 5))",
                "vvt_genmatch_tt"))


class VVJEstimationTT(VVEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(VVEstimation, self).__init__(
            name="VVJ",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv2")

    def get_cuts(self):
        return Cuts(
            Cut("!((gen_match_1==5)&&(gen_match_2==5))", "vvj_genmatch_tt"))


class QCDEstimationET(SStoOSEstimationMethod):
    def __init__(self,
                 era,
                 directory,
                 channel,
                 bg_processes,
                 data_process,
                 friend_directory=None,
                 extrapolation_factor=1.0):
        super(QCDEstimationET, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process,
            extrapolation_factor=extrapolation_factor)


class QCDEstimationMT(QCDEstimationET):
    pass


class QCDEstimationTT(ABCDEstimationMethod):
    def __init__(self,
                 era,
                 directory,
                 channel,
                 bg_processes,
                 data_process,
                 friend_directory=None):
        super(QCDEstimationTT, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process,
            AC_cut_names=
            [  # cuts to be removed to include region for shape derivation
                "tau_2_iso"
            ],
            BD_cuts=
            [  # cuts to be applied to restrict to region for shape derivation
                Cut("byTightIsolationMVArun2v1DBoldDMwLT_2<0.5", "tau_2_iso"),
                Cut("byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5",
                    "tau_2_iso_loose")
            ],
            AB_cut_names=
            [  # cuts to be removed to include region for the determination of the extrapolation derivation
                "os"
            ],
            CD_cuts=
            [  # cuts to be applied to restrict to region for the determination of the extrapolation derivation
                Cut("q_1*q_2>0", "ss")
            ])


class WEstimationWithQCD(EstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process,
                 w_process, qcd_ss_to_os_extrapolation_factor, friend_directory=None):
        super(WEstimationWithQCD, self).__init__(
            name="WJets",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign=None)
        self._bg_processes = bg_processes
        self._data_process = data_process
        self._w_process = w_process
        self._qcd_ss_to_os_extrapolation_factor = qcd_ss_to_os_extrapolation_factor

    def create_root_objects(self, systematic):

        # create category for MC WJets shape estimation in the signal region
        signal_region = copy.deepcopy(systematic.category)
        signal_region.name = (signal_region.name +
                              "_for_wjets_mc").lstrip(self.channel.name + "_")

        # create control regions for W yield estimation
        high_mt_ss_control_region = copy.deepcopy(systematic.category)
        high_mt_ss_control_region.name = "wjets_high_mt_ss_cr"
        high_mt_ss_control_region._variable = None

        high_mt_ss_control_region.cuts.remove("mt")
        high_mt_ss_control_region.cuts.remove("os")

        high_mt_ss_control_region.cuts.add(Cut("mt_1>70", "mt"))
        high_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W high mt to low mt extrapolation factor
        high_mt_os_control_region = copy.deepcopy(
            systematic.category)  # this one also used for W yield estimation
        high_mt_os_control_region.name = "wjets_high_mt_os_cr"
        high_mt_os_control_region._variable = None

        high_mt_os_control_region.cuts.remove("mt")

        high_mt_os_control_region.cuts.add(Cut("mt_1>70", "mt"))

        low_mt_os_control_region = copy.deepcopy(systematic.category)
        low_mt_os_control_region.name = "wjets_low_mt_os_cr"
        low_mt_os_control_region._variable = None

        low_mt_ss_control_region = copy.deepcopy(systematic.category)
        low_mt_ss_control_region.name = "wjets_low_mt_ss_cr"
        low_mt_ss_control_region._variable = None

        low_mt_ss_control_region.cuts.remove("os")

        low_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W ss to os extrapolation factor
        inclusive_os_control_region = copy.deepcopy(systematic.category)
        inclusive_os_control_region.name = "wjets_os_cr"
        inclusive_os_control_region._variable = None

        inclusive_os_control_region.cuts.remove("mt")

        inclusive_ss_control_region = copy.deepcopy(systematic.category)
        inclusive_ss_control_region.name = "wjets_ss_cr"
        inclusive_ss_control_region._variable = None

        inclusive_ss_control_region.cuts.remove("mt")
        inclusive_ss_control_region.cuts.remove("os")

        inclusive_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # initialize root objects and systematics
        root_objects = []
        systematic._WandQCD_systematics = []

        # for extrapolation factors: only W MC is needed
        for category in [
                high_mt_os_control_region, low_mt_os_control_region,
                high_mt_ss_control_region, low_mt_ss_control_region,
                inclusive_os_control_region, inclusive_ss_control_region
        ]:
            s = Systematic(
                category=category,
                process=self._w_process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125)
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects

        # for yields in high mt control regions: data and other bg processes needed
        for process in [self._data_process] + self._bg_processes:
            for category in [
                    high_mt_os_control_region, high_mt_ss_control_region
            ]:
                s = Systematic(
                    category=category,
                    process=process,
                    analysis=systematic.analysis,
                    era=self.era,
                    variation=systematic.variation,
                    mass=125)
                systematic._WandQCD_systematics.append(s)
                s.create_root_objects()
                root_objects += s.root_objects

        # for signal region shape
        s = Systematic(
            category=signal_region,
            process=self._w_process,
            analysis=systematic.analysis,
            era=self.era,
            variation=systematic.variation,
            mass=125)
        systematic._WandQCD_systematics.append(s)
        s.create_root_objects()
        root_objects += s.root_objects

        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_WandQCD_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _WandQCD_systematics needed for WandQCD estimation.",
                systematic.name)
            raise Exception

        # Sort shapes and counts
        wjets_mc_shape = None
        wjets_high_mt_ss_cr_counts = {}
        wjets_high_mt_os_cr_counts = {}
        wjets_low_mt_os_cr_count = None
        wjets_low_mt_ss_cr_count = None
        wjets_os_cr_count = None
        wjets_ss_cr_count = None
        for s in systematic._WandQCD_systematics:
            s.do_estimation()
            if s.category.name.endswith("for_wjets_mc"):
                wjets_mc_shape = s.shape
            elif s.category.name.endswith("wjets_high_mt_ss_cr"):
                wjets_high_mt_ss_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_os_cr"):
                wjets_high_mt_os_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_low_mt_os_cr"):
                wjets_low_mt_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_low_mt_ss_cr"):
                wjets_low_mt_ss_cr_count = s.shape
            elif s.category.name.endswith("wjets_os_cr"):
                wjets_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_ss_cr"):
                wjets_ss_cr_count = s.shape

        # Determine extrapolation factors
        R_ss_to_os = wjets_os_cr_count.result / wjets_ss_cr_count.result

        wjets_integral_low_mt_os = wjets_low_mt_os_cr_count.result
        wjets_integral_high_mt_os = wjets_high_mt_os_cr_counts.pop(
            self._w_process.name).result
        logger.debug("Integral of WJets MC in low mt OS region: %s",
                     str(wjets_integral_low_mt_os))
        logger.debug("Integral of WJets MC in high mt OS region: %s",
                     str(wjets_integral_high_mt_os))

        R_high_to_low_mt_os = wjets_integral_low_mt_os / wjets_integral_high_mt_os
        R_high_to_low_mt_ss = wjets_low_mt_ss_cr_count.result / wjets_high_mt_ss_cr_counts.pop(
            self._w_process.name).result
        logger.debug("WJets SS to OS extrapolation factor: %s",
                     str(R_ss_to_os))
        logger.debug("WJets high to low mt os extrapolation factor: %s",
                     str(R_high_to_low_mt_os))
        logger.debug("WJets high to low mt ss extrapolation factor: %s",
                     str(R_high_to_low_mt_ss))

        # Determine yields in wjets CRs
        logger.debug(
            "Data yield in ss high mt region: %s",
            str(wjets_high_mt_ss_cr_counts[self._data_process.name].result))
        high_mt_ss_yield = wjets_high_mt_ss_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_ss_cr_counts.values()])
        sum_mc = sum([s.result for s in wjets_high_mt_ss_cr_counts.values()])
        logger.debug("MC yield to be subtracted: %s", str(sum_mc))
        for name, s in wjets_high_mt_ss_cr_counts.items():
            logger.debug(name + " : " + str(s.result / sum_mc))

        logger.debug(
            "Data yield in os high mt region: %s",
            str(wjets_high_mt_os_cr_counts[self._data_process.name].result))
        high_mt_os_yield = wjets_high_mt_os_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_os_cr_counts.values()])
        sum_mc = sum([s.result for s in wjets_high_mt_os_cr_counts.values()])
        logger.debug("MC yield to be subtracted: %s", str(sum_mc))
        for name, s in wjets_high_mt_os_cr_counts.items():
            logger.debug(name + " : " + str(s.result / sum_mc))

        logger.debug("WJets + QCD yield in ss high mt region: %s",
                     str(high_mt_ss_yield))
        logger.debug("WJets + QCD yield in os high mt region: %s",
                     str(high_mt_os_yield))

        # Derive and normalize final shape
        logger.debug("WJets MC yield in signal region: %s",
                     str(wjets_integral_low_mt_os))
        sf = R_ss_to_os * (
            high_mt_os_yield -
            self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor
            ) / wjets_integral_high_mt_os
        estimated_yield = R_high_to_low_mt_os * R_ss_to_os * (
            high_mt_os_yield -
            self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor)
        logger.debug("WJets Estimated yield in signal region: %s",
                     str(estimated_yield))
        logger.debug("Scale WJets by %s", str(sf))
        wjets_shape = copy.deepcopy(wjets_mc_shape)
        wjets_shape.result.Scale(sf)

        # Rename root object accordingly
        wjets_shape.name = systematic.name

        # Replace negative entries by zeros and renormalize shape
        wjets_shape.replace_negative_entries_and_renormalize(tolerance=100.5)

        return wjets_shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError


class QCDEstimationWithW(EstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process,
                 w_process, qcd_ss_to_os_extrapolation_factor, friend_directory=None):
        super(QCDEstimationWithW, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign=None)
        self._bg_processes = bg_processes
        self._data_process = data_process
        self._w_process = w_process
        self._qcd_ss_to_os_extrapolation_factor = qcd_ss_to_os_extrapolation_factor

    def create_root_objects(self, systematic):

        # create category for WJets and QCD shape estimation in the qcd control region
        qcd_control_region = copy.deepcopy(systematic.category)
        qcd_control_region.name = (qcd_control_region.name + "_ss_for_qcd"
                                   ).lstrip(self.channel.name + "_")

        qcd_control_region.cuts.remove("os")

        qcd_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W yield estimation
        high_mt_ss_control_region = copy.deepcopy(systematic.category)
        high_mt_ss_control_region.name = "wjets_high_mt_ss_cr"
        high_mt_ss_control_region._variable = None

        high_mt_ss_control_region.cuts.remove("mt")
        high_mt_ss_control_region.cuts.remove("os")

        high_mt_ss_control_region.cuts.add(Cut("mt_1>70", "mt"))
        high_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W high mt to low mt extrapolation factor
        high_mt_os_control_region = copy.deepcopy(
            systematic.category)  # this one also used for W yield estimation
        high_mt_os_control_region.name = "wjets_high_mt_os_cr"
        high_mt_os_control_region._variable = None

        high_mt_os_control_region.cuts.remove("mt")

        high_mt_os_control_region.cuts.add(Cut("mt_1>70", "mt"))

        low_mt_os_control_region = copy.deepcopy(systematic.category)
        low_mt_os_control_region.name = "wjets_low_mt_os_cr"
        low_mt_os_control_region._variable = None

        low_mt_ss_control_region = copy.deepcopy(systematic.category)
        low_mt_ss_control_region.name = "wjets_low_mt_ss_cr"
        low_mt_ss_control_region._variable = None

        low_mt_ss_control_region.cuts.remove("os")

        low_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W ss to os extrapolation factor
        inclusive_os_control_region = copy.deepcopy(systematic.category)
        inclusive_os_control_region.name = "wjets_os_cr"
        inclusive_os_control_region._variable = None

        inclusive_os_control_region.cuts.remove("mt")

        inclusive_ss_control_region = copy.deepcopy(systematic.category)
        inclusive_ss_control_region.name = "wjets_ss_cr"
        inclusive_ss_control_region._variable = None

        inclusive_ss_control_region.cuts.remove("mt")
        inclusive_ss_control_region.cuts.remove("os")

        inclusive_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # initialize root objects and systematics
        root_objects = []
        systematic._WandQCD_systematics = []

        # for extrapolation factors: only W MC is needed
        for category in [
                high_mt_os_control_region, low_mt_os_control_region,
                high_mt_ss_control_region, low_mt_ss_control_region,
                inclusive_os_control_region, inclusive_ss_control_region
        ]:
            s = Systematic(
                category=category,
                process=self._w_process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125)
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects

        # for yields in high mt control regions: data and other bg processes needed
        for process in [self._data_process] + self._bg_processes:
            for category in [
                    high_mt_os_control_region, high_mt_ss_control_region
            ]:
                s = Systematic(
                    category=category,
                    process=process,
                    analysis=systematic.analysis,
                    era=self.era,
                    variation=systematic.variation,
                    mass=125)
                systematic._WandQCD_systematics.append(s)
                s.create_root_objects()
                root_objects += s.root_objects

        # for Wjets and QCD shape
        for process in [self._data_process, self._w_process
                        ] + self._bg_processes:
            s = Systematic(
                category=qcd_control_region,
                process=process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125)
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects

        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_WandQCD_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _WandQCD_systematics needed for WandQCD estimation.",
                systematic.name)
            raise Exception

        # Sort shapes and counts
        qcd_control_region_shapes = {}
        wjets_high_mt_ss_cr_counts = {}
        wjets_high_mt_os_cr_counts = {}
        wjets_low_mt_os_cr_count = None
        wjets_low_mt_ss_cr_count = None
        wjets_os_cr_count = None
        wjets_ss_cr_count = None
        for s in systematic._WandQCD_systematics:
            s.do_estimation()
            if s.category.name.endswith("ss_for_qcd"):
                qcd_control_region_shapes[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_ss_cr"):
                wjets_high_mt_ss_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_os_cr"):
                wjets_high_mt_os_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_low_mt_os_cr"):
                wjets_low_mt_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_low_mt_ss_cr"):
                wjets_low_mt_ss_cr_count = s.shape
            elif s.category.name.endswith("wjets_os_cr"):
                wjets_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_ss_cr"):
                wjets_ss_cr_count = s.shape

        # Determine extrapolation factors
        R_ss_to_os = wjets_os_cr_count.result / wjets_ss_cr_count.result

        wjets_integral_low_mt_ss = wjets_low_mt_ss_cr_count.result
        wjets_integral_high_mt_ss = wjets_high_mt_ss_cr_counts.pop(
            self._w_process.name).result

        R_high_to_low_mt_os = wjets_low_mt_os_cr_count.result / wjets_high_mt_os_cr_counts.pop(
            self._w_process.name).result
        R_high_to_low_mt_ss = wjets_integral_low_mt_ss / wjets_integral_high_mt_ss

        # Determine yields in wjets CRs
        high_mt_ss_yield = wjets_high_mt_ss_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_ss_cr_counts.values()])

        high_mt_os_yield = wjets_high_mt_os_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_os_cr_counts.values()])

        # Derive and normalize final shape for QCD
        wjets_shape = qcd_control_region_shapes.pop(self._w_process.name)
        logger.debug("WJets MC yield in qcd control region: %s",
                     str(wjets_integral_low_mt_ss))
        sf = (high_mt_os_yield -
              self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                  R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor
              ) / wjets_integral_high_mt_ss
        estimated_yield = R_high_to_low_mt_ss * (
            high_mt_os_yield -
            self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor)
        logger.debug("WJets Estimated yield in qcd control region: %s",
                     str(estimated_yield))
        logger.debug("Scale WJets by %s", str(sf))
        wjets_shape.result.Scale(sf)
        wjets_shape._result.Write()

        qcd_shape = copy.deepcopy(
            qcd_control_region_shapes.pop(self._data_process.name))
        qcd_shape.result.Add(wjets_shape.result, -1.0)
        for sh in qcd_control_region_shapes.values():
            qcd_shape.result.Add(sh.result, -1.0)
        # Saving QCD shape in ss control region
        qcd_ss_shape = copy.deepcopy(qcd_shape)
        ss_category_name = ""
        for s in systematic._WandQCD_systematics:
            if s.category.name.endswith("ss_for_qcd"):
                ss_category_name = s.category._name
        qcd_ss_shape.name = systematic.name.replace(systematic.category._name,
                                                    ss_category_name)
        qcd_ss_shape._result.Write()

        # Rescale QCD shape for signal region
        qcd_shape.result.Scale(self._qcd_ss_to_os_extrapolation_factor)

        # Rename root object accordingly
        qcd_shape.name = systematic.name

        # Replace negative entries by zeros and renormalize shape
        qcd_shape.replace_negative_entries_and_renormalize(tolerance=100.5)

        return qcd_shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError
