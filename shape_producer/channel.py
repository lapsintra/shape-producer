# -*- coding: utf-8 -*-

from cutstring import Cuts, Cut
import logging
logger = logging.getLogger(__name__)
"""
"""


class Channel(object):
    pass


class MT(Channel):
    def __init__(self):
        self._name = "mt"

    @property
    def cuts(self):
        return Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("againstMuonTight3_2>0.5", "againstMuonTight"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstElectronVLooseMVA6_2>0.5", "againstElectronVeto"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"),
            Cut("q_1*q_2<0", "os"), Cut("trg_singlemuon==1", "trg_singlemuon"))

    @property
    def name(self):
        return self._name


class ET(MT):
    def __init__(self):
        self._name = "et"

    @property
    def cuts(self):
        return Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonTight"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstElectronTightMVA6_2>0.5", "againstElectronVeto"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_iso"),
            Cut("iso_1<0.1", "ele_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("trg_singleelectron==1", "trg_singleelectron"))

    @property
    def name(self):
        return self._name


# collection of channels an analysis can be ran on
class Channels(object):
    def __init__(self, name):
        self._name = name
        self._channels = []

    def add(self, channel):
        self._channels.append(channel)

    @property
    def name(self):
        return self._name
