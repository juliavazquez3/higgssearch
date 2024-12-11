import sys
import ROOT
import os
from ROOT import *
import json
import argparse
import numpy as np
from os.path import isfile, join, isdir

#if not sys.flags.interactive: ROOT.EnableImplicitMT()

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--stack", action="store_true", default=False,
                    help="Stack simulation or not")
parser.add_argument("--norm", action="store_true", default=False,
                    help="normalising test")
parser.add_argument("--nodata", action="store_true", default=False,
                    help="Do not plot data")
parser.add_argument("--wcs", action="store_true", default=False,
                    help="MC classififcation")
parser.add_argument("--mudhad", action="store_true", default=False,
                    help="MC classififcation fot charm signal sample")
parser.add_argument("--year", type=string, default="2016",
                    help="Select year of process to run")
parser.add_argument("--channel", type=string, default="btagMM",
                    help="Select year of process to run")
parser.add_argument("--sumEM", action="store_true", default=False,
                    help="union of M and E channels")

# Use like:
# python higgssearch/fromJF/hist_fromJFwqq_syst_total.py --stack --ratio --png --norm --year="all" --channel="btagMM_chitest_sl"

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')


if args.channel == "btagMM_chitest": term_path = ""
elif args.channel == "btagMM_chitest_ctag": term_path = "/ctag"
elif args.channel == "btagMM_chitest_noctag": term_path = "/noctag"
elif args.channel == "btagMM_chitest_sl": term_path = "/sl"
elif args.channel == "btagMM_chitest_antisl": term_path = "/antisl"
elif args.channel == "btagMM_chitest_slss": term_path = "/sl/ss"
elif args.channel == "btagMM_chitest_slos": term_path = "/sl/os"
elif args.channel == "btagMM_chitest_slssos": term_path = "/sl/ssos"
else: raise NameError('Incorrect data option')

sl_channel = ["btagMM_chitest_sl","lepton50_chitest_sl","btagMM_chitest_slss","lepton50_chitest_slss",
      "btagMM_chitest_slos","lepton50_chitest_slos","btagMM_chitest_slssos","lepton50_chitest_slssos"]

plotdir = '/nfs/cms/vazqueze/higgssearch/plotspng/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

c_rat = 1.5
c_rat2 = 0.5

norm_factorM = {}
norm_factorE = {}

norm_factorM["all"] = {}; norm_factorM["2016"] = {}; norm_factorM["2016B"] = {}; norm_factorM["2017"] = {}; norm_factorM["2018"] = {};
norm_factorE["all"] = {}; norm_factorE["2016"] = {}; norm_factorE["2016B"] = {}; norm_factorE["2017"] = {}; norm_factorE["2018"] = {};

### btagMM
norm_factorM["all"]["btagMM"] = 0.95;
norm_factorE["all"]["btagMM"] = 0.94;

### btagMM smeared, chitest
norm_factorM["all"]["btagMM_chitest"] = 0.915;
norm_factorE["all"]["btagMM_chitest"] = 0.894;

norm_factorM["2016"]["btagMM_chitest"] = 0.89;norm_factorE["2016"]["btagMM_chitest"] = 0.87;
norm_factorM["2016B"]["btagMM_chitest"] = 0.92;norm_factorE["2016B"]["btagMM_chitest"] = 0.94;
norm_factorM["2017"]["btagMM_chitest"] = 0.92;norm_factorE["2017"]["btagMM_chitest"] = 0.88;
norm_factorM["2018"]["btagMM_chitest"] = 0.91;norm_factorE["2018"]["btagMM_chitest"] = 0.90;

### btagMM smeared, chitest, SL
norm_factorM["all"]["btagMM_chitest_sl"] = 0.92; norm_factorM["all"]["btagMM_chitest_slss"] = 0.92; 
norm_factorM["all"]["btagMM_chitest_slos"] = 0.92; norm_factorM["all"]["btagMM_chitest_slssos"] = 0.92;
norm_factorE["all"]["btagMM_chitest_sl"] = 0.90; norm_factorE["all"]["btagMM_chitest_slss"] = 0.90; 
norm_factorE["all"]["btagMM_chitest_slos"] = 0.90; norm_factorE["all"]["btagMM_chitest_slssos"] = 0.90;

norm_factorM["all"]["btagMM_chitest_ctag"] = 0.92; norm_factorM["all"]["btagMM_chitest_noctag"] = 0.92; 
norm_factorE["all"]["btagMM_chitest_ctag"] = 0.90; norm_factorE["all"]["btagMM_chitest_noctag"] = 0.90; 
norm_factorM["all"]["btagMM_chitest_antisl"] = 0.92; norm_factorE["all"]["btagMM_chitest_antisl"] = 0.90; 

observable_names = ["InvM_2jets","jet_1_pt", "jet_1_nmu", "jet_1_eta", "jet_2_pt", "jet_2_eta", "jet_2_mass", "jet_2_qgl","jet_2_nmu","jet_1_qgl",
   "lepton_pt", "lepton_eta", "lepton_pt_detail", "lepton_eta_thick", "InvM_bot_closer", "InvM_bot_farther",
   "deltaR_jet1_jet2", "deltaphi_jet1_jet2", "deltaeta_jet1_jet2", "MET_pt_aux", "MET_sig", "MET_my_sig",
   "transverse_mass", "tracks_jet1", "tracks_jet2", "deltaphi_MET_jets_1", "deltaphi_MET_jets_2", "pT_Wlep",
   "jet_1_btag", "jet_2_btag", "deltaphi_MET_lep", "jet_bot1_btag", "jet_bot2_btag", "jet_bot1_pt", "jet_bot1_eta", "jet_bot2_pt", "jet_bot2_eta",
   "jet_bot1_btag_thick", "jet_bot2_btag_thick", "jet_1_btag_thick", "jet_2_btag_thick",
   "jet_bot1_btagnumber", "jet_bot2_btagnumber", "jet_1_btagnumber", "jet_2_btagnumber",
   "jet_1_cvltag_csv", "jet_2_cvltag_csv", "jet_1_cvltag", "jet_2_cvltag","InvM30","InvM31","InvMl0","InvMl1","chi2_test0","chi2_test1",
   "InvM3_good","InvM3_bad","InvMl_good","InvMl_bad","chi2_test_good","chi2_test_bad","jet_max_cvltag","jet_min_cvltag",
   "jet_1_cvbtag_csv", "jet_2_cvbtag_csv", "jet_1_cvbtag", "jet_2_cvbtag", "jet_max_cvbtag", "jet_min_cvbtag",
   "jet_1_eta_thick","jet_2_eta_thick","jet_bot1_eta_thick","jet_bot2_eta_thick",
   "InvM_2jets_thick","InvM_2jets_short","bot1_muons","bot2_muons","muon_bot1_eta","muon_bot2_eta","muon_bot1_pt","muon_bot2_pt","nJetGood",
   "jet_bot1_tracks","jet_bot2_tracks","tau_discr_jet1","tau_discr_jet2","tau_discr_jetbot1","tau_discr_jetbot2"]

if "sl" in str(args.channel): observable_names = observable_names + ["muon_jet_pt","muon_jet_z","muon_jet_eta","muon_jet_pt_rel","muon_jet_iso",
         "muon_jet_iso_log","muon_jet_z_short","InvM3_good_short","muon_jet_sigr","muon_jet_sigxy","muon_jet_sigdz","muon_jet_r","deltaR_jet1_muon",
         "muon_jet_z2","muon_jet_z3"]

not_rebin = ["nJetGood","InvM3_good","InvM3_bad","InvMl_good","InvMl_bad","lepton_eta_thick","jet_bot1_btagnumber", "jet_bot2_btagnumber", 
      "jet_1_btagnumber", "jet_2_btagnumber","muon_jet_pt","muon_jet_relpt","muon_jet_eta","tau_discr","transverse_mass",
      "jet_1_flavourP", "jet_2_flavourP", "jet_bot1_flavourP", "jet_bot2_flavourP","lepton_pt", "muon_jet_z","tau_discr_jet1","tau_discr_jet2","tau_discr_jetbot1",
      "tau_discr_jetbot2","muon_jet_iso","InvM3_good_short"]

observable_names = ["InvM_2jets_short","InvM3_good"]
#observable_names = ["jet_bot1_btagnumber"]

datayears = ["2016","2016B","2017","2018"]
#datayears = ["2018","2016","2016B"]

samplesHT = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_light","ttbar_sl_bottom","ttbar_dl","ttbar_dh","zz","wz",
        "st_1","st_2","st_3","st_4", "ttbar_sl_else", "ttbar_sl_charmgluon","ttbar_sl_bottomgluon"]

if args.wcs:
   samplesHT = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_nocharm","ttbar_dl","ttbar_dh","zz","wz",
        "st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm",
        "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else"]

samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "ttbar_sl","ttbar_dl","ttbar_dh","zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6",
        "zjets_7","zjets_8","st_1","st_2","st_3","st_4","zz","wz"]

samples_d = ["2016","2016B","2017","2018"]

## lumi info
lumi = {}
xsecs = {}
nevents = {}

for data_op in samples_d:
        files = json.load(open("/nfs/cms/vazqueze/higgssearch/mcinfo"+data_op+".json"))
        lumi[data_op] = {}
        for p in samples:
                #num_events = files[p]["events"] # Number of events
                num_files = files[p]["files"] # Number of files
                luminosity = files[p]["lumi"] # Luminosity
                #print(files[p]["type"])
                lumi[data_op][p] = luminosity
        lumi[data_op]["ttbar_sl_charm"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_charm_mudhad"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_nocharm"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_light"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_bottom"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_else"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_charmgluon"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_bottomgluon"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["st_1_charm"] = lumi[data_op]["st_1"]
        lumi[data_op]["st_2_charm"] = lumi[data_op]["st_2"]
        lumi[data_op]["st_3_charm"] = lumi[data_op]["st_3"]
        lumi[data_op]["st_4_charm"] = lumi[data_op]["st_4"]
        lumi[data_op]["st_1_nocharm"] = lumi[data_op]["st_1"]
        lumi[data_op]["st_2_nocharm"] = lumi[data_op]["st_2"]
        lumi[data_op]["st_3_nocharm"] = lumi[data_op]["st_3"]
        lumi[data_op]["st_4_nocharm"] = lumi[data_op]["st_4"]
        lumi[data_op]["st_1_else"] = lumi[data_op]["st_1"]
        lumi[data_op]["st_2_else"] = lumi[data_op]["st_2"]
        lumi[data_op]["st_3_else"] = lumi[data_op]["st_3"]
        lumi[data_op]["st_4_else"] = lumi[data_op]["st_4"]

listsampl = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "ttbar_sl","ttbar_dl","ttbar_dh","zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6",
        "zjets_7","zjets_8","st_1","st_2","st_3","st_4","zz","wz", "ttbar_sl_charm", "ttbar_sl_charmgluon",
        "ttbar_sl_bottom", "ttbar_sl_bottomgluon", "ttbar_sl_else", "ttbar_sl_light", "ttbar_sl_nocharm",
        "st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm",
        "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else","ttbar_sl_charm_mudhad"]

for s in listsampl:
   lumi["2016B"][s] = lumi["2016"][s]

lumi_d = {}
lumi_d["2016"] = 19.5
lumi_d["2016B"] = 16.8
lumi_d["2017"] = 41.5
lumi_d["2018"] = 59.8

#######################################################
########### Start of plot creation ####################
#######################################################

if args.year == "all": datayears = ["2016","2016B","2017","2018"]
elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): datayears = [str(args.year)]
else: raise NameError('Incorrect year')

histFile = {}
histFileDM = {}
histFileDE = {}

myfile = {}
#path_histfile = "/nfs/cms/vazqueze/higgssearch/datacards/"
path_histfile = "/nfs/cms/vazqueze/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit/data/tutorials/longexercise/jec_syst/"

for name in observable_names:
  #### Systematics to take into account
  #list_syst = ["btaglightup", "btaglightdown", "btagheavyup", "btagheavydown", "seclepup", "seclepdown", "lepidup", "lepiddown", 
  #    "lepisoup", "lepisodown", "leptrigup", "leptrigdown", "notoppt","puwup","puwdown"] 
  list_syst = ["jecup","jecdown"]
  #list_syst = ["ctagup","ctagdown"] 
  ## Open hists files
  filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/chi_test"+term_path+"/"
  if args.wcs: filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/chi_test/wcs_classes_jec"+term_path+"/"
  #if args.wcs: filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/chi_test/wcs_classes"+term_path+"/"
  term = "hist_wqqfromJF_"
  end_term = ".root"
  ## mc files
  histFile[name] = {}
  histFileDM[name] = {}
  histFileDE[name] = {}
  for data_op in datayears:
       if isfile(filePath + term+"MC_"+data_op+"_"+name+end_term):
           histFile[name][data_op] = TFile.Open(filePath + term+"MC_"+data_op+"_"+name+end_term,"READ")
       # data files
       if not args.nodata: histFileDM[name][data_op] = TFile.Open(filePath + term+"dataM_"+data_op+"_"+name+".root","READ")
       if not args.nodata: histFileDE[name][data_op] = TFile.Open(filePath + term+"dataE_"+data_op+"_"+name+".root","READ")
  #print(data_op)
  #print(histFile[data_op].keys())
  #print(histFileDM[data_op].keys())

  samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_light","ttbar_sl_bottom","ttbar_dl","ttbar_dh","zz","wz",
        "st_1","st_2","st_3","st_4","ttbar_sl_else","ttbar_sl_charmgluon","ttbar_sl_bottomgluon"]

  if args.wcs:
    if args.mudhad:
       samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
            "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
            "ttbar_sl_charm_mudhad","ttbar_sl_nocharm","ttbar_dl","ttbar_dh","zz","wz",
            "st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm",
            "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else"]
    else:
       samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
            "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
            "ttbar_sl_charm","ttbar_sl_nocharm","ttbar_dl","ttbar_dh","zz","wz",
            "st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm",
            "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else"]

  ## HISTS
  samples_foryear = {}
  hist_nom_M = {}
  hist_nom_E = {}
  hist_syst_M = {}
  hist_syst_E = {}
  for syst_name in list_syst:
     hist_syst_M[syst_name] = {}
     hist_syst_E[syst_name] = {}
  histdata_M = {}
  histdata_E = {}
  for data_op in datayears:
    hist_nom_M[data_op] = {}
    hist_nom_E[data_op] = {}
    for syst_name in list_syst:
       hist_syst_M[syst_name][data_op] = {}
       hist_syst_E[syst_name][data_op] = {}
    data_term = data_op
    #print(data_op+"M_"+name+"_M")
    #print(histFile[data_op][name].ls())
    for s in samples:
      if s[0:8] == "ttbar_sl": 
         s_term = s[0:8]+data_term+s[8:] 
      elif s[0:2] == "st":
         s_term = s[0:4]+data_term+s[4:]
      else:
         s_term = s+data_term
      hist_nom_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M")
      hist_nom_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E")
      for syst_name in list_syst:
         hist_syst_M[syst_name][data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_"+str(syst_name))
         hist_syst_E[syst_name][data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_"+str(syst_name))
    if not args.nodata: histdata_M[data_op] = histFileDM[name][data_op].Get("data"+data_op+"M_"+name+"_M")
    if not args.nodata: histdata_E[data_op] = histFileDE[name][data_op].Get("data"+data_op+"E_"+name+"_E")

  samples_st = {}
  samples_wjets = {}
  samples_zjets = {}
  ## Scaling to lumi
  #print(name) 
  for data_op in datayears:
    lumi_data = lumi_d[data_op]
    for s in samples:
      #print(s)
      #print(data_op)
      hist_nom_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_nom_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
      for syst_name in list_syst:
         #print(syst_name)
         hist_syst_M[syst_name][data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_syst_E[syst_name][data_op][s].Scale(lumi_data/lumi[data_op][s])
      if args.norm: 
         hist_nom_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
         hist_nom_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
         for syst_name in list_syst:
            hist_syst_M[syst_name][data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
            hist_syst_E[syst_name][data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])

    ## Fixing single top
    #print(samples_foryear[data_op]) 
    #### List of summing samples:
    list_st = ["st_1","st_2","st_3","st_4"]
    list_wjets = ["wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8"]
    list_zjets = ["zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8"]
    list_vv = ["ww","wz","zz"]
    if args.wcs:
       list_st_aux = ["_charm","_nocharm","_else"]
    else:
       list_st_aux = [""]
    for apx in list_st_aux:
      hist_nom_M[data_op]["st"+apx] = hist_nom_M[data_op][list_st[0]+apx]
      hist_nom_E[data_op]["st"+apx] = hist_nom_E[data_op][list_st[0]+apx]
      for syst_name in list_syst:
         hist_syst_M[syst_name][data_op]["st"+apx] = hist_syst_M[syst_name][data_op][list_st[0]+apx]
         hist_syst_E[syst_name][data_op]["st"+apx] = hist_syst_E[syst_name][data_op][list_st[0]+apx]
      for l in list_st[1:]:
        hist_nom_M[data_op]["st"+apx].Add(hist_nom_M[data_op][l+apx])
        hist_nom_E[data_op]["st"+apx].Add(hist_nom_E[data_op][l+apx])
        for syst_name in list_syst:
           hist_syst_M[syst_name][data_op]["st"+apx].Add(hist_syst_M[syst_name][data_op][l+apx])
           hist_syst_E[syst_name][data_op]["st"+apx].Add(hist_syst_E[syst_name][data_op][l+apx])
    hist_nom_M[data_op]["wjets"] = hist_nom_M[data_op][list_wjets[0]]
    hist_nom_E[data_op]["wjets"] = hist_nom_E[data_op][list_wjets[0]]
    for syst_name in list_syst:
       hist_syst_M[syst_name][data_op]["wjets"] = hist_syst_M[syst_name][data_op][list_wjets[0]]
       hist_syst_E[syst_name][data_op]["wjets"] = hist_syst_E[syst_name][data_op][list_wjets[0]]
    for l in list_wjets[1:]:
      hist_nom_M[data_op]["wjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["wjets"].Add(hist_nom_E[data_op][l])
      for syst_name in list_syst:
         hist_syst_M[syst_name][data_op]["wjets"].Add(hist_syst_M[syst_name][data_op][l])
         hist_syst_E[syst_name][data_op]["wjets"].Add(hist_syst_E[syst_name][data_op][l])
    hist_nom_M[data_op]["zjets"] = hist_nom_M[data_op][list_zjets[0]]
    hist_nom_E[data_op]["zjets"] = hist_nom_E[data_op][list_zjets[0]]
    for syst_name in list_syst:
       hist_syst_M[syst_name][data_op]["zjets"] = hist_syst_M[syst_name][data_op][list_zjets[0]]
       hist_syst_E[syst_name][data_op]["zjets"] = hist_syst_E[syst_name][data_op][list_zjets[0]]
    for l in list_zjets[1:]:
      hist_nom_M[data_op]["zjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["zjets"].Add(hist_nom_E[data_op][l])
      for syst_name in list_syst:
         hist_syst_M[syst_name][data_op]["zjets"].Add(hist_syst_M[syst_name][data_op][l])
         hist_syst_E[syst_name][data_op]["zjets"].Add(hist_syst_E[syst_name][data_op][l])
    hist_nom_M[data_op]["vv"] = hist_nom_M[data_op][list_vv[0]]
    hist_nom_E[data_op]["vv"] = hist_nom_E[data_op][list_vv[0]]
    for syst_name in list_syst:
       hist_syst_M[syst_name][data_op]["vv"] = hist_syst_M[syst_name][data_op][list_vv[0]]
       hist_syst_E[syst_name][data_op]["vv"] = hist_syst_E[syst_name][data_op][list_vv[0]]
    for l in list_vv[1:]:
      hist_nom_M[data_op]["vv"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["vv"].Add(hist_nom_E[data_op][l])
      for syst_name in list_syst:
         hist_syst_M[syst_name][data_op]["vv"].Add(hist_syst_M[syst_name][data_op][l])
         hist_syst_E[syst_name][data_op]["vv"].Add(hist_syst_E[syst_name][data_op][l])

  samples = ["ttbar_sl_bottom","ttbar_sl_charm","ttbar_sl_else","ttbar_sl_light","ttbar_dl","ttbar_dh","zjets","vv","st","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon"]
  if args.wcs: 
     if args.mudhad:
        samples = ["ttbar_sl_nocharm","ttbar_sl_charm_mudhad","ttbar_dl","ttbar_dh","zjets","vv","st_nocharm","st_charm","st_else","wjets"]
     else:
        samples = ["ttbar_sl_nocharm","ttbar_sl_charm","ttbar_dl","ttbar_dh","zjets","vv","st_nocharm","st_charm","st_else","wjets"]

  ## sumamos todos los anos para todos los histogramas
  histT_nom_M = {}
  histT_nom_E = {}
  histT_syst_M = {}
  histT_syst_E = {}
  for syst_name in list_syst:
     histT_syst_M[syst_name] = {}
     histT_syst_E[syst_name] = {}
  for s in samples:
       histT_nom_M[s] = hist_nom_M[datayears[0]][s]
       histT_nom_E[s] = hist_nom_E[datayears[0]][s]
       for syst_name in list_syst:
         histT_syst_M[syst_name][s] = hist_syst_M[syst_name][datayears[0]][s]
         histT_syst_E[syst_name][s] = hist_syst_E[syst_name][datayears[0]][s]
       for d in datayears[1:]:
          histT_nom_M[s].Add(hist_nom_M[d][s])
          histT_nom_E[s].Add(hist_nom_E[d][s])
          for syst_name in list_syst:
             histT_syst_M[syst_name][s].Add(hist_syst_M[syst_name][d][s])
             histT_syst_E[syst_name][s].Add(hist_syst_E[syst_name][d][s])
       #if (args.channel in sl_channel) and (name not in not_rebin):
       #   histT_nom_M[s].Rebin(2)
       #   histT_nom_E[s].Rebin(2)
       #   for syst_name in list_syst:
       #      histT_syst_M[syst_name][s].Rebin(2)
       #      histT_syst_E[syst_name][s].Rebin(2)

  if not args.nodata:
    histD_M = histdata_M[datayears[0]]
    histD_E = histdata_E[datayears[0]]
    for d in datayears[1:]:
       histD_M.Add(histdata_M[d])
       histD_E.Add(histdata_E[d])
    #if (args.channel in sl_channel) and (name not in not_rebin): 
    #   histD_M.Rebin(2)
    #   histD_E.Rebin(2)

  env_syst_T = {}
  env_syst_M = {}
  env_syst_E = {}
  for syst_name in list_syst:
     env_syst_T[syst_name] = {}
     env_syst_M[syst_name] = {}
     env_syst_E[syst_name] = {}
  if args.sumEM:
    for s in samples:
       for syst_name in list_syst:
          env_syst_T[syst_name][s] = histT_syst_M[syst_name][s]
          env_syst_T[syst_name][s].Add(histT_syst_E[syst_name][s])
  else:
    for s in samples:
       for syst_name in list_syst:
         env_syst_M[syst_name][s] = histT_syst_M[syst_name][s]
         env_syst_E[syst_name][s] = histT_syst_E[syst_name][s]

  ## Stack creation
  samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon","ttbar_sl_else","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]
  if args.wcs:
     if args.mudhad: 
        samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm_mudhad"]
     else:
        samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]

  if args.sumEM:
    histT_nom_T = {}
    for s in samples:
      histT_nom_T[s] = histT_nom_M[s]
      histT_nom_T[s].Add(histT_nom_E[s])
    if not args.nodata:
      histD_T = histD_M
      histD_T.Add(histD_E)

  sysName = {};
  for syst in list_syst+["mcstat"]:
     sysName[syst] = syst;
  sysName["btaglightup"] = "Bot_mistagUp"; sysName["btaglightdown"] = "Bot_mistagDown"; sysName["btagheavyup"] = "Bot_tagUp"; sysName["btagheavydown"] = "Bot_tagDown";
  sysName["puwup"] = "PU_weightUp"; sysName["puwdown"] = "PU_weightDown"; sysName["seclepup"] = "muonBjet_SFUp"; sysName["seclepdown"] = "muonBjet_SFDown"; 
  sysName["jecup"] = "jecUp";sysName["jecdown"] = "jecDown";
  dataCname = {};
  if args.wcs:
      dataCname["vv"] = "bck";dataCname["wjets"] = "vjets";dataCname["ttbar_dl"] = "ttdl";dataCname["ttbar_sl_nocharm"] = "ttWuq";dataCname["ttbar_sl_charm"] = "ttWcq"; 
      dataCname["st_else"] = "stnW";dataCname["st_nocharm"] = "stWuq";dataCname["st_charm"] = "stWcq"; 
  else:
      for s in samples: 
         dataCname[s] = s; 
  if args.sumEM:
    ###############################
    ###### sumEM channel ##########
    ###############################
    #### syst from MC stat creation
    env_syst_T["mcstat"] = {}
    error_tot = {}
    for s in samples:
        env_syst_T["mcstat"][s] = histT_nom_T[s].Clone(str(s)+"mcstaterr")
        error_tot[s] = 0
        for bin in range(histT_nom_T[s].GetNbinsX()):
            error_tot[s] = error_tot[s] + histT_nom_T[s].GetBinError(bin+1)**2
            #env_syst_T["mcstat"][s].SetBinContent(bin+1,env_syst_T["mcstat"][s].GetBinContent(bin+1)+histT_nom_T[s].GetBinError(bin+1))
        error_tot[s] = sqrt(error_tot[s])
        env_syst_T["mcstat"][s].Scale((histT_nom_T[s].Integral()+error_tot[s])/histT_nom_T[s].Integral())
        print("MC stat error is %s for %s sample with %s percentual effect"%(round(error_tot[s],2),str(s),round(error_tot[s]/histT_nom_T[s].Integral(),3)))

    list_syst = list_syst+["mcstat"]
    ### hist renamings
    #data_obs = histD_T.Clone("data_obs")
    data_obs = histT_nom_T[samples[0]].Clone("data_obs")
    for	s in samples[1:]:
      data_obs.Add(histT_nom_T[s])
    mchists = {}
    mchists_syst = {}
    for syst_name in list_syst:
       mchists_syst[syst_name] = {}
    for s in samples:
      if (s == "ttbar_sl_charm" or s == "ttbar_sl_charm_mudhad"):
         mchists[s] = histT_nom_T[s].Clone("signal")
         for syst_name in list_syst:
            mchists_syst[syst_name][s] = env_syst_T[syst_name][s].Clone("signal_"+syst_name)
      else:
         mchists[s] = histT_nom_T[s].Clone(str(s))
         for syst_name in list_syst:
            mchists_syst[syst_name][s] = env_syst_T[syst_name][s].Clone(s+"_"+syst_name)

    print("Integral of data is "+str(data_obs.Integral()))
    if args.wcs:
      if args.mudhad:
         aux_list = ["vv","ttbar_dl","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm_mudhad"]
      else:
         aux_list = ["vv","ttbar_dl","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]
      for s in aux_list:
        print("-------------------------------------------------------")
        if s == "vv":
          number_nom = mchists["vv"].Integral()+mchists["ttbar_dh"].Integral()+mchists["zjets"].Integral()+mchists["wjets"].Integral()+mchists["st_else"].Integral()
          print("Integral of backgrounds is %s" %str(number_nom))
          for syst_name in list_syst:
             number = mchists_syst[syst_name]["vv"].Integral()+mchists_syst[syst_name]["ttbar_dh"].Integral()+mchists_syst[syst_name]["zjets"].Integral()+mchists_syst[syst_name]["wjets"].Integral()+mchists_syst[syst_name]["st_else"].Integral()
             print(syst_name+" effect is %s with effect of %s percent" %(str(round(number,2)),str(round(100*abs(number-number_nom)/number_nom,2)))) 
        else:
          number_nom = mchists[s].Integral()
          print("Integral of "+str(s)+" MC process is "+str(mchists[s].Integral()))
          for syst_name	in list_syst:
             number = mchists_syst[syst_name][s].Integral()
             print(syst_name+" effect is %s with effect of %s percent" %(str(round(number,2)),str(round(100*abs(number-number_nom)/number_nom,2)))) 
      print("-------------------------------------------------------")
      print("-------------------------------------------------------")
      quant = 0
      quant_syst = {}
      for syst_name in list_syst:
         quant_syst[syst_name] = 0
      for s in samples:
        quant = quant + mchists[s].Integral()
        for syst_name in list_syst:
           quant_syst[syst_name] = quant_syst[syst_name] + mchists_syst[syst_name][s].Integral()
      print("Integral of all MC processes is "+str(round(quant,2)))
      for syst_name in list_syst:
        print(syst_name+" effect is %s with effect of %s percent" %(str(quant_syst[syst_name]),str(round(100*abs(quant-quant_syst[syst_name])/quant,2)))) 
    else:
      for s in samples:
        print("Integral of "+str(s)+" MC process is "+str(mchists[s].Integral()))

    ### file saving
    if args.wcs:
      termfile = "mydatacard_"+str(args.channel)+"_"+str(name)+"_wcs.root"
    else:
      termfile = "mydatacard_"+str(args.channel)+"_"+str(name)+".root"
    myfile[name] = TFile( path_histfile+termfile, 'RECREATE' )
    data_obs.Write()
    for s in samples:
      mchists[s].Write()
      for syst_name in list_syst:
         mchists_syst[syst_name][s].Write()
    myfile[name].Close()

  ## Samples reminder
  ## samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottomgluon",
  ## "ttbar_sl_charmgluon","ttbar_sl_else","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]
  ## if wcs: samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]

  else:
    ###############################
    ###### M channel ##############
    ###############################
    #### syst from MC stat creation
    env_syst_M["mcstat"] = {}
    error_tot = {}
    for s in samples:
        env_syst_M["mcstat"][s] = histT_nom_M[s].Clone(str(s)+"mcstaterr")
        error_tot[s] = 0
        for bin in range(histT_nom_M[s].GetNbinsX()):
            error_tot[s] = error_tot[s] + histT_nom_M[s].GetBinError(bin+1)**2
            #env_syst_M["mcstat"][s].SetBinContent(bin+1,env_syst_M["mcstat"][s].GetBinContent(bin+1)+histT_nom_M[s].GetBinError(bin+1))
        error_tot[s] = sqrt(error_tot[s])
        env_syst_M["mcstat"][s].Scale((histT_nom_M[s].Integral()+error_tot[s])/histT_nom_M[s].Integral())
        print("MC stat error is %s for %s sample with %s percentual effect"%(round(error_tot[s],2),str(s),round(error_tot[s]/histT_nom_M[s].Integral(),3)))

    list_syst = list_syst+["mcstat"]
    print(list_syst)
    ### hist renamings
    data_obs = histD_M.Clone("data_obs")
    mchists = {}
    mchists_syst = {}
    for syst_name in list_syst:
       mchists_syst[syst_name] = {}
    for s in samples:
       if s == "vv":
         mchists["bck"] = histT_nom_M["vv"].Clone(dataCname[s])
         mchists["bck"].Add(histT_nom_M["ttbar_dh"])
         for syst_name in list_syst:
            mchists_syst[syst_name]["bck"] = env_syst_M[syst_name]["vv"].Clone(dataCname[s]+"_"+sysName[syst_name])
            mchists_syst[syst_name]["bck"].Add(env_syst_M[syst_name]["ttbar_dh"])
       elif s == "wjets":
         mchists["vjets"] = histT_nom_M["wjets"].Clone(dataCname[s])
       	 mchists["vjets"].Add(histT_nom_M["zjets"])
         for syst_name in list_syst:
            mchists_syst[syst_name]["vjets"] = env_syst_M[syst_name]["wjets"].Clone(dataCname[s]+"_"+sysName[syst_name])
            mchists_syst[syst_name]["vjets"].Add(env_syst_M[syst_name]["zjets"])
       elif s in ["vv","ttbar_dl","wjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]:
         mchists[s] = histT_nom_M[s].Clone(dataCname[s])
         for syst_name in list_syst:
            mchists_syst[syst_name][s] = env_syst_M[syst_name][s].Clone(dataCname[s]+"_"+sysName[syst_name])
    print("-------------------------------------------------------------------------------------------------------------------------")
    print("Integral of M data is "+str(data_obs.Integral()))
    print("-------------------------------------------------------------------------------------------------------------------------")
    if args.wcs:
      if args.mudhad:
         aux_list = ["vv","ttbar_dl","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm_mudhad"]
      else:
         #aux_list = ["vv","ttbar_dl","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]
         #aux_list = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]
         aux_list = ["bck","ttbar_dl","vjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]
      for s in aux_list:
        print("-------------------------------------------------------")
        #if s == "vv":
        if False:
          number_nom = mchists["vv"].Integral()+mchists["ttbar_dh"].Integral()+mchists["zjets"].Integral()+mchists["wjets"].Integral()+mchists["st_else"].Integral()
          print("Integral of M channel backgrounds is %s" %str(round(number_nom,2)))
          for syst_name in list_syst:
             number = mchists_syst[syst_name]["vv"].Integral()+mchists_syst[syst_name]["ttbar_dh"].Integral()+mchists_syst[syst_name]["zjets"].Integral()+mchists_syst[syst_name]["wjets"].Integral()+mchists_syst[syst_name]["st_else"].Integral()
             print(syst_name+" effect is %s with effect of %s percent for M channel" %(str(round(number,2)),str(round(100*abs(number-number_nom)/number_nom,2))))
        else:
          number_nom = mchists[s].Integral()
          print("Integral of "+str(s)+" MC in M channel process is "+str(round(number_nom,2)))
          for syst_name in list_syst:
             number = mchists_syst[syst_name][s].Integral()
             print(syst_name+" effect is %s with effect of %s percent for M channel" %(str(round(number,2)),str(round(100*abs(number-number_nom)/number_nom,2))))
      print("-------------------------------------------------------")
      print("-------------------------------------------------------")
      quant = 0
      quant_syst = {}
      for syst_name in list_syst:
         quant_syst[syst_name] = 0
      for s in aux_list:
        quant = quant + mchists[s].Integral()
        for syst_name in list_syst:
           quant_syst[syst_name] = quant_syst[syst_name] + mchists_syst[syst_name][s].Integral()
      print("Integral of all MC processes in M channel is "+str(round(quant,2)) + " for M cahnnel")
      for syst_name in list_syst:
        print(syst_name+" effect is %s with effect of %s percent for M channel" %(str(quant_syst[syst_name]),str(round(100*abs(quant-quant_syst[syst_name])/quant,2))))
    else:
      for s in samples:
        print("Integral of "+str(s)+" MC process is "+str(mchists[s].Integral()) + " for M channel")

    ####### FILE SAVING
    if args.wcs:
      termfile = "mydatacard_leptonM_"+str(args.channel)+"_"+str(name)+"_wcs.root"
    else:
      termfile = "mydatacard_leptonM_"+str(args.channel)+"_"+str(name)+".root"
    myfile[name] = TFile( path_histfile+termfile, 'RECREATE' )
    data_obs.Write()
    for s in aux_list:
       mchists[s].Write()
       for syst_name in ["jecup", "jecdown"]:
         if syst_name in list_syst:
           mchists_syst[syst_name][s].Write()
    myfile[name].Close()
    ###############################
    ###### E channel ##############
    ###############################
    #### syst from MC stat creation
    env_syst_E["mcstat"] = {}
    error_tot = {}
    for s in samples:
        env_syst_E["mcstat"][s] = histT_nom_E[s].Clone(str(s)+"mcstaterr")
        error_tot[s] = 0
        for bin in range(histT_nom_E[s].GetNbinsX()):
            error_tot[s] = error_tot[s] + histT_nom_E[s].GetBinError(bin+1)**2
            #env_syst_E["mcstat"][s].SetBinContent(bin+1,env_syst_E["mcstat"][s].GetBinContent(bin+1)+histT_nom_E[s].GetBinError(bin+1))
        error_tot[s] = sqrt(error_tot[s])
        env_syst_E["mcstat"][s].Scale((histT_nom_E[s].Integral()+error_tot[s])/histT_nom_E[s].Integral())
        print("MC stat error is %s for %s sample with %s percentual effect"%(round(error_tot[s],2),str(s),round(error_tot[s]/histT_nom_E[s].Integral(),3)))

    #list_syst = list_syst+["mcstat"]
    ### hist renamings
    data_obs = histD_E.Clone("data_obs")
    mchists = {}
    mchists_syst = {}
    for syst_name in list_syst:
       mchists_syst[syst_name] = {}
    for s in samples:
       if s == "vv":
         mchists["bck"] = histT_nom_E["vv"].Clone(dataCname[s])
       	 mchists["bck"].Add(histT_nom_E["ttbar_dh"])
         for syst_name in list_syst:
            mchists_syst[syst_name]["bck"] = env_syst_E[syst_name]["vv"].Clone(dataCname[s]+"_"+sysName[syst_name])
            mchists_syst[syst_name]["bck"].Add(env_syst_E[syst_name]["ttbar_dh"])
       elif s == "wjets":
         mchists["vjets"] = histT_nom_E["wjets"].Clone(dataCname[s])
         mchists["vjets"].Add(histT_nom_E["zjets"])
         for syst_name in list_syst:
            mchists_syst[syst_name]["vjets"] = env_syst_E[syst_name]["wjets"].Clone(dataCname[s]+"_"+sysName[syst_name])
            mchists_syst[syst_name]["vjets"].Add(env_syst_E[syst_name]["zjets"])
       elif s in ["vv","ttbar_dl","wjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]:
         mchists[s] = histT_nom_E[s].Clone(dataCname[s])
         for syst_name in list_syst:
            mchists_syst[syst_name][s] = env_syst_E[syst_name][s].Clone(dataCname[s]+"_"+sysName[syst_name])

    print("-------------------------------------------------------------------------------------------------------------------------")
    print("Integral of E data is "+str(data_obs.Integral()))
    print("-------------------------------------------------------------------------------------------------------------------------")
    if args.wcs:
      if args.mudhad:
         aux_list = ["vv","ttbar_dl","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm_mudhad"]
      else:
         #aux_list = ["vv","ttbar_dl","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]
         #aux_list = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]
         aux_list = ["bck","ttbar_dl","vjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]
      for s in aux_list:
        print("-------------------------------------------------------")
        #if s == "vv":
        if False:
          number_nom = mchists["vv"].Integral()+mchists["ttbar_dh"].Integral()+mchists["zjets"].Integral()+mchists["wjets"].Integral()+mchists["st_else"].Integral()
          print("Integral of E channel backgrounds is %s" %str(round(number_nom,2)))
          for syst_name in list_syst:
             number = mchists_syst[syst_name]["vv"].Integral()+mchists_syst[syst_name]["ttbar_dh"].Integral()+mchists_syst[syst_name]["zjets"].Integral()+mchists_syst[syst_name]["wjets"].Integral()+mchists_syst[syst_name]["st_else"].Integral()
             print(syst_name+" effect is %s with effect of %s percent for E channel" %(str(round(number,2)),str(round(100*abs(number-number_nom)/number_nom,2))))
        else:
          number_nom = mchists[s].Integral()
          print("Integral of "+str(s)+" MC in E channel process is "+str(round(number_nom,2)))
          for syst_name in list_syst:
             number = mchists_syst[syst_name][s].Integral()
             print(syst_name+" effect is %s with effect of %s percent for E channel" %(str(round(number,2)),str(round(100*abs(number-number_nom)/number_nom,2))))
      print("-------------------------------------------------------")
      print("-------------------------------------------------------")
      quant = 0
      quant_syst = {}
      for syst_name in list_syst:
         quant_syst[syst_name] = 0
      for s in aux_list:
        quant = quant + mchists[s].Integral()
        for syst_name in list_syst:
           quant_syst[syst_name] = quant_syst[syst_name] + mchists_syst[syst_name][s].Integral()
      print("Integral of all MC processes in E channel is "+str(round(quant,2)))
      for syst_name in list_syst:
        print(syst_name+" effect is %s with effect of %s percent for E channel" %(str(quant_syst[syst_name]),str(round(100*abs(quant-quant_syst[syst_name])/quant,2))))
    else:
      for s in samples:
        print("Integral of "+str(s)+" MC process is "+str(mchists[s].Integral()))

    ####### FILE SAVING
    if args.wcs:
      termfile = "mydatacard_leptonE_"+str(args.channel)+"_"+str(name)+"_wcs.root"
    else:
      termfile = "mydatacard_leptonE_"+str(args.channel)+"_"+str(name)+".root"
    myfile[name] = TFile( path_histfile+termfile, 'RECREATE' )
    data_obs.Write()
    for s in aux_list:
         mchists[s].Write()
         for syst_name in ["jecup", "jecdown"]:
            if syst_name in list_syst:
               mchists_syst[syst_name][s].Write()
    myfile[name].Close()

  for data_op in datayears:
         histFile[name][data_op].Close()
         if not args.nodata:
            histFileDM[name][data_op].Close()
            histFileDE[name][data_op].Close()



