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
parser.add_argument("--png", action="store_true", default=False,
                    help="png format")
parser.add_argument("--nodata", action="store_true", default=False,
                    help="Do not plot data")
parser.add_argument("--wcs", action="store_true", default=False,
                    help="classification for ttbar and st")
parser.add_argument("--norm", action="store_true", default=False,
                    help="normalisation")
parser.add_argument("--year", type=string, default="2016",
                    help="Select year of process to run")
parser.add_argument("--syst", type=string, default="jec",
                    help="systematic to plot")
parser.add_argument("--channel", type=string, default="btagMM",
                    help="Select year of process to run")
parser.add_argument("--nosyst", action="store_true", default=False,
                    help="systematics inclusion")
parser.add_argument("--sumEM", action="store_true", default=False,
                    help="union of M and E channels")

# Use like:
# python higgssearch/fromJF/hist_fromJFwqq_syst_total.py --stack --ratio --png --norm --year="all" --channel="btagMM_chitest_sl"

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

if args.channel == "btagMM_chitest": term_path = ""
elif args.channel == "btagMM_chitest_sl": term_path = "/sl"
elif args.channel == "btagMM_chitest_slss": term_path = "/sl/ss"
elif args.channel == "btagMM_chitest_slos": term_path = "/sl/os"
elif args.channel == "btagMM_chitest_slssos": term_path = "/sl/ssos"
else: raise NameError('Incorrect data option')

sl_channel = ["btagMM_chitest_sl","lepton50_chitest_sl","btagMM_chitest_slss","lepton50_chitest_slss",
      "btagMM_chitest_slos","lepton50_chitest_slos","btagMM_chitest_slssos","lepton50_chitest_slssos"]

plotdir = '/nfs/cms/vazqueze/higgssearch/plotspng/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

c_rat = 1.2
c_rat2 = 0.8
nrebin = 2

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
         "muon_jet_z2_v2","muon_jet_z3","muon_jet_iso_abs"]

not_rebin = ["nJetGood","InvM3_bad","InvMl_good","InvMl_bad","lepton_eta_thick","jet_bot1_btagnumber", "jet_bot2_btagnumber", 
      "jet_1_btagnumber", "jet_2_btagnumber","muon_jet_pt","muon_jet_relpt","muon_jet_eta","tau_discr","transverse_mass",
      "jet_1_flavourP", "jet_2_flavourP", "jet_bot1_flavourP", "jet_bot2_flavourP","lepton_pt", "muon_jet_z","tau_discr_jet1","tau_discr_jet2","tau_discr_jetbot1",
      "tau_discr_jetbot2","muon_jet_iso","muon_jet_iso_abs","InvM3_good_short"]

#observable_names = ["jet_max_cvltag","InvM_2jets"]

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
        "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else"]

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

for name in observable_names:
  ## Open hists files
  filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/chi_test"+term_path+"/"
  if args.wcs: 
     #filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/chi_test/wcs_classes_aux/test_aux"+term_path+"/"
     filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/chi_test/wcs_classes"+term_path+"/"
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

  if args.png: c1 = TCanvas("c1","",1200,800)
  else: c1 = TCanvas("c1","",600,400)

  samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_light","ttbar_sl_bottom","ttbar_dl","ttbar_dh","zz","wz",
        "st_1","st_2","st_3","st_4","ttbar_sl_else","ttbar_sl_charmgluon","ttbar_sl_bottomgluon"]

  if args.wcs: 
    samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
            "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
            "ttbar_sl_charm","ttbar_sl_nocharm","ttbar_dl","ttbar_dh","zz","wz",
            "st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm",
            "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else"]

  ## HISTS
  samples_foryear = {}
  hist_nom_M = {}
  hist_nom_E = {}
  hist_jecup_M = {}
  hist_jecup_E = {}
  hist_jecdown_M = {}
  hist_jecdown_E = {}
  hist_smearup_M = {}
  hist_smearup_E = {}
  hist_smeardown_M = {}
  hist_smeardown_E = {}
  histdata_M = {}
  histdata_E = {}
  for data_op in datayears:
    hist_nom_M[data_op] = {}
    hist_nom_E[data_op] = {}
    hist_jecup_M[data_op] = {}
    hist_jecup_E[data_op] = {}
    hist_jecdown_M[data_op] = {}
    hist_jecdown_E[data_op] = {}
    hist_smearup_M[data_op] = {}
    hist_smearup_E[data_op] = {}
    hist_smeardown_M[data_op] = {}
    hist_smeardown_E[data_op] = {}
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
      if not args.nosyst:
         hist_jecup_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_jecup")
         hist_jecup_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_jecup")
         hist_jecdown_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_jecdown")
         hist_jecdown_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_jecdown")
         hist_smearup_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_smearup")
         hist_smearup_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_smearup")
         hist_smeardown_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_smeardown")
         hist_smeardown_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_smeardown")
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
      if not args.nosyst:
         hist_jecup_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_jecup_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_jecdown_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_jecdown_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_smearup_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_smearup_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_smeardown_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_smeardown_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
      if args.norm: hist_nom_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
      if args.norm: hist_nom_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
      if not args.nosyst:
         if args.norm: hist_jecup_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
         if args.norm: hist_jecup_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
         if args.norm: hist_jecdown_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
         if args.norm: hist_jecdown_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
         if args.norm: hist_smearup_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
         if args.norm: hist_smearup_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
         if args.norm: hist_smeardown_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
         if args.norm: hist_smeardown_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
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
      if not args.nosyst:
        hist_jecup_M[data_op]["st"+apx] = hist_jecup_M[data_op][list_st[0]+apx]
        hist_jecup_E[data_op]["st"+apx] = hist_jecup_E[data_op][list_st[0]+apx]
        hist_jecdown_M[data_op]["st"+apx] = hist_jecdown_M[data_op][list_st[0]+apx]
        hist_jecdown_E[data_op]["st"+apx] = hist_jecdown_E[data_op][list_st[0]+apx]
        hist_smearup_M[data_op]["st"+apx] = hist_smearup_M[data_op][list_st[0]+apx]
        hist_smearup_E[data_op]["st"+apx] = hist_smearup_E[data_op][list_st[0]+apx]
        hist_smeardown_M[data_op]["st"+apx] = hist_smeardown_M[data_op][list_st[0]+apx]
        hist_smeardown_E[data_op]["st"+apx] = hist_smeardown_E[data_op][list_st[0]+apx]
      for l in list_st[1:]:
        hist_nom_M[data_op]["st"+apx].Add(hist_nom_M[data_op][l+apx])
        hist_nom_E[data_op]["st"+apx].Add(hist_nom_E[data_op][l+apx])
        if not args.nosyst:
          hist_jecup_M[data_op]["st"+apx].Add(hist_jecup_M[data_op][l+apx])
          hist_jecup_E[data_op]["st"+apx].Add(hist_jecup_E[data_op][l+apx])
          hist_jecdown_M[data_op]["st"+apx].Add(hist_jecdown_M[data_op][l+apx])
          hist_jecdown_E[data_op]["st"+apx].Add(hist_jecdown_E[data_op][l+apx])
          hist_smearup_M[data_op]["st"+apx].Add(hist_smearup_M[data_op][l+apx])
          hist_smearup_E[data_op]["st"+apx].Add(hist_smearup_E[data_op][l+apx])
          hist_smeardown_M[data_op]["st"+apx].Add(hist_smeardown_M[data_op][l+apx])
          hist_smeardown_E[data_op]["st"+apx].Add(hist_smeardown_E[data_op][l+apx])
    hist_nom_M[data_op]["wjets"] = hist_nom_M[data_op][list_wjets[0]]
    hist_nom_E[data_op]["wjets"] = hist_nom_E[data_op][list_wjets[0]]
    if not args.nosyst:
      hist_jecup_M[data_op]["wjets"] = hist_jecup_M[data_op][list_wjets[0]]
      hist_jecup_E[data_op]["wjets"] = hist_jecup_E[data_op][list_wjets[0]]
      hist_jecdown_M[data_op]["wjets"] = hist_jecdown_M[data_op][list_wjets[0]]
      hist_jecdown_E[data_op]["wjets"] = hist_jecdown_E[data_op][list_wjets[0]]
      hist_smearup_M[data_op]["wjets"] = hist_smearup_M[data_op][list_wjets[0]]
      hist_smearup_E[data_op]["wjets"] = hist_smearup_E[data_op][list_wjets[0]]
      hist_smeardown_M[data_op]["wjets"] = hist_smeardown_M[data_op][list_wjets[0]]
      hist_smeardown_E[data_op]["wjets"] = hist_smeardown_E[data_op][list_wjets[0]]
    for l in list_wjets[1:]:
      hist_nom_M[data_op]["wjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["wjets"].Add(hist_nom_E[data_op][l])
      if not args.nosyst:
        hist_jecup_M[data_op]["wjets"].Add(hist_jecup_M[data_op][l])
        hist_jecup_E[data_op]["wjets"].Add(hist_jecup_E[data_op][l])
        hist_jecdown_M[data_op]["wjets"].Add(hist_jecdown_M[data_op][l])
        hist_jecdown_E[data_op]["wjets"].Add(hist_jecdown_E[data_op][l])
        hist_smearup_M[data_op]["wjets"].Add(hist_smearup_M[data_op][l])
        hist_smearup_E[data_op]["wjets"].Add(hist_smearup_E[data_op][l])
        hist_smeardown_M[data_op]["wjets"].Add(hist_smeardown_M[data_op][l])
        hist_smeardown_E[data_op]["wjets"].Add(hist_smeardown_E[data_op][l])
    hist_nom_M[data_op]["zjets"] = hist_nom_M[data_op][list_zjets[0]]
    hist_nom_E[data_op]["zjets"] = hist_nom_E[data_op][list_zjets[0]]
    if not args.nosyst:
      hist_jecup_M[data_op]["zjets"] = hist_jecup_M[data_op][list_zjets[0]]
      hist_jecup_E[data_op]["zjets"] = hist_jecup_E[data_op][list_zjets[0]]
      hist_jecdown_M[data_op]["zjets"] = hist_jecdown_M[data_op][list_zjets[0]]
      hist_jecdown_E[data_op]["zjets"] = hist_jecdown_E[data_op][list_zjets[0]]
      hist_smearup_M[data_op]["zjets"] = hist_smearup_M[data_op][list_zjets[0]]
      hist_smearup_E[data_op]["zjets"] = hist_smearup_E[data_op][list_zjets[0]]
      hist_smeardown_M[data_op]["zjets"] = hist_smeardown_M[data_op][list_zjets[0]]
      hist_smeardown_E[data_op]["zjets"] = hist_smeardown_E[data_op][list_zjets[0]]
    for l in list_zjets[1:]:
      hist_nom_M[data_op]["zjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["zjets"].Add(hist_nom_E[data_op][l])
      if not args.nosyst:
        hist_jecup_M[data_op]["zjets"].Add(hist_jecup_M[data_op][l])
        hist_jecup_E[data_op]["zjets"].Add(hist_jecup_E[data_op][l])
        hist_jecdown_M[data_op]["zjets"].Add(hist_jecdown_M[data_op][l])
        hist_jecdown_E[data_op]["zjets"].Add(hist_jecdown_E[data_op][l])
        hist_smearup_M[data_op]["zjets"].Add(hist_smearup_M[data_op][l])
        hist_smearup_E[data_op]["zjets"].Add(hist_smearup_E[data_op][l])
        hist_smeardown_M[data_op]["zjets"].Add(hist_smeardown_M[data_op][l])
        hist_smeardown_E[data_op]["zjets"].Add(hist_smeardown_E[data_op][l])
    hist_nom_M[data_op]["vv"] = hist_nom_M[data_op][list_vv[0]]
    hist_nom_E[data_op]["vv"] = hist_nom_E[data_op][list_vv[0]]
    if not args.nosyst:
      hist_jecup_M[data_op]["vv"] = hist_jecup_M[data_op][list_vv[0]]
      hist_jecup_E[data_op]["vv"] = hist_jecup_E[data_op][list_vv[0]]
      hist_jecdown_M[data_op]["vv"] = hist_jecdown_M[data_op][list_vv[0]]
      hist_jecdown_E[data_op]["vv"] = hist_jecdown_E[data_op][list_vv[0]]
      hist_smearup_M[data_op]["vv"] = hist_smearup_M[data_op][list_vv[0]]
      hist_smearup_E[data_op]["vv"] = hist_smearup_E[data_op][list_vv[0]]
      hist_smeardown_M[data_op]["vv"] = hist_smeardown_M[data_op][list_vv[0]]
      hist_smeardown_E[data_op]["vv"] = hist_smeardown_E[data_op][list_vv[0]]
    for l in list_vv[1:]:
      hist_nom_M[data_op]["vv"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["vv"].Add(hist_nom_E[data_op][l])
      if not args.nosyst:
        hist_jecup_M[data_op]["vv"].Add(hist_jecup_M[data_op][l])
        hist_jecup_E[data_op]["vv"].Add(hist_jecup_E[data_op][l])
        hist_jecdown_M[data_op]["vv"].Add(hist_jecdown_M[data_op][l])
        hist_jecdown_E[data_op]["vv"].Add(hist_jecdown_E[data_op][l])
        hist_smearup_M[data_op]["vv"].Add(hist_smearup_M[data_op][l])
        hist_smearup_E[data_op]["vv"].Add(hist_smearup_E[data_op][l])
        hist_smeardown_M[data_op]["vv"].Add(hist_smeardown_M[data_op][l])
        hist_smeardown_E[data_op]["vv"].Add(hist_smeardown_E[data_op][l])

  samples = ["ttbar_sl_bottom","ttbar_sl_charm","ttbar_sl_else","ttbar_sl_light","ttbar_dl","ttbar_dh","zjets","vv","st","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon"]
  if args.wcs: samples = ["ttbar_sl_nocharm","ttbar_sl_charm","ttbar_dl","ttbar_dh","zjets","vv","st_nocharm","st_charm","st_else","wjets"]

  ## sumamos todos los anos para todos los histogramas
  histT_nom_M = {}
  histT_nom_E = {}
  histT_jecup_M = {}
  histT_jecup_E = {}
  histT_jecdown_M = {}
  histT_jecdown_E = {}
  histT_smearup_M = {}
  histT_smearup_E = {}
  histT_smeardown_M = {}
  histT_smeardown_E = {}
  for s in samples:
       histT_nom_M[s] = hist_nom_M[datayears[0]][s]
       histT_nom_E[s] = hist_nom_E[datayears[0]][s]
       if not args.nosyst:
         histT_jecup_M[s] = hist_jecup_M[datayears[0]][s]
         histT_jecup_E[s] = hist_jecup_E[datayears[0]][s]
         histT_jecdown_M[s] = hist_jecdown_M[datayears[0]][s]
         histT_jecdown_E[s] = hist_jecdown_E[datayears[0]][s]
         histT_smearup_M[s] = hist_smearup_M[datayears[0]][s]
         histT_smearup_E[s] = hist_smearup_E[datayears[0]][s]
         histT_smeardown_M[s] = hist_smeardown_M[datayears[0]][s]
         histT_smeardown_E[s] = hist_smeardown_E[datayears[0]][s]
       for d in datayears[1:]:
          histT_nom_M[s].Add(hist_nom_M[d][s])
          histT_nom_E[s].Add(hist_nom_E[d][s])
          if not args.nosyst:
            histT_jecup_M[s].Add(hist_jecup_M[d][s])
            histT_jecup_E[s].Add(hist_jecup_E[d][s])
            histT_jecdown_M[s].Add(hist_jecdown_M[d][s])
            histT_jecdown_E[s].Add(hist_jecdown_E[d][s])
            histT_smearup_M[s].Add(hist_smearup_M[d][s])
            histT_smearup_E[s].Add(hist_smearup_E[d][s])
            histT_smeardown_M[s].Add(hist_smeardown_M[d][s])
            histT_smeardown_E[s].Add(hist_smeardown_E[d][s])
       if (args.channel in sl_channel) and (name not in not_rebin):
          histT_nom_M[s].Rebin(nrebin)
          histT_nom_E[s].Rebin(nrebin)
          if not args.nosyst:
            histT_jecup_M[s].Rebin(nrebin)
            histT_jecup_E[s].Rebin(nrebin)
            histT_jecdown_M[s].Rebin(nrebin)
            histT_jecdown_E[s].Rebin(nrebin)
            histT_smearup_M[s].Rebin(nrebin)
            histT_smearup_E[s].Rebin(nrebin)
            histT_smeardown_M[s].Rebin(nrebin)
            histT_smeardown_E[s].Rebin(nrebin)

  if not args.nodata:
    histD_M = histdata_M[datayears[0]]
    histD_E = histdata_E[datayears[0]]
    for d in datayears[1:]:
       histD_M.Add(histdata_M[d])
       histD_E.Add(histdata_E[d])
    if (args.channel in sl_channel) and (name not in not_rebin): 
       histD_M.Rebin(nrebin)
       histD_E.Rebin(nrebin)

  if not args.nosyst:
    ### Sumamos todas las contribuciones
    histT_sT_nom_M = histT_nom_M[samples[0]]
    histT_sT_nom_E = histT_nom_E[samples[0]]
    histT_sT_jecup_M = histT_jecup_M[samples[0]]
    histT_sT_jecup_E = histT_jecup_E[samples[0]]
    histT_sT_jecdown_M = histT_jecdown_M[samples[0]]
    histT_sT_jecdown_E = histT_jecdown_E[samples[0]]
    histT_sT_smearup_M = histT_smearup_M[samples[0]]
    histT_sT_smearup_E = histT_smearup_E[samples[0]]
    histT_sT_smeardown_M = histT_smeardown_M[samples[0]]
    histT_sT_smeardown_E = histT_smeardown_E[samples[0]]
    for s in samples[1:]:
      histT_sT_nom_M.Add(histT_nom_M[s])
      histT_sT_nom_E.Add(histT_nom_E[s])
      histT_sT_jecup_M.Add(histT_jecup_M[s])
      histT_sT_jecup_E.Add(histT_jecup_E[s])
      histT_sT_jecdown_M.Add(histT_jecdown_M[s])
      histT_sT_jecdown_E.Add(histT_jecdown_E[s])
      histT_sT_smearup_M.Add(histT_smearup_M[s])
      histT_sT_smearup_E.Add(histT_smearup_E[s])
      histT_sT_smeardown_M.Add(histT_smeardown_M[s])
      histT_sT_smeardown_E.Add(histT_smeardown_E[s])

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  colors = {}
  colors["ttbar_dl"] = (222,90,106)
  colors["wjets"] = (155,152,204)
  colors["ttbar_sl_bottomgluon"] = (236,164,207)
  colors["ttbar_sl_charm"] = (204,255,153)
  colors["ttbar_sl_light"] = (120,154,86)
  colors["ttbar_sl_nocharm"] = (120,154,86)
  colors["ttbar_sl_bottom"] = (201,79,152)
  colors["ttbar_sl_else"] = (242,193,121)
  colors["ttbar_sl_charmgluon"] = (222,212,74)
  colors["vv"] = (255,180,85)
  colors["ttbar_dh"] = (204,204,0)
  colors["zjets"] = (113,209,223)
  colors["st"] = (153,51,255)
  colors["st_charm"] = (102,0,204)
  colors["st_nocharm"] = (198,101,222)
  colors["st_else"] = (207,176,235)

  if args.sumEM:
      histT_sT_nom_M.Add(histT_sT_nom_E)
      histT_sT_jecup_M.Add(histT_sT_jecup_E)
      histT_sT_jecdown_M.Add(histT_sT_jecdown_E)
      histT_sT_smearup_M.Add(histT_sT_smearup_E)
      histT_sT_smeardown_M.Add(histT_sT_smeardown_E)

      histT_sT_nom_M.Scale(1/histT_sT_nom_M.Integral())
      histT_sT_jecup_M.Scale(1/histT_sT_jecup_M.Integral())
      histT_sT_jecdown_M.Scale(1/histT_sT_jecdown_M.Integral())
      histT_sT_smearup_M.Scale(1/histT_sT_smearup_M.Integral())
      histT_sT_smeardown_M.Scale(1/histT_sT_smeardown_M.Integral())

      histT_sT_nom_M.SetLineWidth(2)
      histT_sT_nom_M.SetLineColor(kBlack)
      histT_sT_nom_M.GetYaxis().SetTitle("Number of events")
      histT_sT_nom_M.GetXaxis().SetTitle(name)
      histT_sT_jecup_M.SetLineWidth(2)
      histT_sT_jecup_M.SetLineColor(ROOT.TColor.GetColor(*(222,90,106)))
      histT_sT_jecup_M.GetYaxis().SetTitle("Number of events")
      histT_sT_jecup_M.GetXaxis().SetTitle(name)
      histT_sT_jecdown_M.SetLineWidth(2)
      histT_sT_jecdown_M.SetLineColor(ROOT.TColor.GetColor(*(120,154,86)))
      histT_sT_jecdown_M.GetYaxis().SetTitle("Number of events")
      histT_sT_jecdown_M.GetXaxis().SetTitle(name)
      histT_sT_smearup_M.SetLineWidth(2)
      histT_sT_smearup_M.SetLineColor(ROOT.TColor.GetColor(*(222,90,106)))
      histT_sT_smearup_M.GetYaxis().SetTitle("Number of events")
      histT_sT_smearup_M.GetXaxis().SetTitle(name)
      histT_sT_smeardown_M.SetLineWidth(2)
      histT_sT_smeardown_M.SetLineColor(ROOT.TColor.GetColor(*(120,154,86)))
      histT_sT_smeardown_M.GetYaxis().SetTitle("Number of events")
      histT_sT_smeardown_M.GetXaxis().SetTitle(name)
  else:
      histT_sT_nom_M.Scale(1/histT_sT_nom_M.Integral())
      histT_sT_jecup_M.Scale(1/histT_sT_jecup_M.Integral())
      histT_sT_jecdown_M.Scale(1/histT_sT_jecdown_M.Integral())
      histT_sT_smearup_M.Scale(1/histT_sT_smearup_M.Integral())
      histT_sT_smeardown_M.Scale(1/histT_sT_smeardown_M.Integral())

      histT_sT_nom_M.SetLineWidth(2)
      histT_sT_nom_M.SetLineColor(kBlack)
      histT_sT_nom_M.GetYaxis().SetTitle("Number of events")
      histT_sT_nom_M.GetXaxis().SetTitle(name)
      histT_sT_jecup_M.SetLineWidth(2)
      histT_sT_jecup_M.SetLineColor(ROOT.TColor.GetColor(*(222,90,106)))
      histT_sT_jecup_M.GetYaxis().SetTitle("Number of events")
      histT_sT_jecup_M.GetXaxis().SetTitle(name)
      histT_sT_jecdown_M.SetLineWidth(2)
      histT_sT_jecdown_M.SetLineColor(ROOT.TColor.GetColor(*(120,154,86)))
      histT_sT_jecdown_M.GetYaxis().SetTitle("Number of events")
      histT_sT_jecdown_M.GetXaxis().SetTitle(name)
      histT_sT_smearup_M.SetLineWidth(2)
      histT_sT_smearup_M.SetLineColor(ROOT.TColor.GetColor(*(222,90,106)))
      histT_sT_smearup_M.GetYaxis().SetTitle("Number of events")
      histT_sT_smearup_M.GetXaxis().SetTitle(name)
      histT_sT_smeardown_M.SetLineWidth(2)
      histT_sT_smeardown_M.SetLineColor(ROOT.TColor.GetColor(*(120,154,86)))
      histT_sT_smeardown_M.GetYaxis().SetTitle("Number of events")
      histT_sT_smeardown_M.GetXaxis().SetTitle(name)

      histT_sT_nom_E.Scale(1/histT_sT_nom_E.Integral())
      histT_sT_jecup_E.Scale(1/histT_sT_jecup_E.Integral())
      histT_sT_jecdown_E.Scale(1/histT_sT_jecdown_E.Integral())
      histT_sT_smearup_E.Scale(1/histT_sT_smearup_E.Integral())
      histT_sT_smeardown_E.Scale(1/histT_sT_smeardown_E.Integral())

      histT_sT_nom_E.SetLineWidth(2)
      histT_sT_nom_E.SetLineColor(kBlack)
      histT_sT_nom_E.GetYaxis().SetTitle("Number of events")
      histT_sT_nom_E.GetXaxis().SetTitle(name)
      histT_sT_jecup_E.SetLineWidth(2)
      histT_sT_jecup_E.SetLineColor(ROOT.TColor.GetColor(*(222,90,106)))
      histT_sT_jecup_E.GetYaxis().SetTitle("Number of events")
      histT_sT_jecup_E.GetXaxis().SetTitle(name)
      histT_sT_jecdown_E.SetLineWidth(2)
      histT_sT_jecdown_E.SetLineColor(ROOT.TColor.GetColor(*(120,154,86)))
      histT_sT_jecdown_E.GetYaxis().SetTitle("Number of events")
      histT_sT_jecdown_E.GetXaxis().SetTitle(name)
      histT_sT_smearup_E.SetLineWidth(2)
      histT_sT_smearup_E.SetLineColor(ROOT.TColor.GetColor(*(222,90,106)))
      histT_sT_smearup_E.GetYaxis().SetTitle("Number of events")
      histT_sT_smearup_E.GetXaxis().SetTitle(name)
      histT_sT_smeardown_E.SetLineWidth(2)
      histT_sT_smeardown_E.SetLineColor(ROOT.TColor.GetColor(*(120,154,86)))
      histT_sT_smeardown_E.GetYaxis().SetTitle("Number of events")
      histT_sT_smeardown_E.GetXaxis().SetTitle(name)

  ## Stack creation
  samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon","ttbar_sl_else","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]
  if args.wcs: samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]

  if args.sumEM:
      #####################################
      ########## sumEM channel ############
      #####################################

      leg = TLegend(0.8,0.75,0.95,0.95)
      leg.SetBorderSize(1)
      if args.syst == "jec":
         histT_sT_nom_M.Draw("HIST C")
         histT_sT_jecup_M.Draw("HIST SAME C")
         histT_sT_jecdown_M.Draw("HIST SAME C")
         leg.AddEntry(histT_sT_nom_M,"nominal","f")
         leg.AddEntry(histT_sT_jecup_M,"JEC up","f")
         leg.AddEntry(histT_sT_jecdown_M,"JEC down","f")
      if args.syst == "smear":
         histT_sT_nom_M.Draw("HIST C")
         histT_sT_smearup_M.Draw("HIST SAME C")
         histT_sT_smeardown_M.Draw("HIST SAME C")
         leg.AddEntry(histT_sT_nom_M,"nominal","f")
         leg.AddEntry(histT_sT_smearup_M,"smearing up","f")
         leg.AddEntry(histT_sT_smeardown_M,"smearing down","f")
      leg.Draw()
      termp= "totalHT_wqq"

      if args.png: c1.Print(plotdir+termp+"_normed_"+name + ".png")
      else: c1.Print(plotdir+termp+"_normed_"+name + ".root")

      if c1: 
         c1.Close(); gSystem.ProcessEvents();

  else:
      #####################################
      ############ M channel ##############
      #####################################

      leg = TLegend(0.8,0.75,0.95,0.95)
      leg.SetBorderSize(1)
      if args.syst == "jec":
         histT_sT_nom_M.Draw("HIST C")
         histT_sT_jecup_M.Draw("HIST SAME C")
         histT_sT_jecdown_M.Draw("HIST SAME C")
         leg.AddEntry(histT_sT_nom_M,"nominal","f")
         leg.AddEntry(histT_sT_jecup_M,"JEC up","f")
         leg.AddEntry(histT_sT_jecdown_M,"JEC down","f")
      if args.syst == "smear":
         histT_sT_nom_M.Draw("HIST C")
         histT_sT_smearup_M.Draw("HIST SAME C")
         histT_sT_smeardown_M.Draw("HIST SAME C")
         leg.AddEntry(histT_sT_nom_M,"nominal","f")
         leg.AddEntry(histT_sT_smearup_M,"smearing up","f")
         leg.AddEntry(histT_sT_smeardown_M,"smearing down","f")
      leg.Draw()
      termp= "totalHT_wqq"

      if args.png: c1.Print(plotdir+termp+"_normed_"+name + "_M.png")
      else: c1.Print(plotdir+termp+"_normed_"+name + "_M.root")

      if c1:
         c1.Close(); gSystem.ProcessEvents();

      #####################################
      ############ E channel ##############
      #####################################

      leg = TLegend(0.8,0.75,0.95,0.95)
      leg.SetBorderSize(1)
      if args.syst == "jec":
         histT_sT_nom_E.Draw("HIST C")
         histT_sT_jecup_E.Draw("HIST SAME C")
         histT_sT_jecdown_E.Draw("HIST SAME C")
         leg.AddEntry(histT_sT_nom_E,"nominal","f")
         leg.AddEntry(histT_sT_jecup_E,"JEC up","f")
         leg.AddEntry(histT_sT_jecdown_E,"JEC down","f")
      if args.syst == "smear":
         histT_sT_nom_E.Draw("HIST C")
         histT_sT_smearup_E.Draw("HIST SAME C")
         histT_sT_smeardown_E.Draw("HIST SAME C")
         leg.AddEntry(histT_sT_nom_E,"nominal","f")
         leg.AddEntry(histT_sT_smearup_E,"smearing up","f")
         leg.AddEntry(histT_sT_smeardown_E,"smearing down","f")
      leg.Draw()
      termp= "totalHT_wqq"

      if args.png: c1.Print(plotdir+termp+"_normed_"+name + ".png")
      else: c1.Print(plotdir+termp+"_normed_"+name + ".root")

      if c1:
         c1.Close(); gSystem.ProcessEvents();

  for data_op in datayears:
      histFile[name][data_op].Close()
      if not args.nodata:
            histFileDM[name][data_op].Close()
            histFileDE[name][data_op].Close()


