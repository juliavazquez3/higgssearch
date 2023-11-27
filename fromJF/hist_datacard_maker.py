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
parser.add_argument("--ratio", action="store_true", default=False,
                    help="Plot ratio or not")
parser.add_argument("--linear", action="store_true", default=False,
                    help="Plot linearly")
parser.add_argument("--png", action="store_true", default=False,
                    help="png format")
parser.add_argument("--norm", action="store_true", default=False,
                    help="normalising test")
parser.add_argument("--nodata", action="store_true", default=False,
                    help="Do not plot data")
parser.add_argument("--year", type=string, default="2016",
                    help="Select year of process to run")
parser.add_argument("--channel", type=string, default="btagMM",
                    help="Select year of process to run")

# Use like:
# python higgssearch/fromJF/hist_fromJFwqq_syst_total.py --stack --ratio --png --norm --year="all" --channel="btagMM_chitest_sl"

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

if (args.channel == "lepton50" or args.channel == "btagMM" or args.channel == "nobtag"): term_path = str(args.channel) 
elif args.channel == "lepton50_chitest": term_path = "lepton50/chi_test"
elif args.channel == "btagMM_chitest": term_path = "btagMM/chi_test"
elif args.channel == "lepton50_chitest_sl": term_path = "lepton50/chi_test/sl"
elif args.channel == "btagMM_chitest_sl": term_path = "btagMM/chi_test/sl"
elif args.channel == "lepton50_chitest_slss": term_path = "lepton50/chi_test/sl/ss"
elif args.channel == "btagMM_chitest_slss": term_path = "btagMM/chi_test/sl/ss"
elif args.channel == "lepton50_chitest_slos": term_path = "lepton50/chi_test/sl/os"
elif args.channel == "btagMM_chitest_slos": term_path = "btagMM/chi_test/sl/os"
elif args.channel == "lepton50_chitest_slssos": term_path = "lepton50/chi_test/sl/ssos"
elif args.channel == "btagMM_chitest_slssos": term_path = "btagMM/chi_test/sl/ssos"
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

observable_names = ["chi2_test_good"]

datayears = ["2016","2016B","2017","2018"]
#datayears = ["2018","2016","2016B"]

samplesHT = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_light","ttbar_sl_bottom","ttbar_dl","ttbar_dh","zz","wz",
        "st_1","st_2","st_3","st_4", "ttbar_sl_else", "ttbar_sl_charmgluon","ttbar_sl_bottomgluon"]

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
        lumi[data_op]["ttbar_sl_light"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_bottom"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_else"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_charmgluon"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_bottomgluon"] = lumi[data_op]["ttbar_sl"]

listsampl = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "ttbar_sl","ttbar_dl","ttbar_dh","zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6",
        "zjets_7","zjets_8","st_1","st_2","st_3","st_4","zz","wz", "ttbar_sl_charm", "ttbar_sl_charmgluon",
        "ttbar_sl_bottom", "ttbar_sl_bottomgluon", "ttbar_sl_else", "ttbar_sl_light"]

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
path_histfile = "/nfs/cms/vazqueze/higgssearch/datacards/"

for name in observable_names:
  ## Open hists files
  filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/"+term_path+"/"
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

  ## HISTS
  samples_foryear = {}
  hist_nom_M = {}
  hist_nom_E = {}
  hist_btaglightup_M = {}
  hist_btaglightup_E = {}
  hist_btaglightdown_M = {}
  hist_btaglightdown_E = {}
  hist_btagheavyup_M = {}
  hist_btagheavyup_E = {}
  hist_btagheavydown_M = {}
  hist_btagheavydown_E = {}
  histdata_M = {}
  histdata_E = {}
  for data_op in datayears:
    hist_nom_M[data_op] = {}
    hist_nom_E[data_op] = {}
    hist_btaglightup_M[data_op] = {}
    hist_btaglightup_E[data_op] = {}
    hist_btaglightdown_M[data_op] = {}
    hist_btaglightdown_E[data_op] = {}
    hist_btagheavyup_M[data_op] = {}
    hist_btagheavyup_E[data_op] = {}
    hist_btagheavydown_M[data_op] = {}
    hist_btagheavydown_E[data_op] = {}
    data_term = data_op
    #print(data_op+"M_"+name+"_M")
    #print(histFile[data_op][name].ls())
    for s in samples:
      if s[0:8] == "ttbar_sl": 
         s_term = s[0:8]+data_term+s[8:] 
      else:
         s_term = s+data_term
      hist_nom_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M")
      hist_nom_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E")
      hist_btaglightup_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_btaglightup")
      hist_btaglightup_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_btaglightup")
      hist_btaglightdown_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_btaglightdown")
      hist_btaglightdown_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_btaglightdown")
      hist_btagheavyup_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_btagheavyup")
      hist_btagheavyup_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_btagheavyup")
      hist_btagheavydown_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_btagheavydown")
      hist_btagheavydown_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_btagheavydown")
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
      hist_btaglightup_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_btaglightup_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_btaglightdown_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_btaglightdown_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_btagheavyup_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_btagheavyup_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_btagheavydown_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_btagheavydown_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
      if args.norm: hist_nom_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
      if args.norm: hist_nom_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
      if args.norm: hist_btaglightup_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
      if args.norm: hist_btaglightup_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
      if args.norm: hist_btaglightdown_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
      if args.norm: hist_btaglightdown_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
      if args.norm: hist_btagheavyup_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
      if args.norm: hist_btagheavyup_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
      if args.norm: hist_btagheavydown_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
      if args.norm: hist_btagheavydown_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
    ## Fixing single top
    #print(samples_foryear[data_op]) 
    #### List of summing samples:
    list_st = ["st_1","st_2","st_3","st_4"]
    list_wjets = ["wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8"]
    list_zjets = ["zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8"]
    list_vv = ["ww","wz","zz"]
    hist_nom_M[data_op]["st"] = hist_nom_M[data_op][list_st[0]]
    hist_nom_E[data_op]["st"] = hist_nom_E[data_op][list_st[0]]
    hist_btaglightup_M[data_op]["st"] = hist_btaglightup_M[data_op][list_st[0]]
    hist_btaglightup_E[data_op]["st"] = hist_btaglightup_E[data_op][list_st[0]]
    hist_btaglightdown_M[data_op]["st"] = hist_btaglightdown_M[data_op][list_st[0]]
    hist_btaglightdown_E[data_op]["st"] = hist_btaglightdown_E[data_op][list_st[0]]
    hist_btagheavyup_M[data_op]["st"] = hist_btagheavyup_M[data_op][list_st[0]]
    hist_btagheavyup_E[data_op]["st"] = hist_btagheavyup_E[data_op][list_st[0]]
    hist_btagheavydown_M[data_op]["st"] = hist_btagheavydown_M[data_op][list_st[0]]
    hist_btagheavydown_E[data_op]["st"] = hist_btagheavydown_E[data_op][list_st[0]]
    for l in list_st[1:]:
      hist_nom_M[data_op]["st"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["st"].Add(hist_nom_E[data_op][l])
      hist_btaglightup_M[data_op]["st"].Add(hist_btaglightup_M[data_op][l])
      hist_btaglightup_E[data_op]["st"].Add(hist_btaglightup_E[data_op][l])
      hist_btaglightdown_M[data_op]["st"].Add(hist_btaglightdown_M[data_op][l])
      hist_btaglightdown_E[data_op]["st"].Add(hist_btaglightdown_E[data_op][l])
      hist_btagheavyup_M[data_op]["st"].Add(hist_btagheavyup_M[data_op][l])
      hist_btagheavyup_E[data_op]["st"].Add(hist_btagheavyup_E[data_op][l])
      hist_btagheavydown_M[data_op]["st"].Add(hist_btagheavydown_M[data_op][l])
      hist_btagheavydown_E[data_op]["st"].Add(hist_btagheavydown_E[data_op][l])
    hist_nom_M[data_op]["wjets"] = hist_nom_M[data_op][list_wjets[0]]
    hist_nom_E[data_op]["wjets"] = hist_nom_E[data_op][list_wjets[0]]
    hist_btaglightup_M[data_op]["wjets"] = hist_btaglightup_M[data_op][list_wjets[0]]
    hist_btaglightup_E[data_op]["wjets"] = hist_btaglightup_E[data_op][list_wjets[0]]
    hist_btaglightdown_M[data_op]["wjets"] = hist_btaglightdown_M[data_op][list_wjets[0]]
    hist_btaglightdown_E[data_op]["wjets"] = hist_btaglightdown_E[data_op][list_wjets[0]]
    hist_btagheavyup_M[data_op]["wjets"] = hist_btagheavyup_M[data_op][list_wjets[0]]
    hist_btagheavyup_E[data_op]["wjets"] = hist_btagheavyup_E[data_op][list_wjets[0]]
    hist_btagheavydown_M[data_op]["wjets"] = hist_btagheavydown_M[data_op][list_wjets[0]]
    hist_btagheavydown_E[data_op]["wjets"] = hist_btagheavydown_E[data_op][list_wjets[0]]
    for l in list_wjets[1:]:
      hist_nom_M[data_op]["wjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["wjets"].Add(hist_nom_E[data_op][l])
      hist_btaglightup_M[data_op]["wjets"].Add(hist_btaglightup_M[data_op][l])
      hist_btaglightup_E[data_op]["wjets"].Add(hist_btaglightup_E[data_op][l])
      hist_btaglightdown_M[data_op]["wjets"].Add(hist_btaglightdown_M[data_op][l])
      hist_btaglightdown_E[data_op]["wjets"].Add(hist_btaglightdown_E[data_op][l])
      hist_btagheavyup_M[data_op]["wjets"].Add(hist_btagheavyup_M[data_op][l])
      hist_btagheavyup_E[data_op]["wjets"].Add(hist_btagheavyup_E[data_op][l])
      hist_btagheavydown_M[data_op]["wjets"].Add(hist_btagheavydown_M[data_op][l])
      hist_btagheavydown_E[data_op]["wjets"].Add(hist_btagheavydown_E[data_op][l])
    hist_nom_M[data_op]["zjets"] = hist_nom_M[data_op][list_zjets[0]]
    hist_nom_E[data_op]["zjets"] = hist_nom_E[data_op][list_zjets[0]]
    hist_btaglightup_M[data_op]["zjets"] = hist_btaglightup_M[data_op][list_zjets[0]]
    hist_btaglightup_E[data_op]["zjets"] = hist_btaglightup_E[data_op][list_zjets[0]]
    hist_btaglightdown_M[data_op]["zjets"] = hist_btaglightdown_M[data_op][list_zjets[0]]
    hist_btaglightdown_E[data_op]["zjets"] = hist_btaglightdown_E[data_op][list_zjets[0]]
    hist_btagheavyup_M[data_op]["zjets"] = hist_btagheavyup_M[data_op][list_zjets[0]]
    hist_btagheavyup_E[data_op]["zjets"] = hist_btagheavyup_E[data_op][list_zjets[0]]
    hist_btagheavydown_M[data_op]["zjets"] = hist_btagheavydown_M[data_op][list_zjets[0]]
    hist_btagheavydown_E[data_op]["zjets"] = hist_btagheavydown_E[data_op][list_zjets[0]]
    for l in list_zjets[1:]:
      hist_nom_M[data_op]["zjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["zjets"].Add(hist_nom_E[data_op][l])
      hist_btaglightup_M[data_op]["zjets"].Add(hist_btaglightup_M[data_op][l])
      hist_btaglightup_E[data_op]["zjets"].Add(hist_btaglightup_E[data_op][l])
      hist_btaglightdown_M[data_op]["zjets"].Add(hist_btaglightdown_M[data_op][l])
      hist_btaglightdown_E[data_op]["zjets"].Add(hist_btaglightdown_E[data_op][l])
      hist_btagheavyup_M[data_op]["zjets"].Add(hist_btagheavyup_M[data_op][l])
      hist_btagheavyup_E[data_op]["zjets"].Add(hist_btagheavyup_E[data_op][l])
      hist_btagheavydown_M[data_op]["zjets"].Add(hist_btagheavydown_M[data_op][l])
      hist_btagheavydown_E[data_op]["zjets"].Add(hist_btagheavydown_E[data_op][l])
    hist_nom_M[data_op]["vv"] = hist_nom_M[data_op][list_vv[0]]
    hist_nom_E[data_op]["vv"] = hist_nom_E[data_op][list_vv[0]]
    hist_btaglightup_M[data_op]["vv"] = hist_btaglightup_M[data_op][list_vv[0]]
    hist_btaglightup_E[data_op]["vv"] = hist_btaglightup_E[data_op][list_vv[0]]
    hist_btaglightdown_M[data_op]["vv"] = hist_btaglightdown_M[data_op][list_vv[0]]
    hist_btaglightdown_E[data_op]["vv"] = hist_btaglightdown_E[data_op][list_vv[0]]
    hist_btagheavyup_M[data_op]["vv"] = hist_btagheavyup_M[data_op][list_vv[0]]
    hist_btagheavyup_E[data_op]["vv"] = hist_btagheavyup_E[data_op][list_vv[0]]
    hist_btagheavydown_M[data_op]["vv"] = hist_btagheavydown_M[data_op][list_vv[0]]
    hist_btagheavydown_E[data_op]["vv"] = hist_btagheavydown_E[data_op][list_vv[0]]
    for l in list_vv[1:]:
      hist_nom_M[data_op]["vv"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["vv"].Add(hist_nom_E[data_op][l])
      hist_btaglightup_M[data_op]["vv"].Add(hist_btaglightup_M[data_op][l])
      hist_btaglightup_E[data_op]["vv"].Add(hist_btaglightup_E[data_op][l])
      hist_btaglightdown_M[data_op]["vv"].Add(hist_btaglightdown_M[data_op][l])
      hist_btaglightdown_E[data_op]["vv"].Add(hist_btaglightdown_E[data_op][l])
      hist_btagheavyup_M[data_op]["vv"].Add(hist_btagheavyup_M[data_op][l])
      hist_btagheavyup_E[data_op]["vv"].Add(hist_btagheavyup_E[data_op][l])
      hist_btagheavydown_M[data_op]["vv"].Add(hist_btagheavydown_M[data_op][l])
      hist_btagheavydown_E[data_op]["vv"].Add(hist_btagheavydown_E[data_op][l])

  samples = ["ttbar_sl_bottom","ttbar_sl_charm","ttbar_sl_else","ttbar_sl_light","ttbar_dl","ttbar_dh","zjets","vv","st","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon"]

  ## sumamos todos los anos para todos los histogramas
  histT_nom_M = {}
  histT_nom_E = {}
  histT_btaglightup_M = {}
  histT_btaglightup_E = {}
  histT_btaglightdown_M = {}
  histT_btaglightdown_E = {}
  histT_btagheavyup_M = {}
  histT_btagheavyup_E = {}
  histT_btagheavydown_M = {}
  histT_btagheavydown_E = {}
  for s in samples:
       histT_nom_M[s] = hist_nom_M[datayears[0]][s]
       histT_nom_E[s] = hist_nom_E[datayears[0]][s]
       histT_btaglightup_M[s] = hist_btaglightup_M[datayears[0]][s]
       histT_btaglightup_E[s] = hist_btaglightup_E[datayears[0]][s]
       histT_btaglightdown_M[s] = hist_btaglightdown_M[datayears[0]][s]
       histT_btaglightdown_E[s] = hist_btaglightdown_E[datayears[0]][s]
       histT_btagheavyup_M[s] = hist_btagheavyup_M[datayears[0]][s]
       histT_btagheavyup_E[s] = hist_btagheavyup_E[datayears[0]][s]
       histT_btagheavydown_M[s] = hist_btagheavydown_M[datayears[0]][s]
       histT_btagheavydown_E[s] = hist_btagheavydown_E[datayears[0]][s]
       for d in datayears[1:]:
          histT_nom_M[s].Add(hist_nom_M[d][s])
          histT_nom_E[s].Add(hist_nom_E[d][s])
          histT_btaglightup_M[s].Add(hist_btaglightup_M[d][s])
          histT_btaglightup_E[s].Add(hist_btaglightup_E[d][s])
          histT_btaglightdown_M[s].Add(hist_btaglightdown_M[d][s])
          histT_btaglightdown_E[s].Add(hist_btaglightdown_E[d][s])
          histT_btagheavyup_M[s].Add(hist_btagheavyup_M[d][s])
          histT_btagheavyup_E[s].Add(hist_btagheavyup_E[d][s])
          histT_btagheavydown_M[s].Add(hist_btagheavydown_M[d][s])
          histT_btagheavydown_E[s].Add(hist_btagheavydown_E[d][s])
       if (args.channel in sl_channel) and (name not in not_rebin):
          histT_nom_M[s].Rebin(2)
          histT_nom_E[s].Rebin(2)
          histT_btaglightup_M[s].Rebin(2)
          histT_btaglightup_E[s].Rebin(2)
          histT_btaglightdown_M[s].Rebin(2)
          histT_btaglightdown_E[s].Rebin(2)
          histT_btagheavyup_M[s].Rebin(2)
          histT_btagheavyup_E[s].Rebin(2)
          histT_btagheavydown_M[s].Rebin(2)
          histT_btagheavydown_E[s].Rebin(2)

  if not args.nodata:
    histD_M = histdata_M[datayears[0]]
    histD_E = histdata_E[datayears[0]]
    for d in datayears[1:]:
       histD_M.Add(histdata_M[d])
       histD_E.Add(histdata_E[d])
    if (args.channel in sl_channel) and (name not in not_rebin): 
       histD_M.Rebin(2)
       histD_E.Rebin(2)

  envHi_btaglight_T = {}
  envLo_btaglight_T = {}
  envHi_btagheavy_T = {}
  envLo_btagheavy_T = {}
  for s in samples:
      envHi_btaglight_T[s] = histT_btaglightup_M[s]
      envHi_btaglight_T[s].Add(histT_btaglightup_E[s])
      envLo_btaglight_T[s] = histT_btaglightdown_M[s] 
      envLo_btaglight_T[s].Add(histT_btaglightdown_E[s])
      envHi_btagheavy_T[s] = histT_btagheavyup_M[s] 
      envHi_btagheavy_T[s].Add(histT_btagheavyup_E[s])
      envLo_btagheavy_T[s] = histT_btagheavydown_M[s] 
      envLo_btagheavy_T[s].Add(histT_btagheavydown_E[s])

  ## Stack creation
  samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon","ttbar_sl_else","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]

  histT_nom_T = {}
  for s in samples:
    histT_nom_T[s] = histT_nom_M[s]
    histT_nom_T[s].Add(histT_nom_E[s])
  stack_T = ROOT.THStack()
  for s in samples:
    stack_T.Add(histT_nom_T[s])
  last_T = stack_T.GetStack().Last()

  if not args.nodata:
    histD_T = histD_M
    histD_T.Add(histD_E)

  ## Samples reminder
  ## samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottomgluon",
  ## "ttbar_sl_charmgluon","ttbar_sl_else","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]

  ### hist renamings
  data_obs = histD_T.Clone("data_obs")
  mchists = {}
  mchists_lightup = {}
  mchists_lightdown = {}
  mchists_heavyup = {}
  mchists_heavydown = {}
  for s in samples:
      if s == "ttbar_sl_charm":
         mchists[s] = histT_nom_T[s].Clone("signal")
         mchists_lightup[s] = envHi_btaglight_T[s].Clone("signal_lightUp")
         mchists_lightdown[s] = envLo_btaglight_T[s].Clone("signal_lightDown")
         mchists_heavyup[s] = envHi_btagheavy_T[s].Clone("signal_heavyUp")
         mchists_heavydown[s] = envLo_btagheavy_T[s].Clone("signal_heavyDown")
      else:
         mchists[s] = histT_nom_T[s].Clone(str(s))
         mchists_lightup[s] = envHi_btaglight_T[s].Clone(s+"_lightUp")
         mchists_lightdown[s] = envLo_btaglight_T[s].Clone(s+"_lightDown")
         mchists_heavyup[s] = envHi_btagheavy_T[s].Clone(s+"_heavyUp")
         mchists_heavydown[s] = envLo_btagheavy_T[s].Clone(s+"_heavyDown")

  print("Integral of data is "+str(data_obs.Integral()))
  for s in samples:
     print("Integral of "+str(s)+" MC process is "+str(mchists[s].Integral()))

  ### file saving
  termfile = "mydatacard_"+str(args.channel)+"_"+str(name)+".root"
  myfile[name] = TFile( path_histfile+termfile, 'RECREATE' )
  data_obs.Write()
  for s in samples:
      mchists[s].Write()
      mchists_lightup[s].Write()
      mchists_lightdown[s].Write()
      mchists_heavyup[s].Write()
      mchists_heavydown[s].Write()
  myfile[name].Close()

  for data_op in datayears:
         histFile[name][data_op].Close()
         if not args.nodata:
            histFileDM[name][data_op].Close()
            histFileDE[name][data_op].Close()


