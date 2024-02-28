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
parser.add_argument("--wcs", action="store_true", default=False,
                    help="classification for ttbar and st")
parser.add_argument("--year", type=string, default="2016",
                    help="Select year of process to run")
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
elif args.channel == "btagMM_chitest_ctag": term_path = "/ctag"
elif args.channel == "btagMM_chitest_noctag": term_path = "/noctag"
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

norm_factorM["all"]["btagMM_chitest_ctag"] = 0.92; norm_factorE["all"]["btagMM_chitest_ctag"] = 0.9; 
norm_factorM["all"]["btagMM_chitest_noctag"] = 0.92; norm_factorE["all"]["btagMM_chitest_noctag"] = 0.9; 

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

  if args.ratio and not args.nodata:
    ## In case of ratio plot
    upper_pad = ROOT.TPad("upper_pad", "", 0, 0.35, 1, 1)
    lower_pad = ROOT.TPad("lower_pad", "", 0, 0, 1, 0.35)
    for p in [upper_pad, lower_pad]:
        p.SetLeftMargin(0.14)
        p.SetRightMargin(0.05)
        p.SetTickx(False)
        p.SetTicky(False)
    upper_pad.SetBottomMargin(0)
    lower_pad.SetTopMargin(0)
    lower_pad.SetBottomMargin(0.3)
 
    if not args.linear: upper_pad.SetLogy()
    upper_pad.Draw()
    lower_pad.Draw()

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
  hist_btaglightup_M = {}
  hist_btaglightup_E = {}
  hist_btaglightdown_M = {}
  hist_btaglightdown_E = {}
  hist_btagheavyup_M = {}
  hist_btagheavyup_E = {}
  hist_btagheavydown_M = {}
  hist_btagheavydown_E = {}
  hist_seclepup_M = {}
  hist_seclepup_E = {}
  hist_seclepdown_M = {}
  hist_seclepdown_E = {}
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
    hist_seclepup_M[data_op] = {}
    hist_seclepup_E[data_op] = {}
    hist_seclepdown_M[data_op] = {}
    hist_seclepdown_E[data_op] = {}
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
         hist_btaglightup_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_btaglightup")
         hist_btaglightup_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_btaglightup")
         hist_btaglightdown_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_btaglightdown")
         hist_btaglightdown_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_btaglightdown")
         hist_btagheavyup_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_btagheavyup")
         hist_btagheavyup_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_btagheavyup")
         hist_btagheavydown_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_btagheavydown")
         hist_btagheavydown_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_btagheavydown")
         hist_seclepup_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_jecup")
         hist_seclepup_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_jecup")
         hist_seclepdown_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_jecdown")
         hist_seclepdown_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_jecdown")
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
         hist_btaglightup_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_btaglightup_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_btaglightdown_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_btaglightdown_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_btagheavyup_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_btagheavyup_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_btagheavydown_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_btagheavydown_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_seclepup_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_seclepup_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_seclepdown_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_seclepdown_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
      if args.norm: hist_nom_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
      if args.norm: hist_nom_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
      if not args.nosyst:
         if args.norm: hist_btaglightup_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
         if args.norm: hist_btaglightup_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
         if args.norm: hist_btaglightdown_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
         if args.norm: hist_btaglightdown_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
         if args.norm: hist_btagheavyup_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
         if args.norm: hist_btagheavyup_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
         if args.norm: hist_btagheavydown_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
         if args.norm: hist_btagheavydown_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
         if args.norm: hist_seclepup_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
         if args.norm: hist_seclepup_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
         if args.norm: hist_seclepdown_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
         if args.norm: hist_seclepdown_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
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
        hist_btaglightup_M[data_op]["st"+apx] = hist_btaglightup_M[data_op][list_st[0]+apx]
        hist_btaglightup_E[data_op]["st"+apx] = hist_btaglightup_E[data_op][list_st[0]+apx]
        hist_btaglightdown_M[data_op]["st"+apx] = hist_btaglightdown_M[data_op][list_st[0]+apx]
        hist_btaglightdown_E[data_op]["st"+apx] = hist_btaglightdown_E[data_op][list_st[0]+apx]
        hist_btagheavyup_M[data_op]["st"+apx] = hist_btagheavyup_M[data_op][list_st[0]+apx]
        hist_btagheavyup_E[data_op]["st"+apx] = hist_btagheavyup_E[data_op][list_st[0]+apx]
        hist_btagheavydown_M[data_op]["st"+apx] = hist_btagheavydown_M[data_op][list_st[0]+apx]
        hist_btagheavydown_E[data_op]["st"+apx] = hist_btagheavydown_E[data_op][list_st[0]+apx]
        hist_seclepup_M[data_op]["st"+apx] = hist_seclepup_M[data_op][list_st[0]+apx]
        hist_seclepup_E[data_op]["st"+apx] = hist_seclepup_E[data_op][list_st[0]+apx]
        hist_seclepdown_M[data_op]["st"+apx] = hist_seclepdown_M[data_op][list_st[0]+apx]
        hist_seclepdown_E[data_op]["st"+apx] = hist_seclepdown_E[data_op][list_st[0]+apx]
      for l in list_st[1:]:
        hist_nom_M[data_op]["st"+apx].Add(hist_nom_M[data_op][l+apx])
        hist_nom_E[data_op]["st"+apx].Add(hist_nom_E[data_op][l+apx])
        if not args.nosyst:
          hist_btaglightup_M[data_op]["st"+apx].Add(hist_btaglightup_M[data_op][l+apx])
          hist_btaglightup_E[data_op]["st"+apx].Add(hist_btaglightup_E[data_op][l+apx])
          hist_btaglightdown_M[data_op]["st"+apx].Add(hist_btaglightdown_M[data_op][l+apx])
          hist_btaglightdown_E[data_op]["st"+apx].Add(hist_btaglightdown_E[data_op][l+apx])
          hist_btagheavyup_M[data_op]["st"+apx].Add(hist_btagheavyup_M[data_op][l+apx])
          hist_btagheavyup_E[data_op]["st"+apx].Add(hist_btagheavyup_E[data_op][l+apx])
          hist_btagheavydown_M[data_op]["st"+apx].Add(hist_btagheavydown_M[data_op][l+apx])
          hist_btagheavydown_E[data_op]["st"+apx].Add(hist_btagheavydown_E[data_op][l+apx])
          hist_seclepup_M[data_op]["st"+apx].Add(hist_seclepup_M[data_op][l+apx])
          hist_seclepup_E[data_op]["st"+apx].Add(hist_seclepup_E[data_op][l+apx])
          hist_seclepdown_M[data_op]["st"+apx].Add(hist_seclepdown_M[data_op][l+apx])
          hist_seclepdown_E[data_op]["st"+apx].Add(hist_seclepdown_E[data_op][l+apx])
    hist_nom_M[data_op]["wjets"] = hist_nom_M[data_op][list_wjets[0]]
    hist_nom_E[data_op]["wjets"] = hist_nom_E[data_op][list_wjets[0]]
    if not args.nosyst:
      hist_btaglightup_M[data_op]["wjets"] = hist_btaglightup_M[data_op][list_wjets[0]]
      hist_btaglightup_E[data_op]["wjets"] = hist_btaglightup_E[data_op][list_wjets[0]]
      hist_btaglightdown_M[data_op]["wjets"] = hist_btaglightdown_M[data_op][list_wjets[0]]
      hist_btaglightdown_E[data_op]["wjets"] = hist_btaglightdown_E[data_op][list_wjets[0]]
      hist_btagheavyup_M[data_op]["wjets"] = hist_btagheavyup_M[data_op][list_wjets[0]]
      hist_btagheavyup_E[data_op]["wjets"] = hist_btagheavyup_E[data_op][list_wjets[0]]
      hist_btagheavydown_M[data_op]["wjets"] = hist_btagheavydown_M[data_op][list_wjets[0]]
      hist_btagheavydown_E[data_op]["wjets"] = hist_btagheavydown_E[data_op][list_wjets[0]]
      hist_seclepup_M[data_op]["wjets"] = hist_seclepup_M[data_op][list_wjets[0]]
      hist_seclepup_E[data_op]["wjets"] = hist_seclepup_E[data_op][list_wjets[0]]
      hist_seclepdown_M[data_op]["wjets"] = hist_seclepdown_M[data_op][list_wjets[0]]
      hist_seclepdown_E[data_op]["wjets"] = hist_seclepdown_E[data_op][list_wjets[0]]
    for l in list_wjets[1:]:
      hist_nom_M[data_op]["wjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["wjets"].Add(hist_nom_E[data_op][l])
      if not args.nosyst:
        hist_btaglightup_M[data_op]["wjets"].Add(hist_btaglightup_M[data_op][l])
        hist_btaglightup_E[data_op]["wjets"].Add(hist_btaglightup_E[data_op][l])
        hist_btaglightdown_M[data_op]["wjets"].Add(hist_btaglightdown_M[data_op][l])
        hist_btaglightdown_E[data_op]["wjets"].Add(hist_btaglightdown_E[data_op][l])
        hist_btagheavyup_M[data_op]["wjets"].Add(hist_btagheavyup_M[data_op][l])
        hist_btagheavyup_E[data_op]["wjets"].Add(hist_btagheavyup_E[data_op][l])
        hist_btagheavydown_M[data_op]["wjets"].Add(hist_btagheavydown_M[data_op][l])
        hist_btagheavydown_E[data_op]["wjets"].Add(hist_btagheavydown_E[data_op][l])
        hist_seclepup_M[data_op]["wjets"].Add(hist_seclepup_M[data_op][l])
        hist_seclepup_E[data_op]["wjets"].Add(hist_seclepup_E[data_op][l])
        hist_seclepdown_M[data_op]["wjets"].Add(hist_seclepdown_M[data_op][l])
        hist_seclepdown_E[data_op]["wjets"].Add(hist_seclepdown_E[data_op][l])
    hist_nom_M[data_op]["zjets"] = hist_nom_M[data_op][list_zjets[0]]
    hist_nom_E[data_op]["zjets"] = hist_nom_E[data_op][list_zjets[0]]
    if not args.nosyst:
      hist_btaglightup_M[data_op]["zjets"] = hist_btaglightup_M[data_op][list_zjets[0]]
      hist_btaglightup_E[data_op]["zjets"] = hist_btaglightup_E[data_op][list_zjets[0]]
      hist_btaglightdown_M[data_op]["zjets"] = hist_btaglightdown_M[data_op][list_zjets[0]]
      hist_btaglightdown_E[data_op]["zjets"] = hist_btaglightdown_E[data_op][list_zjets[0]]
      hist_btagheavyup_M[data_op]["zjets"] = hist_btagheavyup_M[data_op][list_zjets[0]]
      hist_btagheavyup_E[data_op]["zjets"] = hist_btagheavyup_E[data_op][list_zjets[0]]
      hist_btagheavydown_M[data_op]["zjets"] = hist_btagheavydown_M[data_op][list_zjets[0]]
      hist_btagheavydown_E[data_op]["zjets"] = hist_btagheavydown_E[data_op][list_zjets[0]]
      hist_seclepup_M[data_op]["zjets"] = hist_seclepup_M[data_op][list_zjets[0]]
      hist_seclepup_E[data_op]["zjets"] = hist_seclepup_E[data_op][list_zjets[0]]
      hist_seclepdown_M[data_op]["zjets"] = hist_seclepdown_M[data_op][list_zjets[0]]
      hist_seclepdown_E[data_op]["zjets"] = hist_seclepdown_E[data_op][list_zjets[0]]
    for l in list_zjets[1:]:
      hist_nom_M[data_op]["zjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["zjets"].Add(hist_nom_E[data_op][l])
      if not args.nosyst:
        hist_btaglightup_M[data_op]["zjets"].Add(hist_btaglightup_M[data_op][l])
        hist_btaglightup_E[data_op]["zjets"].Add(hist_btaglightup_E[data_op][l])
        hist_btaglightdown_M[data_op]["zjets"].Add(hist_btaglightdown_M[data_op][l])
        hist_btaglightdown_E[data_op]["zjets"].Add(hist_btaglightdown_E[data_op][l])
        hist_btagheavyup_M[data_op]["zjets"].Add(hist_btagheavyup_M[data_op][l])
        hist_btagheavyup_E[data_op]["zjets"].Add(hist_btagheavyup_E[data_op][l])
        hist_btagheavydown_M[data_op]["zjets"].Add(hist_btagheavydown_M[data_op][l])
        hist_btagheavydown_E[data_op]["zjets"].Add(hist_btagheavydown_E[data_op][l])
        hist_seclepup_M[data_op]["zjets"].Add(hist_seclepup_M[data_op][l])
        hist_seclepup_E[data_op]["zjets"].Add(hist_seclepup_E[data_op][l])
        hist_seclepdown_M[data_op]["zjets"].Add(hist_seclepdown_M[data_op][l])
        hist_seclepdown_E[data_op]["zjets"].Add(hist_seclepdown_E[data_op][l])
    hist_nom_M[data_op]["vv"] = hist_nom_M[data_op][list_vv[0]]
    hist_nom_E[data_op]["vv"] = hist_nom_E[data_op][list_vv[0]]
    if not args.nosyst:
      hist_btaglightup_M[data_op]["vv"] = hist_btaglightup_M[data_op][list_vv[0]]
      hist_btaglightup_E[data_op]["vv"] = hist_btaglightup_E[data_op][list_vv[0]]
      hist_btaglightdown_M[data_op]["vv"] = hist_btaglightdown_M[data_op][list_vv[0]]
      hist_btaglightdown_E[data_op]["vv"] = hist_btaglightdown_E[data_op][list_vv[0]]
      hist_btagheavyup_M[data_op]["vv"] = hist_btagheavyup_M[data_op][list_vv[0]]
      hist_btagheavyup_E[data_op]["vv"] = hist_btagheavyup_E[data_op][list_vv[0]]
      hist_btagheavydown_M[data_op]["vv"] = hist_btagheavydown_M[data_op][list_vv[0]]
      hist_btagheavydown_E[data_op]["vv"] = hist_btagheavydown_E[data_op][list_vv[0]]
      hist_seclepup_M[data_op]["vv"] = hist_seclepup_M[data_op][list_vv[0]]
      hist_seclepup_E[data_op]["vv"] = hist_seclepup_E[data_op][list_vv[0]]
      hist_seclepdown_M[data_op]["vv"] = hist_seclepdown_M[data_op][list_vv[0]]
      hist_seclepdown_E[data_op]["vv"] = hist_seclepdown_E[data_op][list_vv[0]]
    for l in list_vv[1:]:
      hist_nom_M[data_op]["vv"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["vv"].Add(hist_nom_E[data_op][l])
      if not args.nosyst:
        hist_btaglightup_M[data_op]["vv"].Add(hist_btaglightup_M[data_op][l])
        hist_btaglightup_E[data_op]["vv"].Add(hist_btaglightup_E[data_op][l])
        hist_btaglightdown_M[data_op]["vv"].Add(hist_btaglightdown_M[data_op][l])
        hist_btaglightdown_E[data_op]["vv"].Add(hist_btaglightdown_E[data_op][l])
        hist_btagheavyup_M[data_op]["vv"].Add(hist_btagheavyup_M[data_op][l])
        hist_btagheavyup_E[data_op]["vv"].Add(hist_btagheavyup_E[data_op][l])
        hist_btagheavydown_M[data_op]["vv"].Add(hist_btagheavydown_M[data_op][l])
        hist_btagheavydown_E[data_op]["vv"].Add(hist_btagheavydown_E[data_op][l])
        hist_seclepup_M[data_op]["vv"].Add(hist_seclepup_M[data_op][l])
        hist_seclepup_E[data_op]["vv"].Add(hist_seclepup_E[data_op][l])
        hist_seclepdown_M[data_op]["vv"].Add(hist_seclepdown_M[data_op][l])
        hist_seclepdown_E[data_op]["vv"].Add(hist_seclepdown_E[data_op][l])

  samples = ["ttbar_sl_bottom","ttbar_sl_charm","ttbar_sl_else","ttbar_sl_light","ttbar_dl","ttbar_dh","zjets","vv","st","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon"]
  if args.wcs: samples = ["ttbar_sl_nocharm","ttbar_sl_charm","ttbar_dl","ttbar_dh","zjets","vv","st_nocharm","st_charm","st_else","wjets"]

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
  histT_seclepup_M = {}
  histT_seclepup_E = {}
  histT_seclepdown_M = {}
  histT_seclepdown_E = {}
  for s in samples:
       histT_nom_M[s] = hist_nom_M[datayears[0]][s]
       histT_nom_E[s] = hist_nom_E[datayears[0]][s]
       if not args.nosyst:
         histT_btaglightup_M[s] = hist_btaglightup_M[datayears[0]][s]
         histT_btaglightup_E[s] = hist_btaglightup_E[datayears[0]][s]
         histT_btaglightdown_M[s] = hist_btaglightdown_M[datayears[0]][s]
         histT_btaglightdown_E[s] = hist_btaglightdown_E[datayears[0]][s]
         histT_btagheavyup_M[s] = hist_btagheavyup_M[datayears[0]][s]
         histT_btagheavyup_E[s] = hist_btagheavyup_E[datayears[0]][s]
         histT_btagheavydown_M[s] = hist_btagheavydown_M[datayears[0]][s]
         histT_btagheavydown_E[s] = hist_btagheavydown_E[datayears[0]][s]
         histT_seclepup_M[s] = hist_seclepup_M[datayears[0]][s]
         histT_seclepup_E[s] = hist_seclepup_E[datayears[0]][s]
         histT_seclepdown_M[s] = hist_seclepdown_M[datayears[0]][s]
         histT_seclepdown_E[s] = hist_seclepdown_E[datayears[0]][s]
       for d in datayears[1:]:
          histT_nom_M[s].Add(hist_nom_M[d][s])
          histT_nom_E[s].Add(hist_nom_E[d][s])
          if not args.nosyst:
            histT_btaglightup_M[s].Add(hist_btaglightup_M[d][s])
            histT_btaglightup_E[s].Add(hist_btaglightup_E[d][s])
            histT_btaglightdown_M[s].Add(hist_btaglightdown_M[d][s])
            histT_btaglightdown_E[s].Add(hist_btaglightdown_E[d][s])
            histT_btagheavyup_M[s].Add(hist_btagheavyup_M[d][s])
            histT_btagheavyup_E[s].Add(hist_btagheavyup_E[d][s])
            histT_btagheavydown_M[s].Add(hist_btagheavydown_M[d][s])
            histT_btagheavydown_E[s].Add(hist_btagheavydown_E[d][s])
            histT_seclepup_M[s].Add(hist_seclepup_M[d][s])
            histT_seclepup_E[s].Add(hist_seclepup_E[d][s])
            histT_seclepdown_M[s].Add(hist_seclepdown_M[d][s])
            histT_seclepdown_E[s].Add(hist_seclepdown_E[d][s])
       if (args.channel in sl_channel) and (name not in not_rebin):
          histT_nom_M[s].Rebin(nrebin)
          histT_nom_E[s].Rebin(nrebin)
          if not args.nosyst:
            histT_btaglightup_M[s].Rebin(nrebin)
            histT_btaglightup_E[s].Rebin(nrebin)
            histT_btaglightdown_M[s].Rebin(nrebin)
            histT_btaglightdown_E[s].Rebin(nrebin)
            histT_btagheavyup_M[s].Rebin(nrebin)
            histT_btagheavyup_E[s].Rebin(nrebin)
            histT_btagheavydown_M[s].Rebin(nrebin)
            histT_btagheavydown_E[s].Rebin(nrebin)
            histT_seclepup_M[s].Rebin(nrebin)
            histT_seclepup_E[s].Rebin(nrebin)
            histT_seclepdown_M[s].Rebin(nrebin)
            histT_seclepdown_E[s].Rebin(nrebin)

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
    ### Para los histogramas de sistemática up y down sumamos todas las contribuciones
    histT_sT_btaglightup_M = histT_btaglightup_M[samples[0]]
    histT_sT_btaglightup_E = histT_btaglightup_E[samples[0]]
    histT_sT_btaglightdown_M = histT_btaglightdown_M[samples[0]]
    histT_sT_btaglightdown_E = histT_btaglightdown_E[samples[0]]
    histT_sT_btagheavyup_M = histT_btagheavyup_M[samples[0]]
    histT_sT_btagheavyup_E = histT_btagheavyup_E[samples[0]]
    histT_sT_btagheavydown_M = histT_btagheavydown_M[samples[0]]
    histT_sT_btagheavydown_E = histT_btagheavydown_E[samples[0]]
    histT_sT_seclepup_M = histT_seclepup_M[samples[0]]
    histT_sT_seclepup_E = histT_seclepup_E[samples[0]]
    histT_sT_seclepdown_M = histT_seclepdown_M[samples[0]]
    histT_sT_seclepdown_E = histT_seclepdown_E[samples[0]]
    for s in samples[1:]:
      histT_sT_btaglightup_M.Add(histT_btaglightup_M[s])
      histT_sT_btaglightup_E.Add(histT_btaglightup_E[s])
      histT_sT_btaglightdown_M.Add(histT_btaglightdown_M[s])
      histT_sT_btaglightdown_E.Add(histT_btaglightdown_E[s])
      histT_sT_btagheavyup_M.Add(histT_btagheavyup_M[s])
      histT_sT_btagheavyup_E.Add(histT_btagheavyup_E[s])
      histT_sT_btagheavydown_M.Add(histT_btagheavydown_M[s])
      histT_sT_btagheavydown_E.Add(histT_btagheavydown_E[s])
      histT_sT_seclepup_M.Add(histT_seclepup_M[s])
      histT_sT_seclepup_E.Add(histT_seclepup_E[s])
      histT_sT_seclepdown_M.Add(histT_seclepdown_M[s])
      histT_sT_seclepdown_E.Add(histT_seclepdown_E[s])

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

  if args.stack:
    ymax_M = 0
    ymin_M = 0
    ymax_E = 0
    ymin_E = 0
    for s in samples:
      histT_nom_M[s].SetLineWidth(1)
      histT_nom_M[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_nom_M[s].GetYaxis().SetTitle("Number of events")
      histT_nom_M[s].GetXaxis().SetTitle(name)
      histT_nom_E[s].SetLineWidth(1)
      histT_nom_E[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_nom_E[s].GetYaxis().SetTitle("Number of events")
      histT_nom_E[s].GetXaxis().SetTitle(name)

      y = histT_nom_M[s].GetMaximum()
      ym = histT_nom_M[s].GetMinimum()
      if y>ymax_M: ymax_M=y
      if ymin_M>ym: ymin_M=ym
      if not args.nodata:
        y = histD_M.GetMaximum()
        if y>ymax_M: ymax_M=y

      y = histT_nom_E[s].GetMaximum()
      ym = histT_nom_E[s].GetMinimum()
      if y>ymax_E: ymax_E=y
      if ymin_E>ym: ymin_E=ym
      if not args.nodata:
        y = histD_E.GetMaximum()
        if y>ymax_E: ymax_E=y

      #histT_M["vv"].SetMinimum(1000.)
      #histT_M["vv"].SetMaximum(3*ymax_M)
      if args.linear: histT_nom_M["vv"].SetMaximum(1.3*ymax_M)
      if args.linear: histT_nom_M["vv"].SetMinimum(1.3*ymin_M)

      #histT_E["vv"].SetMinimum(1000.)
      #histT_E["vv"].SetMaximum(3*ymax_E)
      if args.linear: histT_nom_E["vv"].SetMaximum(1.3*ymax_E)
      if args.linear: histT_nom_E["vv"].SetMinimum(1.3*ymin_E)

    ## Stack creation
    samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon","ttbar_sl_else","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]
    if args.wcs: samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]

    if args.ratio and not args.nodata: upper_pad.cd()

    if args.sumEM:
       histT_nom_T = {}
       for s in samples:
         histT_nom_T[s] = histT_nom_M[s]
         histT_nom_T[s].Add(histT_nom_E[s])
       stack_T = ROOT.THStack()
       for s in samples:
         stack_T.Add(histT_nom_T[s])
       last_T = stack_T.GetStack().Last()

    else:
       stack_M = ROOT.THStack()
       stack_E = ROOT.THStack()
       for s in samples:
         stack_M.Add(histT_nom_M[s])
         stack_E.Add(histT_nom_E[s])

       last_M = stack_M.GetStack().Last()
       last_E = stack_E.GetStack().Last()

    if args.sumEM:
      if not args.nosyst:
        ### Creamos el tgraph de sistematicos
        envHi_btaglight_T = histT_sT_btaglightup_M;
        envLo_btaglight_T = histT_sT_btaglightdown_M;
        envHi_btagheavy_T = histT_sT_btagheavyup_M;
        envLo_btagheavy_T = histT_sT_btagheavydown_M;
        envHi_seclep_T = histT_sT_seclepup_M;
        envLo_seclep_T = histT_sT_seclepdown_M;
        envHi_btaglight_T.Add(histT_sT_btaglightup_E);
        envLo_btaglight_T.Add(histT_sT_btaglightdown_E);
        envHi_btagheavy_T.Add(histT_sT_btagheavyup_E);
        envLo_btagheavy_T.Add(histT_sT_btagheavydown_E);
        envHi_seclep_T.Add(histT_sT_seclepup_E);
        envLo_seclep_T.Add(histT_sT_seclepdown_E);

        graph_err_T = TGraphAsymmErrors();
        ratio_graph_err_T = TGraphAsymmErrors();

        for bin in range(last_T.GetNbinsX()):
          bin_nom_T = float(last_T.GetBinContent(bin+1));
          bin_btaglightup_T = float(envHi_btaglight_T.GetBinContent(bin+1));
          bin_btaglightdown_T = float(envLo_btaglight_T.GetBinContent(bin+1));
          bin_btagheavyup_T = float(envHi_btagheavy_T.GetBinContent(bin+1));
          bin_btagheavydown_T = float(envLo_btagheavy_T.GetBinContent(bin+1));
          bin_seclepup_T = float(envHi_seclep_T.GetBinContent(bin+1));
          bin_seclepdown_T = float(envLo_seclep_T.GetBinContent(bin+1));

          sysup_btaglight_T = abs(bin_nom_T-bin_btaglightup_T);
          sysdown_btaglight_T = abs(bin_nom_T-bin_btaglightdown_T);
          sysup_btagheavy_T = abs(bin_nom_T-bin_btagheavyup_T);
          sysdown_btagheavy_T = abs(bin_nom_T-bin_btagheavydown_T);
          sysup_seclep_T = abs(bin_nom_T-bin_seclepup_T);
          sysdown_seclep_T = abs(bin_nom_T-bin_seclepdown_T);
          sysup_total_T = sqrt(sysup_btaglight_T**2 + sysup_btagheavy_T**2);
          sysdown_total_T = sqrt(sysdown_btaglight_T**2 + sysdown_btagheavy_T**2);

          graph_err_T.SetPointEXhigh(bin,last_T.GetBinWidth(bin+1)/2);
          graph_err_T.SetPointEXlow(bin,last_T.GetBinWidth(bin+1)/2);
          graph_err_T.SetPoint(bin,last_T.GetBinCenter(bin+1),last_T.GetBinContent(bin+1));
          #graph_err_T.SetPointEYhigh(bin,sysup_total_T);
          #graph_err_T.SetPointEYlow(bin,sysdown_total_T);
          graph_err_T.SetPointEYhigh(bin,sysup_seclep_T);
          graph_err_T.SetPointEYlow(bin,sysdown_seclep_T);

          ratio_graph_err_T.SetPointEXhigh(bin,last_T.GetBinWidth(bin+1)/2);
          ratio_graph_err_T.SetPointEXlow(bin,last_T.GetBinWidth(bin+1)/2);
          ratio_graph_err_T.SetPoint(bin,last_T.GetBinCenter(bin+1),1.);
          if (last_T.GetBinContent(bin+1)>0.01):
             #ratio_graph_err_T.SetPointEYhigh(bin,abs(sysup_total_T/bin_nom_T));
             #ratio_graph_err_T.SetPointEYlow(bin,abs(sysdown_total_T/bin_nom_T));
             ratio_graph_err_T.SetPointEYhigh(bin,abs(sysup_seclep_T/bin_nom_T));
             ratio_graph_err_T.SetPointEYlow(bin,abs(sysdown_seclep_T/bin_nom_T));
          else:
             ratio_graph_err_T.SetPointEYhigh(bin,0.);
             ratio_graph_err_T.SetPointEYlow(bin,0.);

        graph_err_T.SetFillColorAlpha(kGray+2,0.7);
        graph_err_T.SetMarkerColor(kGray);
        graph_err_T.SetMarkerSize(0);
        graph_err_T.SetFillStyle(3001);
        ratio_graph_err_T.SetFillColorAlpha(kGray+2,0.7);
        ratio_graph_err_T.SetMarkerColor(kGray);
        ratio_graph_err_T.SetMarkerSize(0);
        ratio_graph_err_T.SetFillStyle(3001);
    else:
      if not args.nosyst:
        ### Creamos el tgraph de sistematicos
        envHi_btaglight_M = histT_sT_btaglightup_M;
        envHi_btaglight_E = histT_sT_btaglightup_E;
        envLo_btaglight_M = histT_sT_btaglightdown_M;
        envLo_btaglight_E = histT_sT_btaglightdown_E;
        envHi_btagheavy_M = histT_sT_btagheavyup_M;
        envHi_btagheavy_E = histT_sT_btagheavyup_E;
        envLo_btagheavy_M = histT_sT_btagheavydown_M;
        envLo_btagheavy_E = histT_sT_btagheavydown_E;

        graph_err_M = TGraphAsymmErrors();
        graph_err_E = TGraphAsymmErrors();
        ratio_graph_err_M = TGraphAsymmErrors();
        ratio_graph_err_E = TGraphAsymmErrors();

        for bin in range(last_M.GetNbinsX()):
          #sysup = std::sqrt((nom->GetBinContent(bin+1)-envHi->GetBinContent(bin+1))*(nom->GetBinContent(bin+1)-envHi->GetBinContent(bin+1)));
          #sysdown = std::sqrt((nom->GetBinContent(bin+1)-envLo->GetBinContent(bin+1))*(nom->GetBinContent(bin+1)-envLo->GetBinContent(bin+1)));
          bin_nom_M = float(last_M.GetBinContent(bin+1));
          bin_nom_E = float(last_E.GetBinContent(bin+1));
          bin_btaglightup_M = float(envHi_btaglight_M.GetBinContent(bin+1));
          bin_btaglightup_E = float(envHi_btaglight_E.GetBinContent(bin+1));
          bin_btaglightdown_M = float(envLo_btaglight_M.GetBinContent(bin+1));
          bin_btaglightdown_E = float(envLo_btaglight_E.GetBinContent(bin+1));
          bin_btagheavyup_M = float(envHi_btagheavy_M.GetBinContent(bin+1));
          bin_btagheavyup_E = float(envHi_btagheavy_E.GetBinContent(bin+1));
          bin_btagheavydown_M = float(envLo_btagheavy_M.GetBinContent(bin+1));
          bin_btagheavydown_E = float(envLo_btagheavy_E.GetBinContent(bin+1));

          sysup_btaglight_M = abs(bin_nom_M-bin_btaglightup_M);
          sysup_btaglight_E = abs(bin_nom_E-bin_btaglightup_E);
          sysdown_btaglight_M = abs(bin_nom_M-bin_btaglightdown_M);
          sysdown_btaglight_E = abs(bin_nom_E-bin_btaglightdown_E);

          sysup_btagheavy_M = abs(bin_nom_M-bin_btagheavyup_M);
          sysup_btagheavy_E = abs(bin_nom_E-bin_btagheavyup_E);
          sysdown_btagheavy_M = abs(bin_nom_M-bin_btagheavydown_M);
          sysdown_btagheavy_E = abs(bin_nom_E-bin_btagheavydown_E);

          sysup_total_M = sqrt(sysup_btaglight_M**2 + sysup_btagheavy_M**2);
          sysup_total_E = sqrt(sysup_btaglight_E**2 + sysup_btagheavy_E**2);
          sysdown_total_M = sqrt(sysdown_btaglight_M**2 + sysdown_btagheavy_M**2);
          sysdown_total_E = sqrt(sysdown_btaglight_E**2 + sysdown_btagheavy_E**2);

          graph_err_M.SetPointEXhigh(bin,last_M.GetBinWidth(bin+1)/2);
          graph_err_M.SetPointEXlow(bin,last_M.GetBinWidth(bin+1)/2);
          graph_err_M.SetPoint(bin,last_M.GetBinCenter(bin+1),last_M.GetBinContent(bin+1));
          graph_err_M.SetPointEYhigh(bin,sysup_total_M);
          graph_err_M.SetPointEYlow(bin,sysdown_total_M);

          graph_err_E.SetPointEXhigh(bin,last_E.GetBinWidth(bin+1)/2);
          graph_err_E.SetPointEXlow(bin,last_E.GetBinWidth(bin+1)/2);
          graph_err_E.SetPoint(bin,last_E.GetBinCenter(bin+1),last_E.GetBinContent(bin+1));
          graph_err_E.SetPointEYhigh(bin,sysup_total_E);
          graph_err_E.SetPointEYlow(bin,sysdown_total_E);

          ratio_graph_err_M.SetPointEXhigh(bin,last_M.GetBinWidth(bin+1)/2);
          ratio_graph_err_M.SetPointEXlow(bin,last_M.GetBinWidth(bin+1)/2);
          ratio_graph_err_M.SetPoint(bin,last_M.GetBinCenter(bin+1),1.);
          if (last_M.GetBinContent(bin+1)>0.01):
             ratio_graph_err_M.SetPointEYhigh(bin,abs(sysup_total_M/bin_nom_M));
             ratio_graph_err_M.SetPointEYlow(bin,abs(sysdown_total_M/bin_nom_M));
          else:
             ratio_graph_err_M.SetPointEYhigh(bin,0.);
             ratio_graph_err_M.SetPointEYlow(bin,0.);

          ratio_graph_err_E.SetPointEXhigh(bin,last_E.GetBinWidth(bin+1)/2);
          ratio_graph_err_E.SetPointEXlow(bin,last_E.GetBinWidth(bin+1)/2);
          ratio_graph_err_E.SetPoint(bin,last_E.GetBinCenter(bin+1),1.);
          if (last_E.GetBinContent(bin+1)>0.01):
             ratio_graph_err_E.SetPointEYhigh(bin,abs(sysup_total_E/bin_nom_E));
             ratio_graph_err_E.SetPointEYlow(bin,abs(sysdown_total_E/bin_nom_E));
          else:
             ratio_graph_err_E.SetPointEYhigh(bin,0.);
             ratio_graph_err_E.SetPointEYlow(bin,0.);

        graph_err_M.SetFillColorAlpha(kGray+2,0.7);
        graph_err_M.SetMarkerColor(kGray);
        graph_err_M.SetMarkerSize(0);
        graph_err_M.SetFillStyle(3001);
        graph_err_E.SetFillColorAlpha(kGray+2,0.7);
        graph_err_E.SetMarkerColor(kGray);
        graph_err_E.SetMarkerSize(0);
        graph_err_E.SetFillStyle(3001);
        ratio_graph_err_M.SetFillColorAlpha(kGray+2,0.7);
        ratio_graph_err_M.SetMarkerColor(kGray);
        ratio_graph_err_M.SetMarkerSize(0);
        ratio_graph_err_M.SetFillStyle(3001);
        ratio_graph_err_E.SetFillColorAlpha(kGray+2,0.7);
        ratio_graph_err_E.SetMarkerColor(kGray);
        ratio_graph_err_E.SetMarkerSize(0);
        ratio_graph_err_E.SetFillStyle(3001);

    if args.sumEM:
      #####################################
      ########## sumEM channel ############
      #####################################
      ymax_T = ymax_E+ymax_M
      y_T = stack_T.GetMaximum()
      if y_T>ymax_T: ymax_T=y_T
      #stack_M.SetMinimum(1000.)
      #stack_M.SetMaximum(3*ymax_T)
      if args.linear: stack_T.SetMaximum(1.3*ymax_T)
      if args.linear: stack_T.SetMinimum(1.)

      if not args.nodata:
        histD_T = histD_M
        histD_T.Add(histD_E)

      stack_T.Draw("HIST")
      if not args.nodata:
        histD_T.SetMarkerStyle(20)
        histD_T.SetMarkerSize(0.4)
        histD_T.SetLineWidth(1)
        histD_T.SetLineColor(ROOT.kBlack)
        histD_T.Draw("E SAME")
      if not args.nosyst: graph_err_T.Draw("SAME 2")

      if args.ratio and not args.nodata:
        lower_pad.cd()
        ratio = histD_T.Clone("ratio")
        ratio.SetLineColor(kBlack)
        ratio.SetMarkerStyle(21)
        ratio.SetTitle("")
        ratio.SetMinimum(c_rat2)
        ratio.SetMaximum(c_rat)
        ratio.GetYaxis().SetTitle("Data/MC")
        ratio.GetXaxis().SetTitle(name)
        ratio.GetXaxis().SetLabelSize(0.08)
        ratio.GetXaxis().SetTitleSize(0.12)
        ratio.GetXaxis().SetTitleOffset(1.0)
        ratio.GetYaxis().SetLabelSize(0.05)
        ratio.GetYaxis().SetTitleSize(0.09)
        ratio.GetYaxis().CenterTitle()
        ratio.GetYaxis().SetTitleOffset(0.5)
        # Set up plot for markers and errors
        ratio.Sumw2()
        ratio.SetStats(0)
        hTotal = histT_nom_T["vv"].Clone('hTotal')
        for s in samples[1:]:
          hTotal.Add(histT_nom_T[s])
        ratio.Divide(hTotal)
        ratio.Draw("ep")
        if not args.nosyst: ratio_graph_err_T.Draw("same 2")

      if name == "jet_max_cvltag" and (not args.nodata):
        print("charm tag")
        print("Integral of data is "+str(histD_T.Integral()))
        print("Integral of MC is "+str(hTotal.Integral()))
        for s in samples:
            print("Integral of "+str(s)+" MC is "+str(histT_nom_T[s].Integral()))
        print("Integral of data from 0.3 is "+str(histD_T.Integral(16,51)))
        print("Integral of MC from 0.3 is "+str(hTotal.Integral(16,51)))
       	for s in samples:
            print("Integral of "+str(s)+" MC from 0.3 is "+str(histT_nom_T[s].Integral(16,51)))

      if name == "InvM_2jets" and (not args.nodata):
        print("Integral of data is "+str(histD_T.Integral()))
        print("Integral of MC is "+str(hTotal.Integral()))
        print("Ratio is "+str(histD_T.Integral()/hTotal.Integral()))
        error_tot = {}
        error_fulltot = 0
        for s in samples:
            print("Integral of "+str(s)+" MC process is "+str(histT_nom_T[s].Integral()))
            error_tot[s] = 0
            for bin in range(histT_nom_T[s].GetNbinsX()):
                error_tot[s] = error_tot[s] + histT_nom_T[s].GetBinError(bin+1)**2
            error_tot[s] = sqrt(error_tot[s])
            print("error is "+str(error_tot[s]))
            error_fulltot = error_fulltot + error_tot[s]**2
        print("Final MC error is "+str(sqrt(error_fulltot)))

      ## Legends
      if args.ratio and not args.nodata: upper_pad.cd()
      leg = TLegend(0.7,0.38,0.95,0.95)
      leg.SetBorderSize(1)
      #leg.AddEntry(histT_nom_T["vv"],"VV","f")
      if args.wcs:
        leg.AddEntry(histT_nom_T["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_nocharm"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_T["st_charm"],"Single top tW cq","f")
        leg.AddEntry(histT_nom_T["st_nocharm"],"Single top tW uq","f")
        leg.AddEntry(histT_nom_T["st_else"],"Single top s/t channel","f")
      else:
        leg.AddEntry(histT_nom_T["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_light"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_T["st"],"Single top","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_bottom"],"t#bar{t} bq","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_else"],"t#bar{t} qg,gg","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_charmgluon"],"t#bar{t} cg","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_bottomgluon"],"t#bar{t} bg","f")
      leg.AddEntry(histT_nom_T["zjets"],"Z+jets","f")
      leg.AddEntry(histT_nom_T["wjets"],"W+jets","f")
      leg.AddEntry(histT_nom_T["ttbar_dl"],"Dileptonic t#bar{t}","f")
      #leg.AddEntry(histT_nom_T["ttbar_dh"],"Hadronic t#bar{t}","f")
      if args.stack and not args.nodata: leg.AddEntry(histD_T, "Data" ,"lep")
      #leg.Draw()
      termp= "totalHT_wqq"
      if args.ratio:
        notation = "_ratio_"
        if args.linear:
          notation = "_linratio_"
      else:
        notation = "_normed_"

      if args.year == "all": term_d = ""
      elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): term_d = str(args.year)

      if args.png: c1.Print(plotdir+termp+notation+ term_d+name + ".png")
      else: c1.Print(plotdir+termp+notation + term_d+name + ".root")

      if c1: 
         c1.Close(); gSystem.ProcessEvents();

    else:
      #####################################
      ############ M channel ##############
      #####################################
      y_M = stack_M.GetMaximum()
      if y_M>ymax_M: ymax_M=y_M
      #stack_M.SetMinimum(1000.)
      #stack_M.SetMaximum(3*ymax_M)
      if args.linear: stack_M.SetMaximum(1.3*ymax_M)
      if args.linear: stack_M.SetMinimum(1.)
      y_E = stack_E.GetMaximum()
      if y_E>ymax_E: ymax_E=y_E
      #stack_E.SetMinimum(1000.)
      #stack_E.SetMaximum(3*ymax_E)
      if args.linear: stack_E.SetMaximum(1.3*ymax_E)
      if args.linear: stack_M.SetMinimum(1.)

      stack_M.Draw("HIST")
      if not args.nodata:
        histD_M.SetMarkerStyle(20)
        histD_M.SetMarkerSize(0.3)
        histD_M.SetLineWidth(1)
        histD_M.SetLineColor(ROOT.kBlack)
        histD_M.Draw("E SAME")
      if not args.nosyst: graph_err_M.Draw("SAME 2")

      if args.ratio and not args.nodata:
        lower_pad.cd()
        ratio = histD_M.Clone("ratio")
        ratio.SetLineColor(kBlack)
        ratio.SetMarkerStyle(21)
        ratio.SetTitle("")
        ratio.SetMinimum(c_rat2)
        ratio.SetMaximum(c_rat)
        ratio.GetYaxis().SetTitle("Data/MC")
        ratio.GetXaxis().SetTitle(name)
        ratio.GetXaxis().SetLabelSize(0.08)
        ratio.GetXaxis().SetTitleSize(0.12)
        ratio.GetXaxis().SetTitleOffset(1.0)
        ratio.GetYaxis().SetLabelSize(0.05)
        ratio.GetYaxis().SetTitleSize(0.09)
        ratio.GetYaxis().CenterTitle()
        ratio.GetYaxis().SetTitleOffset(0.5)
        # Set up plot for markers and errors
        ratio.Sumw2()
        ratio.SetStats(0)
        hTotal = histT_nom_M["vv"].Clone('hTotal')
        for s in samples[1:]:
          hTotal.Add(histT_nom_M[s])
        ratio.Divide(hTotal)
        ratio.Draw("ep")
        if not args.nosyst: ratio_graph_err_M.Draw("same 2")

      if name == "InvM_2jets" and (not args.nodata):
        print("Integral of M data is "+str(histD_M.Integral()))
        print("Integral of M MC is "+str(hTotal.Integral()))
        print("Ratio is "+str(histD_M.Integral()/hTotal.Integral()))
        error_tot = {}
        error_fulltot = 0
        for s in samples:
            print("Integral of "+str(s)+" MC process is "+str(histT_nom_M[s].Integral()))
            error_tot[s] = 0
            for bin in range(histT_nom_M[s].GetNbinsX()):
                error_tot[s] = error_tot[s] + histT_nom_M[s].GetBinError(bin+1)**2
            error_tot[s] = sqrt(error_tot[s])
            print("error is "+str(error_tot[s]))
            error_fulltot = error_fulltot + error_tot[s]**2
        print("Final MC error is "+str(sqrt(error_fulltot)))

      ## Legends
      if args.ratio and not args.nodata: upper_pad.cd()
      leg = TLegend(0.78,0.4,0.95,0.95)
      leg.SetBorderSize(1)
      #leg.AddEntry(histT_nom_M["vv"],"VV","f")
      if args.wcs:
        leg.AddEntry(histT_nom_M["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_nocharm"],"t#bar{t} qq'","f")
        leg.AddEntry(histT_nom_M["st_charm"],"Single top charm","f")
        leg.AddEntry(histT_nom_M["st_nocharm"],"Single top light","f")
        leg.AddEntry(histT_nom_M["st_else"],"Single top no had W","f")
      else:
        leg.AddEntry(histT_nom_M["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_light"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_M["st"],"Single top","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_bottom"],"t#bar{t} bq","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_else"],"t#bar{t} qg,gg","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_charmgluon"],"t#bar{t} cg","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_bottomgluon"],"t#bar{t} bg","f")
      leg.AddEntry(histT_nom_M["zjets"],"Z+jets","f")
      leg.AddEntry(histT_nom_M["wjets"],"W+jets","f")
      leg.AddEntry(histT_nom_M["ttbar_dl"],"Dileptonic t#bar{t}","f")
      #leg.AddEntry(histT_nom_M["ttbar_dh"],"Hadronic t#bar{t}","f")
      if args.stack and not args.nodata: leg.AddEntry(histD_M, "Data" ,"lep")
      #leg.Draw()
      termp= "totalHT_wqq_M"
      if args.ratio: 
        notation = "_ratio_"
        if args.linear:
          notation = "_linratio_"
      else: 
        notation = "_normed_"

      if args.year == "all": term_d = ""
      elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): term_d = str(args.year)

      if args.png: c1.Print(plotdir+termp+notation+ term_d+name + ".png")
      else: c1.Print(plotdir+termp+notation + term_d+name + ".root")

      #####################################
      ############ E channel ##############
      #####################################
      stack_E.Draw("HIST")
      if not args.nodata:
        histD_E.SetMarkerStyle(20)
        histD_E.SetMarkerSize(0.3)
        histD_E.SetLineWidth(1)
        histD_E.SetLineColor(ROOT.kBlack)
        histD_E.Draw("E SAME")
      if not args.nosyst: graph_err_E.Draw("SAME 2")

      if args.ratio and not args.nodata:
        lower_pad.cd()
        ratio = histD_E.Clone("ratio")
        ratio.SetLineColor(kBlack)
        ratio.SetMarkerStyle(21)
        ratio.SetTitle("")
        ratio.SetMinimum(c_rat2)
        ratio.SetMaximum(c_rat)
        ratio.GetYaxis().SetTitle("Data/MC")
        ratio.GetXaxis().SetTitle(name)
        ratio.GetXaxis().SetLabelSize(0.08)
        ratio.GetXaxis().SetTitleSize(0.12)
        ratio.GetXaxis().SetTitleOffset(1.0)
        ratio.GetYaxis().SetLabelSize(0.05)
        ratio.GetYaxis().SetTitleSize(0.09)
        ratio.GetYaxis().CenterTitle()
        ratio.GetYaxis().SetTitleOffset(0.5)
        # Set up plot for markers and errors
        ratio.Sumw2()
        ratio.SetStats(0)
        hTotal = histT_nom_E["vv"].Clone('hTotal')
        for s in samples[1:]:
          hTotal.Add(histT_nom_E[s])
        ratio.Divide(hTotal)
        ratio.Draw("ep")
        if not args.nosyst: ratio_graph_err_E.Draw("same 2")

      if name == "InvM_2jets" and (not args.nodata):
        print("Integral of E data is "+str(histD_E.Integral()))
        print("Integral of E MC is "+str(hTotal.Integral()))
        print("Ratio is "+str(histD_E.Integral()/hTotal.Integral()))
        for s in samples:
            print("Integral of "+str(s)+" MC process is "+str(histT_nom_E[s].Integral()))
            error_tot[s] = 0
            for bin in range(histT_nom_E[s].GetNbinsX()):
                error_tot[s] = error_tot[s] + histT_nom_E[s].GetBinError(bin+1)**2
            error_tot[s] = sqrt(error_tot[s])
            print("error is "+str(error_tot[s]))
            error_fulltot = error_fulltot + error_tot[s]**2
        print("Final MC error is "+str(sqrt(error_fulltot)))

      ## Legends
      if args.ratio and not args.nodata: upper_pad.cd()
      leg = TLegend(0.7,0.6,0.89,0.89)
      leg.SetBorderSize(1)
      #leg.AddEntry(histT_nom_E["vv"],"VV","f")
      if args.wcs:
        leg.AddEntry(histT_nom_E["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_nocharm"],"t#bar{t} qq'","f")
        leg.AddEntry(histT_nom_E["st_charm"],"Single top charm","f")
        leg.AddEntry(histT_nom_E["st_nocharm"],"Single top light","f")
        leg.AddEntry(histT_nom_E["st_else"],"Single top no had W","f")
      else:
        leg.AddEntry(histT_nom_E["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_light"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_E["st"],"Single top","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_bottom"],"t#bar{t} bq","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_else"],"t#bar{t} qg,gg","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_charmgluon"],"t#bar{t} cg","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_bottomgluon"],"t#bar{t} bg","f")
      leg.AddEntry(histT_nom_E["ttbar_dl"],"Dileptonic top antitop","f")
      leg.AddEntry(histT_nom_E["zjets"],"Z plus jets","f")
      leg.AddEntry(histT_nom_E["wjets"],"W plus jets","f")
      if args.stack and not args.nodata: leg.AddEntry(histD_E, "Data" ,"lep")
      #leg.Draw()
      termp= "totalHT_wqq_E"
      if args.ratio:
        notation = "_ratio_"
        if args.linear:
          notation = "_linratio_"
      else:
        notation = "_normed_"

      if args.year == "all": term_d = ""
      elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): term_d = str(args.year)

      if args.png: c1.Print(plotdir+termp+notation+ term_d+name + ".png")
      else: c1.Print(plotdir+termp+notation + term_d+name + ".root")

      if c1: 
         c1.Close(); gSystem.ProcessEvents();

    for data_op in datayears:
                        histFile[name][data_op].Close()
                        if not args.nodata:
                          histFileDM[name][data_op].Close()
                          histFileDE[name][data_op].Close()


