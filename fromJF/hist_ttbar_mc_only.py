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

# Use like:
# python arg.py --data="No"
# python hist_plot.py --data="No" --stack --ratio

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

plotdir = '/nfs/cms/vazqueze/higgssearch/plotspng/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

c_rat = 1.2
c_rat2 = 0.8

norm_factorM = {}
norm_factorE = {}

### btagMM smeared, ptjets30
norm_factorM["all"] = 0.9520
norm_factorM["2016"] = 0.9374
norm_factorM["2016B"] = 0.9967
norm_factorM["2017"] = 0.9479
norm_factorM["2018"] = 0.9481
norm_factorE["all"] = 0.9382
norm_factorE["2016"] = 0.9340
norm_factorE["2016B"] = 0.9764 
norm_factorE["2017"] = 0.9153
norm_factorE["2018"] = 0.9459

### btagMM smeared
norm_factorM["all"] = 0.9614
norm_factorM["2016"] = 0.9447
norm_factorM["2016B"] = 1.0033
norm_factorM["2017"] = 0.9534
norm_factorM["2018"] = 0.9613
norm_factorE["all"] = 0.9474
norm_factorE["2016"] = 0.9390
norm_factorE["2016B"] = 0.9837
norm_factorE["2017"] = 0.9212
norm_factorE["2018"] = 0.9592

### lepton50 smeared
norm_factorM["all"] = 0.9593
norm_factorM["2016"] = 0.9464
norm_factorM["2016B"] = 1.0044
norm_factorM["2017"] = 0.9283
norm_factorM["2018"] = 0.9747
norm_factorE["all"] = 0.9678
norm_factorE["2016"] = 0.9614
norm_factorE["2016B"] = 1.0128
norm_factorE["2017"] = 0.9174
norm_factorE["2018"] = 0.9961

### lepton50 smeared, ptjets30
norm_factorM["all"] = 0.9522
norm_factorM["2016"] = 0.9357
norm_factorM["2016B"] = 0.9912
norm_factorM["2017"] = 0.9257
norm_factorM["2018"] = 0.9665
norm_factorE["all"] = 0.9561
norm_factorE["2016"] = 0.9460
norm_factorE["2016B"] = 0.9946
norm_factorE["2017"] = 0.9100
norm_factorE["2018"] = 0.9832

## Open hists files

filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/"

term = "hist_wqqfromJF_"

observable_names = ["InvM_2jets", "nJetGood", "jet_1_pt", "jet_1_nmu", "jet_1_eta", "jet_2_pt", "jet_2_eta", "jet_2_mass", "jet_2_qgl","jet_2_nmu","jet_1_qgl",
   "lepton_pt", "lepton_eta", "lepton_pt_detail", "lepton_eta_thick", "InvM_bot_closer", "InvM_bot_farther",
   "deltaR_jet1_jet2", "deltaphi_jet1_jet2", "deltaeta_jet1_jet2", "MET_pt", "MET_sig", "MET_my_sig",
   "transverse_mass", "tracks_jet1", "tracks_jet2", "deltaphi_MET_jets_1", "deltaphi_MET_jets_2", "pT_Wlep",
   "jet_1_btag", "jet_2_btag", "deltaphi_MET_lep", "jet_bot1_btag", "jet_bot2_btag", "jet_bot1_pt", "jet_bot1_eta", "jet_bot2_pt", "jet_bot2_eta",
   "jet_bot1_btag_thick", "jet_bot2_btag_thick", "jet_1_btag_thick", "jet_2_btag_thick",
   "jet_bot1_btagnumber", "jet_bot2_btagnumber", "jet_1_btagnumber", "jet_2_btagnumber"]

datayears = ["2016","2016B","2017","2018"]
#datayears = ["2018","2016","2016B"]

samplesHT = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_light","ttbar_sl_bottom","ttbar_dl","ttbar_dh","zz","wz",
        "st_1","st_2","st_3","st_4", "ttbar_sl_else", "ttbar_sl_charmgluon","ttbar_sl_bottomgluon"]

## Adding QCD

histFile = {}
histFileDM = {}
histFileDE = {}

for data_op in datayears:
        if data_op == "2016B": 
            end_term = "2.root"
        else:
            end_term = ".root"
        ## mc files
        histFile[data_op] = {}
        histFileDM[data_op] = {}
        histFileDE[data_op] = {}
        for name in observable_names:
                if isfile(filePath + term+"MC_"+data_op+"_"+name+end_term):
                        histFile[data_op][name] = TFile.Open(filePath + term+"MC_"+data_op+"_"+name+end_term,"READ")
                # data files
                histFileDM[data_op][name] = TFile.Open(filePath + term+"dataM_"+data_op+"_"+name+".root","READ")
                histFileDE[data_op][name] = TFile.Open(filePath + term+"dataE_"+data_op+"_"+name+".root","READ")
        #print(data_op)
        #print(histFile[data_op].keys())
        #print(histFileDM[data_op].keys())


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

for name in observable_names:

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
    if data_op == "2016B":
      data_term = "2016"
    else:
      data_term = data_op
    for s in samples:
      if s[0:8] == "ttbar_sl": 
         s_term = s[0:8]+data_term+s[8:] 
      else:
         s_term = s+data_term
      hist_nom_M[data_op][s] = histFile[data_op][name].Get(s_term+"_"+name+"_M")
      hist_nom_E[data_op][s] = histFile[data_op][name].Get(s_term+"_"+name+"_E")
      hist_btaglightup_M[data_op][s] = histFile[data_op][name].Get(s_term+"_"+name+"_M_btaglightup")
      hist_btaglightup_E[data_op][s] = histFile[data_op][name].Get(s_term+"_"+name+"_E_btaglightup")
      hist_btaglightdown_M[data_op][s] = histFile[data_op][name].Get(s_term+"_"+name+"_M_btaglightdown")
      hist_btaglightdown_E[data_op][s] = histFile[data_op][name].Get(s_term+"_"+name+"_E_btaglightdown")
      hist_btagheavyup_M[data_op][s] = histFile[data_op][name].Get(s_term+"_"+name+"_M_btagheavyup")
      hist_btagheavyup_E[data_op][s] = histFile[data_op][name].Get(s_term+"_"+name+"_E_btagheavyup")
      hist_btagheavydown_M[data_op][s] = histFile[data_op][name].Get(s_term+"_"+name+"_M_btagheavydown")
      hist_btagheavydown_E[data_op][s] = histFile[data_op][name].Get(s_term+"_"+name+"_E_btagheavydown")
    #print(data_op+"M_"+name+"_M")
    #print(histFileDM[data_op][name].ls())
    if not args.nodata: histdata_M[data_op] = histFileDM[data_op][name].Get("data"+data_op+"M_"+name+"_M")
    if not args.nodata: histdata_E[data_op] = histFileDE[data_op][name].Get("data"+data_op+"E_"+name+"_E")

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
      if args.norm: hist_nom_M[data_op][s].Scale(norm_factorM[str(args.year)])
      if args.norm: hist_nom_E[data_op][s].Scale(norm_factorE[str(args.year)])
      if args.norm: hist_btaglightup_M[data_op][s].Scale(norm_factorM[str(args.year)])
      if args.norm: hist_btaglightup_E[data_op][s].Scale(norm_factorE[str(args.year)])
      if args.norm: hist_btaglightdown_M[data_op][s].Scale(norm_factorM[str(args.year)])
      if args.norm: hist_btaglightdown_E[data_op][s].Scale(norm_factorE[str(args.year)])
      if args.norm: hist_btagheavyup_M[data_op][s].Scale(norm_factorM[str(args.year)])
      if args.norm: hist_btagheavyup_E[data_op][s].Scale(norm_factorE[str(args.year)])
      if args.norm: hist_btagheavydown_M[data_op][s].Scale(norm_factorM[str(args.year)])
      if args.norm: hist_btagheavydown_E[data_op][s].Scale(norm_factorE[str(args.year)])
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

  ### sumamos todos los años para todos los histogramas
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

  if not args.nodata:
    histD_M = histdata_M[datayears[0]]
    histD_E = histdata_E[datayears[0]]
    for d in datayears[1:]:
       histD_M.Add(histdata_M[d])
       histD_E.Add(histdata_E[d])

  ### Para los histogramas de sistemática up y down sumamos todas las contribuciones
  histT_sT_btaglightup_M = histT_btaglightup_M[samples[0]]
  histT_sT_btaglightup_E = histT_btaglightup_E[samples[0]]
  histT_sT_btaglightdown_M = histT_btaglightdown_M[samples[0]]
  histT_sT_btaglightdown_E = histT_btaglightdown_E[samples[0]]
  histT_sT_btagheavyup_M = histT_btagheavyup_M[samples[0]]
  histT_sT_btagheavyup_E = histT_btagheavyup_E[samples[0]]
  histT_sT_btagheavydown_M = histT_btagheavydown_M[samples[0]]
  histT_sT_btagheavydown_E = histT_btagheavydown_E[samples[0]]
  for s in samples[1:]:
    histT_sT_btaglightup_M.Add(histT_btaglightup_M[s])
    histT_sT_btaglightup_E.Add(histT_btaglightup_E[s])
    histT_sT_btaglightdown_M.Add(histT_btaglightdown_M[s])
    histT_sT_btaglightdown_E.Add(histT_btaglightdown_E[s])
    histT_sT_btagheavyup_M.Add(histT_btagheavyup_M[s])
    histT_sT_btagheavyup_E.Add(histT_btagheavyup_E[s])
    histT_sT_btagheavydown_M.Add(histT_btagheavydown_M[s])
    histT_sT_btagheavydown_E.Add(histT_btagheavydown_E[s])

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  colors = {}
  colors["ttbar_dl"] = (222,90,106)
  colors["wjets"] = (155,152,204)
  colors["ttbar_sl_bottomgluon"] = (236,164,207)
  colors["ttbar_sl_charm"] = (204,255,153)
  colors["ttbar_sl_light"] = (120,154,86)
  colors["ttbar_sl_bottom"] = (201,79,152)
  colors["ttbar_sl_else"] = (242,193,121)
  colors["ttbar_sl_charmgluon"] = (222,212,74)
  colors["vv"] = (255,180,85)
  colors["ttbar_dh"] = (204,204,0)
  colors["zjets"] = (113,209,223)
  colors["st"] = (153,51,255)

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

    ## Stack creation
    samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon","ttbar_sl_else","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]
    samples = ["ttbar_sl_light", "ttbar_sl_charm"]

    if args.ratio and not args.nodata: upper_pad.cd()
    stack_M = ROOT.THStack()
    stack_E = ROOT.THStack()
    for s in samples:
      stack_M.Add(histT_nom_M[s])
      stack_E.Add(histT_nom_E[s])

    last_M = stack_M.GetStack().Last()
    last_E = stack_E.GetStack().Last()

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
    #graph_err_M.Draw("SAME 2")

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
      ratio_graph_err_M.Draw("same 2")

    #if name == "InvM_2jets":
      #print("Integral of M data is "+str(histD_M.Integral()))
      #print("Integral of M MC is "+str(hTotal.Integral()))
      #print("Ratio is "+str(histD_M.Integral()/hTotal.Integral()))

    ## Legends
    if args.ratio and not args.nodata: upper_pad.cd()
    leg = TLegend(0.7,0.6,0.89,0.89)
    leg.SetBorderSize(1)
    leg.AddEntry(histT_nom_M["vv"],"VV","f")
    leg.AddEntry(histT_nom_M["ttbar_dl"],"Dileptonic top antitop","f")
    leg.AddEntry(histT_nom_M["zjets"],"Z plus jets","f")
    leg.AddEntry(histT_nom_M["wjets"],"W plus jets","f")
    leg.AddEntry(histT_nom_M["st"],"Single top","f")
    leg.AddEntry(histT_nom_M["ttbar_dh"],"Hadronic top antitop","f")
    leg.AddEntry(histT_nom_M["ttbar_sl_bottom"],"Top antitop, bottom","f")
    leg.AddEntry(histT_nom_M["ttbar_sl_light"],"Top antitop, light","f")
    leg.AddEntry(histT_nom_M["ttbar_sl_charm"],"Top antitop, charm","f")
    leg.AddEntry(histT_nom_M["ttbar_sl_bottomgluon"],"Top antitop, bottom gluon","f")
    leg.AddEntry(histT_nom_M["ttbar_sl_charmgluon"],"Top antitop, charm gluon","f")
    leg.AddEntry(histT_nom_M["ttbar_sl_else"],"Top antitop, gluon gluon","f")
    if args.stack and not args.nodata: leg.AddEntry(histD_M, "Data" ,"lep")
    leg.Draw()
    term= "totalHT_wqq_M"
    if args.ratio: 
      notation = "_ratio_"
      if args.linear:
        notation = "_linratio_"
    else: 
      notation = "_normed_"

    if args.year == "all": term_d = ""
    elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): term_d = str(args.year)

    if args.png: c1.Print(plotdir+term+notation+ term_d+name + ".png")
    else: c1.Print(plotdir+term+notation + term_d+name + ".root")

    stack_E.Draw("HIST")
    #graph_err_E.Draw("SAME 2")

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
      ratio_graph_err_E.Draw("same 2")

    #if name == "InvM_2jets":
      #print("Integral of E data is "+str(histD_E.Integral()))
      #print("Integral of E MC is "+str(hTotal.Integral()))
      #print("Ratio is "+str(histD_E.Integral()/hTotal.Integral()))

    ## Legends
    if args.ratio and not args.nodata: upper_pad.cd()
    leg = TLegend(0.7,0.6,0.89,0.89)
    leg.SetBorderSize(1)
    leg.AddEntry(histT_nom_E["vv"],"VV","f")
    leg.AddEntry(histT_nom_E["ttbar_dl"],"Dileptonic top antitop","f")
    leg.AddEntry(histT_nom_E["zjets"],"Z plus jets","f")
    leg.AddEntry(histT_nom_E["wjets"],"W plus jets","f")
    leg.AddEntry(histT_nom_E["st"],"Single top","f")
    leg.AddEntry(histT_nom_E["ttbar_dh"],"Hadronic top antitop","f")
    leg.AddEntry(histT_nom_E["ttbar_sl_bottom"],"Top antitop, bottom","f")
    leg.AddEntry(histT_nom_E["ttbar_sl_light"],"Top antitop, light","f")
    leg.AddEntry(histT_nom_E["ttbar_sl_charm"],"Top antitop, charm","f")
    leg.AddEntry(histT_nom_E["ttbar_sl_bottomgluon"],"Top antitop, bottom gluon","f")
    leg.AddEntry(histT_nom_E["ttbar_sl_charmgluon"],"Top antitop, charm gluon","f")
    leg.AddEntry(histT_nom_E["ttbar_sl_else"],"Top antitop, gluon gluon","f")
    if args.stack and not args.nodata: leg.AddEntry(histD_E, "Data" ,"lep")
    leg.Draw()
    term= "totalHT_wqq_E"
    if args.ratio:
      notation = "_ratio_"
      if args.linear:
        notation = "_linratio_"
    else:
      notation = "_normed_"

    if args.year == "all": term_d = ""
    elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): term_d = str(args.year)

    if args.png: c1.Print(plotdir+term+notation+ term_d+name + ".png")
    else: c1.Print(plotdir+term+notation + term_d+name + ".root")



for s in observable_names:
        for data_op in datayears:
                        histFile[data_op][s].Close()
                        histFileDM[data_op][s].Close()
                        histFileDE[data_op][s].Close()


