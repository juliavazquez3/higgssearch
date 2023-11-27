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
parser.add_argument("--channel", type=string, default="bot1_chi",
                    help="Select year of process to run")
parser.add_argument("--etabin", type=string, default="none",
                    help="eta muon requirement")
parser.add_argument("--nosyst", action="store_true", default=False,
                    help="systematics inclusion")

# Use like:
# python higgssearch/fromJF/hist_fromJFwqq_syst_total.py --stack --ratio --png --norm --year="all" --channel="btagMM_chitest_sl"

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

plotdir = '/nfs/cms/vazqueze/higgssearch/plotspng/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

c_rat = 1.5
c_rat2 = 0.5

if args.etabin == "none": eta_bin = "none"
elif args.etabin == "one": eta_bin = "one"
elif args.etabin == "two": eta_bin = "two"
elif args.etabin == "three": eta_bin = "three"
elif args.etabin == "four": eta_bin = "four"
elif args.etabin == "five": eta_bin = "five"
else: raise NameError('Incorrect channel')
print(eta_bin)

if args.channel == "one": channel = "jet_bot1_nochi"
elif args.channel == "two": channel = "jet_bot2_nochi"
elif args.channel == "one_chi": channel = "jet_bot1_chi"
elif args.channel == "two_chi": channel = "jet_bot2_chi"
elif args.channel == "both_chi": channel = "jet_both_chi"
else: raise NameError('Incorrect channel')
print(channel)

norm_factor = {}

norm_factor["all"] = {}; norm_factor["2016"] = {}; norm_factor["2016B"] = {}; norm_factor["2017"] = {}; norm_factor["2018"] = {};

### btagMM smeared, chitest or no
if (channel == "jet_bot1_chi" or channel == "jet_bot2_chi" or channel == "jet_both_chi"): norm_factor["all"] = 0.91; ## btagMM_chitest
elif (channel == "jet_bot1_nochi" or channel == "jet_bot2_nochi" or channel == "jet_both_nochi"): norm_factor["all"] = 0.95; ## btagMM

observable_names = ["InvM_2jets","jet_1_pt", "jet_1_nmu", "jet_1_eta", "jet_2_pt", "jet_2_eta", "jet_2_mass", "jet_2_qgl","jet_2_nmu","jet_1_qgl",
   "lepton_pt", "lepton_eta", "lepton_pt_detail", "lepton_eta_thick", "InvM_bot_closer", "InvM_bot_farther",
   "deltaR_jet1_jet2", "deltaphi_jet1_jet2", "deltaeta_jet1_jet2", "MET_pt_aux", "MET_sig", "MET_my_sig",
   "transverse_mass", "tracks_jet1", "tracks_jet2", "deltaphi_MET_jets_1", "deltaphi_MET_jets_2", "pT_Wlep",
   "jet_1_btag", "jet_2_btag", "deltaphi_MET_lep", "jet_bot1_btag", "jet_bot2_btag", "jet_bot1_pt", "jet_bot1_eta", "jet_bot2_pt", "jet_bot2_eta",
   "jet_bot1_btag_thick", "jet_bot2_btag_thick", "jet_1_btag_thick", "jet_2_btag_thick",
   "jet_bot1_btagnumber", "jet_bot2_btagnumber", "jet_1_btagnumber", "jet_2_btagnumber",
   "jet_1_cvltag_csv", "jet_2_cvltag_csv", "jet_1_cvltag", "jet_2_cvltag","InvM30","InvM31","InvMl0","InvMl1","chi2_test0","chi2_test1",
   "InvM3_good","InvM3_bad","InvMl_good","InvMl_bad","chi2_test_good","chi2_test_bad","jet_max_cvltag","jet_min_cvltag",
   "jet_1_cvbtag_csv", "jet_2_cvbtag_csv", "jet_1_cvbtag", "jet_2_cvbtag", "jet_max_cvbtag", "jet_min_cvbtag","tau_discr",
   "jet_1_eta_thick","jet_2_eta_thick","jet_bot1_eta_thick","jet_bot2_eta_thick",
   "InvM_2jets_thick","InvM_2jets_short","bot1_muons","bot2_muons","muon_both_eta","muon_both_z","muon_both_pt","muon_both_iso","muon_both_iso_log",
   "muon_both_pt_rel","muon_both_z_short","jet_bot1_tracks", "jet_bot2_tracks","nJetGood","deltaR_muon_jet","muon_both_z2","muon_both_z3","muon_both_z2_v2","muon_both_iso_abs"]

not_rebin = ["nJetGood","InvM3_good","InvM3_bad","InvMl_good","InvMl_bad","lepton_eta_thick","jet_bot1_btagnumber", "jet_bot2_btagnumber", 
      "jet_1_btagnumber", "jet_2_btagnumber","muon_jet_pt","muon_jet_relpt","muon_jet_eta","tau_discr","transverse_mass",
      "jet_1_flavourP", "jet_2_flavourP", "jet_bot1_flavourP", "jet_bot2_flavourP"]

#observable_names = ["jet_1_flavourP", "jet_2_flavourP", "jet_bot1_flavourP", "jet_bot2_flavourP", "btag_sf", "lep_id_sf", "lep_trig_sf", "lep_iso_sf",
#     "puWeight", "PUjetID_SF", "top_weight","Frag_weight_sl","Br_weight_sl","muon_bot1_mother","muon_bot2_mother","muon_bot1_mother_mine","muon_bot2_mother_mine",
#     "muon_bot1_mother_detail","muon_bot2_mother_detail","muon_from_bot_sf_z"]

#observable_names = ["muon_bot1_iso_abs","muon_bot2_iso_abs"]

datayears = ["2016","2016B","2017","2018"]
#datayears = ["2018","2016","2016B"]

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
        lumi[data_op]["ttbar_sl_nomuon"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_muon_charm"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_muon_bottom"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_muon_else"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_dl_nomuon"] = lumi[data_op]["ttbar_dl"]
        lumi[data_op]["ttbar_dl_muon_charm"] = lumi[data_op]["ttbar_dl"]
        lumi[data_op]["ttbar_dl_muon_bottom"] = lumi[data_op]["ttbar_dl"]
        lumi[data_op]["ttbar_dl_muon_else"] = lumi[data_op]["ttbar_dl"]

listsampl = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "ttbar_sl","ttbar_dl","ttbar_dh","zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6",
        "zjets_7","zjets_8","st_1","st_2","st_3","st_4","zz","wz", "ttbar_sl_nomuon", "ttbar_sl_muon_charm", 
        "ttbar_sl_muon_bottom", "ttbar_sl_muon_else","ttbar_dl_nomuon", "ttbar_dl_muon_charm",
        "ttbar_dl_muon_bottom", "ttbar_dl_muon_else"]

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

if channel == "jet_bot1_chi": term1 = 'botjets_muons_corr/muon_bot1/'
elif channel == "jet_bot2_chi": term1 = "botjets_muons_corr/muon_bot2/"
elif channel == "jet_both_chi": term1 = "botjets_muons_corr/muon_both/"
elif channel == "jet_bot1_nochi": term1 = "botjets_muons_corr/muon_bot1/nochitest/"
elif channel == "jet_bot2_nochi": term1 = "botjets_muons_corr/muon_bot2/nochitest/"
else: raise NameError('Incorrect channel option')

if eta_bin == "one":
   if channel == "jet_bot1_chi":
      term1 = "botjets_muons/eta_bins/muon_bot1/chitest/eta1/"
   if channel == "jet_bot2_chi":
      term1 = "botjets_muons/eta_bins/muon_bot2/chitest/eta1/"
   if channel == "jet_bot1_nochi":
      term1 = "botjets_muons/eta_bins/muon_bot1/nochitest/eta1/"
   if channel == "jet_bot2_nochi":
      term1 = "botjets_muons/eta_bins/muon_bot2/nochitest/eta1/"
if eta_bin == "two":
   if channel == "jet_bot1_chi":
      term1 = "botjets_muons/eta_bins/muon_bot1/chitest/eta2/"
   if channel == "jet_bot2_chi":
      term1 = "botjets_muons/eta_bins/muon_bot2/chitest/eta2/"
   if channel == "jet_bot1_nochi":
      term1 = "botjets_muons/eta_bins/muon_bot1/nochitest/eta2/"
   if channel == "jet_bot2_nochi":
      term1 = "botjets_muons/eta_bins/muon_bot2/nochitest/eta2/"
if eta_bin == "three":
   if channel == "jet_bot1_chi":
      term1 = "botjets_muons/eta_bins/muon_bot1/chitest/eta3/"
   if channel == "jet_bot2_chi":
      term1 = "botjets_muons/eta_bins/muon_bot2/chitest/eta3/"
   if channel == "jet_bot1_nochi":
      term1 = "botjets_muons/eta_bins/muon_bot1/nochitest/eta3/"
   if channel == "jet_bot2_nochi":
      term1 = "botjets_muons/eta_bins/muon_bot2/nochitest/eta3/"
if eta_bin == "four":
   if channel == "jet_bot1_chi":
      term1 = "botjets_muons/eta_bins/muon_bot1/chitest/eta4/"
   if channel == "jet_bot2_chi":
      term1 = "botjets_muons/eta_bins/muon_bot2/chitest/eta4/"
   if channel == "jet_bot1_nochi":
      term1 = "botjets_muons/eta_bins/muon_bot1/nochitest/eta4/"
   if channel == "jet_bot2_nochi":
      term1 = "botjets_muons/eta_bins/muon_bot2/nochitest/eta4/"
if eta_bin == "five":
   if channel == "jet_bot1_chi":
      term1 = "botjets_muons/eta_bins/muon_bot1/chitest/eta5/"
   if channel == "jet_bot2_chi":
      term1 = "botjets_muons/eta_bins/muon_bot2/chitest/eta5/"
   if channel == "jet_bot1_nochi":
      term1 = "botjets_muons/eta_bins/muon_bot1/nochitest/eta5/"
   if channel == "jet_bot2_nochi":
      term1 = "botjets_muons/eta_bins/muon_bot2/nochitest/eta5/"

histFile = {}
histFileDM = {}
histFileDE = {}

myratiofile = TFile('/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/etabins/pt5_test/ratios_muons_'+str(args.year)+'_'+str(args.channel)+'_'+str(eta_bin)+'.root', 'RECREATE' )

for name in observable_names:
  ## Open hists files
  filePath = '/nfs/cms/vazqueze/new_hists/fromJF/wqq/'+term1
  term = 'hist_wqqfromJF_'
  end_term = '.root'
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
       #print(histFile[name].keys())
       #print(histFileDM[name].keys())

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
        "ttbar_sl_muon_charm","ttbar_sl_muon_else","ttbar_sl_muon_bottom","ttbar_dh","zz","wz",
        "ttbar_dl_muon_charm","ttbar_dl_muon_else","ttbar_dl_muon_bottom","ttbar_dl_nomuon",
        "st_1","st_2","st_3","st_4","ttbar_sl_nomuon"]

  ## HISTS
  samples_foryear = {}
  hist_nom = {}
  hist_btaglightup = {}
  hist_btaglightdown = {}
  hist_btagheavyup = {}
  hist_btagheavydown = {}
  histdata_M = {}
  histdata_E = {}
  histdata = {}
  for data_op in datayears:
    hist_nom[data_op] = {}
    hist_btaglightup[data_op] = {}
    hist_btaglightdown[data_op] = {}
    hist_btagheavyup[data_op] = {}
    hist_btagheavydown[data_op] = {}
    data_term = data_op
    #print(name)
    for s in samples:
      if (s[0:8] == "ttbar_sl" or s[0:8] == "ttbar_dl"): 
         s_term = s[0:8]+data_term+s[8:] 
      else:
         s_term = s+data_term
      hist_nom[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name)
      if not args.nosyst:
         hist_btaglightup[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_btaglightup")
         hist_btaglightdown[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_btaglightdown")
         hist_btagheavyup[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_btagheavyup")
         hist_btagheavydown[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_btagheavydown")
    #print(histFileDM[name][data_op].ls())
    if not args.nodata: histdata[data_op] = histFileDM[name][data_op].Get("data"+data_op+"M_"+name)
    if not args.nodata: histdata_E[data_op] = histFileDE[name][data_op].Get("data"+data_op+"E_"+name)
    if not args.nodata: histdata[data_op].Add(histdata_E[data_op])

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
      hist_nom[data_op][s].Scale(lumi_data/lumi[data_op][s])
      if not args.nosyst:
         hist_btaglightup[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_btaglightdown[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_btagheavyup[data_op][s].Scale(lumi_data/lumi[data_op][s])
         hist_btagheavydown[data_op][s].Scale(lumi_data/lumi[data_op][s])
      if args.norm: hist_nom[data_op][s].Scale(norm_factor[str(args.year)])
      if not args.nosyst:
         if args.norm: hist_btaglightup[data_op][s].Scale(norm_factor[str(args.year)])
         if args.norm: hist_btaglightdown[data_op][s].Scale(norm_factor[str(args.year)])
         if args.norm: hist_btagheavyup[data_op][s].Scale(norm_factor[str(args.year)])
         if args.norm: hist_btagheavydown[data_op][s].Scale(norm_factor[str(args.year)])
    ## Fixing single top
    #print(samples_foryear[data_op]) 
    #### List of summing samples:
    list_st = ["st_1","st_2","st_3","st_4"]
    list_wjets = ["wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8"]
    list_zjets = ["zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8"]
    list_vv = ["ww","wz","zz"]
    list_ttbar_nomuon = ["ttbar_sl_nomuon","ttbar_dl_nomuon","ttbar_dh"]
    list_ttbar_muon_charm = ["ttbar_sl_muon_charm","ttbar_dl_muon_charm"]
    list_ttbar_muon_bottom = ["ttbar_sl_muon_bottom","ttbar_dl_muon_bottom"]
    list_ttbar_muon_else = ["ttbar_sl_muon_else","ttbar_dl_muon_else"]
    hist_nom[data_op]["ttbar_nomuon"] = hist_nom[data_op][list_ttbar_nomuon[0]]
    if not args.nosyst:
      hist_btaglightup[data_op]["ttbar_nomuon"] = hist_btaglightup[data_op][list_ttbar_nomuon[0]]
      hist_btaglightdown[data_op]["ttbar_nomuon"] = hist_btaglightdown[data_op][list_ttbar_nomuon[0]]
      hist_btagheavyup[data_op]["ttbar_nomuon"] = hist_btagheavyup[data_op][list_ttbar_nomuon[0]]
      hist_btagheavydown[data_op]["ttbar_nomuon"] = hist_btagheavydown[data_op][list_ttbar_nomuon[0]]
    for l in list_ttbar_nomuon[1:]:
      hist_nom[data_op]["ttbar_nomuon"].Add(hist_nom[data_op][l])
      if not args.nosyst:
        hist_btaglightup[data_op]["ttbar_nomuon"].Add(hist_btaglightup[data_op][l])
        hist_btaglightdown[data_op]["ttbar_nomuon"].Add(hist_btaglightdown[data_op][l])
        hist_btagheavyup[data_op]["ttbar_nomuon"].Add(hist_btagheavyup[data_op][l])
        hist_btagheavydown[data_op]["ttbar_nomuon"].Add(hist_btagheavydown[data_op][l])
    hist_nom[data_op]["ttbar_muon_charm"] = hist_nom[data_op][list_ttbar_muon_charm[0]]
    if not args.nosyst:
      hist_btaglightup[data_op]["ttbar_muon_charm"] = hist_btaglightup[data_op][list_ttbar_muon_charm[0]]
      hist_btaglightdown[data_op]["ttbar_muon_charm"] = hist_btaglightdown[data_op][list_ttbar_muon_charm[0]]
      hist_btagheavyup[data_op]["ttbar_muon_charm"] = hist_btagheavyup[data_op][list_ttbar_muon_charm[0]]
      hist_btagheavydown[data_op]["ttbar_muon_charm"] = hist_btagheavydown[data_op][list_ttbar_muon_charm[0]]
    for l in list_ttbar_muon_charm[1:]:
      hist_nom[data_op]["ttbar_muon_charm"].Add(hist_nom[data_op][l])
      if not args.nosyst:
        hist_btaglightup[data_op]["ttbar_muon_charm"].Add(hist_btaglightup[data_op][l])
        hist_btaglightdown[data_op]["ttbar_muon_charm"].Add(hist_btaglightdown[data_op][l])
        hist_btagheavyup[data_op]["ttbar_muon_charm"].Add(hist_btagheavyup[data_op][l])
        hist_btagheavydown[data_op]["ttbar_muon_charm"].Add(hist_btagheavydown[data_op][l])
    hist_nom[data_op]["ttbar_muon_bottom"] = hist_nom[data_op][list_ttbar_muon_bottom[0]]
    if not args.nosyst:
      hist_btaglightup[data_op]["ttbar_muon_bottom"] = hist_btaglightup[data_op][list_ttbar_muon_bottom[0]]
      hist_btaglightdown[data_op]["ttbar_muon_bottom"] = hist_btaglightdown[data_op][list_ttbar_muon_bottom[0]]
      hist_btagheavyup[data_op]["ttbar_muon_bottom"] = hist_btagheavyup[data_op][list_ttbar_muon_bottom[0]]
      hist_btagheavydown[data_op]["ttbar_muon_bottom"] = hist_btagheavydown[data_op][list_ttbar_muon_bottom[0]]
    for l in list_ttbar_muon_bottom[1:]:
      hist_nom[data_op]["ttbar_muon_bottom"].Add(hist_nom[data_op][l])
      if not args.nosyst:
        hist_btaglightup[data_op]["ttbar_muon_bottom"].Add(hist_btaglightup[data_op][l])
        hist_btaglightdown[data_op]["ttbar_muon_bottom"].Add(hist_btaglightdown[data_op][l])
        hist_btagheavyup[data_op]["ttbar_muon_bottom"].Add(hist_btagheavyup[data_op][l])
        hist_btagheavydown[data_op]["ttbar_muon_bottom"].Add(hist_btagheavydown[data_op][l])
    hist_nom[data_op]["ttbar_muon_else"] = hist_nom[data_op][list_ttbar_muon_else[0]]
    if not args.nosyst:
      hist_btaglightup[data_op]["ttbar_muon_else"] = hist_btaglightup[data_op][list_ttbar_muon_else[0]]
      hist_btaglightdown[data_op]["ttbar_muon_else"] = hist_btaglightdown[data_op][list_ttbar_muon_else[0]]
      hist_btagheavyup[data_op]["ttbar_muon_else"] = hist_btagheavyup[data_op][list_ttbar_muon_else[0]]
      hist_btagheavydown[data_op]["ttbar_muon_else"] = hist_btagheavydown[data_op][list_ttbar_muon_else[0]]
    for l in list_ttbar_muon_else[1:]:
      hist_nom[data_op]["ttbar_muon_else"].Add(hist_nom[data_op][l])
      if not args.nosyst:
        hist_btaglightup[data_op]["ttbar_muon_else"].Add(hist_btaglightup[data_op][l])
        hist_btaglightdown[data_op]["ttbar_muon_else"].Add(hist_btaglightdown[data_op][l])
        hist_btagheavyup[data_op]["ttbar_muon_else"].Add(hist_btagheavyup[data_op][l])
        hist_btagheavydown[data_op]["ttbar_muon_else"].Add(hist_btagheavydown[data_op][l])
    hist_nom[data_op]["st"] = hist_nom[data_op][list_st[0]]
    if not args.nosyst:
      hist_btaglightup[data_op]["st"] = hist_btaglightup[data_op][list_st[0]]
      hist_btaglightdown[data_op]["st"] = hist_btaglightdown[data_op][list_st[0]]
      hist_btagheavyup[data_op]["st"] = hist_btagheavyup[data_op][list_st[0]]
      hist_btagheavydown[data_op]["st"] = hist_btagheavydown[data_op][list_st[0]]
    for l in list_st[1:]:
      hist_nom[data_op]["st"].Add(hist_nom[data_op][l])
      if not args.nosyst:
        hist_btaglightup[data_op]["st"].Add(hist_btaglightup[data_op][l])
        hist_btaglightdown[data_op]["st"].Add(hist_btaglightdown[data_op][l])
        hist_btagheavyup[data_op]["st"].Add(hist_btagheavyup[data_op][l])
        hist_btagheavydown[data_op]["st"].Add(hist_btagheavydown[data_op][l])
    hist_nom[data_op]["wjets"] = hist_nom[data_op][list_wjets[0]]
    if not args.nosyst:
      hist_btaglightup[data_op]["wjets"] = hist_btaglightup[data_op][list_wjets[0]]
      hist_btaglightdown[data_op]["wjets"] = hist_btaglightdown[data_op][list_wjets[0]]
      hist_btagheavyup[data_op]["wjets"] = hist_btagheavyup[data_op][list_wjets[0]]
      hist_btagheavydown[data_op]["wjets"] = hist_btagheavydown[data_op][list_wjets[0]]
    for l in list_wjets[1:]:
      hist_nom[data_op]["wjets"].Add(hist_nom[data_op][l])
      if not args.nosyst:
        hist_btaglightup[data_op]["wjets"].Add(hist_btaglightup[data_op][l])
        hist_btaglightdown[data_op]["wjets"].Add(hist_btaglightdown[data_op][l])
        hist_btagheavyup[data_op]["wjets"].Add(hist_btagheavyup[data_op][l])
        hist_btagheavydown[data_op]["wjets"].Add(hist_btagheavydown[data_op][l])
    hist_nom[data_op]["zjets"] = hist_nom[data_op][list_zjets[0]]
    if not args.nosyst:
      hist_btaglightup[data_op]["zjets"] = hist_btaglightup[data_op][list_zjets[0]]
      hist_btaglightdown[data_op]["zjets"] = hist_btaglightdown[data_op][list_zjets[0]]
      hist_btagheavyup[data_op]["zjets"] = hist_btagheavyup[data_op][list_zjets[0]]
      hist_btagheavydown[data_op]["zjets"] = hist_btagheavydown[data_op][list_zjets[0]]
    for l in list_zjets[1:]:
      hist_nom[data_op]["zjets"].Add(hist_nom[data_op][l])
      if not args.nosyst:
        hist_btaglightup[data_op]["zjets"].Add(hist_btaglightup[data_op][l])
        hist_btaglightdown[data_op]["zjets"].Add(hist_btaglightdown[data_op][l])
        hist_btagheavyup[data_op]["zjets"].Add(hist_btagheavyup[data_op][l])
        hist_btagheavydown[data_op]["zjets"].Add(hist_btagheavydown[data_op][l])
    hist_nom[data_op]["vv"] = hist_nom[data_op][list_vv[0]]
    if not args.nosyst:
      hist_btaglightup[data_op]["vv"] = hist_btaglightup[data_op][list_vv[0]]
      hist_btaglightdown[data_op]["vv"] = hist_btaglightdown[data_op][list_vv[0]]
      hist_btagheavyup[data_op]["vv"] = hist_btagheavyup[data_op][list_vv[0]]
      hist_btagheavydown[data_op]["vv"] = hist_btagheavydown[data_op][list_vv[0]]
    for l in list_vv[1:]:
      hist_nom[data_op]["vv"].Add(hist_nom[data_op][l])
      if not args.nosyst:
        hist_btaglightup[data_op]["vv"].Add(hist_btaglightup[data_op][l])
        hist_btaglightdown[data_op]["vv"].Add(hist_btaglightdown[data_op][l])
        hist_btagheavyup[data_op]["vv"].Add(hist_btagheavyup[data_op][l])
        hist_btagheavydown[data_op]["vv"].Add(hist_btagheavydown[data_op][l])

  samples = ["ttbar_muon_bottom","ttbar_muon_charm","ttbar_muon_else","ttbar_nomuon","zjets","vv","st","wjets"]

  ## sumamos todos los anos para todos los histogramas
  histT_nom = {}
  histT_btaglightup = {}
  histT_btaglightdown = {}
  histT_btagheavyup = {}
  histT_btagheavydown = {}
  for s in samples:
       histT_nom[s] = hist_nom[datayears[0]][s]
       if not args.nosyst:
         histT_btaglightup[s] = hist_btaglightup[datayears[0]][s]
         histT_btaglightdown[s] = hist_btaglightdown[datayears[0]][s]
         histT_btagheavyup[s] = hist_btagheavyup[datayears[0]][s]
         histT_btagheavydown[s] = hist_btagheavydown[datayears[0]][s]
       for d in datayears[1:]:
          histT_nom[s].Add(hist_nom[d][s])
          if not args.nosyst:
            histT_btaglightup[s].Add(hist_btaglightup[d][s])
            histT_btaglightdown[s].Add(hist_btaglightdown[d][s])
            histT_btagheavyup[s].Add(hist_btagheavyup[d][s])
            histT_btagheavydown[s].Add(hist_btagheavydown[d][s])

  if not args.nodata:
    histD = histdata[datayears[0]]
    for d in datayears[1:]:
       histD.Add(histdata[d])

  if not args.nosyst:
    ### Para los histogramas de sistemática up y down sumamos todas las contribuciones
    histT_sT_btaglightup = histT_btaglightup[samples[0]]
    histT_sT_btaglightdown = histT_btaglightdown[samples[0]]
    histT_sT_btagheavyup = histT_btagheavyup[samples[0]]
    histT_sT_btagheavydown = histT_btagheavydown[samples[0]]
    for s in samples[1:]:
      histT_sT_btaglightup.Add(histT_btaglightup[s])
      histT_sT_btaglightdown.Add(histT_btaglightdown[s])
      histT_sT_btagheavyup.Add(histT_btagheavyup[s])
      histT_sT_btagheavydown.Add(histT_btagheavydown[s])

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  colors = {}
  colors["ttbar_dl_muon_charm"] = (240,76,76)
  colors["wjets"] = (45,66,172)
  colors["ttbar_muon_charm"] = (229,167,204)
  colors["ttbar_muon_bottom"] = (167,229,214)
  colors["ttbar_muon_else"] = (222,212,74)
  colors["ttbar_nomuon"] = (242,193,121)
  colors["vv"] = (255,180,85)
  colors["ttbar_dh"] = (204,204,0)
  colors["zjets"] = (26,117,35)
  colors["st"] = (153,51,255)

  if args.stack:
    ymax = 0
    ymin = 0
    for s in samples:
      histT_nom[s].SetLineWidth(1)
      histT_nom[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_nom[s].GetYaxis().SetTitle("Number of events")
      histT_nom[s].GetXaxis().SetTitle(name)

      y = histT_nom[s].GetMaximum()
      ym = histT_nom[s].GetMinimum()
      if y>ymax: ymax=y
      if ymin>ym: ymin=ym
      if not args.nodata:
        y = histD.GetMaximum()
        if y>ymax: ymax=y

      #histT_M["vv"].SetMinimum(1000.)
      #histT_M["vv"].SetMaximum(3*ymax_M)
      if args.linear: histT_nom["vv"].SetMaximum(1.3*ymax)
      if args.linear: histT_nom["vv"].SetMinimum(1.3*ymin)

    ## Stack creation
    samples = ["vv","zjets","wjets","ttbar_nomuon","st","ttbar_muon_else","ttbar_muon_charm","ttbar_muon_bottom"]

    if args.ratio and not args.nodata: upper_pad.cd()
    stack = ROOT.THStack()
    for s in samples:
      stack.Add(histT_nom[s])

    last = stack.GetStack().Last()

    if not args.nosyst:
      ### Creamos el tgraph de sistematicos
      envHi_btaglight = histT_sT_btaglightup;
      envLo_btaglight = histT_sT_btaglightdown;
      envHi_btagheavy = histT_sT_btagheavyup;
      envLo_btagheavy = histT_sT_btagheavydown;

      graph_err = TGraphAsymmErrors();
      ratio_graph_err = TGraphAsymmErrors();

      for bin in range(last.GetNbinsX()):
        #sysup = std::sqrt((nom->GetBinContent(bin+1)-envHi->GetBinContent(bin+1))*(nom->GetBinContent(bin+1)-envHi->GetBinContent(bin+1)));
        #sysdown = std::sqrt((nom->GetBinContent(bin+1)-envLo->GetBinContent(bin+1))*(nom->GetBinContent(bin+1)-envLo->GetBinContent(bin+1)));
        bin_nom = float(last.GetBinContent(bin+1));
        bin_btaglightup = float(envHi_btaglight.GetBinContent(bin+1));
        bin_btaglightdown = float(envLo_btaglight.GetBinContent(bin+1));
        bin_btagheavyup = float(envHi_btagheavy.GetBinContent(bin+1));
        bin_btagheavydown = float(envLo_btagheavy.GetBinContent(bin+1));

        sysup_btaglight = abs(bin_nom-bin_btaglightup);
        sysdown_btaglight = abs(bin_nom-bin_btaglightdown);

        sysup_btagheavy = abs(bin_nom-bin_btagheavyup);
        sysdown_btagheavy = abs(bin_nom-bin_btagheavydown);

        sysup_total = sqrt(sysup_btaglight**2 + sysup_btagheavy**2);
        sysdown_total = sqrt(sysdown_btaglight**2 + sysdown_btagheavy**2);

        graph_err.SetPointEXhigh(bin,last.GetBinWidth(bin+1)/2);
        graph_err.SetPointEXlow(bin,last.GetBinWidth(bin+1)/2);
        graph_err.SetPoint(bin,last.GetBinCenter(bin+1),last.GetBinContent(bin+1));
        graph_err.SetPointEYhigh(bin,sysup_total);
        graph_err.SetPointEYlow(bin,sysdown_total);

        ratio_graph_err.SetPointEXhigh(bin,last.GetBinWidth(bin+1)/2);
        ratio_graph_err.SetPointEXlow(bin,last.GetBinWidth(bin+1)/2);
        ratio_graph_err.SetPoint(bin,last.GetBinCenter(bin+1),1.);
        if (last.GetBinContent(bin+1)>0.01):
           ratio_graph_err.SetPointEYhigh(bin,abs(sysup_total/bin_nom));
           ratio_graph_err.SetPointEYlow(bin,abs(sysdown_total/bin_nom));
        else:
           ratio_graph_err.SetPointEYhigh(bin,0.);
           ratio_graph_err.SetPointEYlow(bin,0.);

      graph_err.SetFillColorAlpha(kGray+2,0.7);
      graph_err.SetMarkerColor(kGray);
      graph_err.SetMarkerSize(0);
      graph_err.SetFillStyle(3001);
      ratio_graph_err.SetFillColorAlpha(kGray+2,0.7);
      ratio_graph_err.SetMarkerColor(kGray);
      ratio_graph_err.SetMarkerSize(0);
      ratio_graph_err.SetFillStyle(3001);

    y = stack.GetMaximum()
    if y>ymax: ymax=y
    #stack.SetMinimum(1000.)
    #stack.SetMaximum(3*ymax)
    if args.linear: stack.SetMaximum(1.3*ymax)
    if args.linear: stack.SetMinimum(1.)

    stack.Draw("HIST")
    if not args.nodata:
      histD.SetMarkerStyle(20)
      histD.SetMarkerSize(0.4)
      histD.SetLineWidth(1)
      histD.SetLineColor(ROOT.kBlack)
      histD.Draw("E SAME")
    if not args.nosyst: graph_err.Draw("SAME 2")

    if args.ratio and not args.nodata:
      lower_pad.cd()
      ratio = histD.Clone("ratio")
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
      hTotal = histT_nom["vv"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histT_nom[s])
      if name in ["muon_bot1_iso","muon_bot2_iso","muon_bot1_iso_log","muon_bot2_iso_log","muon_bot1_z","muon_bot2_z","muon_bot1_iso_abs","muon_bot2_iso_abs"]:
        ratio_num_aux = ratio.Clone("ratio_num_"+str(name))
        ratio_den_aux = hTotal.Clone("ratio_den_"+str(name))
        myratiofile.WriteObject(ratio_num_aux,"ratio_num_"+str(name))
        myratiofile.WriteObject(ratio_den_aux,"ratio_den_"+str(name))
      ratio.Divide(hTotal)
      ratio.Draw("ep")
      if not args.nosyst: ratio_graph_err.Draw("same 2")

    if name == "InvM_2jets" and (not args.nodata):
      print("Integral of data is "+str(histD.Integral()))
      print("Integral of MC is "+str(hTotal.Integral()))
      print("Ratio is "+str(histD.Integral()/hTotal.Integral()))
      error_tot = {}
      error_fulltot = 0
      for s in samples:
          print("Integral of "+str(s)+" MC process is "+str(histT_nom[s].Integral()))
          error_tot[s] = 0
          for bin in range(histT_nom[s].GetNbinsX()):
              error_tot[s] = error_tot[s] + histT_nom[s].GetBinError(bin+1)**2
          error_tot[s] = sqrt(error_tot[s])
          print("error is "+str(error_tot[s]))
          error_fulltot = error_fulltot + error_tot[s]**2
      print("Final MC error is "+str(sqrt(error_fulltot)))

    ## Legends
    if args.ratio and not args.nodata: upper_pad.cd()
    leg = TLegend(0.68,0.59,0.95,0.95)
    leg.SetBorderSize(1)
    #leg.AddEntry(histT_nom["vv"],"VV","f")
    leg.AddEntry(histT_nom["ttbar_muon_charm"],"t#bar{t}, muon from D hadron","f")
    leg.AddEntry(histT_nom["ttbar_muon_bottom"],"t#bar{t}, muon from B hadron","f")
    leg.AddEntry(histT_nom["ttbar_muon_else"],"t#bar{t}, else","f")
    #leg.AddEntry(histT_nom["ttbar_nomuon"],"ttbar, no muon","f")
    leg.AddEntry(histT_nom["st"],"Single top","f")
    leg.AddEntry(histT_nom["zjets"],"Z+jets","f")
    leg.AddEntry(histT_nom["wjets"],"W+jets","f")
    #leg.AddEntry(histT_nom["ttbar_dh"],"Hadronic t#bar{t}","f")
    if args.stack and not args.nodata: leg.AddEntry(histD, "Data" ,"lep")
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

  for data_op in datayears:
                        histFile[name][data_op].Close()
                        if not args.nodata:
                          histFileDM[name][data_op].Close()
                          histFileDE[name][data_op].Close()

myratiofile.Close()

