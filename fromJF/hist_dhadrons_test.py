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
parser.add_argument("--nosyst", action="store_true", default=False,
                    help="systematics inclusion")

# Use like:
# python higgssearch/fromJF/hist_fromJFwqq_syst_total.py --stack --ratio --png --norm --year="all" --channel="btagMM_chitest_sl"

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

sl_channel = ["btagMM_chitest_sl","lepton50_chitest_sl","btagMM_chitest_slss","lepton50_chitest_slss",
      "btagMM_chitest_slos","lepton50_chitest_slos","btagMM_chitest_slssos","lepton50_chitest_slssos"]

plotdir = '/nfs/cms/vazqueze/higgssearch/plotspng/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

c_rat = 1.2
c_rat2 = 0.8

filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/dhadrons_test/"

term = "hist_wqqfromJF_"

observable_names = ["nJetGood", "jet_1_pt", "counts_dplus","counts_dzero","counts_dpluss","counts_dlambda","counts_drest","counts_charm","indx_dplus",
      "indx_dzero","indx_dpluss","indx_dlambda","indx_drest", "indx_sum"]

if "sl" in str(args.channel): observable_names = observable_names + ["muon_jet_pt","muon_jet_relpt","muon_jet_eta"]

#observable_names = [""]

not_rebin = ["nJetGood","InvM3_good","InvM3_bad","InvMl_good","InvMl_bad","lepton_eta_thick","jet_bot1_btagnumber", "jet_bot2_btagnumber", 
      "jet_1_btagnumber", "jet_2_btagnumber","muon_jet_pt","muon_jet_relpt","muon_jet_eta","tau_discr","transverse_mass",
      "jet_1_flavourP", "jet_2_flavourP", "jet_bot1_flavourP", "jet_bot2_flavourP"]

#observable_names = ["InvM_2jets","chi2_test0","chi2_test1","chi2_test_good","chi2_test_bad"]
#observable_names = ["jet_1_flavourP", "jet_2_flavourP", "jet_bot1_flavourP", "jet_bot2_flavourP"]

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
        end_term = ".root"
        ## mc files
        histFile[data_op] = {}
        for name in observable_names:
                if isfile(filePath + term+"MC_"+data_op+"_"+name+end_term):
                        histFile[data_op][name] = TFile.Open(filePath + term+"MC_"+data_op+"_"+name+end_term,"READ")
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
  hist_nom = {}
  for data_op in datayears:
    hist_nom[data_op] = {}
    data_term = data_op
    #print(data_op+"M_"+name+"_M")
    #print(histFile[data_op][name].ls())
    for s in samples:
      if s[0:8] == "ttbar_sl": 
         s_term = s[0:8]+data_term+s[8:] 
      else:
         s_term = s+data_term
      hist_nom[data_op][s] = histFile[data_op][name].Get(s_term+"_"+name)

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
    ## Fixing single top
    #print(samples_foryear[data_op]) 
    #### List of summing samples:
    list_st = ["st_1","st_2","st_3","st_4"]
    list_wjets = ["wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8"]
    list_zjets = ["zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8"]
    list_vv = ["ww","wz","zz"]
    hist_nom[data_op]["st"] = hist_nom[data_op][list_st[0]]
    for l in list_st[1:]:
      hist_nom[data_op]["st"].Add(hist_nom[data_op][l])
    hist_nom[data_op]["wjets"] = hist_nom[data_op][list_wjets[0]]
    for l in list_wjets[1:]:
      hist_nom[data_op]["wjets"].Add(hist_nom[data_op][l])
    hist_nom[data_op]["zjets"] = hist_nom[data_op][list_zjets[0]]
    for l in list_zjets[1:]:
      hist_nom[data_op]["zjets"].Add(hist_nom[data_op][l])
    hist_nom[data_op]["vv"] = hist_nom[data_op][list_vv[0]]
    for l in list_vv[1:]:
      hist_nom[data_op]["vv"].Add(hist_nom[data_op][l])

  samples = ["ttbar_sl_bottom","ttbar_sl_charm","ttbar_sl_else","ttbar_sl_light","ttbar_dl","ttbar_dh","zjets","vv","st","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon"]

  ## sumamos todos los anos para todos los histogramas
  histT_nom = {}
  for s in samples:
       histT_nom[s] = hist_nom[datayears[0]][s]
       for d in datayears[1:]:
          histT_nom[s].Add(hist_nom[d][s])

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

      #histT_M["vv"].SetMinimum(1000.)
      #histT_M["vv"].SetMaximum(3*ymax_M)
      if args.linear: histT_nom["vv"].SetMaximum(1.3*ymax)
      if args.linear: histT_nom["vv"].SetMinimum(1.3*ymin)

    ## Stack creation
    samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon","ttbar_sl_else","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]

    if args.ratio and not args.nodata: upper_pad.cd()
    stack = ROOT.THStack()
    for s in samples:
      stack.Add(histT_nom[s])

    last = stack.GetStack().Last()

    y = stack.GetMaximum()
    if y>ymax: ymax=y
    if args.linear: stack.SetMaximum(1.3*ymax)
    if args.linear: stack.SetMinimum(1.)

    stack.Draw("HIST")

    stack.GetXaxis().SetTitle(name);gPad.Modified();gPad.Update();
    if args.ratio and not args.nodata: upper_pad.cd()
    leg = TLegend(0.7,0.6,0.89,0.89)
    leg.SetBorderSize(1)
    #leg.AddEntry(histT_nom["vv"],"VV","f")
    leg.AddEntry(histT_nom["ttbar_sl_charm"],"t#bar{t} cq","f")
    leg.AddEntry(histT_nom["ttbar_sl_light"],"t#bar{t} qq'","f")
    leg.AddEntry(histT_nom["st"],"Single top","f")
    leg.AddEntry(histT_nom["ttbar_sl_bottom"],"t#bar{t} bq","f")
    leg.AddEntry(histT_nom["ttbar_sl_else"],"t#bar{t} qg,gg","f")
    leg.AddEntry(histT_nom["ttbar_sl_charmgluon"],"t#bar{t} cg","f")
    leg.AddEntry(histT_nom["ttbar_sl_bottomgluon"],"t#bar{t} bg","f")
    leg.AddEntry(histT_nom["zjets"],"Z+jets","f")
    leg.AddEntry(histT_nom["wjets"],"W+jets","f")
    leg.AddEntry(histT_nom["ttbar_dl"],"Dileptonic t#bar{t}","f")
    #leg.AddEntry(histT_nom_M["ttbar_dh"],"Hadronic t#bar{t}","f")
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

    for data_op in datayears:
                        histFile[data_op][name].Close()

