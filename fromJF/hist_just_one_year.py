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
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Are these ssos plots?")
parser.add_argument("--ratio", action="store_true", default=False,
                    help="Plot ratio or not")
parser.add_argument("--linear", action="store_true", default=False,
                    help="Plot linearly")
parser.add_argument("--png", action="store_true", default=False,
                    help="png format")
parser.add_argument("--qcd", action="store_true", default=False,
                    help="include qcd samples")
parser.add_argument("--nodata", action="store_true", default=False,
                    help="Do not plot data")
parser.add_argument("--year", type=string, default="2016",
                    help="Select year of process to run")
parser.add_argument("--type", type=string, default="sl",
                    help="Selec type of channel to run")

# Use like:
# python arg.py --data="No"
# python hist_plot.py --data="No" --stack --ratio

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

plotdir = '/nfs/cms/vazqueze/higgssearch/plotspng/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

## Open hists files

if args.type == "sv": filePath = "/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/sv/"
else:  filePath = "/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/"

if args.type == "sv": term = "histstt_SV_fromJF"
else: term = "histstt_SL_fromJF"

datayears = ["2016","2016B","2017","2018"]
#datayears = ["2018","2016","2016B"]

samplesHT = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_light","ttbar_sl_bottom","ttbar_dl","ttbar_dh","zz","wz",
        "st_1","st_2","st_3","st_4"]

## Adding QCD

histFile = {}

for data_op in datayears:
        ## mc files
        histFile[data_op] = {}
        for s in samplesHT:
                if s[0:8] == "ttbar_sl" and isfile(filePath + term+"_"+s[0:8]+data_op+s[8:]+".root"):
                        histFile[data_op][s] = TFile.Open(filePath + term+"_"+s[0:8]+data_op+s[8:]+".root","READ")
                elif isfile(filePath + term+"_"+s+data_op+".root"):
                        histFile[data_op][s] = TFile.Open(filePath + term+"_"+s+data_op+".root","READ")
        #print(histFile[data_op].keys())

histFileD = {}

for data_op in datayears:
	histFileD[data_op] = {}
	# data files
	histFileD[data_op]["M"] = TFile.Open(filePath + term+"_"+data_op+"M.root","READ")
	histFileD[data_op]["E"] = TFile.Open(filePath + term+"_"+data_op+"E.root","READ")

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

lumi_d = {}
lumi_d["2016"] = 19.5
lumi_d["2016B"] = 16.8
lumi_d["2017"] = 41.5
lumi_d["2018"] = 59.8

histNames = ["InvM_2jets"]

#######################################################
########### Start of plot creation ####################
#######################################################

for name in histNames:

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
        "st_1","st_2","st_3","st_4"]

  ## HISTS
  samples_foryear = {}
  hist_M_ss = {}
  hist_M_os = {}
  hist_E_ss = {}
  hist_E_os = {}
  histdata_M_ss = {}
  histdata_M_os = {}
  histdata_E_ss = {}
  histdata_E_os = {}
  for data_op in datayears:
    samples_foryear[data_op] = [s for s in samples if s in histFile[data_op].keys()]
    hist_M_ss[data_op] = {}
    hist_M_os[data_op] = {}
    hist_E_ss[data_op] = {}
    hist_E_os[data_op] = {}
    for s in samples_foryear[data_op]:
      hist_M_ss[data_op][s] = histFile[data_op][s].Get(name+"_M_ss")
      hist_M_os[data_op][s] = histFile[data_op][s].Get(name+"_M_os")
      hist_E_ss[data_op][s] = histFile[data_op][s].Get(name+"_E_ss")
      hist_E_os[data_op][s] = histFile[data_op][s].Get(name+"_E_os")
    if not args.nodata: histdata_M_ss[data_op] = histFileD[data_op]["M"].Get(name+"_M_ss")
    if not args.nodata: histdata_M_os[data_op] = histFileD[data_op]["M"].Get(name+"_M_os")
    if not args.nodata: histdata_E_ss[data_op] = histFileD[data_op]["E"].Get(name+"_E_ss")
    if not args.nodata: histdata_E_os[data_op] = histFileD[data_op]["E"].Get(name+"_E_os")
  
  samples_st = {}
  samples_wjets = {}
  samples_zjets = {}
  ## Scaling to lumi
  #print(name) 
  for data_op in datayears:
    lumi_data = lumi_d[data_op]
    for s in samples_foryear[data_op]:
      #print(s)
      #print(data_op)
      hist_M_ss[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_M_os[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_E_ss[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_E_os[data_op][s].Scale(lumi_data/lumi[data_op][s])
      ## Fixing single top
    #print(samples_stc[data_op],samples_stnc[data_op],samples_wjetsc[data_op],samples_wjetsdc[data_op],samples_wjetsl[data_op],samples_wjetsb[data_op],samples_zjets[data_op]) 
    hist_M_ss[data_op]["st"] = hist_M_ss[data_op]["st_1"]
    hist_M_os[data_op]["st"] = hist_M_os[data_op]["st_1"]
    hist_E_ss[data_op]["st"] = hist_E_ss[data_op]["st_1"]
    hist_E_os[data_op]["st"] = hist_E_os[data_op]["st_1"]
    hist_M_ss[data_op]["st"].Add(hist_M_ss[data_op]["st_2"])
    hist_M_os[data_op]["st"].Add(hist_M_os[data_op]["st_2"])
    hist_E_ss[data_op]["st"].Add(hist_E_ss[data_op]["st_2"])
    hist_E_os[data_op]["st"].Add(hist_E_os[data_op]["st_2"])
    hist_M_ss[data_op]["st"].Add(hist_M_ss[data_op]["st_3"])
    hist_M_os[data_op]["st"].Add(hist_M_os[data_op]["st_3"])
    hist_E_ss[data_op]["st"].Add(hist_E_ss[data_op]["st_3"])
    hist_E_os[data_op]["st"].Add(hist_E_os[data_op]["st_3"])
    hist_M_ss[data_op]["st"].Add(hist_M_ss[data_op]["st_4"])
    hist_M_os[data_op]["st"].Add(hist_M_os[data_op]["st_4"])
    hist_E_ss[data_op]["st"].Add(hist_E_ss[data_op]["st_4"])
    hist_E_os[data_op]["st"].Add(hist_E_os[data_op]["st_4"])
    hist_M_ss[data_op]["wjets"] = hist_M_ss[data_op]["wjets_1"]
    hist_M_os[data_op]["wjets"] = hist_M_os[data_op]["wjets_1"]
    hist_E_ss[data_op]["wjets"] = hist_E_ss[data_op]["wjets_1"]
    hist_E_os[data_op]["wjets"] = hist_E_os[data_op]["wjets_1"]
    hist_M_ss[data_op]["wjets"].Add(hist_M_ss[data_op]["wjets_2"])
    hist_M_os[data_op]["wjets"].Add(hist_M_os[data_op]["wjets_2"])
    hist_E_ss[data_op]["wjets"].Add(hist_E_ss[data_op]["wjets_2"])
    hist_E_os[data_op]["wjets"].Add(hist_E_os[data_op]["wjets_2"])
    hist_M_ss[data_op]["wjets"].Add(hist_M_ss[data_op]["wjets_3"])
    hist_M_os[data_op]["wjets"].Add(hist_M_os[data_op]["wjets_3"])
    hist_E_ss[data_op]["wjets"].Add(hist_E_ss[data_op]["wjets_3"])
    hist_E_os[data_op]["wjets"].Add(hist_E_os[data_op]["wjets_3"])
    hist_M_ss[data_op]["wjets"].Add(hist_M_ss[data_op]["wjets_4"])
    hist_M_os[data_op]["wjets"].Add(hist_M_os[data_op]["wjets_4"])
    hist_E_ss[data_op]["wjets"].Add(hist_E_ss[data_op]["wjets_4"])
    hist_E_os[data_op]["wjets"].Add(hist_E_os[data_op]["wjets_4"])
    hist_M_ss[data_op]["wjets"].Add(hist_M_ss[data_op]["wjets_5"])
    hist_M_os[data_op]["wjets"].Add(hist_M_os[data_op]["wjets_5"])
    hist_E_ss[data_op]["wjets"].Add(hist_E_ss[data_op]["wjets_5"])
    hist_E_os[data_op]["wjets"].Add(hist_E_os[data_op]["wjets_5"])
    hist_M_ss[data_op]["wjets"].Add(hist_M_ss[data_op]["wjets_6"])
    hist_M_os[data_op]["wjets"].Add(hist_M_os[data_op]["wjets_6"])
    hist_E_ss[data_op]["wjets"].Add(hist_E_ss[data_op]["wjets_6"])
    hist_E_os[data_op]["wjets"].Add(hist_E_os[data_op]["wjets_6"])
    hist_M_ss[data_op]["wjets"].Add(hist_M_ss[data_op]["wjets_7"])
    hist_M_os[data_op]["wjets"].Add(hist_M_os[data_op]["wjets_7"])
    hist_E_ss[data_op]["wjets"].Add(hist_E_ss[data_op]["wjets_7"])
    hist_E_os[data_op]["wjets"].Add(hist_E_os[data_op]["wjets_7"])
    hist_M_ss[data_op]["wjets"].Add(hist_M_ss[data_op]["wjets_8"])
    hist_M_os[data_op]["wjets"].Add(hist_M_os[data_op]["wjets_8"])
    hist_E_ss[data_op]["wjets"].Add(hist_E_ss[data_op]["wjets_8"])
    hist_E_os[data_op]["wjets"].Add(hist_E_os[data_op]["wjets_8"])
    hist_M_ss[data_op]["zjets"] = hist_M_ss[data_op]["zjets_1"]
    hist_M_os[data_op]["zjets"] = hist_M_os[data_op]["zjets_1"]
    hist_E_ss[data_op]["zjets"] = hist_E_ss[data_op]["zjets_1"]
    hist_E_os[data_op]["zjets"] = hist_E_os[data_op]["zjets_1"]
    hist_M_ss[data_op]["zjets"].Add(hist_M_ss[data_op]["zjets_2"])
    hist_M_os[data_op]["zjets"].Add(hist_M_os[data_op]["zjets_2"])
    hist_E_ss[data_op]["zjets"].Add(hist_E_ss[data_op]["zjets_2"])
    hist_E_os[data_op]["zjets"].Add(hist_E_os[data_op]["zjets_2"])
    hist_M_ss[data_op]["zjets"].Add(hist_M_ss[data_op]["zjets_3"])
    hist_M_os[data_op]["zjets"].Add(hist_M_os[data_op]["zjets_3"])
    hist_E_ss[data_op]["zjets"].Add(hist_E_ss[data_op]["zjets_3"])
    hist_E_os[data_op]["zjets"].Add(hist_E_os[data_op]["zjets_3"])
    hist_M_ss[data_op]["zjets"].Add(hist_M_ss[data_op]["zjets_4"])
    hist_M_os[data_op]["zjets"].Add(hist_M_os[data_op]["zjets_4"])
    hist_E_ss[data_op]["zjets"].Add(hist_E_ss[data_op]["zjets_4"])
    hist_E_os[data_op]["zjets"].Add(hist_E_os[data_op]["zjets_4"])
    hist_M_ss[data_op]["zjets"].Add(hist_M_ss[data_op]["zjets_5"])
    hist_M_os[data_op]["zjets"].Add(hist_M_os[data_op]["zjets_5"])
    hist_E_ss[data_op]["zjets"].Add(hist_E_ss[data_op]["zjets_5"])
    hist_E_os[data_op]["zjets"].Add(hist_E_os[data_op]["zjets_5"])
    hist_M_ss[data_op]["zjets"].Add(hist_M_ss[data_op]["zjets_6"])
    hist_M_os[data_op]["zjets"].Add(hist_M_os[data_op]["zjets_6"])
    hist_E_ss[data_op]["zjets"].Add(hist_E_ss[data_op]["zjets_6"])
    hist_E_os[data_op]["zjets"].Add(hist_E_os[data_op]["zjets_6"])
    hist_M_ss[data_op]["zjets"].Add(hist_M_ss[data_op]["zjets_7"])
    hist_M_os[data_op]["zjets"].Add(hist_M_os[data_op]["zjets_7"])
    hist_E_ss[data_op]["zjets"].Add(hist_E_ss[data_op]["zjets_7"])
    hist_E_os[data_op]["zjets"].Add(hist_E_os[data_op]["zjets_7"])
    hist_M_ss[data_op]["zjets"].Add(hist_M_ss[data_op]["zjets_8"])
    hist_M_os[data_op]["zjets"].Add(hist_M_os[data_op]["zjets_8"])
    hist_E_ss[data_op]["zjets"].Add(hist_E_ss[data_op]["zjets_8"])
    hist_E_os[data_op]["zjets"].Add(hist_E_os[data_op]["zjets_8"])
    hist_M_ss[data_op]["vv"] = hist_M_ss[data_op]["ww"]
    hist_M_os[data_op]["vv"] = hist_M_os[data_op]["ww"]
    hist_E_ss[data_op]["vv"] = hist_E_ss[data_op]["ww"]
    hist_E_os[data_op]["vv"] = hist_E_os[data_op]["ww"]
    hist_M_ss[data_op]["vv"].Add(hist_M_ss[data_op]["wz"])
    hist_M_os[data_op]["vv"].Add(hist_M_os[data_op]["wz"])
    hist_E_ss[data_op]["vv"].Add(hist_E_ss[data_op]["wz"])
    hist_E_os[data_op]["vv"].Add(hist_E_os[data_op]["wz"])
    hist_M_ss[data_op]["vv"].Add(hist_M_ss[data_op]["zz"])
    hist_M_os[data_op]["vv"].Add(hist_M_os[data_op]["zz"])
    hist_E_ss[data_op]["vv"].Add(hist_E_ss[data_op]["zz"])
    hist_E_os[data_op]["vv"].Add(hist_E_os[data_op]["zz"])

  samples = ["ttbar_sl_bottom","ttbar_sl_charm","ttbar_sl_light","ttbar_dl","ttbar_dh","zjets","vv","st","wjets"]

  histT_M_ss = {}
  histT_M_os = {}
  histT_E_ss = {}
  histT_E_os = {}
  #for s in samples:
  #     histT_M_ss[s] = hist_M_ss[datayears[0]][s]
  #     histT_M_os[s] = hist_M_os[datayears[0]][s]
  #     histT_E_ss[s] = hist_E_ss[datayears[0]][s]
  #     histT_E_os[s] = hist_E_os[datayears[0]][s]
  #     for d in datayears[1:]:
  #        histT_M_ss[s].Add(hist_M_ss[d][s])
  #        histT_M_os[s].Add(hist_M_os[d][s])
  #        histT_E_ss[s].Add(hist_E_ss[d][s])
  #        histT_E_os[s].Add(hist_E_os[d][s])
  for s in samples:
       histT_M_ss[s] = hist_M_ss[str(args.year)][s]
       histT_M_os[s] = hist_M_os[str(args.year)][s]
       histT_E_ss[s] = hist_E_ss[str(args.year)][s]
       histT_E_os[s] = hist_E_os[str(args.year)][s]

  if not args.nodata:
    #histD_M_ss = histdata_M_ss[datayears[0]]
    #histD_M_os = histdata_M_os[datayears[0]]
    #histD_E_ss = histdata_E_ss[datayears[0]]
    #histD_E_os = histdata_E_os[datayears[0]]
    #for d in datayears[1:]:
    #   histD_M_ss.Add(histdata_M_ss[d])
    #   histD_M_os.Add(histdata_M_os[d])
    #   histD_E_ss.Add(histdata_E_ss[d])
    #   histD_E_os.Add(histdata_E_os[d])
    histD_M_ss = histdata_M_ss[str(args.year)]
    histD_M_os = histdata_M_os[str(args.year)]
    histD_E_ss = histdata_E_ss[str(args.year)]
    histD_E_os = histdata_E_os[str(args.year)]

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  colors = {}
  colors["vv"] = (222,90,106)
  colors["WW_semi_nocharm"] = (246,165,42)
  colors["WW_hadronic"] = (183,2,2)
  colors["WW_leptonic"] = (153,76,0)
  colors["wjets"] = (155,152,204)
  colors["wjets_bottom"] = (255,0,127)
  colors["wjets_light"] = (255,153,204)
  colors["ttbar_sl_charm"] = (204,255,153)
  colors["ttbar_sl_light"] = (120,154,86)
  colors["ttbar_sl_bottom"] = (255,0,127)
  colors["ttbar_dl"] = (255,255,0)
  colors["ttbar_dh"] = (204,204,0)
  colors["zjets"] = (153,255,255)
  colors["wz"] = (153,153,0)
  colors["st"] = (153,51,255)
  colors["st_nocharm"] = (190,153,228)
  colors["QCD"] = (0,153,76)

  if args.stack:
    ymax_M_ss = 0
    ymin_M_ss = 0
    ymax_M_os = 0
    ymin_M_os = 0
    ymax_E_ss = 0
    ymin_E_ss = 0
    ymax_E_os = 0
    ymin_E_os = 0
    for s in samples:
      histT_M_ss[s].SetLineWidth(1)
      histT_M_ss[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_M_ss[s].GetYaxis().SetTitle("Number of events")
      histT_M_ss[s].GetXaxis().SetTitle(name)
      histT_M_os[s].SetLineWidth(1)
      histT_M_os[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_M_os[s].GetYaxis().SetTitle("Number of events")
      histT_M_os[s].GetXaxis().SetTitle(name)
      histT_E_ss[s].SetLineWidth(1)
      histT_E_ss[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_E_ss[s].GetYaxis().SetTitle("Number of events")
      histT_E_ss[s].GetXaxis().SetTitle(name)
      histT_E_os[s].SetLineWidth(1)
      histT_E_os[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_E_os[s].GetYaxis().SetTitle("Number of events")
      histT_E_os[s].GetXaxis().SetTitle(name)

      y = histT_M_ss[s].GetMaximum()
      ym = histT_M_ss[s].GetMinimum()
      if y>ymax_M_ss: ymax_M_ss=y
      if ymin_M_ss>ym: ymin_M_ss=ym
      histD_M_ss.GetMaximum()
      if y>ymax_M_ss: ymax_M_ss=y

      y = histT_M_os[s].GetMaximum()
      ym = histT_M_os[s].GetMinimum()
      if y>ymax_M_os: ymax_M_os=y
      if ymin_M_os>ym: ymin_M_os=ym
      histD_M_os.GetMaximum()
      if y>ymax_M_os: ymax_M_os=y

      y = histT_E_ss[s].GetMaximum()
      ym = histT_E_ss[s].GetMinimum()
      if y>ymax_E_ss: ymax_E_ss=y
      if ymin_E_ss>ym: ymin_E_ss=ym
      histD_E_ss.GetMaximum()
      if y>ymax_E_ss: ymax_E_ss=y

      y = histT_E_os[s].GetMaximum()
      ym = histT_E_os[s].GetMinimum()
      if y>ymax_E_os: ymax_E_os=y
      if ymin_E_os>ym: ymin_E_os=ym
      histD_E_os.GetMaximum()
      if y>ymax_E_os: ymax_E_os=y

      histT_M_ss["vv"].SetMinimum(ymin_M_ss)
      histT_M_ss["vv"].SetMaximum(5*ymax_M_ss)
      if args.linear: histT_M_ss["vv"].SetMaximum(1.3*ymax_M_ss)
      if args.linear: histT_M_ss["vv"].SetMinimum(1.3*ymin_M_ss)

      histT_M_os["vv"].SetMinimum(ymin_M_os)
      histT_M_os["vv"].SetMaximum(5*ymax_M_os)
      if args.linear: histT_M_os["vv"].SetMaximum(1.3*ymax_M_os)
      if args.linear: histT_M_os["vv"].SetMinimum(1.3*ymin_M_os)

      histT_E_ss["vv"].SetMinimum(ymin_E_ss)
      histT_E_ss["vv"].SetMaximum(5*ymax_E_ss)
      if args.linear: histT_E_ss["vv"].SetMaximum(1.3*ymax_E_ss)
      if args.linear: histT_E_ss["vv"].SetMinimum(1.3*ymin_E_ss)

      histT_E_os["vv"].SetMinimum(ymin_E_os)
      histT_E_os["vv"].SetMaximum(5*ymax_E_os)
      if args.linear: histT_E_os["vv"].SetMaximum(1.3*ymax_E_os)
      if args.linear: histT_E_os["vv"].SetMinimum(1.3*ymin_E_os)

    ## Stack creation
    samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]

    if args.ratio and not args.nodata: upper_pad.cd()
    stack_M_ss = ROOT.THStack()
    stack_M_os = ROOT.THStack()
    stack_E_ss = ROOT.THStack()
    stack_E_os = ROOT.THStack()
    for s in samples:
      stack_M_ss.Add(histT_M_ss[s])
      stack_M_os.Add(histT_M_os[s])
      stack_E_ss.Add(histT_E_ss[s])
      stack_E_os.Add(histT_E_os[s])

    y_M_ss = stack_M_ss.GetMaximum()
    if y_M_ss>ymax_M_ss: ymax_M_ss=y_M_ss
    stack_M_ss.SetMinimum(1.)
    stack_M_ss.SetMaximum(5*ymax_M_ss)
    if args.linear: stack_M_ss.SetMaximum(1.3*ymax_M_ss)
    y_M_os = stack_M_os.GetMaximum()
    if y_M_os>ymax_M_os: ymax_M_os=y_M_os
    stack_M_os.SetMinimum(1.)
    stack_M_os.SetMaximum(5*ymax_M_os)
    if args.linear: stack_M_os.SetMaximum(1.3*ymax_M_os)
    y_E_ss = stack_E_ss.GetMaximum()
    if y_E_ss>ymax_E_ss: ymax_E_ss=y_E_ss
    stack_E_ss.SetMinimum(1.)
    stack_E_ss.SetMaximum(5*ymax_E_ss)
    if args.linear: stack_E_ss.SetMaximum(1.3*ymax_E_ss)
    y_E_os = stack_E_os.GetMaximum()
    if y_E_os>ymax_E_os: ymax_E_os=y_E_os
    stack_E_os.SetMinimum(1.)
    stack_E_os.SetMaximum(5*ymax_E_os)
    if args.linear: stack_E_os.SetMaximum(1.3*ymax_E_os)

    stack_M_ss.Draw("HIST")
    histD_M_ss.SetMarkerStyle(20)
    histD_M_ss.SetMarkerSize(0.3)
    histD_M_ss.SetLineWidth(1)
    histD_M_ss.SetLineColor(ROOT.kBlack)
    histD_M_ss.Draw("E SAME")

    if args.ratio and not args.nodata:
      lower_pad.cd()
      ratio = histD_M_ss.Clone("ratio")
      ratio.SetLineColor(kBlack)
      ratio.SetMarkerStyle(21)
      ratio.SetTitle("")
      ratio.SetMinimum(0)
      ratio.SetMaximum(2)
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
      hTotal = histT_M_ss["vv"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histT_M_ss[s])
      ratio.Divide(hTotal)
      ratio.Draw("ep")

    ## Legends
    if args.ratio and not args.nodata: upper_pad.cd()
    leg = TLegend(0.7,0.6,0.89,0.89)
    leg.SetBorderSize(1)
    leg.AddEntry(histT_M_ss["vv"],"VV","f")
    leg.AddEntry(histT_M_ss["ttbar_dl"],"Dileptonic top antitop","f")
    leg.AddEntry(histT_M_ss["zjets"],"Z plus jets","f")
    leg.AddEntry(histT_M_ss["wjets"],"W plus jets","f")
    leg.AddEntry(histT_M_ss["st"],"Single top","f")
    leg.AddEntry(histT_M_ss["ttbar_dh"],"Hadronic top antitop","f")
    leg.AddEntry(histT_M_ss["ttbar_sl_bottom"],"Top antitop, bottom","f")
    leg.AddEntry(histT_M_ss["ttbar_sl_light"],"Top antitop, light","f")
    leg.AddEntry(histT_M_ss["ttbar_sl_charm"],"Top antitop, charm","f")
    if args.stack and not args.nodata: leg.AddEntry(histD_M_ss, "Data" ,"lep")
    leg.Draw()
    term= "totalHT_"+str(args.type)+"_M_ss"+str(args.year)
    if args.ratio: 
      notation = "_ratio_"
      if args.linear:
        notation = "_linratio_"
    else: 
      notation = "_normed_"

    if args.png: c1.Print(plotdir+term+notation+ name + ".png")
    else: c1.Print(plotdir+term+notation + name + ".pdf")

    stack_M_os.Draw("HIST")
    histD_M_os.SetMarkerStyle(20)
    histD_M_os.SetMarkerSize(0.3)
    histD_M_os.SetLineWidth(1)
    histD_M_os.SetLineColor(ROOT.kBlack)
    histD_M_os.Draw("E SAME")

    if args.ratio and not args.nodata:
      lower_pad.cd()
      ratio = histD_M_os.Clone("ratio")
      ratio.SetLineColor(kBlack)
      ratio.SetMarkerStyle(21)
      ratio.SetTitle("")
      ratio.SetMinimum(0)
      ratio.SetMaximum(2)
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
      hTotal = histT_M_os["vv"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histT_M_os[s])
      ratio.Divide(hTotal)
      ratio.Draw("ep")

    ## Legends
    if args.ratio and not args.nodata: upper_pad.cd()
    leg = TLegend(0.7,0.6,0.89,0.89)
    leg.SetBorderSize(1)
    leg.AddEntry(histT_M_os["vv"],"VV","f")
    leg.AddEntry(histT_M_os["ttbar_dl"],"Dileptonic top antitop","f")
    leg.AddEntry(histT_M_os["zjets"],"Z plus jets","f")
    leg.AddEntry(histT_M_os["wjets"],"W plus jets","f")
    leg.AddEntry(histT_M_os["st"],"Single top","f")
    leg.AddEntry(histT_M_os["ttbar_dh"],"Hadronic top antitop","f")
    leg.AddEntry(histT_M_os["ttbar_sl_bottom"],"Top antitop, bottom","f")
    leg.AddEntry(histT_M_os["ttbar_sl_light"],"Top antitop, light","f")
    leg.AddEntry(histT_M_os["ttbar_sl_charm"],"Top antitop, charm","f")
    if args.stack and not args.nodata: leg.AddEntry(histD_M_os, "Data" ,"lep")
    leg.Draw()
    term= "totalHT_"+str(args.type)+"_M_os"+str(args.year)
    if args.ratio:
      notation = "_ratio_"
      if args.linear:
        notation = "_linratio_"
    else:
      notation = "_normed_"

    if args.png: c1.Print(plotdir+term+notation+ name + ".png")
    else: c1.Print(plotdir+term+notation + name + ".pdf")

    stack_E_ss.Draw("HIST")
    histD_E_ss.SetMarkerStyle(20)
    histD_E_ss.SetMarkerSize(0.3)
    histD_E_ss.SetLineWidth(1)
    histD_E_ss.SetLineColor(ROOT.kBlack)
    histD_E_ss.Draw("E SAME")

    if args.ratio and not args.nodata:
      lower_pad.cd()
      ratio = histD_E_ss.Clone("ratio")
      ratio.SetLineColor(kBlack)
      ratio.SetMarkerStyle(21)
      ratio.SetTitle("")
      ratio.SetMinimum(0)
      ratio.SetMaximum(2)
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
      hTotal = histT_E_ss["vv"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histT_E_ss[s])
      ratio.Divide(hTotal)
      ratio.Draw("ep")

    ## Legends
    if args.ratio and not args.nodata: upper_pad.cd()
    leg = TLegend(0.7,0.6,0.89,0.89)
    leg.SetBorderSize(1)
    leg.AddEntry(histT_E_ss["vv"],"VV","f")
    leg.AddEntry(histT_E_ss["ttbar_dl"],"Dileptonic top antitop","f")
    leg.AddEntry(histT_E_ss["zjets"],"Z plus jets","f")
    leg.AddEntry(histT_E_ss["wjets"],"W plus jets","f")
    leg.AddEntry(histT_E_ss["st"],"Single top","f")
    leg.AddEntry(histT_E_ss["ttbar_dh"],"Hadronic top antitop","f")
    leg.AddEntry(histT_E_ss["ttbar_sl_bottom"],"Top antitop, bottom","f")
    leg.AddEntry(histT_E_ss["ttbar_sl_light"],"Top antitop, light","f")
    leg.AddEntry(histT_E_ss["ttbar_sl_charm"],"Top antitop, charm","f")
    if args.stack and not args.nodata: leg.AddEntry(histD_E_ss, "Data" ,"lep")
    leg.Draw()
    term= "totalHT_"+str(args.type)+"_E_ss"+str(args.year)
    if args.ratio:
      notation = "_ratio_"
      if args.linear:
        notation = "_linratio_"
    else:
      notation = "_normed_"

    if args.png: c1.Print(plotdir+term+notation+ name + ".png")
    else: c1.Print(plotdir+term+notation + name + ".pdf")

    stack_E_os.Draw("HIST")
    histD_E_os.SetMarkerStyle(20)
    histD_E_os.SetMarkerSize(0.3)
    histD_E_os.SetLineWidth(1)
    histD_E_os.SetLineColor(ROOT.kBlack)
    histD_E_os.Draw("E SAME")
    if args.ratio and not args.nodata:
      lower_pad.cd()
      ratio = histD_E_os.Clone("ratio")
      ratio.SetLineColor(kBlack)
      ratio.SetMarkerStyle(21)
      ratio.SetTitle("")
      ratio.SetMinimum(0)
      ratio.SetMaximum(2)
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
      hTotal = histT_E_os["vv"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histT_E_os[s])
      ratio.Divide(hTotal)
      ratio.Draw("ep")

    ## Legends
    if args.ratio and not args.nodata: upper_pad.cd()
    leg = TLegend(0.7,0.6,0.89,0.89)
    leg.SetBorderSize(1)
    leg.AddEntry(histT_E_os["vv"],"VV","f")
    leg.AddEntry(histT_E_os["ttbar_dl"],"Dileptonic top antitop","f")
    leg.AddEntry(histT_E_os["zjets"],"Z plus jets","f")
    leg.AddEntry(histT_E_os["wjets"],"W plus jets","f")
    leg.AddEntry(histT_E_os["st"],"Single top","f")
    leg.AddEntry(histT_E_os["ttbar_dh"],"Hadronic top antitop","f")
    leg.AddEntry(histT_E_os["ttbar_sl_bottom"],"Top antitop, bottom","f")
    leg.AddEntry(histT_E_os["ttbar_sl_light"],"Top antitop, light","f")
    leg.AddEntry(histT_E_os["ttbar_sl_charm"],"Top antitop, charm","f")
    if args.stack and not args.nodata: leg.AddEntry(histD_E_os, "Data" ,"lep")
    leg.Draw()
    term= "totalHT_"+str(args.type)+"_E_os"+str(args.year)
    if args.ratio:
      notation = "_ratio_"
      if args.linear:
        notation = "_linratio_"
    else:
      notation = "_normed_"

    if args.png: c1.Print(plotdir+term+notation+ name + ".png")
    else: c1.Print(plotdir+term+notation + name + ".pdf")



for s in samplesHT:
        for data_op in datayears:
                        histFile[data_op][s].Close()

for data_op in datayears:
	histFileD[data_op]["M"].Close()
	histFileD[data_op]["E"].Close()


