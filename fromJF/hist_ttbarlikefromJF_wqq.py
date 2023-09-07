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
parser.add_argument("--norm", action="store_true", default=False,
                    help="normalising test")
parser.add_argument("--qcd", action="store_true", default=False,
                    help="include qcd samples")
parser.add_argument("--nodata", action="store_true", default=False,
                    help="Do not plot data")

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

## no pile up id and jet id, lepton50, smearing
norm_factorM = 0.9583
norm_factorE = 0.9669

## no pile up id and jet id, lepton50, smearing, ptjets30
norm_factorM = 0.9513
norm_factorE = 0.9551

## no pile up id and jet id, btagM, smearing, ptjets30
norm_factorM = 0.9515
norm_factorE = 0.9377

## no pile up id and jet id, btagM, smearing
norm_factorM = 0.9608
norm_factorE = 0.9469

## Open hists files

filePath = "/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/wqq/btagMM/"

if args.ssos: ssos_add = "SSOS"
else: ssos_add = "" 

term = "histstt_wqq_fromJF_"

datayears = ["2016","2016B","2017","2018"]
#datayears = ["2018","2016","2016B"]

samplesHT = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_light","ttbar_sl_bottom","ttbar_dl","ttbar_dh","zz","wz",
        "st_1","st_2","st_3","st_4", "ttbar_sl_else", "ttbar_sl_charmgluon","ttbar_sl_bottomgluon"]

## Adding QCD

histFile = {}

for data_op in datayears:
        if data_op == "2016B": 
            end_term = "2.root"
        else:
            end_term = ".root"
        ## mc files
        histFile[data_op] = {}
        for s in samplesHT:
                if s[0:8] == "ttbar_sl" and isfile(filePath + term+s[0:8]+data_op+s[8:]+end_term):
                        histFile[data_op][s] = TFile.Open(filePath + term+s[0:8]+data_op+s[8:]+end_term,"READ")
                elif isfile(filePath + term+s+data_op+end_term):
                        histFile[data_op][s] = TFile.Open(filePath + term+s+data_op+end_term,"READ")
        #print(histFile[data_op].keys())

histFileD = {}

for data_op in datayears:
	histFileD[data_op] = {}
	# data files
	histFileD[data_op]["M"] = TFile.Open(filePath + term+data_op+"M.root","READ")
	histFileD[data_op]["E"] = TFile.Open(filePath + term+data_op+"E.root","READ")

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

histNames = ["InvM_2jets", "nJetGood", "jet_1_pt", "jet_1_nmu", "jet_1_eta", "jet_2_pt", "jet_2_eta", "jet_2_mass", "jet_2_qgl", "jet_2_nmu",
   "lepton_pt", "lepton_eta", "InvM_bot_closer", "InvM_bot_farther", "deltaR_jet1_jet2", "deltaphi_jet1_jet2", "deltaeta_jet1_jet2", 
   "MET_pt", "tracks_jet1", "tracks_jet2", "EMN_jet1", "EMtotal_jet1", "EMC_jet1", "pT_sum", "pT_product", "deltaR_lep_2jets", 
   "deltaphi_MET_2jets", "deltaphi_lephad", "eta_2jets", "transverse_mass", "pt_2jets", "pt_Wlep", "deltaR_lep_jet1", 
   "deltaR_lep_jet2", "deltaPhi_lep_jet1", "deltaPhi_lep_jet2", "deltaEta_lep_jet1", "deltaEta_lep_jet2", "jet_1_btag", "jet_2_btag", 
   "jet_bot1_btag", "jet_bot2_btag", "jet_bot1_btag_thick", "jet_bot2_btag_thick", "pT_proy", "pT_sum_2J","deltaphi_MET_jets_1","deltaphi_MET_jets_2", 
   "lepton_pt_detail", "jet_1_qgl", "jet_bot1_pt", "jet_bot1_eta", "jet_bot2_pt", "jet_bot2_eta", "lepton_eta_thick", "MET_sig", "deltaphi_MET_lep", "MET_my_sig",
   "jet_1_btag_thick", "jet_2_btag_thick", "jet_bot1_btagnumber", "jet_bot2_btagnumber", "jet_1_btagnumber", "jet_2_btagnumber"]

#histNames = ["jet_bot1_btag_thick", "jet_bot2_btag_thick","jet_1_btag_thick", "jet_2_btag_thick", "jet_1_nmu", "jet_2_nmu",
#     "jet_bot1_btag", "jet_bot2_btag","jet_1_btag", "jet_2_btag", "jet_bot1_btagnumber", "jet_bot2_btagnumber","jet_1_btagnumber", "jet_2_btagnumber"]

if args.nodata: histNames = ["jet_1_flavourP","jet_2_flavourP","jet_bot1_flavourP","jet_bot2_flavourP","btag_sf","lep_id_sf","lep_trig_sf","puWeight","PUjetID_SF","topweight"]


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
        "st_1","st_2","st_3","st_4","ttbar_sl_else","ttbar_sl_charmgluon","ttbar_sl_bottomgluon"]

  ## HISTS
  samples_foryear = {}
  hist_M = {}
  hist_E = {}
  histdata_M = {}
  histdata_E = {}
  for data_op in datayears:
    samples_foryear[data_op] = [s for s in samples if s in histFile[data_op].keys()]
    hist_M[data_op] = {}
    hist_E[data_op] = {}
    for s in samples_foryear[data_op]:
      hist_M[data_op][s] = histFile[data_op][s].Get(name+"_M")
      hist_E[data_op][s] = histFile[data_op][s].Get(name+"_E")
    if not args.nodata: histdata_M[data_op] = histFileD[data_op]["M"].Get(name+"_M")
    if not args.nodata: histdata_E[data_op] = histFileD[data_op]["E"].Get(name+"_E")

  samples_st = {}
  samples_wjets = {}
  samples_zjets = {}
  ## Scaling to lumi
  print(name) 
  for data_op in datayears:
    lumi_data = lumi_d[data_op]
    for s in samples_foryear[data_op]:
      #print(s)
      #print(data_op)
      hist_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
      if args.norm: hist_M[data_op][s].Scale(norm_factorM)
      if args.norm: hist_E[data_op][s].Scale(norm_factorE)
      ## Fixing single top
    #print(samples_foryear[data_op]) 
    hist_M[data_op]["st"] = hist_M[data_op]["st_1"]
    hist_E[data_op]["st"] = hist_E[data_op]["st_1"]
    hist_M[data_op]["st"].Add(hist_M[data_op]["st_2"])
    hist_E[data_op]["st"].Add(hist_E[data_op]["st_2"])
    hist_M[data_op]["st"].Add(hist_M[data_op]["st_3"])
    hist_E[data_op]["st"].Add(hist_E[data_op]["st_3"])
    hist_M[data_op]["st"].Add(hist_M[data_op]["st_4"])
    hist_E[data_op]["st"].Add(hist_E[data_op]["st_4"])
    hist_M[data_op]["wjets"] = hist_M[data_op]["wjets_1"]
    hist_E[data_op]["wjets"] = hist_E[data_op]["wjets_1"]
    hist_M[data_op]["wjets"].Add(hist_M[data_op]["wjets_2"])
    hist_E[data_op]["wjets"].Add(hist_E[data_op]["wjets_2"])
    hist_M[data_op]["wjets"].Add(hist_M[data_op]["wjets_3"])
    hist_E[data_op]["wjets"].Add(hist_E[data_op]["wjets_3"])
    hist_M[data_op]["wjets"].Add(hist_M[data_op]["wjets_4"])
    hist_E[data_op]["wjets"].Add(hist_E[data_op]["wjets_4"])
    hist_M[data_op]["wjets"].Add(hist_M[data_op]["wjets_5"])
    hist_E[data_op]["wjets"].Add(hist_E[data_op]["wjets_5"])
    hist_M[data_op]["wjets"].Add(hist_M[data_op]["wjets_6"])
    hist_E[data_op]["wjets"].Add(hist_E[data_op]["wjets_6"])
    hist_M[data_op]["wjets"].Add(hist_M[data_op]["wjets_7"])
    hist_E[data_op]["wjets"].Add(hist_E[data_op]["wjets_7"])
    hist_M[data_op]["wjets"].Add(hist_M[data_op]["wjets_8"])
    hist_E[data_op]["wjets"].Add(hist_E[data_op]["wjets_8"])
    hist_M[data_op]["zjets"] = hist_M[data_op]["zjets_1"]
    hist_E[data_op]["zjets"] = hist_E[data_op]["zjets_1"]
    hist_M[data_op]["zjets"].Add(hist_M[data_op]["zjets_2"])
    hist_E[data_op]["zjets"].Add(hist_E[data_op]["zjets_2"])
    hist_M[data_op]["zjets"].Add(hist_M[data_op]["zjets_3"])
    hist_E[data_op]["zjets"].Add(hist_E[data_op]["zjets_3"])
    hist_M[data_op]["zjets"].Add(hist_M[data_op]["zjets_4"])
    hist_E[data_op]["zjets"].Add(hist_E[data_op]["zjets_4"])
    hist_M[data_op]["zjets"].Add(hist_M[data_op]["zjets_5"])
    hist_E[data_op]["zjets"].Add(hist_E[data_op]["zjets_5"])
    hist_M[data_op]["zjets"].Add(hist_M[data_op]["zjets_6"])
    hist_E[data_op]["zjets"].Add(hist_E[data_op]["zjets_6"])
    hist_M[data_op]["zjets"].Add(hist_M[data_op]["zjets_7"])
    hist_E[data_op]["zjets"].Add(hist_E[data_op]["zjets_7"])
    hist_M[data_op]["zjets"].Add(hist_M[data_op]["zjets_8"])
    hist_E[data_op]["zjets"].Add(hist_E[data_op]["zjets_8"])
    hist_M[data_op]["vv"] = hist_M[data_op]["ww"]
    hist_E[data_op]["vv"] = hist_E[data_op]["ww"]
    hist_M[data_op]["vv"].Add(hist_M[data_op]["wz"])
    hist_E[data_op]["vv"].Add(hist_E[data_op]["wz"])
    hist_M[data_op]["vv"].Add(hist_M[data_op]["zz"])
    hist_E[data_op]["vv"].Add(hist_E[data_op]["zz"])

  samples = ["ttbar_sl_bottom","ttbar_sl_charm","ttbar_sl_else","ttbar_sl_light","ttbar_dl","ttbar_dh","zjets","vv","st","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon"]

  histT_M = {}
  histT_E = {}
  for s in samples:
       histT_M[s] = hist_M[datayears[0]][s]
       histT_E[s] = hist_E[datayears[0]][s]
       for d in datayears[1:]:
          histT_M[s].Add(hist_M[d][s])
          histT_E[s].Add(hist_E[d][s])

  if not args.nodata:
    histD_M = histdata_M[datayears[0]]
    histD_E = histdata_E[datayears[0]]
    for d in datayears[1:]:
       histD_M.Add(histdata_M[d])
       histD_E.Add(histdata_E[d])

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
      histT_M[s].SetLineWidth(1)
      histT_M[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_M[s].GetYaxis().SetTitle("Number of events")
      histT_M[s].GetXaxis().SetTitle(name)
      histT_E[s].SetLineWidth(1)
      histT_E[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_E[s].GetYaxis().SetTitle("Number of events")
      histT_E[s].GetXaxis().SetTitle(name)

      y = histT_M[s].GetMaximum()
      ym = histT_M[s].GetMinimum()
      if y>ymax_M: ymax_M=y
      if ymin_M>ym: ymin_M=ym
      y = histD_M.GetMaximum()
      if y>ymax_M: ymax_M=y

      y = histT_E[s].GetMaximum()
      ym = histT_E[s].GetMinimum()
      if y>ymax_E: ymax_E=y
      if ymin_E>ym: ymin_E=ym
      y = histD_E.GetMaximum()
      if y>ymax_E: ymax_E=y

      #histT_M["vv"].SetMinimum(1000.)
      #histT_M["vv"].SetMaximum(3*ymax_M)
      if args.linear: histT_M["vv"].SetMaximum(1.3*ymax_M)
      if args.linear: histT_M["vv"].SetMinimum(1.3*ymin_M)

      #histT_E["vv"].SetMinimum(1000.)
      #histT_E["vv"].SetMaximum(3*ymax_E)
      if args.linear: histT_E["vv"].SetMaximum(1.3*ymax_E)
      if args.linear: histT_E["vv"].SetMinimum(1.3*ymin_E)

    ## Stack creation
    samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon","ttbar_sl_else","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]

    if args.ratio and not args.nodata: upper_pad.cd()
    stack_M = ROOT.THStack()
    stack_E = ROOT.THStack()
    for s in samples:
      stack_M.Add(histT_M[s])
      stack_E.Add(histT_E[s])

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
    histD_M.SetMarkerStyle(20)
    histD_M.SetMarkerSize(0.3)
    histD_M.SetLineWidth(1)
    histD_M.SetLineColor(ROOT.kBlack)
    histD_M.Draw("E SAME")

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
      hTotal = histT_M["vv"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histT_M[s])
      ratio.Divide(hTotal)
      ratio.Draw("ep")

    if name == "InvM_2jets":
      print("Integral of M data is "+str(histD_M.Integral()))
      print("Integral of M MC is "+str(hTotal.Integral()))
      print("Ratio is "+str(histD_M.Integral()/hTotal.Integral()))

    ## Legends
    if args.ratio and not args.nodata: upper_pad.cd()
    leg = TLegend(0.7,0.6,0.89,0.89)
    leg.SetBorderSize(1)
    leg.AddEntry(histT_M["vv"],"VV","f")
    leg.AddEntry(histT_M["ttbar_dl"],"Dileptonic top antitop","f")
    leg.AddEntry(histT_M["zjets"],"Z plus jets","f")
    leg.AddEntry(histT_M["wjets"],"W plus jets","f")
    leg.AddEntry(histT_M["st"],"Single top","f")
    leg.AddEntry(histT_M["ttbar_dh"],"Hadronic top antitop","f")
    leg.AddEntry(histT_M["ttbar_sl_bottom"],"Top antitop, bottom","f")
    leg.AddEntry(histT_M["ttbar_sl_light"],"Top antitop, light","f")
    leg.AddEntry(histT_M["ttbar_sl_charm"],"Top antitop, charm","f")
    leg.AddEntry(histT_M["ttbar_sl_bottomgluon"],"Top antitop, bottom gluon","f")
    leg.AddEntry(histT_M["ttbar_sl_charmgluon"],"Top antitop, charm gluon","f")
    leg.AddEntry(histT_M["ttbar_sl_else"],"Top antitop, gluon gluon","f")
    if args.stack and not args.nodata: leg.AddEntry(histD_M, "Data" ,"lep")
    leg.Draw()
    term= "totalHT_wqq_M"
    if args.ratio: 
      notation = "_ratio_"
      if args.linear:
        notation = "_linratio_"
    else: 
      notation = "_normed_"

    if args.png: c1.Print(plotdir+term+notation+ name + ".png")
    else: c1.Print(plotdir+term+notation + name + ".root")

    stack_E.Draw("HIST")
    histD_E.SetMarkerStyle(20)
    histD_E.SetMarkerSize(0.3)
    histD_E.SetLineWidth(1)
    histD_E.SetLineColor(ROOT.kBlack)
    histD_E.Draw("E SAME")
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
      hTotal = histT_E["vv"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histT_E[s])
      ratio.Divide(hTotal)
      ratio.Draw("ep")

    if name == "InvM_2jets":
      print("Integral of E data is "+str(histD_E.Integral()))
      print("Integral of E MC is "+str(hTotal.Integral()))
      print("Ratio is "+str(histD_E.Integral()/hTotal.Integral()))

    ## Legends
    if args.ratio and not args.nodata: upper_pad.cd()
    leg = TLegend(0.7,0.6,0.89,0.89)
    leg.SetBorderSize(1)
    leg.AddEntry(histT_E["vv"],"VV","f")
    leg.AddEntry(histT_E["ttbar_dl"],"Dileptonic top antitop","f")
    leg.AddEntry(histT_E["zjets"],"Z plus jets","f")
    leg.AddEntry(histT_E["wjets"],"W plus jets","f")
    leg.AddEntry(histT_E["st"],"Single top","f")
    leg.AddEntry(histT_E["ttbar_dh"],"Hadronic top antitop","f")
    leg.AddEntry(histT_E["ttbar_sl_bottom"],"Top antitop, bottom","f")
    leg.AddEntry(histT_E["ttbar_sl_light"],"Top antitop, light","f")
    leg.AddEntry(histT_E["ttbar_sl_charm"],"Top antitop, charm","f")
    leg.AddEntry(histT_E["ttbar_sl_bottomgluon"],"Top antitop, bottom gluon","f")
    leg.AddEntry(histT_E["ttbar_sl_charmgluon"],"Top antitop, charm gluon","f")
    leg.AddEntry(histT_E["ttbar_sl_else"],"Top antitop, gluon gluon","f")
    if args.stack and not args.nodata: leg.AddEntry(histD_E, "Data" ,"lep")
    leg.Draw()
    term= "totalHT_wqq_E"
    if args.ratio:
      notation = "_ratio_"
      if args.linear:
        notation = "_linratio_"
    else:
      notation = "_normed_"

    if args.png: c1.Print(plotdir+term+notation+ name + ".png")
    else: c1.Print(plotdir+term+notation + name + ".root")



for s in samplesHT:
        for data_op in datayears:
                        histFile[data_op][s].Close()

for data_op in datayears:
	histFileD[data_op]["M"].Close()
	histFileD[data_op]["E"].Close()


