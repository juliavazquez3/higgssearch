import sys
import ROOT
import os
from ROOT import *
import json
import argparse
import numpy as np
from os.path import isfile, join, isdir
import CMS_lumi, tdrstyle

#if not sys.flags.interactive: ROOT.EnableImplicitMT()

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

#set the tdr style
#tdrstyle.setTDRStyle()
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"

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
parser.add_argument("--nodata", action="store_true", default=False,
                    help="Do not plot data")
parser.add_argument("--var", type=string, default="",
                    help="Variable to plot")

# Use like:
# python higgssearch/fromJF/hist_fromJFwqq_syst_total.py --stack --ratio --png --norm --year="all" --channel="btagMM_chitest_sl"

args = parser.parse_args()

plotdir = '/nfs/cms/vazqueze/higgssearch/plotspng/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

##### some settings ######

iPeriod = 0
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12

c_rat = 1.5
c_rat2 = 0.5
nrebin = 2; nrat = 11;
mrk_size = 5; titY_off = 0.3; tit_size = 0.15; titX_off = 0.6;
if not args.nodata:
   leg1 = 0.4; leg2 = 0.7; leg3 = 0.96; leg4 = 0.92;
else:
   leg1 = 0.4; leg2 = 0.76; leg3 = 0.92; leg4 = 0.9;

observable_names = ["InvM_2jets","jet_1_pt", "jet_2_pt", "jet_1_eta", "jet_1_nmu", "jet_2_eta", "jet_2_mass", "jet_2_qgl","jet_2_nmu","jet_1_qgl",
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
   "jet_bot1_tracks","jet_bot2_tracks","tau_discr_jet1","tau_discr_jet2","tau_discr_jetbot1","tau_discr_jetbot2","muon_jet_pt","muon_jet_z","muon_jet_eta",
   "muon_jet_pt_rel","muon_jet_iso",
   "muon_jet_iso_log","muon_jet_z_short","InvM3_good_short","muon_jet_sigr","muon_jet_sigxy","muon_jet_sigdz","muon_jet_r","deltaR_jet1_muon",
   "muon_jet_z2_v2","muon_jet_z3","muon_jet_iso_abs","jet_muon_pt","jet_muon_eta","jet_muon_btag","jet_muon_btagnumber",
   "jet_notmuon_pt","jet_notmuon_eta","jet_notmuon_btag","jet_notmuon_btagnumber","jet_muon_nmu","jet_notmuon_nmu"]

cosmetic_names = {}

for name in observable_names:
    cosmetic_names[name] = name

cosmetic_names["InvM_2jets_short"] = "m_{jj} [GeV]";cosmetic_names["nJetGood"] = "Number of jets";
cosmetic_names["InvM_2jets"] = "m_{jj} [GeV]"; cosmetic_names["InvM3_good"] = "m_{jjb} [GeV]"; cosmetic_names["InvMl_good"] = "m_{lb} [GeV]";
cosmetic_names["jet_1_pt"] = "p_{T}^{jet} [GeV]";cosmetic_names["jet_2_pt"] = "p_{T}^{jet} [GeV]";cosmetic_names["jet_bot1_pt"] = "p_{T}^{jet} [GeV]";cosmetic_names["jet_bot2_pt"] = "p_{T}^{jet} [GeV]";
cosmetic_names["jet_1_eta"] = "\eta^{ jet}";cosmetic_names["jet_2_eta"] = "\eta^{ jet}";cosmetic_names["jet_bot1_eta"] = "\eta^{ jet}";cosmetic_names["jet_bot2_eta"] = "\eta^{ jet}";
cosmetic_names["jet_max_pt"] = "p_{T}^{jet} [GeV]";cosmetic_names["jet_min_pt"] = "p_{T}^{jet} [GeV]";
cosmetic_names["jet_max_eta"] = "\eta^{ jet}";cosmetic_names["jet_min_eta"] = "\eta^{ jet}";
cosmetic_names["jet_1_flavourP"] = "Parton flavour";cosmetic_names["jet_2_flavourP"] = "Parton flavour";
cosmetic_names["jet_bot1_flavourP"] = "Parton flavour";cosmetic_names["jet_bot2_flavourP"] = "Parton flavour";
cosmetic_names["lepton_pt"] = "p_{T}^{l} [GeV]";cosmetic_names["lepton_eta"] = "\eta^{l}"; cosmetic_names["lepton_eta_thick"] = "\eta^{l}";
cosmetic_names["chi2_test_good"] = "\chi^{2}";cosmetic_names["transverse_mass"] = "m_{T} [GeV]"; cosmetic_names["MET_pt_aux"] = "MET p_{T} [GeV]";
cosmetic_names["jet_bot1_btag"] = "DeepJet (b)"; cosmetic_names["jet_bot2_btag"] = "DeepJet (b)";
cosmetic_names["jet_bot1_btagnumber"] = "DeepJet WP (b)"; cosmetic_names["jet_bot2_btagnumber"] = "DeepJet WP (b)";
cosmetic_names["jet_1_btagnumber"] = "DeepJet WP (b)"; cosmetic_names["jet_2_btagnumber"] = "DeepJet WP (b)";
cosmetic_names["deltaR_jet1_jet2"] = "\Delta R (jet 1, jet 2)"; cosmetic_names["deltaphi_jet1_jet2"] = "\Delta \phi (jet 1, jet 2)";
cosmetic_names["deltaR_jet1_jet2"] = "\Delta R (jet 1, jet 2)"; cosmetic_names["deltaphi_MET_lep"] = "\Delta \phi (l, MET)";
cosmetic_names["pT_Wlep"] = "W^{l} p_{T} [GeV]";
cosmetic_names["muon_jet_pt"] = "#mu p_{T} [GeV]"; cosmetic_names["muon_jet_eta"] = "#mu #eta"; cosmetic_names["muon_jet_iso_abs"] = "I_{PF} [GeV]";
cosmetic_names["muon_jet_z_short"] = "p_{T}^{#mu}/p_{T}^{jet}";cosmetic_names["muon_jet_z"] = "p_{T}^{#mu}/p_{T}^{jet}";
cosmetic_names["jet_1_eta_thick"] = "\eta^{ jet}";cosmetic_names["jet_2_eta_thick"] = "\eta^{ jet}";
cosmetic_names["jet_bot1_eta_thick"] = "\eta^{ jet}";cosmetic_names["jet_bot2_eta_thick"] = "\eta^{ jet}";

#######################################################
########### Start of plot creation ####################
#######################################################

hist_content = {}
hist_content["prefit"] = {}
hist_content["postfit"] = {}
hist_content["prefit"]["ch1"] = {}
hist_content["postfit"]["ch1"] = {}
hist_content["prefit"]["ch2"] = {}
hist_content["postfit"]["ch2"] = {}
hist_content["prefit"]["ch3"] = {}
hist_content["postfit"]["ch3"] = {}
hist_content["prefit"]["ch4"] = {}
hist_content["postfit"]["ch4"] = {}

## Open hists files
#filePath = "/nfs/cms/vazqueze/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit/data/tutorials/longexercise/makePrePostplots/output_histograms.root"
filePath = "/nfs/cms/vazqueze/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit/data/tutorials/longexercise/prefit_plots/output_histograms"+str(args.var)+".root"
## mc files
file_aux = TFile.Open(filePath,"READ")

samples = ["ttWcq","ttWuq","stWcq","stWuq","bck","stnW","vjets","mcss","ttdl","data","TotalProcs"]
samplesMC = ["ttWcq","ttWuq","stWcq","stWuq","bck","stnW","vjets","mcss","ttdl"]
for s in samples:
    if s=="data": 
       termS = "data_obs"
    else: 
       termS = str(s)
    hist_content["prefit"]["ch1"][s] = file_aux.Get("ch1_ch1_prefit/"+termS)
    hist_content["prefit"]["ch2"][s] = file_aux.Get("ch1_ch2_prefit/"+termS)
    hist_content["prefit"]["ch3"][s] = file_aux.Get("ch1_ch3_prefit/"+termS)
    hist_content["prefit"]["ch4"][s] = file_aux.Get("ch1_ch4_prefit/"+termS)
    hist_content["postfit"]["ch1"][s] = file_aux.Get("ch1_ch1_postfit/"+termS)
    hist_content["postfit"]["ch2"][s] = file_aux.Get("ch1_ch2_postfit/"+termS)
    hist_content["postfit"]["ch3"][s] = file_aux.Get("ch1_ch3_postfit/"+termS)
    hist_content["postfit"]["ch4"][s] = file_aux.Get("ch1_ch4_postfit/"+termS)

#### Case of data being TGraph instead of THF1D
#nPoints = hist_content["postfit"]["ch4"]["data"].GetN()
#for step in ["prefit","postfit"]:
#    for channel in ["ch1","ch2","ch3","ch4"]:
#        hist_content[step][channel]["data_aux"] = hist_content[step][channel]["bck"].Clone("data_aux")
#        hist_content[step][channel]["data_aux"].Scale(0.0)
#        for i in range(nPoints):
#            ygr = hist_content[step][channel]["data"].GetPointY(i)
#            hist_content[step][channel]["data_aux"].SetBinContent(i+1,ygr)


if args.png: c1 = TCanvas("c1","",1200,1000)
else: c1 = TCanvas("c1","",600,400)

## In case of ratio plot
upper_pad = ROOT.TPad("upper_pad", "", 0, 0.35, 1, 1)
lower_pad = ROOT.TPad("lower_pad", "", 0, 0, 1, 0.35)
for p in [upper_pad, lower_pad]:
        p.SetLeftMargin(0.095)
        p.SetRightMargin(0.04)
        p.SetTickx(False)
        p.SetTicky(False)
upper_pad.SetBottomMargin(0)
upper_pad.SetTopMargin(0.075)
lower_pad.SetTopMargin(0)
lower_pad.SetBottomMargin(0.25)

if not args.linear: upper_pad.SetLogy()
upper_pad.Draw()
lower_pad.Draw()

gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

colors = {}
colors["ttdl"] = (222,90,106)
colors["vjets"] = (155,152,204)
colors["ttWcq"] = (204,255,153)
colors["ttWuq"] = (120,154,86)
colors["bck"] = (255,180,85)
colors["mcss"] = (204,204,0)
colors["zjets"] = (113,209,223)
colors["stWcq"] = (102,0,204)
colors["stWuq"] = (198,101,222)
colors["stnW"] = (207,176,235)

if args.stack:
  ymax = {}
  for step in ["prefit","postfit"]:
    ymax[step] = {}
    for channel in ["ch1","ch2","ch3","ch4"]:
       ymax[step][channel] = 0
       for s in samplesMC:
         #print(hist_content[step][channel][s].Integral())
         if not (s == "mcss" and (channel == "ch1" or channel=="ch2")):
           hist_content[step][channel][s].SetLineWidth(1)
           hist_content[step][channel][s].SetLineColor(kBlack)
           hist_content[step][channel][s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
           hist_content[step][channel][s].GetYaxis().SetTitle("Number of events")
           hist_content[step][channel][s].GetXaxis().SetTitle(cosmetic_names[str(args.var)])
           y1 = hist_content[step][channel][s].GetMaximum()
           if y1>ymax[step][channel]: ymax[step][channel]=y1
       if args.linear: hist_content[step][channel]["bck"].SetMaximum(1.45*ymax[step][channel])

  ## Stack creation
  samplesMC = ["bck","mcss","ttdl","vjets","stnW","stWuq","ttWuq","stWcq","ttWcq"]

  if args.ratio and not args.nodata: upper_pad.cd()

  stack_T = {}
  for chan in ["ch1","ch2","ch3","ch4"]:
    stack_T[chan] = {}
    for step in ["prefit","postfit"]:
      stack_T[chan][step] = ROOT.THStack("hs", ";;Events")
      for s in samplesMC:
        if not (s == "mcss" and (chan=="ch1" or chan=="ch2")):
          stack_T[chan][step].Add(hist_content[step][chan][s])

  for chan in ["ch1","ch2","ch3","ch4"]:
    for step in ["prefit","postfit"]:
      y_T = stack_T[chan][step].GetMaximum()
      if y_T>ymax[step][chan]: ymax[step][chan]=y_T
      if args.linear: stack_T[chan][step].SetMaximum(1.45*ymax[step][chan])
      stack_T[chan][step].Draw("HIST")
      hist_content[step][chan]["TotalProcs"].SetFillColorAlpha(14, 1)
      hist_content[step][chan]["TotalProcs"].SetFillStyle(3002)
      hist_content[step][chan]["TotalProcs"].Draw("E2SAME")
      hist_content[step][chan]["data"].Draw("E SAME")
      hist_content[step][chan]["data"].SetMarkerStyle(21)

      CMS_lumi.CMS_lumi(upper_pad, iPeriod, iPos)
      upper_pad.cd()
      upper_pad.Update()
      upper_pad.RedrawAxis()
      frame = upper_pad.GetFrame()
      #frame.Draw()

      if args.ratio and not args.nodata:
        lower_pad.cd()
        ratio = hist_content[step][chan]["data"].Clone("ratio")
        ratio.SetLineColor(kBlack)
        ratio.SetMarkerStyle(21)
        ratio.SetTitle("")
        ratio.SetMinimum(c_rat2)
        ratio.SetMaximum(c_rat)
        ratio.GetYaxis().SetTitle("Data/MC")
        ratio.GetXaxis().SetTitle(cosmetic_names[str(args.var)])
        ratio.GetXaxis().SetLabelSize(0.08)
        ratio.GetXaxis().SetTitleSize(tit_size)
        ratio.GetXaxis().SetTitleOffset(titX_off)
        ratio.GetYaxis().SetLabelSize(0.05)
        ratio.GetYaxis().SetTitleSize(0.12)
        ratio.GetYaxis().CenterTitle(False)
        ratio.GetYaxis().ChangeLabel(nrat,-1,-1,-1,-1,-1,"  ");
        ratio.GetYaxis().SetTitleOffset(titY_off)
        # Set up plot for markers and errors
        #ratio.Sumw2()
        ratio.SetStats(0)
        ratio.Draw("ep")
        hTotal = hist_content[step][chan]["bck"].Clone('hTotal')
        for s in samplesMC[1:]:
          if not (s == "mcss" and (chan=="ch1" or chan=="ch2")):
             #print(s)
             hTotal.Add(hist_content[step][chan][s])
        ratio.Divide(hTotal)
        aux_err = hist_content[step][chan]["TotalProcs"].Clone("ratio_err")
        ratio_err = hist_content[step][chan]["ttWcq"].Clone("ratio_err")
        for bin in range(ratio_err.GetNbinsX()):
          ratio_err.SetBinContent(bin,ratio.GetBinContent(bin));
          if aux_err.GetBinContent(bin) > 0:
            ratio_err.SetBinError(bin,aux_err.GetBinError(bin)/aux_err.GetBinContent(bin));
          else:
            ratio_err.SetBinError(bin,0.);
        ratio_err.SetFillColorAlpha(15, 1)
        ratio_err.SetFillStyle(3001)
        ratio_err.Draw("E2SAME")
        ratio.Draw("epSAME")


      ## Legends
      if args.ratio and not args.nodata: upper_pad.cd()
      leg = TLegend(leg1,leg2,leg3,leg4)
      leg.SetBorderSize(0)
      leg.SetNColumns(2)
      leg.AddEntry(hist_content[step][chan]["ttWcq"],"t#bar{t} cq","f")
      leg.AddEntry(hist_content[step][chan]["ttWuq"],"t#bar{t} uq","f")
      leg.AddEntry(hist_content[step][chan]["stWcq"],"Single top charm","f")
      leg.AddEntry(hist_content[step][chan]["stWuq"],"Single top light","f")
      leg.AddEntry(hist_content[step][chan]["stnW"],"Single top no had W","f")
      leg.AddEntry(hist_content[step][chan]["vjets"],"Z+jets","f")
      leg.AddEntry(hist_content[step][chan]["ttdl"],"dileptonic t#bar{t}","f")
      if (chan=="ch3" or chan=="ch4"):leg.AddEntry(hist_content[step][chan]["mcss"],"Same sign contribution","f")
      leg.AddEntry(hist_content[step][chan]["bck"],"vv plus hadronic t#bar{t}","f")
      leg.AddEntry(hist_content[step][chan]["data"], "Data" ,"lep")
      leg.AddEntry(hist_content[step][chan]["TotalProcs"], "Uncertainty" ,"f")
      leg.Draw()
      termp= "fit_distribution_"
      notation = str(step)+"_"+str(chan)
      if args.png: c1.Print(plotdir+termp+notation+str(args.var)+".pdf")
      else: c1.Print(plotdir+termp+notation+str(args.var)+".root")
  if c1: 
    c1.Close(); gSystem.ProcessEvents();

file_aux.Close()


