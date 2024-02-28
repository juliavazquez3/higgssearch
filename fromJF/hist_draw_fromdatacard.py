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
parser.add_argument("--nodata", action="store_true", default=False,
                    help="Do not plot data")

# Use like:
# python higgssearch/fromJF/hist_fromJFwqq_syst_total.py --stack --ratio --png --norm --year="all" --channel="btagMM_chitest_sl"

args = parser.parse_args()

plotdir = '/nfs/cms/vazqueze/higgssearch/plotspng/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

c_rat = 1.2
c_rat2 = 0.8

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
filePath = "/nfs/cms/vazqueze/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit/data/tutorials/longexercise/fitDiagnostics.postfit.root"
## mc files
file_aux = TFile.Open(filePath,"READ")

samples = ["ttWcq","ttWuq","stWcq","stWuq","bck","stnW","vjets","mcss","ttdl","data","total"]
samplesMC = ["ttWcq","ttWuq","stWcq","stWuq","bck","stnW","vjets","mcss","ttdl"]
for s in samples:
    hist_content["prefit"]["ch1"][s] = file_aux.Get("shapes_prefit/ch1_ch1/"+str(s))
    hist_content["prefit"]["ch2"][s] = file_aux.Get("shapes_prefit/ch1_ch2/"+str(s))
    hist_content["prefit"]["ch3"][s] = file_aux.Get("shapes_prefit/ch1_ch3/"+str(s))
    hist_content["prefit"]["ch4"][s] = file_aux.Get("shapes_prefit/ch1_ch4/"+str(s))
    hist_content["postfit"]["ch1"][s] = file_aux.Get("shapes_fit_s/ch1_ch1/"+str(s))
    hist_content["postfit"]["ch2"][s] = file_aux.Get("shapes_fit_s/ch1_ch2/"+str(s))
    hist_content["postfit"]["ch3"][s] = file_aux.Get("shapes_fit_s/ch1_ch3/"+str(s))
    hist_content["postfit"]["ch4"][s] = file_aux.Get("shapes_fit_s/ch1_ch4/"+str(s))

nPoints = hist_content["postfit"]["ch4"]["data"].GetN()
for step in ["prefit","postfit"]:
    for channel in ["ch1","ch2","ch3","ch4"]:
        hist_content[step][channel]["data_aux"] = hist_content[step][channel]["bck"].Clone("data_aux")
        hist_content[step][channel]["data_aux"].Scale(0.0)
        for i in range(nPoints):
            ygr = hist_content[step][channel]["data"].GetPointY(i)
            hist_content[step][channel]["data_aux"].SetBinContent(i+1,ygr)


if args.png: c1 = TCanvas("c1","",1200,800)
else: c1 = TCanvas("c1","",600,400)

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
           hist_content[step][channel][s].GetXaxis().SetTitle("chi2 value")
           y1 = hist_content[step][channel][s].GetMaximum()
           if y1>ymax[step][channel]: ymax[step][channel]=y1
       if args.linear: hist_content[step][channel]["bck"].SetMaximum(1.3*ymax[step][channel])

  ## Stack creation
  samplesMC = ["bck","mcss","ttdl","vjets","stnW","stWuq","ttWuq","stWcq","ttWcq"]

  if args.ratio and not args.nodata: upper_pad.cd()

  stack_T = {}
  for chan in ["ch1","ch2","ch3","ch4"]:
    stack_T[chan] = {}
    for step in ["prefit","postfit"]:
      stack_T[chan][step] = ROOT.THStack()
      for s in samplesMC:
        if not (s == "mcss" and (chan=="ch1" or chan=="ch2")):
          stack_T[chan][step].Add(hist_content[step][chan][s])

  for chan in ["ch1","ch2","ch3","ch4"]:
    for step in ["prefit","postfit"]:
      y_T = stack_T[chan][step].GetMaximum()
      if y_T>ymax[step][chan]: ymax[step][chan]=y_T
      if args.linear: stack_T[chan][step].SetMaximum(1.3*ymax[step][chan])
      stack_T[chan][step].Draw("HIST")
      hist_content[step][chan]["data"].Draw("PSAME")
      hist_content[step][chan]["total"].SetFillColorAlpha(15, 0.02)
      hist_content[step][chan]["total"].SetFillStyle(3001)
      hist_content[step][chan]["total"].Draw("E2SAME")

      if args.ratio and not args.nodata:
        lower_pad.cd()
        ratio = hist_content[step][chan]["data_aux"].Clone("ratio")
        ratio.SetLineColor(kBlack)
        ratio.SetMarkerStyle(21)
        ratio.SetTitle("")
        ratio.SetMinimum(c_rat2)
        ratio.SetMaximum(c_rat)
        ratio.GetYaxis().SetTitle("Data/MC")
        ratio.GetXaxis().SetTitle("chi2_value")
        ratio.GetXaxis().SetLabelSize(0.08)
        ratio.GetXaxis().SetTitleSize(0.12)
        ratio.GetXaxis().SetTitleOffset(1.0)
        ratio.GetYaxis().SetLabelSize(0.05)
        ratio.GetYaxis().SetTitleSize(0.09)
        ratio.GetYaxis().CenterTitle()
        ratio.GetYaxis().SetTitleOffset(0.5)
        # Set up plot for markers and errors
        #ratio.Sumw2()
        ratio.SetStats(0)
        hTotal = hist_content[step][chan]["bck"].Clone('hTotal')
        for s in samplesMC[1:]:
          if not (s == "mcss" and (chan=="ch1" or chan=="ch2")):
             #print(s)
             hTotal.Add(hist_content[step][chan][s])
        ratio.Divide(hTotal)
        ratio.Draw("ep")

      ## Legends
      if args.ratio and not args.nodata: upper_pad.cd()
      leg = TLegend(0.78,0.4,0.95,0.95)
      leg.SetBorderSize(1)
      leg.AddEntry(hist_content[step][chan]["ttWcq"],"t#bar{t} cq","f")
      leg.AddEntry(hist_content[step][chan]["ttWuq"],"t#bar{t} uq","f")
      leg.AddEntry(hist_content[step][chan]["stWcq"],"Single top charm","f")
      leg.AddEntry(hist_content[step][chan]["stWuq"],"Single top light","f")
      leg.AddEntry(hist_content[step][chan]["stnW"],"Single top no had W","f")
      leg.AddEntry(hist_content[step][chan]["vjets"],"Z+jets","f")
      leg.AddEntry(hist_content[step][chan]["ttdl"],"dileptonic t#bar{t}","f")
      leg.AddEntry(hist_content[step][chan]["mcss"],"Same sign contribution","f")
      leg.AddEntry(hist_content[step][chan]["bck"],"vv plus hadronic t#bar{t}","f")
      leg.AddEntry(hist_content[step][chan]["data"], "Data" ,"lep")
      leg.Draw()
      termp= "fit_distribution_"
      notation = str(step)+"_"+str(chan)
      if args.png: c1.Print(plotdir+termp+notation+".png")
      else: c1.Print(plotdir+termp+notation+".root")
  if c1: 
    c1.Close(); gSystem.ProcessEvents();

file_aux.Close()


