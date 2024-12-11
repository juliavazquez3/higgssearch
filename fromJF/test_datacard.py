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

## Open hists files
filePath = "/nfs/cms/vazqueze/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit/data/tutorials/longexercise/makePrePostplots/output_histograms.root"
## mc files
file_aux = TFile.Open(filePath,"READ")

dir = file_aux.GetDirectory("ch1_ch1_prefit");
dir.cd()

c1 = TCanvas("c1","",1200,800)

TotalProcs.SetFillColorAlpha(15, 0.02)
TotalProcs.SetFillStyle(3001)
TotalProcs.Draw("E2")


c1.Print("test_datacard.png")
c1.Close(); gSystem.ProcessEvents();

file_aux.Close()


