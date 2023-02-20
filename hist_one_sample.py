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
if args.ssos: filePath = "/nfs/cms/vazqueze/higgssearch/hists/ssos/"
else: filePath = "/nfs/cms/vazqueze/higgssearch/hists/"

if args.ssos: ssos_add = "SSOS"
else: ssos_add = "" 

term = "hists_higgs"

#datayears = ["2016","2016B","2017","2018"]
datayears = ["2016"]

samplesHT = ["ww","wjets_2_light","wjets_2_bottom","wjets_1_charm","wjets_1_doublecharm","wjets_2_light","wjets_2_bottom","wjets_2_charm","wjets_2_doublecharm",
        "wjets_3_light","wjets_3_bottom","wjets_3_charm","wjets_3_doublecharm","wjets_4_light","wjets_4_bottom","wjets_4_charm","wjets_4_doublecharm",
        "wjets_5_light","wjets_5_bottom","wjets_5_charm","wjets_5_doublecharm","wjets_6_light","wjets_6_bottom","wjets_6_charm","wjets_6_doublecharm",
        "wjets_7_light","wjets_7_bottom","wjets_7_charm","wjets_7_doublecharm","wjets_8_light","wjets_8_bottom","wjets_8_charm","wjets_8_doublecharm",
        "ttbar_sl_charm","ttbar_sl_nocharm","ttbar_dl_charm","ttbar_dl_nocharm","ttbar_dh_charm","ttbar_dh_nocharm","zjets_1",
        "zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8","wz","zz","st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm",
        "st_2_nocharm","st_3_nocharm","st_4_nocharm", "higgs_tocs_m80", "higgs_tocs_m90", "higgs_tocs_m100", "higgs_tocs_m120", "higgs_tocs_m140",
        "higgs_tocs_m150", "higgs_tocs_m155", "higgs_tocs_m160"]

## Adding QCD

histFile = {}

for data_op in datayears:
	## mc files
	histFile[data_op] = {}
	for s in samplesHT:
		if (s[0:6] == "wjets_") and isfile(filePath + term+ssos_add+"_"+s[0:7]+data_op+s[7:]+".root"):
			histFile[data_op][s] = TFile.Open(filePath + term+ssos_add+"_"+s[0:7]+data_op+s[7:]+".root","READ")
		elif s[0:6] == "ttbar_" and isfile(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root"):
			histFile[data_op][s] = TFile.Open(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root","READ")
		elif s[0:2] == "st" and isfile(filePath + term+ssos_add+"_"+s[0:4]+data_op+s[4:]+".root"):
			histFile[data_op][s] = TFile.Open(filePath + term+ssos_add+"_"+s[0:4]+data_op+s[4:]+".root","READ")
		elif isfile(filePath + term+ssos_add+"_"+s+data_op+".root"):
			histFile[data_op][s] = TFile.Open(filePath + term+ssos_add+"_"+s+data_op+".root","READ")
	#print(histFile[data_op].keys())
#print(histFile.keys())

histNames = []

histNames.append("nGenJetGood_M")
histNames.append("nGenJetGood_E")
histNames.append("nLooseLepton_M")
histNames.append("nLooseLepton_E")
histNames.append("second_muon_pt_M")
histNames.append("second_muon_pt_E")
histNames.append("second_electron_pt_M")
histNames.append("second_electron_pt_E")
histNames.append("second_muon_eta_M")
histNames.append("second_muon_eta_E")
histNames.append("second_electron_eta_M")
histNames.append("second_electron_eta_E")
histNames.append("jet_muon_pt_M")
histNames.append("jet_muon_mass_M")
histNames.append("jet_not_muon_pt_M")
histNames.append("jet_muon_eta_M")
histNames.append("jet_not_muon_eta_M")
histNames.append("jet_notmuon_mass_M")
histNames.append("jet_muon_pt_E")
histNames.append("jet_muon_mass_E")
histNames.append("jet_not_muon_pt_E")
histNames.append("jet_muon_eta_E")
histNames.append("jet_not_muon_eta_E")
histNames.append("jet_notmuon_mass_E")
histNames.append("jet_bot1_pt_M")
histNames.append("jet_bot1_mass_M")
histNames.append("jet_bot2_pt_M")
histNames.append("jet_bot1_eta_M")
histNames.append("jet_bot2_eta_M")
histNames.append("jet_bot2_mass_M")
histNames.append("jet_bot1_pt_E")
histNames.append("jet_bot1_mass_E")
histNames.append("jet_bot2_pt_E")
histNames.append("jet_bot1_eta_E")
histNames.append("jet_bot2_eta_E")
histNames.append("jet_bot2_mass_E")
histNames.append("lepton_pt_M")
histNames.append("lepton_eta_M")
histNames.append("lepton_pt_E")
histNames.append("lepton_eta_E")
histNames.append("InvM_2jets_M")
histNames.append("InvM_2jets_E")
histNames.append("transverse_massM")
histNames.append("transverse_massE")

#histNames = ["nMuon_nobot_M","nMuon_nobot_E"]
#histNames = ["jet_muon_coef_M","jet_muon_coef_E","jet_muon_ptresol_E","jet_muon_ptresol_M",
#	"jet_muon_scalefactor_M","jet_muon_scalefactor_E"]

#histNames = [
#       "jet_muon_flavourH_M","jet_muon_flavourH_E","jet_notmuon_flavourH_M","jet_notmuon_flavourH_E",
#       "jet_muon_flavourP_M","jet_muon_flavourP_E","jet_notmuon_flavourP_M","jet_notmuon_flavourP_E"
#]

not_rebin = ["nJetGood_M","nJetGood_E","nMuoninJet_M","nMuoninJet_E","jet_muon_nmu_M","jet_muon_nmu_E","SSOS_M","SSOS_E",
       "nLooseLepton_M","nLooseLepton_E", "muon_jet_tight_M","muon_jet_tight_E",
       "jet_muon_flavourH_M","jet_muon_flavourH_E","jet_notmuon_flavourH_M","jet_notmuon_flavourH_E",
       "jet_muon_flavourP_M","jet_muon_flavourP_E","jet_notmuon_flavourP_M","jet_notmuon_flavourP_E","nMuon_nobot_M","nMuon_nobot_E",
       "second_muon_pt_M","second_muon_pt_E","second_electron_pt_M","second_electron_pt_E"]

samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8","ttbar_sl","ttbar_dl","ttbar_dh","zjets_1",
        "zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8","wz","zz","st_1","st_2","st_3","st_4"]

#samples_d = ["2018"]
samples_d = ["2016","2016B","2017","2018"]
#samples = ["WW","Wjets","ttbar"]

samples_mc = []

for s in samples:
  for y in samples_d:
    samples_mc.append(s+y)

samples_year = []

for y in samples_d:
  samples_year.append(y+"M")
  samples_year.append(y+"E")

## lumi info

lumi = {}
xsecs = {}
nevents = {}

for data_op in samples_d:
	files = json.load(open("/nfs/cms/vazqueze/ttbaranalisis/mcinfo"+data_op+".json"))
	lumi[data_op] = {}
	for p in samples:
		#num_events = files[p]["events"] # Number of events
		num_files = files[p]["files"] # Number of files
		luminosity = files[p]["lumi"] # Luminosity
		#print(files[p]["type"])
		lumi[data_op][p] = luminosity

for data_op in samples_d:
	lumi[data_op]["ttbar_sl_charm"] = lumi[data_op]["ttbar_sl"]
	lumi[data_op]["ttbar_sl_nocharm"] = lumi[data_op]["ttbar_sl"]
	lumi[data_op]["ttbar_dl_charm"] = lumi[data_op]["ttbar_dl"]
	lumi[data_op]["ttbar_dl_nocharm"] = lumi[data_op]["ttbar_dl"]
	lumi[data_op]["ttbar_dh_charm"] = lumi[data_op]["ttbar_dh"]
	lumi[data_op]["ttbar_dh_nocharm"] = lumi[data_op]["ttbar_dh"]

	for i in np.arange(8)+1:
		lumi[data_op]["wjets_"+str(i)+"_charm"] = lumi[data_op]["wjets_"+str(i)]
		lumi[data_op]["wjets_"+str(i)+"_doublecharm"] = lumi[data_op]["wjets_"+str(i)]
		lumi[data_op]["wjets_"+str(i)+"_bottom"] = lumi[data_op]["wjets_"+str(i)]
		lumi[data_op]["wjets_"+str(i)+"_light"] = lumi[data_op]["wjets_"+str(i)]

	for i in np.arange(4)+1:
		lumi[data_op]["st_"+str(i)+"_charm"] = lumi[data_op]["st_"+str(i)]
		lumi[data_op]["st_"+str(i)+"_nocharm"] = lumi[data_op]["st_"+str(i)]

	lumi[data_op]["higgs_tocs_m80"] = 1.
	lumi[data_op]["higgs_tocs_m90"] = 1.
	lumi[data_op]["higgs_tocs_m100"] = 1.
	lumi[data_op]["higgs_tocs_m120"] = 1.
	lumi[data_op]["higgs_tocs_m140"] = 1.
	lumi[data_op]["higgs_tocs_m150"] = 1. 
	lumi[data_op]["higgs_tocs_m155"] = 1.
	lumi[data_op]["higgs_tocs_m160"] = 1.

#print("mc lumis",lumi)

files_d = json.load(open("/nfs/cms/vazqueze/ttbaranalisis/datainfo.json"))

lumi_d = {}
xsecs_d = {}
nevents_d = {}

for p in samples_year:
    num_events = files_d[p]["events"] # Number of events
    num_files = files_d[p]["files"] # Number of files
    luminosity = files_d[p]["lumi"] # Luminosity
    #print(len(list_files))
    lumi_d[p[:-1]] = luminosity
    nevents_d[p[:-1]] = num_events

#print("data lumis are",lumi_d)
#print(nevents_d)

#histNames = [h for h in histNames if (h[-1]=="M")]

#######################################################
########### Start of plot creation ####################
#######################################################

for name in histNames:

  samples = ["higgs_tocs_m80", "higgs_tocs_m90", "higgs_tocs_m100", "higgs_tocs_m120", "higgs_tocs_m140",
        "higgs_tocs_m150", "higgs_tocs_m155", "higgs_tocs_m160"]

  ## HISTS
  samples_foryear = {}
  histss = {}
  hdataM = {}
  hdataE = {}
  for data_op in datayears:
    samples_foryear[data_op] = [s for s in samples if s in histFile[data_op].keys()]
    histss[data_op] = {}
    for s in samples_foryear[data_op]:
      histss[data_op][s] = histFile[data_op][s].Get(name)
  
  #for data_op in datayears:
  #  lumi_data = lumi_d[data_op]
  #  for s in samples_foryear[data_op]:
  #    #print(s)
  #    if s != "QCD": histss[data_op][s].Scale(lumi_data/lumi[data_op][s])

  histT = {}
  for s in samples:
    if s in samples_foryear[datayears[0]]:
       histT[s] = histss[datayears[0]][s]
    for d in datayears[1:]:
       if s in samples_foryear[d]: histT[s].Add(histss[d][s])

  if (name=="nGenJetGood_E" or name=="nGenJetGood_M"):
    for s in samples:
      print(s)
      print("Number of events for "+name[-1]+" channel in the sample "+s+" is "+str(histT[s].Integral()))

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  canv = {}

  for s in samples:
    gStyle.SetOptTitle(0)
    t = TPaveText(0.0, 0.9, 0.3, 1.0, str(s))
    canv[s] = TCanvas(s,"",1200,800)
    histT[s].SetLineWidth(1)
    histT[s].GetYaxis().SetTitle("Number of events")
    histT[s].GetXaxis().SetTitle(name)
    ymax = histT[s].GetMaximum()
    histT[s].SetMaximum(1.3*ymax)
    histT[s].Draw("HIST")
    term = str(s)
    notation = "_one_sample_"
    canv[s].Print(plotdir+ssos_add+term+notation+ name + ".png")

samplesHT = ["higgs_tocs_m80", "higgs_tocs_m90", "higgs_tocs_m100", "higgs_tocs_m120", "higgs_tocs_m140",
        "higgs_tocs_m150", "higgs_tocs_m155", "higgs_tocs_m160"]

for s in samplesHT:
        for data_op in datayears:
                        histFile[s][data_op].Close()

for data_op in datayears:
	histFileD[data_op]["M"].Close()
	histFileD[data_op]["E"].Close()

