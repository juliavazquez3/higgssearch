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
parser.add_argument("--higgsmass", type=string, default="80",
                    help="Select type of process to run")

args = parser.parse_args()

## Open hists files

filePath = "/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/"
term = "histstt_deeptag_fromJF"

datayears = ["2016","2016B","2017","2018"]
#datayears = ["2018","2016","2016B"]

samplesHT=["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "ttbar_sl_charm","ttbar_sl_bottom","ttbar_sl_light","ttbar_dl","ttbar_dh",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6",
        "zjets_7","zjets_8","st_1","st_2","st_3","st_4","zz","wz"]

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

histFile["2016"]["higgs"] = TFile.Open(filePath + term+"_"+"higgs_tocs_m"+str(args.higgsmass)+"2016"+".root","READ")
#print(histFile.keys())

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

n_eventsH = {}
n_eventsH["80"] = 978770 
n_eventsH["90"]	= 994040
n_eventsH["100"] = 997702
n_eventsH["120"] = 1000000
n_eventsH["140"] = 995092
n_eventsH["150"] = 994950
n_eventsH["155"] = 954536
n_eventsH["160"] = 996144

xsec_higgs_E = 831*0.326*0.02
xsec_higgs_M = 831*0.326*0.02
lumi_higgs_E = n_eventsH[str(args.higgsmass)]/(xsec_higgs_E*1000)
lumi_higgs_M = n_eventsH[str(args.higgsmass)]/(xsec_higgs_M*1000)
lumi_d = {"2016": 19.5, "2016B": 16.8,"2017": 41.5,"2018": 59.8}

print("Lumi for higgs process is ",lumi_higgs_E)

print("data lumis are",lumi_d)
#print(nevents_d)
print("the channel is wqq")
print("the higgs mass is ",args.higgsmass)

#######################################################
########### Start of plot creation ####################
#######################################################

samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "ttbar_sl_charm","ttbar_sl_bottom","ttbar_sl_light","ttbar_dl","ttbar_dh","zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6",
        "zjets_7","zjets_8","st_1","st_2","st_3","st_4","zz","wz"]

## HISTS
hists_M = {}
hists_E = {}
for data_op in datayears:
    hists_M[data_op] = {}
    hists_E[data_op] = {}
    for s in samples:
      hists_M[data_op][s] = histFile[data_op][s].Get("InvM_2jets_M")
      hists_E[data_op][s] = histFile[data_op][s].Get("InvM_2jets_E")
    hists_M[data_op]["data"] = histFileD[data_op]["M"].Get("InvM_2jets_M")
    hists_E[data_op]["data"] = histFileD[data_op]["E"].Get("InvM_2jets_E")

hists_M["2016"]["higgs"] = histFile["2016"]["higgs"].Get("InvM_2jets_M")
hists_E["2016"]["higgs"] = histFile["2016"]["higgs"].Get("InvM_2jets_E")
  
#print(samples)

## Scaling to lumi
#print(name) 
for data_op in datayears:
    lumi_data = lumi_d[data_op]
    for s in samples:
      #print(s)
      #print(data_op)
      hists_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hists_E[data_op][s].Scale(lumi_data/lumi[data_op][s])

hists_M["2016"]["higgs"].Scale((lumi_d["2016"]/lumi_higgs_M)+(lumi_d["2016B"]/lumi_higgs_M)+(lumi_d["2017"]/lumi_higgs_M)+(lumi_d["2018"]/lumi_higgs_M))
hists_E["2016"]["higgs"].Scale((lumi_d["2016"]/lumi_higgs_E)+(lumi_d["2016B"]/lumi_higgs_E)+(lumi_d["2017"]/lumi_higgs_E)+(lumi_d["2018"]/lumi_higgs_E))

##########################################
############# Adding all years ###########
##########################################

histT_M = {}
histT_E = {}
for s in samples:
   histT_M[s] = hists_M[datayears[0]][s]
   histT_E[s] = hists_E[datayears[0]][s]
   for d in datayears[1:]:
        histT_M[s].Add(hists_M[d][s])
        histT_E[s].Add(hists_E[d][s])

histT_M["data"] = hists_M[datayears[0]]["data"]
histT_E["data"] = hists_E[datayears[0]]["data"]
for d in datayears[1:]:
    histT_M["data"].Add(hists_M[d]["data"])
    histT_E["data"].Add(hists_E[d]["data"])

histT_M["higgs"] = hists_M["2016"]["higgs"]
histT_E["higgs"] = hists_E["2016"]["higgs"]

#################################
### Joining certain processes ###
#################################

histT_M["vv"] = histT_M["ww"]
histT_M["vv"].Add(histT_M["wz"])
histT_M["vv"].Add(histT_M["zz"])
histT_M["wjets"] = histT_M["wjets_1"]
histT_M["wjets"].Add(histT_M["wjets_2"])
histT_M["wjets"].Add(histT_M["wjets_3"])
histT_M["wjets"].Add(histT_M["wjets_4"])
histT_M["wjets"].Add(histT_M["wjets_5"])
histT_M["wjets"].Add(histT_M["wjets_6"])
histT_M["wjets"].Add(histT_M["wjets_7"])
histT_M["wjets"].Add(histT_M["wjets_8"])
histT_M["zjets"] = histT_M["zjets_1"]
histT_M["zjets"].Add(histT_M["zjets_2"])
histT_M["zjets"].Add(histT_M["zjets_3"])
histT_M["zjets"].Add(histT_M["zjets_4"])
histT_M["zjets"].Add(histT_M["zjets_5"])
histT_M["zjets"].Add(histT_M["zjets_6"])
histT_M["zjets"].Add(histT_M["zjets_7"])
histT_M["zjets"].Add(histT_M["zjets_8"])
histT_M["st"] = histT_M["st_1"]
histT_M["st"].Add(histT_M["st_2"])
histT_M["st"].Add(histT_M["st_3"])
histT_M["st"].Add(histT_M["st_4"])
histT_M["ttbar_rest"] = histT_M["ttbar_dl"]
histT_M["ttbar_rest"].Add(histT_M["ttbar_dh"])

histT_E["vv"] = histT_E["ww"]
histT_E["vv"].Add(histT_E["wz"])
histT_E["vv"].Add(histT_E["zz"])
histT_E["wjets"] = histT_E["wjets_1"]
histT_E["wjets"].Add(histT_E["wjets_2"])
histT_E["wjets"].Add(histT_E["wjets_3"])
histT_E["wjets"].Add(histT_E["wjets_4"])
histT_E["wjets"].Add(histT_E["wjets_5"])
histT_E["wjets"].Add(histT_E["wjets_6"])
histT_E["wjets"].Add(histT_E["wjets_7"])
histT_E["wjets"].Add(histT_E["wjets_8"])
histT_E["zjets"] = histT_E["zjets_1"]
histT_E["zjets"].Add(histT_E["zjets_2"])
histT_E["zjets"].Add(histT_E["zjets_3"])
histT_E["zjets"].Add(histT_E["zjets_4"])
histT_E["zjets"].Add(histT_E["zjets_5"])
histT_E["zjets"].Add(histT_E["zjets_6"])
histT_E["zjets"].Add(histT_E["zjets_7"])
histT_E["zjets"].Add(histT_E["zjets_8"])
histT_E["st"] = histT_E["st_1"]
histT_E["st"].Add(histT_E["st_2"])
histT_E["st"].Add(histT_E["st_3"])
histT_E["st"].Add(histT_E["st_4"])
histT_E["ttbar_rest"] = histT_E["ttbar_dl"]
histT_E["ttbar_rest"].Add(histT_E["ttbar_dh"])

##########################################
############# Change of names ############
##########################################

histT_E["data"].SetName("mjj4jdeeptag"+"_1")
histT_M["data"].SetName("mjj4jdeeptag"+"_2")
histT_E["ttbar_sl_charm"].SetName("mjj4jdeeptag"+"_3")
histT_M["ttbar_sl_charm"].SetName("mjj4jdeeptag"+"_4")
histT_E["ttbar_rest"].SetName("mjj4jdeeptag"+"_5")
histT_M["ttbar_rest"].SetName("mjj4jdeeptag"+"_6")
histT_E["st"].SetName("mjj4jdeeptag"+"_7")
histT_M["st"].SetName("mjj4jdeeptag"+"_8")
histT_E["wjets"].SetName("mjj4jdeeptag"+"_9")
histT_M["wjets"].SetName("mjj4jdeeptag"+"_10")
histT_E["zjets"].SetName("mjj4jdeeptag"+"_11")
histT_M["zjets"].SetName("mjj4jdeeptag"+"_12")
histT_E["vv"].SetName("mjj4jdeeptag"+"_13")
histT_M["vv"].SetName("mjj4jdeeptag"+"_14")
histT_E["higgs"].SetName("mjj4jdeeptag"+"_15")
histT_M["higgs"].SetName("mjj4jdeeptag"+"_16")
histT_E["ttbar_sl_light"].SetName("mjj4jdeeptag"+"_17")
histT_M["ttbar_sl_light"].SetName("mjj4jdeeptag"+"_18")
histT_E["ttbar_sl_bottom"].SetName("mjj4jdeeptag"+"_19")
histT_M["ttbar_sl_bottom"].SetName("mjj4jdeeptag"+"_20")

samplesF = ["data","ttbar_sl_charm","ttbar_rest","st","wjets","zjets","vv","higgs","ttbar_sl_light","ttbar_sl_bottom"]
 
path = "/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/deeptag/hm"+str(args.higgsmass)+"/mjj4jdeeptag.root"

myfile = TFile( path , 'RECREATE')
for s in samplesF:
  histT_E[s].Write()
  histT_M[s].Write()
myfile.Close()

