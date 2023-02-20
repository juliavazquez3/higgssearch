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
parser.add_argument("--type", type=string, default="sl",
                    help="Selec type of data to run")

args = parser.parse_args()

## Open hists files

if args.type == "sl":
  filePath = "/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/"
  term = "histstt_SL_fromJF"
else:
  filePath = "/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/sv/"
  term = "histstt_SV_fromJF"

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
print("the channel is ",args.type)
print("the higgs mass is ",args.higgsmass)

sv_sel = {}
sv_sel["ss"] = 0.30
sv_sel["os"] = 0.30

#######################################################
########### Start of plot creation ####################
#######################################################

samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "ttbar_sl_charm","ttbar_sl_bottom","ttbar_sl_light","ttbar_dl","ttbar_dh","zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6",
        "zjets_7","zjets_8","st_1","st_2","st_3","st_4","zz","wz"]

## HISTS
hists_M_ss = {}
hists_M_os = {}
hists_E_ss = {}
hists_E_os = {}
for data_op in datayears:
    hists_M_ss[data_op] = {}
    hists_M_os[data_op] = {}
    hists_E_ss[data_op] = {}
    hists_E_os[data_op] = {}
    for s in samples:
      hists_M_ss[data_op][s] = histFile[data_op][s].Get("InvM_2jets_M_ss")
      hists_M_os[data_op][s] = histFile[data_op][s].Get("InvM_2jets_M_os")
      hists_E_ss[data_op][s] = histFile[data_op][s].Get("InvM_2jets_E_ss")
      hists_E_os[data_op][s] = histFile[data_op][s].Get("InvM_2jets_E_os")
    hists_M_ss[data_op]["data"] = histFileD[data_op]["M"].Get("InvM_2jets_M_ss")
    hists_M_os[data_op]["data"] = histFileD[data_op]["M"].Get("InvM_2jets_M_os")
    hists_E_ss[data_op]["data"] = histFileD[data_op]["E"].Get("InvM_2jets_E_ss")
    hists_E_os[data_op]["data"] = histFileD[data_op]["E"].Get("InvM_2jets_E_os")

hists_M_ss["2016"]["higgs"] = histFile["2016"]["higgs"].Get("InvM_2jets_M_ss")
hists_M_os["2016"]["higgs"] = histFile["2016"]["higgs"].Get("InvM_2jets_M_os")
hists_E_ss["2016"]["higgs"] = histFile["2016"]["higgs"].Get("InvM_2jets_E_ss")
hists_E_os["2016"]["higgs"] = histFile["2016"]["higgs"].Get("InvM_2jets_E_os")
  
#print(samples)

## Scaling to lumi
#print(name) 
for data_op in datayears:
    lumi_data = lumi_d[data_op]
    for s in samples:
      #print(s)
      #print(data_op)
      hists_M_ss[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hists_M_os[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hists_E_ss[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hists_E_os[data_op][s].Scale(lumi_data/lumi[data_op][s])

hists_M_ss["2016"]["higgs"].Scale((lumi_d["2016"]/lumi_higgs_M)+(lumi_d["2016B"]/lumi_higgs_M)+(lumi_d["2017"]/lumi_higgs_M)+(lumi_d["2018"]/lumi_higgs_M))
hists_M_os["2016"]["higgs"].Scale((lumi_d["2016"]/lumi_higgs_M)+(lumi_d["2016B"]/lumi_higgs_M)+(lumi_d["2017"]/lumi_higgs_M)+(lumi_d["2018"]/lumi_higgs_M))
hists_E_ss["2016"]["higgs"].Scale((lumi_d["2016"]/lumi_higgs_E)+(lumi_d["2016B"]/lumi_higgs_E)+(lumi_d["2017"]/lumi_higgs_E)+(lumi_d["2018"]/lumi_higgs_E))
hists_E_os["2016"]["higgs"].Scale((lumi_d["2016"]/lumi_higgs_E)+(lumi_d["2016B"]/lumi_higgs_E)+(lumi_d["2017"]/lumi_higgs_E)+(lumi_d["2018"]/lumi_higgs_E))

if args.type == "sv":
   hists_M_ss["2016"]["higgs"].Scale(sv_sel["ss"]*0.33)
   hists_M_os["2016"]["higgs"].Scale(sv_sel["os"]*0.66) 
   hists_E_ss["2016"]["higgs"].Scale(sv_sel["ss"]*0.33)
   hists_E_os["2016"]["higgs"].Scale(sv_sel["os"]*0.66)

##########################################
############# Adding all years ###########
##########################################

histT_M_ss = {}
histT_M_os = {}
histT_E_ss = {}
histT_E_os = {}
for s in samples:
   histT_M_ss[s] = hists_M_ss[datayears[0]][s]
   histT_M_os[s] = hists_M_os[datayears[0]][s]
   histT_E_ss[s] = hists_E_ss[datayears[0]][s]
   histT_E_os[s] = hists_E_os[datayears[0]][s]
   for d in datayears[1:]:
        histT_M_ss[s].Add(hists_M_ss[d][s])
        histT_M_os[s].Add(hists_M_os[d][s])
        histT_E_ss[s].Add(hists_E_ss[d][s])
        histT_E_os[s].Add(hists_E_os[d][s])

histT_M_ss["data"] = hists_M_ss[datayears[0]]["data"]
histT_M_os["data"] = hists_M_os[datayears[0]]["data"]
histT_E_ss["data"] = hists_E_ss[datayears[0]]["data"]
histT_E_os["data"] = hists_E_os[datayears[0]]["data"]
for d in datayears[1:]:
    histT_M_ss["data"].Add(hists_M_ss[d]["data"])
    histT_M_os["data"].Add(hists_M_os[d]["data"])
    histT_E_ss["data"].Add(hists_E_ss[d]["data"])
    histT_E_os["data"].Add(hists_E_os[d]["data"])

histT_M_ss["higgs"] = hists_M_ss["2016"]["higgs"]
histT_M_os["higgs"] = hists_M_os["2016"]["higgs"]
histT_E_ss["higgs"] = hists_E_ss["2016"]["higgs"]
histT_E_os["higgs"] = hists_E_os["2016"]["higgs"]

#################################
### Joining certain processes ###
#################################

histT_M_ss["vv"] = histT_M_ss["ww"]
histT_M_ss["vv"].Add(histT_M_ss["wz"])
histT_M_ss["vv"].Add(histT_M_ss["zz"])
histT_M_ss["wjets"] = histT_M_ss["wjets_1"]
histT_M_ss["wjets"].Add(histT_M_ss["wjets_2"])
histT_M_ss["wjets"].Add(histT_M_ss["wjets_3"])
histT_M_ss["wjets"].Add(histT_M_ss["wjets_4"])
histT_M_ss["wjets"].Add(histT_M_ss["wjets_5"])
histT_M_ss["wjets"].Add(histT_M_ss["wjets_6"])
histT_M_ss["wjets"].Add(histT_M_ss["wjets_7"])
histT_M_ss["wjets"].Add(histT_M_ss["wjets_8"])
histT_M_ss["zjets"] = histT_M_ss["zjets_1"]
histT_M_ss["zjets"].Add(histT_M_ss["zjets_2"])
histT_M_ss["zjets"].Add(histT_M_ss["zjets_3"])
histT_M_ss["zjets"].Add(histT_M_ss["zjets_4"])
histT_M_ss["zjets"].Add(histT_M_ss["zjets_5"])
histT_M_ss["zjets"].Add(histT_M_ss["zjets_6"])
histT_M_ss["zjets"].Add(histT_M_ss["zjets_7"])
histT_M_ss["zjets"].Add(histT_M_ss["zjets_8"])
histT_M_ss["st"] = histT_M_ss["st_1"]
histT_M_ss["st"].Add(histT_M_ss["st_2"])
histT_M_ss["st"].Add(histT_M_ss["st_3"])
histT_M_ss["st"].Add(histT_M_ss["st_4"])
histT_M_ss["ttbar_rest"] = histT_M_ss["ttbar_dl"]
histT_M_ss["ttbar_rest"].Add(histT_M_ss["ttbar_dh"])

histT_M_os["vv"] = histT_M_os["ww"]
histT_M_os["vv"].Add(histT_M_os["wz"])
histT_M_os["vv"].Add(histT_M_os["zz"])
histT_M_os["wjets"] = histT_M_os["wjets_1"]
histT_M_os["wjets"].Add(histT_M_os["wjets_2"])
histT_M_os["wjets"].Add(histT_M_os["wjets_3"])
histT_M_os["wjets"].Add(histT_M_os["wjets_4"])
histT_M_os["wjets"].Add(histT_M_os["wjets_5"])
histT_M_os["wjets"].Add(histT_M_os["wjets_6"])
histT_M_os["wjets"].Add(histT_M_os["wjets_7"])
histT_M_os["wjets"].Add(histT_M_os["wjets_8"])
histT_M_os["zjets"] = histT_M_os["zjets_1"]
histT_M_os["zjets"].Add(histT_M_os["zjets_2"])
histT_M_os["zjets"].Add(histT_M_os["zjets_3"])
histT_M_os["zjets"].Add(histT_M_os["zjets_4"])
histT_M_os["zjets"].Add(histT_M_os["zjets_5"])
histT_M_os["zjets"].Add(histT_M_os["zjets_6"])
histT_M_os["zjets"].Add(histT_M_os["zjets_7"])
histT_M_os["zjets"].Add(histT_M_os["zjets_8"])
histT_M_os["st"] = histT_M_os["st_1"]
histT_M_os["st"].Add(histT_M_os["st_2"])
histT_M_os["st"].Add(histT_M_os["st_3"])
histT_M_os["st"].Add(histT_M_os["st_4"])
histT_M_os["ttbar_rest"] = histT_M_os["ttbar_dl"]
histT_M_os["ttbar_rest"].Add(histT_M_os["ttbar_dh"])

histT_E_ss["vv"] = histT_E_ss["ww"]
histT_E_ss["vv"].Add(histT_E_ss["wz"])
histT_E_ss["vv"].Add(histT_E_ss["zz"])
histT_E_ss["wjets"] = histT_E_ss["wjets_1"]
histT_E_ss["wjets"].Add(histT_E_ss["wjets_2"])
histT_E_ss["wjets"].Add(histT_E_ss["wjets_3"])
histT_E_ss["wjets"].Add(histT_E_ss["wjets_4"])
histT_E_ss["wjets"].Add(histT_E_ss["wjets_5"])
histT_E_ss["wjets"].Add(histT_E_ss["wjets_6"])
histT_E_ss["wjets"].Add(histT_E_ss["wjets_7"])
histT_E_ss["wjets"].Add(histT_E_ss["wjets_8"])
histT_E_ss["zjets"] = histT_E_ss["zjets_1"]
histT_E_ss["zjets"].Add(histT_E_ss["zjets_2"])
histT_E_ss["zjets"].Add(histT_E_ss["zjets_3"])
histT_E_ss["zjets"].Add(histT_E_ss["zjets_4"])
histT_E_ss["zjets"].Add(histT_E_ss["zjets_5"])
histT_E_ss["zjets"].Add(histT_E_ss["zjets_6"])
histT_E_ss["zjets"].Add(histT_E_ss["zjets_7"])
histT_E_ss["zjets"].Add(histT_E_ss["zjets_8"])
histT_E_ss["st"] = histT_E_ss["st_1"]
histT_E_ss["st"].Add(histT_E_ss["st_2"])
histT_E_ss["st"].Add(histT_E_ss["st_3"])
histT_E_ss["st"].Add(histT_E_ss["st_4"])
histT_E_ss["ttbar_rest"] = histT_E_ss["ttbar_dl"]
histT_E_ss["ttbar_rest"].Add(histT_E_ss["ttbar_dh"])

histT_E_os["vv"] = histT_E_os["ww"]
histT_E_os["vv"].Add(histT_E_os["wz"])
histT_E_os["vv"].Add(histT_E_os["zz"])
histT_E_os["wjets"] = histT_E_os["wjets_1"]
histT_E_os["wjets"].Add(histT_E_os["wjets_2"])
histT_E_os["wjets"].Add(histT_E_os["wjets_3"])
histT_E_os["wjets"].Add(histT_E_os["wjets_4"])
histT_E_os["wjets"].Add(histT_E_os["wjets_5"])
histT_E_os["wjets"].Add(histT_E_os["wjets_6"])
histT_E_os["wjets"].Add(histT_E_os["wjets_7"])
histT_E_os["wjets"].Add(histT_E_os["wjets_8"])
histT_E_os["zjets"] = histT_E_os["zjets_1"]
histT_E_os["zjets"].Add(histT_E_os["zjets_2"])
histT_E_os["zjets"].Add(histT_E_os["zjets_3"])
histT_E_os["zjets"].Add(histT_E_os["zjets_4"])
histT_E_os["zjets"].Add(histT_E_os["zjets_5"])
histT_E_os["zjets"].Add(histT_E_os["zjets_6"])
histT_E_os["zjets"].Add(histT_E_os["zjets_7"])
histT_E_os["zjets"].Add(histT_E_os["zjets_8"])
histT_E_os["st"] = histT_E_os["st_1"]
histT_E_os["st"].Add(histT_E_os["st_2"])
histT_E_os["st"].Add(histT_E_os["st_3"])
histT_E_os["st"].Add(histT_E_os["st_4"])
histT_E_os["ttbar_rest"] = histT_E_os["ttbar_dl"]
histT_E_os["ttbar_rest"].Add(histT_E_os["ttbar_dh"])

##########################################
############# Change of names ############
##########################################

histT_E_ss["data"].SetName("mjj4j"+str(args.type)+"ss_1")
histT_M_ss["data"].SetName("mjj4j"+str(args.type)+"ss_2")
histT_E_ss["ttbar_sl_charm"].SetName("mjj4j"+str(args.type)+"ss_3")
histT_M_ss["ttbar_sl_charm"].SetName("mjj4j"+str(args.type)+"ss_4")
histT_E_ss["ttbar_rest"].SetName("mjj4j"+str(args.type)+"ss_5")
histT_M_ss["ttbar_rest"].SetName("mjj4j"+str(args.type)+"ss_6")
histT_E_ss["st"].SetName("mjj4j"+str(args.type)+"ss_7")
histT_M_ss["st"].SetName("mjj4j"+str(args.type)+"ss_8")
histT_E_ss["wjets"].SetName("mjj4j"+str(args.type)+"ss_9")
histT_M_ss["wjets"].SetName("mjj4j"+str(args.type)+"ss_10")
histT_E_ss["zjets"].SetName("mjj4j"+str(args.type)+"ss_11")
histT_M_ss["zjets"].SetName("mjj4j"+str(args.type)+"ss_12")
histT_E_ss["vv"].SetName("mjj4j"+str(args.type)+"ss_13")
histT_M_ss["vv"].SetName("mjj4j"+str(args.type)+"ss_14")
histT_E_ss["higgs"].SetName("mjj4j"+str(args.type)+"ss_15")
histT_M_ss["higgs"].SetName("mjj4j"+str(args.type)+"ss_16")
histT_E_ss["ttbar_sl_light"].SetName("mjj4j"+str(args.type)+"ss_17")
histT_M_ss["ttbar_sl_light"].SetName("mjj4j"+str(args.type)+"ss_18")
histT_E_ss["ttbar_sl_bottom"].SetName("mjj4j"+str(args.type)+"ss_19")
histT_M_ss["ttbar_sl_bottom"].SetName("mjj4j"+str(args.type)+"ss_20")

histT_E_os["data"].SetName("mjj4j"+str(args.type)+"os_1")
histT_M_os["data"].SetName("mjj4j"+str(args.type)+"os_2")
histT_E_os["ttbar_sl_charm"].SetName("mjj4j"+str(args.type)+"os_3")
histT_M_os["ttbar_sl_charm"].SetName("mjj4j"+str(args.type)+"os_4")
histT_E_os["ttbar_rest"].SetName("mjj4j"+str(args.type)+"os_5")
histT_M_os["ttbar_rest"].SetName("mjj4j"+str(args.type)+"os_6")
histT_E_os["st"].SetName("mjj4j"+str(args.type)+"os_7")
histT_M_os["st"].SetName("mjj4j"+str(args.type)+"os_8")
histT_E_os["wjets"].SetName("mjj4j"+str(args.type)+"os_9")
histT_M_os["wjets"].SetName("mjj4j"+str(args.type)+"os_10")
histT_E_os["zjets"].SetName("mjj4j"+str(args.type)+"os_11")
histT_M_os["zjets"].SetName("mjj4j"+str(args.type)+"os_12")
histT_E_os["vv"].SetName("mjj4j"+str(args.type)+"os_13")
histT_M_os["vv"].SetName("mjj4j"+str(args.type)+"os_14")
histT_E_os["higgs"].SetName("mjj4j"+str(args.type)+"os_15")
histT_M_os["higgs"].SetName("mjj4j"+str(args.type)+"os_16")
histT_E_os["ttbar_sl_light"].SetName("mjj4j"+str(args.type)+"os_17")
histT_M_os["ttbar_sl_light"].SetName("mjj4j"+str(args.type)+"os_18")
histT_E_os["ttbar_sl_bottom"].SetName("mjj4j"+str(args.type)+"os_19")
histT_M_os["ttbar_sl_bottom"].SetName("mjj4j"+str(args.type)+"os_20")

samplesF = ["data","ttbar_sl_charm","ttbar_rest","st","wjets","zjets","vv","higgs","ttbar_sl_light","ttbar_sl_bottom"]
 
if args.type == "sl":
  path_ss = "/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/hm"+str(args.higgsmass)+"/mjj4jslss.root"
  path_os = "/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/hm"+str(args.higgsmass)+"/mjj4jslos.root"
else:
  path_ss = "/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/hm"+str(args.higgsmass)+"/mjj4jsvss.root"
  path_os = "/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/hm"+str(args.higgsmass)+"/mjj4jsvos.root"

myfile_ss = TFile( path_ss , 'RECREATE')
for s in samplesF:
  histT_E_ss[s].Write()
  histT_M_ss[s].Write()
myfile_ss.Close()

myfile_os = TFile( path_os , 'RECREATE')
for s in samplesF:
  histT_E_os[s].Write()
  histT_M_os[s].Write()
myfile_os.Close()
