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
parser.add_argument("--norm", action="store_true", default=False,
                    help="normalising test")
parser.add_argument("--nodata", action="store_true", default=False,
                    help="Do not plot data")
parser.add_argument("--wcs", action="store_true", default=False,
                    help="classification for ttbar and st")
parser.add_argument("--year", type=string, default="2016",
                    help="Select year of process to run")
parser.add_argument("--channel", type=string, default="btagMM",
                    help="Select year of process to run")
parser.add_argument("--nosyst", action="store_true", default=False,
                    help="systematics inclusion")
parser.add_argument("--sumEM", action="store_true", default=False,
                    help="union of M and E channels")

# Use like:
# python higgssearch/fromJF/hist_fromJFwqq_syst_total.py --stack --ratio --png --norm --year="all" --channel="btagMM_chitest_sl"

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

if args.channel == "btagMM_chitest": term_path = ""
elif args.channel == "btagMM_nochitest": term_path = "/nochitest"
elif args.channel == "btagMM_chitest_auxctag": term_path = "/aux_ctag"
elif args.channel == "btagMM_chitest_ctag": term_path = "/ctag"
elif args.channel == "btagMM_chitest_noctag": term_path = "/noctag"
elif args.channel == "btagMM_chitest_antisl": term_path = "/antisl"
elif args.channel == "btagMM_chitest_sl": term_path = "/sl"
elif args.channel == "btagMM_chitest_slss": term_path = "/sl/ss"
elif args.channel == "btagMM_chitest_slos": term_path = "/sl/os"
elif args.channel == "btagMM_chitest_slssos": term_path = "/sl/ssos"
else: raise NameError('Incorrect data option')

sl_channel = ["btagMM_chitest_sl","lepton50_chitest_sl","btagMM_chitest_slss","lepton50_chitest_slss",
      "btagMM_chitest_slos","lepton50_chitest_slos","btagMM_chitest_slssos","lepton50_chitest_slssos"]

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
mrk_size = 6; titY_off = 0.3; tit_size = 0.14; titX_off = 0.7; titY_size=0.11;
titY_off2 = 1.0; titY_size2=0.05;
if not args.nodata:
   leg1 = 0.4; leg2 = 0.7; leg3 = 0.96; leg4 = 0.92;
else:
   leg1 = 0.4; leg2 = 0.76; leg3 = 0.92; leg4 = 0.9;

if not args.ratio: titY_size2=0.03; titY_off2 = 1.6;
if args.nodata: titY_size = 0.07;titY_off=0.5;

##########################

norm_factorM = {}
norm_factorE = {}

norm_factorM["all"] = {}; norm_factorM["2016"] = {}; norm_factorM["2016B"] = {}; norm_factorM["2017"] = {}; norm_factorM["2018"] = {};
norm_factorE["all"] = {}; norm_factorE["2016"] = {}; norm_factorE["2016B"] = {}; norm_factorE["2017"] = {}; norm_factorE["2018"] = {};

### btagMM
norm_factorM["all"]["btagMM"] = 0.95;
norm_factorE["all"]["btagMM"] = 0.94;

### btagMM smeared, chitest
norm_factorM["all"]["btagMM_chitest"] = 0.915;
norm_factorE["all"]["btagMM_chitest"] = 0.894;

norm_factorM["2016"]["btagMM_chitest"] = 0.89;norm_factorE["2016"]["btagMM_chitest"] = 0.87;
norm_factorM["2016B"]["btagMM_chitest"] = 0.92;norm_factorE["2016B"]["btagMM_chitest"] = 0.94;
norm_factorM["2017"]["btagMM_chitest"] = 0.92;norm_factorE["2017"]["btagMM_chitest"] = 0.88;
norm_factorM["2018"]["btagMM_chitest"] = 0.91;norm_factorE["2018"]["btagMM_chitest"] = 0.90;

### btagMM smeared, chitest, SL
norm_factorM["all"]["btagMM_chitest_sl"] = 0.92; norm_factorM["all"]["btagMM_chitest_slss"] = 0.92; 
norm_factorM["all"]["btagMM_chitest_slos"] = 0.92; norm_factorM["all"]["btagMM_chitest_slssos"] = 0.92;
norm_factorE["all"]["btagMM_chitest_sl"] = 0.90; norm_factorE["all"]["btagMM_chitest_slss"] = 0.90; 
norm_factorE["all"]["btagMM_chitest_slos"] = 0.90; norm_factorE["all"]["btagMM_chitest_slssos"] = 0.90;

norm_factorM["all"]["btagMM_chitest_ctag"] = 0.92; norm_factorE["all"]["btagMM_chitest_ctag"] = 0.9; 
norm_factorM["all"]["btagMM_chitest_noctag"] = 0.92; norm_factorE["all"]["btagMM_chitest_noctag"] = 0.9; 
norm_factorM["all"]["btagMM_chitest_auxctag"] = 0.92; norm_factorE["all"]["btagMM_chitest_auxctag"] = 0.9; 
norm_factorM["all"]["btagMM_chitest_antisl"] = 0.92; norm_factorE["all"]["btagMM_chitest_antisl"] = 0.9; 

########### pT jet not corrected
norm_factorM["all"]["btagMM_chitest"] = 0.930;
norm_factorE["all"]["btagMM_chitest"] = 0.909;
norm_factorM["2016"]["btagMM_chitest"] = 0.91;norm_factorE["2016"]["btagMM_chitest"] = 0.89;
norm_factorM["2016B"]["btagMM_chitest"] = 0.93;norm_factorE["2016B"]["btagMM_chitest"] = 0.95;
norm_factorM["2017"]["btagMM_chitest"] = 0.94;norm_factorE["2017"]["btagMM_chitest"] = 0.89;
norm_factorM["2018"]["btagMM_chitest"] = 0.93;norm_factorE["2018"]["btagMM_chitest"] = 0.92;

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
   "jet_bot1_tracks","jet_bot2_tracks","tau_discr_jet1","tau_discr_jet2","tau_discr_jetbot1","tau_discr_jetbot2"]


if ("sl" in str(args.channel)) and not("anti" in str(args.channel)): observable_names = observable_names + ["muon_jet_pt","muon_jet_z","muon_jet_eta",
         "muon_jet_pt_rel","muon_jet_iso","muon_jet_xy","muon_jet_dz",
         "muon_jet_iso_log","muon_jet_z_short","InvM3_good_short","muon_jet_sigr","muon_jet_sigxy","muon_jet_sigdz","muon_jet_r","deltaR_jet1_muon",
         "muon_jet_z2_v2","muon_jet_z3","muon_jet_iso_abs","jet_muon_pt","jet_muon_eta","jet_muon_btag","jet_muon_btagnumber",
         "jet_notmuon_pt","jet_notmuon_eta","jet_notmuon_btag","jet_notmuon_btagnumber","jet_muon_nmu","jet_notmuon_nmu"]

if ("ctag" in str(args.channel)): observable_names = observable_names + ["jet_max_pt","jet_max_eta","jet_min_pt","jet_min_eta"]
else: observable_names = observable_names + ["second_muon_pt","second_el_pt"]

observable_names = ["jet_bot1_btag", "jet_bot2_btag"]

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
cosmetic_names["muon_jet_iso"] = "I_{PF}/p_{T}^{#mu}";

### Ptmiss, MT, pt y eta de los bjets, btag discriminant cvl for the b-jets, pt y eta
### de los dos jets del W, deltaR(jet1W, jet2W), number of b-tag jets, deltaphi (lepton, MET), pt
### leptonic W

not_rebin = ["nJetGood","InvM3_bad","InvMl_good","InvMl_bad","lepton_eta_thick","jet_bot1_btagnumber", "jet_bot2_btagnumber", 
      "jet_1_btagnumber", "jet_2_btagnumber","muon_jet_relpt","muon_jet_eta","tau_discr","transverse_mass","muon_jet_pt",
      "jet_1_flavourP", "jet_2_flavourP", "jet_bot1_flavourP", "jet_bot2_flavourP","lepton_pt", "muon_jet_z","tau_discr_jet1","tau_discr_jet2","tau_discr_jetbot1",
      "tau_discr_jetbot2","muon_jet_iso","muon_jet_iso_abs","InvM3_good_short","jet_muon_btagnumber","jet_notmuon_btagnumber","MET_pt_aux","MET_my_sig",
      "jet_1_nmu","jet_2_nmu","jet_muon_nmu","jet_notmuon_nmu","ttsl_lepflav","ttdl_lepflav"]

#observable_names = ["jet_max_cvltag","jet_max_cvbtag"]
#observable_names = ["jet_1_btagnumber", "jet_2_btagnumber"]
#observable_names = ["nJetGood","chi2_test_good"]

datayears = ["2016","2016B","2017","2018"]
#datayears = ["2018","2016","2016B"]

samplesHT = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_light","ttbar_sl_bottom","ttbar_dl","ttbar_dh","zz","wz",
        "st_1","st_2","st_3","st_4", "ttbar_sl_else", "ttbar_sl_charmgluon","ttbar_sl_bottomgluon"]

if args.wcs:
   samplesHT = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_nocharm","ttbar_dl","ttbar_dh","zz","wz",
        "st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm",
        "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else"]

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
        lumi[data_op]["ttbar_sl_nocharm"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_light"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_bottom"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_else"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_charmgluon"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["ttbar_sl_bottomgluon"] = lumi[data_op]["ttbar_sl"]
        lumi[data_op]["st_1_charm"] = lumi[data_op]["st_1"]
        lumi[data_op]["st_2_charm"] = lumi[data_op]["st_2"]
        lumi[data_op]["st_3_charm"] = lumi[data_op]["st_3"]
        lumi[data_op]["st_4_charm"] = lumi[data_op]["st_4"]
        lumi[data_op]["st_1_nocharm"] = lumi[data_op]["st_1"]
        lumi[data_op]["st_2_nocharm"] = lumi[data_op]["st_2"]
        lumi[data_op]["st_3_nocharm"] = lumi[data_op]["st_3"]
        lumi[data_op]["st_4_nocharm"] = lumi[data_op]["st_4"]
        lumi[data_op]["st_1_else"] = lumi[data_op]["st_1"]
        lumi[data_op]["st_2_else"] = lumi[data_op]["st_2"]
        lumi[data_op]["st_3_else"] = lumi[data_op]["st_3"]
        lumi[data_op]["st_4_else"] = lumi[data_op]["st_4"]

listsampl = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "ttbar_sl","ttbar_dl","ttbar_dh","zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6",
        "zjets_7","zjets_8","st_1","st_2","st_3","st_4","zz","wz", "ttbar_sl_charm", "ttbar_sl_charmgluon",
        "ttbar_sl_bottom", "ttbar_sl_bottomgluon", "ttbar_sl_else", "ttbar_sl_light", "ttbar_sl_nocharm",
        "st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm",
        "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else"]

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

histFile = {}
histFileDM = {}
histFileDE = {}

for name in observable_names:
  ## Open hists files
  filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/chi_test"+term_path+"/"
  if args.wcs: 
     #filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/chi_test/wcs_classes_aux/test_aux"+term_path+"/"
     #filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/chi_test/wcs_classes_ver2"+term_path+"/"
     filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/chi_test/wcs_classes"+term_path+"/"
     #filePath = "/nfs/cms/vazqueze/new_hists/fromJF/wqq/btagMM/chi_test/wcs_classes/test_pt27"+term_path+"/"
  term = "hist_wqqfromJF_"
  end_term = ".root"
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
  #print(data_op)
  #print(histFile[data_op].keys())
  #print(histFileDM[data_op].keys())

  if args.png: c1 = TCanvas("c1","",1200,1000)
  else: c1 = TCanvas("c1","",600,400)

  if args.ratio:
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

  samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl_charm","ttbar_sl_light","ttbar_sl_bottom","ttbar_dl","ttbar_dh","zz","wz",
        "st_1","st_2","st_3","st_4","ttbar_sl_else","ttbar_sl_charmgluon","ttbar_sl_bottomgluon"]

  if args.wcs: 
    samples = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
            "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
            "ttbar_sl_charm","ttbar_sl_nocharm","ttbar_dl","ttbar_dh","zz","wz",
            "st_1_charm","st_2_charm","st_3_charm","st_4_charm","st_1_nocharm","st_2_nocharm","st_3_nocharm",
            "st_4_nocharm","st_1_else","st_2_else","st_3_else","st_4_else"]

  ## HISTS
  samples_foryear = {}
  hist_nom_M = {}
  hist_nom_E = {}
  hist_syst_M = {}
  hist_syst_E = {}
  #list_syst = ["btaglightup","btaglightdown","btagheavyup","btagheavydown","seclepup","seclepdown"]
  list_syst = ["ctagup","ctagdown"]
  for syst in list_syst:
      hist_syst_M[syst] = {}
      hist_syst_E[syst] = {}
  histdata_M = {}
  histdata_E = {}
  for data_op in datayears:
    hist_nom_M[data_op] = {}
    hist_nom_E[data_op] = {}
    for syst in list_syst:
       hist_syst_M[syst][data_op] = {}
       hist_syst_E[syst][data_op] = {}
    data_term = data_op
    #print(data_op+"M_"+name+"_M")
    #print(histFile[data_op][name].ls())
    for s in samples:
      if s[0:8] == "ttbar_sl": 
         s_term = s[0:8]+data_term+s[8:] 
      elif s[0:2] == "st":
         s_term = s[0:4]+data_term+s[4:] 
      else:
         s_term = s+data_term
      hist_nom_M[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M")
      hist_nom_E[data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E")
      if not args.nosyst:
         for syst in list_syst:
             hist_syst_M[syst][data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_M_"+str(syst))
             hist_syst_E[syst][data_op][s] = histFile[name][data_op].Get(s_term+"_"+name+"_E_"+str(syst))
    if not args.nodata: histdata_M[data_op] = histFileDM[name][data_op].Get("data"+data_op+"M_"+name+"_M")
    if not args.nodata: histdata_E[data_op] = histFileDE[name][data_op].Get("data"+data_op+"E_"+name+"_E")

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
      hist_nom_M[data_op][s].Scale(lumi_data/lumi[data_op][s])
      hist_nom_E[data_op][s].Scale(lumi_data/lumi[data_op][s])
      if not args.nosyst:
         for syst in list_syst:
           hist_syst_M[syst][data_op][s].Scale(lumi_data/lumi[data_op][s])
           hist_syst_E[syst][data_op][s].Scale(lumi_data/lumi[data_op][s])
      if args.norm: hist_nom_M[data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
      if args.norm: hist_nom_E[data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
      if not args.nosyst:
         for syst in list_syst:
             if args.norm: hist_syst_M[syst][data_op][s].Scale(norm_factorM[str(args.year)][str(args.channel)])
             if args.norm: hist_syst_E[syst][data_op][s].Scale(norm_factorE[str(args.year)][str(args.channel)])
    ## Fixing single top
    #print(samples_foryear[data_op]) 
    #### List of summing samples:
    list_st = ["st_1","st_2","st_3","st_4"]
    list_wjets = ["wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8"]
    list_zjets = ["zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8"]
    list_vv = ["ww","wz","zz"]
    if args.wcs:
       list_st_aux = ["_charm","_nocharm","_else"]
    else:
       list_st_aux = [""]
    for apx in list_st_aux:
      hist_nom_M[data_op]["st"+apx] = hist_nom_M[data_op][list_st[0]+apx]
      hist_nom_E[data_op]["st"+apx] = hist_nom_E[data_op][list_st[0]+apx]
      if not args.nosyst:
         for syst in list_syst:
              hist_syst_M[syst][data_op]["st"+apx] = hist_syst_M[syst][data_op][list_st[0]+apx]
              hist_syst_E[syst][data_op]["st"+apx] = hist_syst_E[syst][data_op][list_st[0]+apx]
      for l in list_st[1:]:
        hist_nom_M[data_op]["st"+apx].Add(hist_nom_M[data_op][l+apx])
        hist_nom_E[data_op]["st"+apx].Add(hist_nom_E[data_op][l+apx])
        if not args.nosyst:
           for syst in list_syst:
              hist_syst_M[syst][data_op]["st"+apx].Add(hist_syst_M[syst][data_op][l+apx])
              hist_syst_E[syst][data_op]["st"+apx].Add(hist_syst_E[syst][data_op][l+apx])
    hist_nom_M[data_op]["wjets"] = hist_nom_M[data_op][list_wjets[0]]
    hist_nom_E[data_op]["wjets"] = hist_nom_E[data_op][list_wjets[0]]
    if not args.nosyst:
       for syst in list_syst:
            hist_syst_M[syst][data_op]["wjets"] = hist_syst_M[syst][data_op][list_wjets[0]]
            hist_syst_E[syst][data_op]["wjets"] = hist_syst_E[syst][data_op][list_wjets[0]]
    for l in list_wjets[1:]:
      hist_nom_M[data_op]["wjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["wjets"].Add(hist_nom_E[data_op][l])
      if not args.nosyst:
         for syst in list_syst:
            hist_syst_M[syst][data_op]["wjets"].Add(hist_syst_M[syst][data_op][l])
            hist_syst_E[syst][data_op]["wjets"].Add(hist_syst_E[syst][data_op][l])
    hist_nom_M[data_op]["zjets"] = hist_nom_M[data_op][list_zjets[0]]
    hist_nom_E[data_op]["zjets"] = hist_nom_E[data_op][list_zjets[0]]
    if not args.nosyst:
       for syst in list_syst:
            hist_syst_M[syst][data_op]["zjets"] = hist_syst_M[syst][data_op][list_zjets[0]]
            hist_syst_E[syst][data_op]["zjets"] = hist_syst_E[syst][data_op][list_zjets[0]]
    for l in list_zjets[1:]:
      hist_nom_M[data_op]["zjets"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["zjets"].Add(hist_nom_E[data_op][l])
      if not args.nosyst:
         for syst in list_syst:
            hist_syst_M[syst][data_op]["zjets"].Add(hist_syst_M[syst][data_op][l])
            hist_syst_E[syst][data_op]["zjets"].Add(hist_syst_E[syst][data_op][l])
    hist_nom_M[data_op]["vv"] = hist_nom_M[data_op][list_vv[0]]
    hist_nom_E[data_op]["vv"] = hist_nom_E[data_op][list_vv[0]]
    if not args.nosyst:
       for syst in list_syst:
            hist_syst_M[syst][data_op]["vv"] = hist_syst_M[syst][data_op][list_vv[0]]
            hist_syst_E[syst][data_op]["vv"] = hist_syst_E[syst][data_op][list_vv[0]]
    for l in list_vv[1:]:
      hist_nom_M[data_op]["vv"].Add(hist_nom_M[data_op][l])
      hist_nom_E[data_op]["vv"].Add(hist_nom_E[data_op][l])
      if not args.nosyst:
         for syst in list_syst:
            hist_syst_M[syst][data_op]["vv"].Add(hist_syst_M[syst][data_op][l])
            hist_syst_E[syst][data_op]["vv"].Add(hist_syst_E[syst][data_op][l])

  samples = ["ttbar_sl_bottom","ttbar_sl_charm","ttbar_sl_else","ttbar_sl_light","ttbar_dl","ttbar_dh","zjets","vv","st","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon"]
  if args.wcs: samples = ["ttbar_sl_nocharm","ttbar_sl_charm","ttbar_dl","ttbar_dh","zjets","vv","st_nocharm","st_charm","st_else","wjets"]

  ## sumamos todos los anos para todos los histogramas
  histT_nom_M = {}
  histT_nom_E = {}
  histT_syst_M = {}
  histT_syst_E = {}
  for syst in list_syst:
      histT_syst_M[syst] = {}
      histT_syst_E[syst] = {}
  for s in samples:
       histT_nom_M[s] = hist_nom_M[datayears[0]][s]
       histT_nom_E[s] = hist_nom_E[datayears[0]][s]
       if not args.nosyst:
          for syst in list_syst:
              histT_syst_M[syst][s] =  hist_syst_M[syst][datayears[0]][s]
              histT_syst_E[syst][s] =  hist_syst_E[syst][datayears[0]][s]
       for d in datayears[1:]:
          histT_nom_M[s].Add(hist_nom_M[d][s])
          histT_nom_E[s].Add(hist_nom_E[d][s])
          if not args.nosyst:
             for syst in list_syst:
                 histT_syst_M[syst][s].Add(hist_syst_M[syst][d][s])
                 histT_syst_E[syst][s].Add(hist_syst_E[syst][d][s])
       if (args.channel in sl_channel) and (name not in not_rebin):
          histT_nom_M[s].Rebin(nrebin)
          histT_nom_E[s].Rebin(nrebin)
          if not args.nosyst:
             for syst in list_syst:
                 histT_syst_M[syst][s].Rebin(nrebin)
                 histT_syst_E[syst][s].Rebin(nrebin)

  if not args.nodata:
    histD_M = histdata_M[datayears[0]]
    histD_E = histdata_E[datayears[0]]
    for d in datayears[1:]:
       histD_M.Add(histdata_M[d])
       histD_E.Add(histdata_E[d])
    if (args.channel in sl_channel) and (name not in not_rebin): 
       histD_M.Rebin(nrebin)
       histD_E.Rebin(nrebin)

  if not args.nosyst:
    histT_sT_syst_M = {}
    histT_sT_syst_E = {}
    ### Para los histogramas de sistemática up y down sumamos todas las contribuciones
    for syst in list_syst:
        histT_sT_syst_M[syst] = histT_syst_M[syst][samples[0]]
        histT_sT_syst_E[syst] = histT_syst_E[syst][samples[0]]
        for s in samples[1:]:
            histT_sT_syst_M[syst].Add(histT_syst_M[syst][s])
            histT_sT_syst_E[syst].Add(histT_syst_E[syst][s])

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  colors = {}
  colors["ttbar_dl"] = (222,90,106)
  colors["wjets"] = (155,152,204)
  colors["ttbar_sl_bottomgluon"] = (236,164,207)
  colors["ttbar_sl_charm"] = (204,255,153)
  colors["ttbar_sl_light"] = (120,154,86)
  colors["ttbar_sl_nocharm"] = (120,154,86)
  colors["ttbar_sl_bottom"] = (201,79,152)
  colors["ttbar_sl_else"] = (242,193,121)
  colors["ttbar_sl_charmgluon"] = (222,212,74)
  colors["vv"] = (255,180,85)
  colors["ttbar_dh"] = (204,204,0)
  colors["zjets"] = (113,209,223)
  colors["st"] = (153,51,255)
  colors["st_charm"] = (102,0,204)
  colors["st_nocharm"] = (198,101,222)
  colors["st_else"] = (207,176,235)

  if args.stack:
    ymax_M = 0
    ymin_M = 0
    ymax_E = 0
    ymin_E = 0
    for s in samples:
      histT_nom_M[s].SetLineWidth(1)
      histT_nom_M[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_nom_M[s].GetYaxis().SetTitle("Events")
      histT_nom_M[s].GetXaxis().SetTitle(name)
      histT_nom_M[s].GetYaxis().SetTitleOffset(titY_off2)
      histT_nom_M[s].GetYaxis().SetTitleSize(titY_size2)
      histT_nom_E[s].SetLineWidth(1)
      histT_nom_E[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      histT_nom_E[s].GetYaxis().SetTitle("Events")
      histT_nom_E[s].GetXaxis().SetTitle(name)
      histT_nom_E[s].GetYaxis().SetTitleOffset(titY_off2)
      histT_nom_E[s].GetYaxis().SetTitleSize(titY_size2)

      y = histT_nom_M[s].GetMaximum()
      ym = histT_nom_M[s].GetMinimum()
      if y>ymax_M: ymax_M=y
      if ymin_M>ym: ymin_M=ym
      if not args.nodata:
        y = histD_M.GetMaximum()
        if y>ymax_M: ymax_M=y

      y = histT_nom_E[s].GetMaximum()
      ym = histT_nom_E[s].GetMinimum()
      if y>ymax_E: ymax_E=y
      if ymin_E>ym: ymin_E=ym
      if not args.nodata:
        y = histD_E.GetMaximum()
        if y>ymax_E: ymax_E=y

      #histT_M["vv"].SetMinimum(1000.)
      #histT_M["vv"].SetMaximum(3*ymax_M)
      if args.linear: histT_nom_M["vv"].SetMaximum(1.3*ymax_M)
      if args.linear: histT_nom_M["vv"].SetMinimum(1.3*ymin_M)

      #histT_E["vv"].SetMinimum(1000.)
      #histT_E["vv"].SetMaximum(3*ymax_E)
      if args.linear: histT_nom_E["vv"].SetMaximum(1.3*ymax_E)
      if args.linear: histT_nom_E["vv"].SetMinimum(1.3*ymin_E)

    ## Stack creation
    samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","ttbar_sl_bottomgluon","ttbar_sl_charmgluon","ttbar_sl_else","ttbar_sl_bottom","st","ttbar_sl_light","ttbar_sl_charm"]
    if args.wcs: samples = ["vv","ttbar_dl","ttbar_dh","zjets","wjets","st_else","st_nocharm","st_charm","ttbar_sl_nocharm","ttbar_sl_charm"]

    if args.ratio: upper_pad.cd()

    if args.sumEM:
       histT_nom_T = {}
       for s in samples:
         histT_nom_T[s] = histT_nom_M[s]
         histT_nom_T[s].Add(histT_nom_E[s])
       #stack_T = ROOT.THStack("hs", ";;Events")
       stack_T = ROOT.THStack()
       for s in samples:
         if ("jet" in name) and ("pt" in name) and not("ssos" in str(args.channel)):
            histT_nom_T[s].SetBinContent(5,0.)
         stack_T.Add(histT_nom_T[s])
       last_T = stack_T.GetStack().Last()
       termfile = "/nfs/cms/vazqueze/higgssearch/fromJF/file_aux_jet_btag_allcorrections.root"
       myfile = TFile(termfile, "RECREATE")
       hist_aux_sum = histT_nom_T[samples[0]].Clone("sumMC")
       for s in samples[1:]:
           hist_aux_sum.Add(histT_nom_T[s])
       hist_aux_sum.Write()
       myfile.Close()           
    else:
       stack_M = ROOT.THStack()
       stack_E = ROOT.THStack()
       for s in samples:
         if ("jet" in name) and ("pt" in name) and not("ssos" in str(args.channel)):
            histT_nom_M[s].SetBinContent(5,0.)
            histT_nom_E[s].SetBinContent(5,0.)
         stack_M.Add(histT_nom_M[s])
         stack_E.Add(histT_nom_E[s])

       last_M = stack_M.GetStack().Last()
       last_E = stack_E.GetStack().Last()

    if args.sumEM:
      if not args.nosyst:
        envHi_T = {}
        ### Creamos el tgraph de sistematicos
        for syst in list_syst:
            envHi_T[syst] =  histT_sT_syst_M[syst]
            envHi_T[syst].Add(histT_sT_syst_E[syst])
        graph_err_T = TGraphAsymmErrors();
        ratio_graph_err_T = TGraphAsymmErrors();

        for bin in range(last_T.GetNbinsX()):
          bin_nom_T = float(last_T.GetBinContent(bin+1));
          bin_T = {}
          sysup_T = {}
          sysdown_T = {}
          for syst in list_syst:
              bin_T[syst] = float(envHi_T[syst].GetBinContent(bin+1));
              if syst[-2:] == "up":
                sysup_T[syst] = abs(bin_nom_T-bin_T[syst]);
              else:
                sysdown_T[syst] = abs(bin_nom_T-bin_T[syst]);
          sysup_total_T = 0;
          sysdown_total_T = 0;
          for syst in list_syst:
          #for syst in ["seclepup","seclepdown"]:
             if syst[-2:] == "up":
                sysup_total_T = sysup_total_T+sysup_T[syst]**2;
             else:
                sysdown_total_T = sysdown_total_T+sysdown_T[syst]**2;
          sysup_total_T = sqrt(sysup_total_T);
          sysdown_total_T = sqrt(sysdown_total_T);

          graph_err_T.SetPointEXhigh(bin,last_T.GetBinWidth(bin+1)/2);
          graph_err_T.SetPointEXlow(bin,last_T.GetBinWidth(bin+1)/2);
          graph_err_T.SetPoint(bin,last_T.GetBinCenter(bin+1),last_T.GetBinContent(bin+1));
          graph_err_T.SetPointEYhigh(bin,sysup_total_T);
          graph_err_T.SetPointEYlow(bin,sysdown_total_T);

          ratio_graph_err_T.SetPointEXhigh(bin,last_T.GetBinWidth(bin+1)/2);
          ratio_graph_err_T.SetPointEXlow(bin,last_T.GetBinWidth(bin+1)/2);
          ratio_graph_err_T.SetPoint(bin,last_T.GetBinCenter(bin+1),1.);
          if (last_T.GetBinContent(bin+1)>0.01):
             ratio_graph_err_T.SetPointEYhigh(bin,abs(sysup_total_T/bin_nom_T));
             ratio_graph_err_T.SetPointEYlow(bin,abs(sysdown_total_T/bin_nom_T));
          else:
             ratio_graph_err_T.SetPointEYhigh(bin,0.);
             ratio_graph_err_T.SetPointEYlow(bin,0.);

        graph_err_T.SetFillColorAlpha(kGray+2,0.4);
        graph_err_T.SetMarkerColor(kGray);
        graph_err_T.SetMarkerSize(0);
        graph_err_T.SetFillStyle(3001);
        ratio_graph_err_T.SetFillColorAlpha(kGray+2,0.4);
        ratio_graph_err_T.SetMarkerColor(kGray);
        ratio_graph_err_T.SetMarkerSize(0);
        ratio_graph_err_T.SetFillStyle(3001);
    else:
      if not args.nosyst:
        envHi_M = {}
        envHi_E = {}
        ### Creamos el tgraph de sistematicos
        for syst in list_syst:
            envHi_M[syst] = histT_sT_syst_M[syst]
            envHi_E[syst] = histT_sT_syst_E[syst]

        graph_err_M = TGraphAsymmErrors();
        graph_err_E = TGraphAsymmErrors();
        ratio_graph_err_M = TGraphAsymmErrors();
        ratio_graph_err_E = TGraphAsymmErrors();

        for bin in range(last_T.GetNbinsX()):
          bin_nom_M = float(last_M.GetBinContent(bin+1));
          bin_nom_E = float(last_E.GetBinContent(bin+1));
          bin_M = {}
          bin_E = {}
          sysup_M = {}
          sysdown_M = {}
          sysup_E = {}
          sysdown_E = {}
          for syst in list_syst:
              bin_M[syst] = float(envHi_M[syst].GetBinContent(bin+1));
              bin_E[syst] = float(envHi_E[syst].GetBinContent(bin+1));
              if syst[-2:] == "up":
                sysup_M[syst] = abs(bin_nom_M-bin_M[syst]);
                sysup_E[syst] = abs(bin_nom_E-bin_E[syst]);
              else:
                sysdown_M[syst] = abs(bin_nom_M-bin_M[syst]);
                sysdown_E[syst] = abs(bin_nom_E-bin_E[syst]);
          sysup_total_M = 0;
          sysdown_total_M = 0;
          sysup_total_E = 0;
          sysdown_total_E = 0;
          #for syst in list_syst:
          for syst in ["seclepup","seclepdown"]:
             if syst[-2:] == "up":
                sysup_total_M = sysup_total_T+sysup_M[syst]**2;
                sysup_total_E = sysup_total_T+sysup_E[syst]**2;
             else:
                sysdown_total_M = sysdown_total_M+sysdown_M[syst]**2;
                sysdown_total_E = sysdown_total_E+sysdown_E[syst]**2;
          sysup_total_M = sqrt(sysup_total_M);
          sysdown_total_E = sqrt(sysdown_total_E);
          sysup_total_M = sqrt(sysup_total_M);
          sysdown_total_E = sqrt(sysdown_total_E);

          graph_err_M.SetPointEXhigh(bin,last_M.GetBinWidth(bin+1)/2);
          graph_err_M.SetPointEXlow(bin,last_M.GetBinWidth(bin+1)/2);
          graph_err_M.SetPoint(bin,last_M.GetBinCenter(bin+1),last_M.GetBinContent(bin+1));
          graph_err_M.SetPointEYhigh(bin,sysup_total_M);
          graph_err_M.SetPointEYlow(bin,sysdown_total_M);

          graph_err_E.SetPointEXhigh(bin,last_E.GetBinWidth(bin+1)/2);
          graph_err_E.SetPointEXlow(bin,last_E.GetBinWidth(bin+1)/2);
          graph_err_E.SetPoint(bin,last_E.GetBinCenter(bin+1),last_E.GetBinContent(bin+1));
          graph_err_E.SetPointEYhigh(bin,sysup_total_E);
          graph_err_E.SetPointEYlow(bin,sysdown_total_E);

          ratio_graph_err_M.SetPointEXhigh(bin,last_M.GetBinWidth(bin+1)/2);
          ratio_graph_err_M.SetPointEXlow(bin,last_M.GetBinWidth(bin+1)/2);
          ratio_graph_err_M.SetPoint(bin,last_M.GetBinCenter(bin+1),1.);
          if (last_M.GetBinContent(bin+1)>0.01):
             ratio_graph_err_M.SetPointEYhigh(bin,abs(sysup_total_M/bin_nom_M));
             ratio_graph_err_M.SetPointEYlow(bin,abs(sysdown_total_M/bin_nom_M));
          else:
             ratio_graph_err_M.SetPointEYhigh(bin,0.);
             ratio_graph_err_M.SetPointEYlow(bin,0.);

          ratio_graph_err_E.SetPointEXhigh(bin,last_E.GetBinWidth(bin+1)/2);
          ratio_graph_err_E.SetPointEXlow(bin,last_E.GetBinWidth(bin+1)/2);
          ratio_graph_err_E.SetPoint(bin,last_E.GetBinCenter(bin+1),1.);
          if (last_E.GetBinContent(bin+1)>0.01):
             ratio_graph_err_E.SetPointEYhigh(bin,abs(sysup_total_E/bin_nom_E));
             ratio_graph_err_E.SetPointEYlow(bin,abs(sysdown_total_E/bin_nom_E));
          else:
             ratio_graph_err_E.SetPointEYhigh(bin,0.);
             ratio_graph_err_E.SetPointEYlow(bin,0.);

        graph_err_M.SetFillColorAlpha(kGray+2,0.7);
        graph_err_M.SetMarkerColor(kGray);
        graph_err_M.SetMarkerSize(0);
        graph_err_M.SetFillStyle(3001);
        graph_err_E.SetFillColorAlpha(kGray+2,0.7);
        graph_err_E.SetMarkerColor(kGray);
        graph_err_E.SetMarkerSize(0);
        graph_err_E.SetFillStyle(3001);
        ratio_graph_err_M.SetFillColorAlpha(kGray+2,0.7);
        ratio_graph_err_M.SetMarkerColor(kGray);
        ratio_graph_err_M.SetMarkerSize(0);
        ratio_graph_err_M.SetFillStyle(3001);
        ratio_graph_err_E.SetFillColorAlpha(kGray+2,0.7);
        ratio_graph_err_E.SetMarkerColor(kGray);
        ratio_graph_err_E.SetMarkerSize(0);
        ratio_graph_err_E.SetFillStyle(3001);

    if args.sumEM:
      #####################################
      ########## sumEM channel ############
      #####################################
      ymax_T = ymax_E+ymax_M
      y_T = stack_T.GetMaximum()
      if y_T>ymax_T: ymax_T=y_T
      #stack_M.SetMinimum(1000.)
      #stack_M.SetMaximum(3*ymax_T)
      if args.linear: stack_T.SetMaximum(1.25*ymax_T)
      if args.linear: stack_T.SetMinimum(1.)

      if not args.nodata:
        histD_T = histD_M
        histD_T.Add(histD_E)

      stack_T.Draw("HIST")
      stack_T.GetYaxis().SetTitle("Events")
      if not args.ratio: 
         stack_T.GetXaxis().SetTitle(cosmetic_names[name])
         stack_T.GetYaxis().SetLabelSize(titY_size2)
      stack_T.GetYaxis().SetTitleSize(titY_size2)
      stack_T.GetYaxis().SetTitleOffset(titY_off2)
      if not args.nodata:
        if args.channel=="btagMM_chitest" or args.channel=="btagMM_chitest_antisl":
           if name == "InvM_2jets_short": histD_T.SetBinContent(25,0.)
           if name == "MET_my_sig": histD_T.SetBinContent(1,0.)
           if name == "deltaR_jet1_jet2": histD_T.SetBinContent(4,0.)
        if "sl" in str(args.channel) and not("anti" in str(args.channel)):
           if name == "jet_2_nmu": histD_T.SetBinContent(1,0.)
        histD_T.SetMarkerStyle(20)
        histD_T.SetMarkerSize(0.5)
        histD_T.SetLineWidth(1)
        histD_T.SetLineColor(ROOT.kBlack)
        histD_T.Draw("E SAME")
      if not args.nosyst: graph_err_T.Draw("SAME 2")

      if args.ratio and not args.nodata:
        CMS_lumi.CMS_lumi(upper_pad, iPeriod, iPos)
        upper_pad.cd()
        upper_pad.Update()
        upper_pad.RedrawAxis()
        frame = upper_pad.GetFrame()
        #frame.Draw()

      if args.ratio:
        lower_pad.cd()
        if not args.nodata:
          ratio = histD_T.Clone("ratio")
          ratio.SetLineColor(kBlack)
          ratio.SetMarkerStyle(21)
          ratio.SetTitle("")
          ratio.SetMinimum(c_rat2)
          ratio.SetMaximum(c_rat)
          ratio.GetYaxis().SetTitle("Data/MC")
          ratio.GetXaxis().SetTitle(cosmetic_names[name])
          ratio.GetXaxis().SetLabelSize(0.08)
          ratio.GetXaxis().SetTitleSize(tit_size)
          ratio.GetXaxis().SetTitleOffset(titX_off)
          ratio.GetYaxis().SetLabelSize(0.05)
          ratio.GetYaxis().SetTitleSize(titY_size)
          ratio.GetYaxis().CenterTitle(False)
          ratio.GetYaxis().ChangeLabel(nrat,-1,-1,-1,-1,-1,"  ");
          ratio.GetYaxis().SetTitleOffset(titY_off)
          # Set up plot for markers and errors
          ratio.Sumw2()
          ratio.SetStats(0)
          hTotal = histT_nom_T["vv"].Clone('hTotal')
          for s in samples[1:]:
            hTotal.Add(histT_nom_T[s])
          ratio.Divide(hTotal)
          ratio.Draw("ep")
          if not args.nosyst: ratio_graph_err_T.Draw("same 2")
        else:
          if not args.nosyst:
             ratio = histT_nom_T["vv"].Clone('ratio')
             ratio.SetLineColor(kBlack)
             ratio.SetMarkerStyle(6)
             ratio.SetTitle("")
             ratio.SetMinimum(c_rat2)
             ratio.SetMaximum(c_rat)
             ratio.GetYaxis().SetTitle("Ctag syst ratio")
             ratio.GetXaxis().SetTitle(cosmetic_names[name])
             ratio.GetXaxis().SetLabelSize(0.08)
             ratio.GetXaxis().SetTitleSize(tit_size)
             ratio.GetXaxis().SetTitleOffset(titX_off)
             ratio.GetYaxis().SetLabelSize(0.05)
             ratio.GetYaxis().SetTitleSize(titY_size)
             ratio.GetYaxis().CenterTitle(False)
             ratio.GetYaxis().ChangeLabel(nrat,-1,-1,-1,-1,-1,"  ");
             ratio.GetYaxis().SetTitleOffset(titY_off)
             # Set up plot for markers and errors
             ratio.Sumw2()
             ratio.SetStats(0)
             ratio.Divide(ratio)
             ratio.Draw("hist p")
             ratio_graph_err_T.Draw("same 2")

      if "flavourP" in str(name) and args.wcs:
         print(name)
         for s in ["ttbar_dl","ttbar_sl_nocharm","ttbar_sl_charm"]:
           print(s)
           qaux1 = histT_nom_T[s].Integral()
           print("Total events are "+str(qaux1))
           for i in [2,3,11,12,28]:
              print("Content of bin "+str(i)+" is "+str(histT_nom_T[s].GetBinContent(i)))
              print("Percentage over total is "+str(histT_nom_T[s].GetBinContent(i)/qaux1))

      if name == "jet_max_cvltag" and (not args.nodata):
          print("charm tag")
          print("Integral of data is "+str(histD_T.Integral()))
          print("Integral of MC is "+str(hTotal.Integral()))
          for s in samples:
              print("Integral of "+str(s)+" MC is "+str(histT_nom_T[s].Integral()))
          print("Integral of data from 0.3 is "+str(histD_T.Integral(16,51)))
          print("Integral of MC from 0.3 is "+str(hTotal.Integral(16,51)))
          for s in samples:
              print("Integral of "+str(s)+" MC from 0.3 is "+str(histT_nom_T[s].Integral(16,51)))

      if name == "InvM_2jets" and (not args.nodata):
        print("Integral of data is "+str(histD_T.Integral()))
        print("Integral of MC is "+str(hTotal.Integral()))
        print("Ratio is "+str(histD_T.Integral()/hTotal.Integral()))
        error_tot = {}
        error_fulltot = 0
        for s in samples:
            print("Integral of "+str(s)+" MC process is "+str(histT_nom_T[s].Integral()))
            error_tot[s] = 0
            for bin in range(histT_nom_T[s].GetNbinsX()):
                error_tot[s] = error_tot[s] + histT_nom_T[s].GetBinError(bin+1)**2
            error_tot[s] = sqrt(error_tot[s])
            print("error is "+str(error_tot[s]))
            error_fulltot = error_fulltot + error_tot[s]**2
        print("Final MC error is "+str(sqrt(error_fulltot)))

      ## Legends
      if args.ratio: upper_pad.cd()
      leg = TLegend(leg1, leg2, leg3, leg4)
      leg.SetBorderSize(0)
      leg.SetNColumns(2)
      #leg.AddEntry(histT_nom_T["vv"],"VV","f")
      if args.wcs:
        leg.AddEntry(histT_nom_T["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_nocharm"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_T["st_charm"],"Single top tW cq","f")
        leg.AddEntry(histT_nom_T["st_nocharm"],"Single top tW uq","f")
        leg.AddEntry(histT_nom_T["st_else"],"Single top s/t channel","f")
        leg.AddEntry(histT_nom_T["ttbar_dl"],"Dileptonic t#bar{t}","f")
        leg.AddEntry(histT_nom_T["zjets"],"Z+jets","f")
        leg.AddEntry(histT_nom_T["wjets"],"W+jets","f")
        if args.stack and not args.nodata: leg.AddEntry(histD_T, "Data" ,"lep")
      else:
        leg.AddEntry(histT_nom_T["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_T["ttbar_dl"],"Dileptonic t#bar{t}","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_light"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_T["st"],"Single top","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_bottom"],"t#bar{t} bq","f")
        leg.AddEntry(histT_nom_T["zjets"],"Z+jets","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_else"],"t#bar{t} qg,gg","f")
        leg.AddEntry(histT_nom_T["wjets"],"W+jets","f")
        leg.AddEntry(histT_nom_T["ttbar_sl_charmgluon"],"t#bar{t} cg","f")
        if args.stack and not args.nodata: leg.AddEntry(histD_T, "Data" ,"lep")
        leg.AddEntry(histT_nom_T["ttbar_sl_bottomgluon"],"t#bar{t} bg","f")
      #leg.AddEntry(histT_nom_T["ttbar_dh"],"Hadronic t#bar{t}","f")
      leg.Draw()
      termp= "totalHT_wqq"
      if args.ratio:
        notation = "_ratio_"
        if args.linear:
          notation = "_linratio_"
      else:
        notation = "_normed_"

      if args.year == "all": term_d = ""
      elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): term_d = str(args.year)

      if args.png: 
         c1.Print(plotdir+termp+notation+ term_d+name + ".pdf");
         c1.Print(plotdir+termp+notation+ term_d+name + ".png");
      else: 
         c1.Print(plotdir+termp+notation + term_d+name + ".root")

      if c1: 
         c1.Close(); gSystem.ProcessEvents();

    else:
      #####################################
      ############ M channel ##############
      #####################################
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
      stack_M.GetYaxis().SetTitle("Events")
      stack_M.GetYaxis().SetTitleSize(titY_size2)
      stack_M.GetYaxis().SetTitleOffset(titY_off2)
      if not args.ratio:
         stack_M.GetXaxis().SetTitle(cosmetic_names[name])
         stack_M.GetYaxis().SetLabelSize(titY_size2)
      if not args.nodata:
        histD_M.SetMarkerStyle(20)
        histD_M.SetMarkerSize(0.5)
        histD_M.SetLineWidth(1)
        histD_M.SetLineColor(ROOT.kBlack)
        histD_M.Draw("E SAME")
      if not args.nosyst: graph_err_M.Draw("SAME 2")

      if args.ratio:
        CMS_lumi.CMS_lumi(upper_pad, iPeriod, iPos)
        upper_pad.cd()
        upper_pad.Update()
        upper_pad.RedrawAxis()
        frame = upper_pad.GetFrame()
        lower_pad.cd()
        if not args.nodata:
          ratio = histD_M.Clone("ratio")
          ratio.SetLineColor(kBlack)
          ratio.SetMarkerStyle(21)
          ratio.SetTitle("")
          ratio.SetMinimum(c_rat2)
          ratio.SetMaximum(c_rat)
          ratio.GetYaxis().SetTitle("Data/MC")
          ratio.GetXaxis().SetTitle(cosmetic_names[name])
          ratio.GetXaxis().SetLabelSize(0.08)
          ratio.GetXaxis().SetTitleSize(tit_size)
          ratio.GetXaxis().SetTitleOffset(titX_off)
          ratio.GetYaxis().SetLabelSize(0.05)
          ratio.GetYaxis().SetTitleSize(titY_size)
          ratio.GetYaxis().CenterTitle()
          ratio.GetYaxis().ChangeLabel(nrat,-1,-1,-1,-1,-1,"  ");
          ratio.GetYaxis().SetTitleOffset(titY_off)
          # Set up plot for markers and errors
          ratio.Sumw2()
          ratio.SetStats(0)
          hTotal = histT_nom_M["vv"].Clone('hTotal')
          for s in samples[1:]:
            hTotal.Add(histT_nom_M[s])
          ratio.Divide(hTotal)
          ratio.Draw("ep")
        if not args.nosyst: ratio_graph_err_M.Draw("same 2")

      if name == "InvM_2jets" and (not args.nodata):
        print("Integral of M data is "+str(histD_M.Integral()))
        print("Integral of M MC is "+str(hTotal.Integral()))
        print("Ratio is "+str(histD_M.Integral()/hTotal.Integral()))
        error_tot = {}
        error_fulltot = 0
        for s in samples:
            print("Integral of "+str(s)+" MC process is "+str(histT_nom_M[s].Integral()))
            error_tot[s] = 0
            for bin in range(histT_nom_M[s].GetNbinsX()):
                error_tot[s] = error_tot[s] + histT_nom_M[s].GetBinError(bin+1)**2
            error_tot[s] = sqrt(error_tot[s])
            print("error is "+str(error_tot[s]))
            error_fulltot = error_fulltot + error_tot[s]**2
        print("Final MC error is "+str(sqrt(error_fulltot)))

      ## Legends
      if args.ratio and not args.nodata: upper_pad.cd()
      leg = TLegend(leg1, leg2, leg3, leg4)
      leg.SetBorderSize(0)
      leg.SetNColumns(2)
      #leg.AddEntry(histT_nom_M["vv"],"VV","f")
      if args.wcs:
        leg.AddEntry(histT_nom_M["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_nocharm"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_M["st_charm"],"Single top tW cq","f")
        leg.AddEntry(histT_nom_M["st_nocharm"],"Single top tW uq","f")
        leg.AddEntry(histT_nom_M["st_else"],"Single top s/t channel","f")
      else:
        leg.AddEntry(histT_nom_M["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_light"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_M["st"],"Single top","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_bottom"],"t#bar{t} bq","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_else"],"t#bar{t} qg,gg","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_charmgluon"],"t#bar{t} cg","f")
        leg.AddEntry(histT_nom_M["ttbar_sl_bottomgluon"],"t#bar{t} bg","f")
      leg.AddEntry(histT_nom_M["zjets"],"Z+jets","f")
      leg.AddEntry(histT_nom_M["wjets"],"W+jets","f")
      leg.AddEntry(histT_nom_M["ttbar_dl"],"Dileptonic t#bar{t}","f")
      #leg.AddEntry(histT_nom_M["ttbar_dh"],"Hadronic t#bar{t}","f")
      if args.stack and not args.nodata: leg.AddEntry(histD_M, "Data" ,"lep")
      leg.Draw()
      termp= "totalHT_wqq_M"
      if args.ratio: 
        notation = "_ratio_"
        if args.linear:
          notation = "_linratio_"
      else: 
        notation = "_normed_"

      if args.year == "all": term_d = ""
      elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): term_d = str(args.year)

      if args.png: 
         c1.Print(plotdir+termp+notation+ term_d+name + ".pdf");
         c1.Print(plotdir+termp+notation+ term_d+name + ".png");
      else: c1.Print(plotdir+termp+notation + term_d+name + ".root")

      #####################################
      ############ E channel ##############
      #####################################
      stack_E.Draw("HIST")
      if not args.ratio:
         stack_E.GetXaxis().SetTitle(cosmetic_names[name])
         stack_E.GetYaxis().SetLabelSize(titY_size2)
      stack_E.GetYaxis().SetTitle("Events")
      stack_E.GetYaxis().SetTitleSize(titY_size2)
      stack_E.GetYaxis().SetTitleOffset(titY_off2)
      if not args.nodata:
        histD_E.SetMarkerStyle(20)
        histD_E.SetMarkerSize(0.5)
        histD_E.SetLineWidth(1)
        histD_E.SetLineColor(ROOT.kBlack)
        histD_E.Draw("E SAME")
      if not args.nosyst: graph_err_E.Draw("SAME 2")

      if args.ratio and not args.nodata:
        CMS_lumi.CMS_lumi(upper_pad, iPeriod, iPos)
        upper_pad.cd()
        upper_pad.Update()
        upper_pad.RedrawAxis()
        frame = upper_pad.GetFrame()
        lower_pad.cd()
        if not args.nodata:
          ratio = histD_E.Clone("ratio")
          ratio.SetLineColor(kBlack)
          ratio.SetMarkerStyle(21)
          ratio.SetTitle("")
          ratio.SetMinimum(c_rat2)
          ratio.SetMaximum(c_rat)
          ratio.GetYaxis().SetTitle("Data/MC")
          ratio.GetXaxis().SetTitle(cosmetic_names[name])
          ratio.GetXaxis().SetLabelSize(0.08)
          ratio.GetXaxis().SetTitleSize(tit_size)
          ratio.GetXaxis().SetTitleOffset(titX_off)
          ratio.GetYaxis().SetLabelSize(0.05)
          ratio.GetYaxis().SetTitleSize(titY_size)
          ratio.GetYaxis().CenterTitle()
          ratio.GetYaxis().ChangeLabel(nrat,-1,-1,-1,-1,-1,"  ");
          ratio.GetYaxis().SetTitleOffset(titY_off)
          # Set up plot for markers and errors
          ratio.Sumw2()
          ratio.SetStats(0)
          hTotal = histT_nom_E["vv"].Clone('hTotal')
          for s in samples[1:]:
            hTotal.Add(histT_nom_E[s])
          ratio.Divide(hTotal)
          ratio.Draw("ep")
        if not args.nosyst: ratio_graph_err_E.Draw("same 2")

      if name == "InvM_2jets" and (not args.nodata):
        print("Integral of E data is "+str(histD_E.Integral()))
        print("Integral of E MC is "+str(hTotal.Integral()))
        print("Ratio is "+str(histD_E.Integral()/hTotal.Integral()))
        for s in samples:
            print("Integral of "+str(s)+" MC process is "+str(histT_nom_E[s].Integral()))
            error_tot[s] = 0
            for bin in range(histT_nom_E[s].GetNbinsX()):
                error_tot[s] = error_tot[s] + histT_nom_E[s].GetBinError(bin+1)**2
            error_tot[s] = sqrt(error_tot[s])
            print("error is "+str(error_tot[s]))
            error_fulltot = error_fulltot + error_tot[s]**2
        print("Final MC error is "+str(sqrt(error_fulltot)))

      ## Legends
      if args.ratio and not args.nodata: upper_pad.cd()
      leg = TLegend(leg1, leg2, leg3, leg4)
      leg.SetBorderSize(0)
      leg.SetNColumns(2)
      #leg.AddEntry(histT_nom_E["vv"],"VV","f")
      if args.wcs:
        leg.AddEntry(histT_nom_E["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_nocharm"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_E["st_charm"],"Single top tW cq","f")
        leg.AddEntry(histT_nom_E["st_nocharm"],"Single top tW uq","f")
        leg.AddEntry(histT_nom_E["st_else"],"Single top s/t channel","f")
      else:
        leg.AddEntry(histT_nom_E["ttbar_sl_charm"],"t#bar{t} cq","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_light"],"t#bar{t} uq","f")
        leg.AddEntry(histT_nom_E["st"],"Single top","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_bottom"],"t#bar{t} bq","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_else"],"t#bar{t} qg,gg","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_charmgluon"],"t#bar{t} cg","f")
        leg.AddEntry(histT_nom_E["ttbar_sl_bottomgluon"],"t#bar{t} bg","f")
      leg.AddEntry(histT_nom_E["zjets"],"Z+jets","f")
      leg.AddEntry(histT_nom_E["wjets"],"W+jets","f")
      leg.AddEntry(histT_nom_E["ttbar_dl"],"Dileptonic t#bar{t}","f")
      if args.stack and not args.nodata: leg.AddEntry(histD_E, "Data" ,"lep")
      leg.Draw()
      termp= "totalHT_wqq_E"
      if args.ratio:
        notation = "_ratio_"
        if args.linear:
          notation = "_linratio_"
      else:
        notation = "_normed_"

      if args.year == "all": term_d = ""
      elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): term_d = str(args.year)

      if args.png: 
         c1.Print(plotdir+termp+notation+ term_d+name + ".pdf");
         c1.Print(plotdir+termp+notation+ term_d+name + ".png");
      else: c1.Print(plotdir+termp+notation + term_d+name + ".root")

      if c1: 
         c1.Close(); gSystem.ProcessEvents();

    for data_op in datayears:
                        histFile[name][data_op].Close()
                        if not args.nodata:
                          histFileDM[name][data_op].Close()
                          histFileDE[name][data_op].Close()

