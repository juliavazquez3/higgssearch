######################################                                        
######################################                                        
#######      GOOD  VERSION     #######
######################################                                     
######################################                                     

print('SL CHANNEL')

## Modified version, meant to be run locally, not sent to condor 

## Selection for noth MC and data samples, created to distinguish between years
## Different histogram files are produced for each situation 
## It includes an option for SSOS substraction

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join, isdir

import json
import argparse

sys.path.append('/nfs/cms/vazqueze/ttbaranalisis/last_corrections/')
sys.path.append('/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/')

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

ROOT.EnableImplicitMT()

years = ["2017"]
proc = ["ttbar_sl"]

samples = []

for p in proc:
  for y in years:
    samples.append(p+y)

samples_test = samples[:]

print("the samples treated are",samples)

# Create a ROOT dataframe for each dataset
# Note that we load the filenames from the external json file placed in the same folder than this script.
# Example "python analisisWW/selection_v2.py --process="WW" --notfull -l 0 50"

files = json.load(open("/nfs/cms/vazqueze/higgssearch/mcinfo"+years[0]+".json"))

processes = files.keys()

df = {}
xsecs = {}
sumws = {}
archives = {}
event_test = {}

for s in samples:
  archives[s]=[]

names_proc = {"M":{"2016":["2016M1","2016M2","2016M3","2016M4","2016M5","2016M6"],"2016B":["2016BM1","2016BM2","2016BM3"],
   "2017":["2017M1","2017M2","2017M3","2017M4","2017M5"],"2018":["2018M1","2018M2","2018M3","2018M4"]},
   "E":{"2016":["2016E1","2016E2","2016E3","2016E4","2016E5","2016E6"],"2016B":["2016BE1","2016BE2","2016BE3"],
   "2017":["2017E1","2017E2","2017E3","2017E4","2017E5"],"2018":["2018E1","2018E2","2018E3","2018E4"]}}

files_data = {"M":{"2016":230,"2016B":157,"2017":436,"2018":393},
   "E":{"2016":335,"2016B":156,"2017":245,"2018":738}}

term1 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_ver2/"
term2 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_down/"
term3 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux/"
#term3 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_nosmearing/"
for p in proc:
    # Construct the dataframes
    folder = term3 + "myconfig"+years[0]+"/{sample}/cat_base/prod_test/" # Folder name
    if years[0] == "2016B": folder = term3 + "myconfig2016/{sample}/cat_base/prod_test/"
    folder = folder.format(sample = p) # Sample name
    num_files = files[p]["files"] # Number of files
    list_files = [f for f in listdir(folder) if isfile(join(folder, f))] # Lista de archivos
    if (num_files == len(list_files)):
         for f in list_files:
               file_r = TFile(join(folder, f))
               if file_r.GetListOfKeys().Contains("Events"):
                  archives[p+years[0]].append(join(folder,f))

df_test = ROOT.RDataFrame("Events",set(archives["ttbar_sl2017"][0:30]))

chain_ttbar_sl = TChain("Events","");
chain_ttbar_sl.Add("/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux/myconfig2017/ttbar_sl/cat_base/prod_test/data_*.root");

## Cuts per year btag

cuts_btag = {}
cuts_btag["2016"]=[0.0480, 0.2489, 0.6377];cuts_btag["2016B"]=[0.0480, 0.2489, 0.6377];
cuts_btag["2017"]=[0.0532, 0.3040, 0.7476];cuts_btag["2018"]=[0.0490,0.2783,0.7100];
cuts_btag["2016"]=[0.0614, 0.3093, 0.7221];cuts_btag["2016B"]=[0.0614, 0.3093, 0.7221];

## Cuts per year ctag

cuts_cvl_csv = {}
cuts_cvl_csv["2016"]=[0.088,0.181,0.417];cuts_cvl_csv["2016B"]=[0.088,0.180,0.407];
cuts_cvl_csv["2017"]=[0.040,0.144,0.730];cuts_cvl_csv["2018"]=[0.064,0.153,0.405];

cuts_cvb_csv = {}
cuts_cvb_csv["2016"]=[0.214,0.228,0.138];cuts_cvb_csv["2016B"]=[0.204,0.221,0.136];
cuts_cvb_csv["2017"]=[0.345,0.290,0.100];cuts_cvb_csv["2018"]=[0.313,0.363,0.288];

## Muon and Electron pT cuts per year, not relevant here, we use them in the preprocessing

muon_pt = {}
muon_pt["2016"]=26;muon_pt["2016B"]=26;muon_pt["2017"]=29;muon_pt["2018"]=26;
el_pt = {}
el_pt["2016"]=30;el_pt["2016B"]=30;el_pt["2017"]=35;el_pt["2018"]=35;

## Triggers per year

muon_trig = {}

muon_trig["2016"]="HLT_IsoMu24 || HLT_IsoTkMu24"
muon_trig["2016B"]="HLT_IsoMu24 || HLT_IsoTkMu24"
muon_trig["2017"]="HLT_IsoMu27"
muon_trig["2018"]="HLT_IsoMu24"

el_trig = {}

el_trig["2016"]="HLT_Ele27_WPTight_Gsf"
el_trig["2016B"]="HLT_Ele27_WPTight_Gsf"
el_trig["2017"]="HLT_Ele32_WPTight_Gsf"
el_trig["2018"]="HLT_Ele32_WPTight_Gsf"

### MET Filters per year

met_filter = {}

met_filter["2016"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_BadPFMuonFilter")
met_filter["2016B"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_BadPFMuonFilter")
met_filter["2017"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_BadPFMuonFilter")
met_filter["2018"] = ("Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalSuperTightHalo2016Filter "
             "&& Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_BadPFMuonFilter")

###################################################
################   DEFINITIONS   ##################
###################################################

gInterpreter.Declare("""
      #include <bitset>
      #include <string>
      #include <iostream>
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using Vbool = const ROOT::RVec<bool>&;
      auto muonmotherverprev(Vint pdg, Vint mother, const int muon_id){
            int mother_id_v1 = fabs(pdg[mother[muon_id]]);
            int mother_id_v2 = -10;
            if(mother_id_v1 == 411){
                 mother_id_v2 = 1;
            } else if (mother_id_v1 == 421) {
                 mother_id_v2 = 2;
            } else if (mother_id_v1 == 431) {
                 mother_id_v2 = 3;
            } else if (mother_id_v1 == 4122) {
                 mother_id_v2 = 4;
            } else if (mother_id_v1 == 511) {
                 mother_id_v2 = 5;
            } else if (mother_id_v1 == 521) {
                 mother_id_v2 = 6;
            } else if (mother_id_v1 == 531) {
                 mother_id_v2 = 7;
            } else if (mother_id_v1 == 5122) {
                 mother_id_v2 = 8;
            } else {
                 mother_id_v2 = 0;
            }
            return mother_id_v2;
      };
      auto muongenreco(UInt_t nPart, Vint pdg, Vint mother, Vfloat eta, Vfloat phi, Vfloat pt, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vbool mu_id) {
            vector<bool> vb;
            vector<int> muind;
            int indM = -1;
            float ptM = -10.;
            float ptM2 = -10.;
            bool cond1 = false;
            bool condgen = false;
            bool condreco = false;
            for (unsigned int i=0; i<nPart; ++i) {
                if (pt[i]>ptM && fabs(pdg[i])==13 && muonmotherverprev(pdg,mother,i)>0) {
                   indM = i;
                   ptM = pt[i];
                }
                if (indM>-1) {
                   condgen = true;
                   muind.push_back(indM);
                }
            }
            if (condgen) {
              for (unsigned int k=0; k<mu_pt.size(); ++k) {
                 cond1 = ROOT::VecOps::DeltaR(mu_eta[k],eta[muind[0]],mu_phi[k],phi[muind[0]]) < 0.4;
                 if (mu_id[k] && cond1 && mu_pt[k]>ptM2) {
                   condreco = true;
                   ptM = mu_pt[k];
                 }
              }
            }
            vb.push_back(condgen);
            vb.push_back(condreco);
            return vb;
      };
      auto muongenrecopt(UInt_t nPart, Vint pdg, Vint mother, Vfloat pt) {
            vector<int> muind;
            int indM = -1;
            float ptM = -10.;
            bool condgen = false;
            for (unsigned int i=0; i<nPart; ++i) {
                if (pt[i]>ptM && fabs(pdg[i])==13 && muonmotherverprev(pdg,mother,i)>0) {
                   indM = i;
                   ptM = pt[i];
                }
       	       	if (indM>-1) {
                   condgen = true;
                   muind.push_back(indM);
                }
            }
            return ptM;
      };
""")


#############################################################
######################## Definitions ########################
#############################################################

gInterpreter.Declare("""
   Double_t xbins[10] = {0.,3.,4.,5.,7.,10.,15.,20.,25.,40.};
""")


df_test = df_test.Define('vec_bool','muongenreco(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, GenPart_pt, Muon_pt, Muon_eta, Muon_phi, Muon_tightId)')
df_test = df_test.Define('muongenpt','muongenrecopt(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt)')
df_gen = df_test.Filter('vec_bool[0]')
df_reco = df_gen.Filter('vec_bool[1]')

histgen1 = df_gen.Histo1D(("histgen1","",9,xbins),"muongenpt")
histgen2 = df_reco.Histo1D(("histgen2","",9,xbins),"muongenpt")

#ratio_hist = histgen2.Clone("ratio_hist")
#ratio_hist.Divide(histgen1)

#############################
####     DATA SAVING     ####
#############################

path_hist = "/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/muon_eff_files.root"
myfile = TFile( path_hist, 'RECREATE' )

histgen1.Write()
histgen2.Write()
#ratio_hist.Write()

myfile.Close()

print('Ended succesfully')
