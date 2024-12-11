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

path_df = "/pnfs/ciemat.es/data/cms/store/user/juvazque/data_further_analysis/btagMM/folder"
df_2 = {}

for dyear in ["2016","2016B","2017","2018"]:
    aux_path = path_df+str(dyear)+"/wcs/dataset_wqq_btagMM_fromJF_tt*.root"
    df_2[dyear] = ROOT.RDataFrame("Events",aux_path)

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


## Funciones para seleccionar JETS
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto InvariantM(const float pt, const float eta, const float phi, const float mass, const float pt1, const float eta1, const float phi1, const float mass1) {
            auto x = pt*std::cos(phi);
            auto x1 = pt1*std::cos(phi1);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto z = pt*std::sinh(eta);
            auto z1 = pt1*std::sinh(eta1);
            auto e = std::sqrt(x*x+y*y+z*z+mass*mass);
            auto e1 = std::sqrt(x1*x1+y1*y1+z1*z1+mass1*mass1);

            auto mJet = std::sqrt((e+e1)*(e+e1)-(x+x1)*(x+x1)-(y+y1)*(y+y1)-(z+z1)*(z+z1));
            return mJet;
      };
      auto InvariantM3(const float pt, const float eta, const float phi, const float pt1, const float eta1, const float phi1, const float pt2, const float eta2, const float phi2) {
            auto x = pt*std::cos(phi);
            auto x1 = pt1*std::cos(phi1);
            auto x2 = pt2*std::cos(phi2);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto y2 = pt2*std::sin(phi2);
            auto z = pt*std::sinh(eta);
            auto z1 = pt1*std::sinh(eta1);
            auto z2 = pt2*std::sinh(eta2);
            auto e = pt*std::cosh(eta);
            auto e1 =  pt1*std::cosh(eta1);
            auto e2 = pt2*std::cosh(eta2);

            auto mJet = std::sqrt((e+e1+e2)*(e+e1+e2)-(x+x1+x2)*(x+x1+x2)-(y+y1+y2)*(y+y1+y2)-(z+z1+z2)*(z+z1+z2));
            return mJet;
      };
""")

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto chi2calcver2(const float mw, const float mt){
        Double_t xx2 = 0.;
        float mwval = 0.;
        float mtopval = 0.;
        float sigw = 0.;
        float sigtop = 0.;
        mwval = 77.5;
        mtopval = 161;
        sigw = 11.5;
        sigtop = 20;
	xx2  = (mw-mwval)*(mw-mwval)/sigw/sigw;
        xx2 += (mt-mtopval)*(mt-mtopval)/sigtop/sigtop;
        xx2 += -0.6*2.*(mw-mwval)*(mt-mtopval)/sigw/sigtop;
        xx2 = xx2/(1-0.6*0.6);

        if (xx2<0) return -999.;

        return xx2;
      };
      auto JetInds(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, Vint jetbotind) {
            vector<int> vb;
            bool condb = true;
            int ind1 = -1;
            int ind2 = -1;
            float ptJ{-10.};
            for (unsigned int j=0; j<4; ++j){
                        if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                        if(pt[good[j]]>ptJ && condb){
                                ind1 = good[j];
                                ptJ = pt[good[j]];
                        }
            }
            if (ind1 > -1) {
                vb.push_back(ind1);
            }
            ptJ = -10.;
            for (unsigned int j=0; j<4; ++j){
                        if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                        if(pt[good[j]]>ptJ && condb && good[j] != ind1){
                                ind2 = good[j];
                                ptJ = pt[good[j]];
                        }
            }
            if (ind2 > -1) {
                vb.push_back(ind2);
            }
            return vb;
      };
      auto InvMhad(Vint good, Vint bgood, Vfloat pt, Vfloat eta, Vfloat phi)  {
           vector<float> vb;
           auto vec = JetInds(4,good,pt,eta,phi,bgood);
           auto jet1 = vec[0];
           auto jet2 = vec[1];
           auto m30 = InvariantM3(pt[good[bgood[0]]],eta[good[bgood[0]]],phi[good[bgood[0]]],pt[jet1],eta[jet1],phi[jet1],pt[jet2],eta[jet2],phi[jet2]);
           auto m31 = InvariantM3(pt[good[bgood[1]]],eta[good[bgood[1]]],phi[good[bgood[1]]],pt[jet1],eta[jet1],phi[jet1],pt[jet2],eta[jet2],phi[jet2]);
           auto mw = InvariantM(pt[jet1],eta[jet1],phi[jet1],0.,pt[jet2],eta[jet2],phi[jet2],0.);
           vb.push_back(m30);
           vb.push_back(m31);
           vb.push_back(mw);
           vb.push_back(chi2calcver2(mw,m30));
           vb.push_back(chi2calcver2(mw,m31));
           return vb;
      };
""")


gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto InvMasslep(Vint jetbot, Vint jetgood, Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi, Vint muongood, Vint elgood, Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass, Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass) {
            vector<float> vb;
            if (muongood.size() > 0) {
               vb.push_back(InvariantM(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],0., muon_pt[muongood[0]], muon_eta[muongood[0]], muon_phi[muongood[0]], muon_mass[muongood[0]]));
               vb.push_back(InvariantM(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],0., muon_pt[muongood[0]], muon_eta[muongood[0]], muon_phi[muongood[0]], muon_mass[muongood[0]]));
            } else {
               vb.push_back(InvariantM(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],0., el_pt[elgood[0]], el_eta[elgood[0]], el_phi[elgood[0]], el_mass[elgood[0]]));
               vb.push_back(InvariantM(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],0., el_pt[elgood[0]], el_eta[elgood[0]], el_phi[elgood[0]], el_mass[elgood[0]]));
            }
            return vb;
      };
      auto bjetschoose(Vint partonF, Vint jetbot, Vint jetgood, Vfloat pt, Vfloat eta, Vfloat phi, Vint muongood, Vint elgood, Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass, Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Vint muon_charge, Vint el_charge) {
            vector<float> vb;
            int indM = -1;
            float ptM = -10.;
            bool condgen = false;
            auto vec_had = InvMhad(jetgood, jetbot, pt,eta,phi);
            auto vec_aux = InvMasslep(jetbot,jetgood,pt,eta,phi,muongood,elgood,muon_pt,muon_eta,muon_phi,muon_mass,el_pt,el_eta,el_phi,el_mass);
            if (muongood.size() > 0) {
               if (partonF[jetgood[jetbot[0]]]*muon_charge[muongood[0]] > 0) {
                   vb.push_back(vec_aux[0]);
                   vb.push_back(vec_aux[1]);
                   vb.push_back(vec_had[4]);
                   vb.push_back(vec_had[3]);
               } else {
                   vb.push_back(vec_aux[1]);
                   vb.push_back(vec_aux[0]);
                   vb.push_back(vec_had[3]);
                   vb.push_back(vec_had[4]);
               }
            } else {
               if (partonF[jetgood[jetbot[0]]]*el_charge[elgood[0]] > 0) {
                   vb.push_back(vec_aux[0]);
                   vb.push_back(vec_aux[1]);
                   vb.push_back(vec_had[4]);
                   vb.push_back(vec_had[3]);
               } else {
                   vb.push_back(vec_aux[1]);
                   vb.push_back(vec_aux[0]);
                   vb.push_back(vec_had[3]);
                   vb.push_back(vec_had[4]);
               }
            }
            return vb;
      };
""")


#############################################################
######################## Definitions ########################
#############################################################

df_test = df_test.Define('vec_aux','bjetschoose(Jet_partonFlavour, JetBotInd, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, MuonGoodInd, ElectronGoodInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, Muon_charge, Electron_charge)')
df_test = df_test.Define('InvMl_good','vec_aux[0]')
df_test = df_test.Define('InvMl_bad','vec_aux[1]')
df_test = df_test.Define('chi2_good','vec_aux[2]')
df_test = df_test.Define('chi2_bad','vec_aux[3]')

hist_InvMl_good = df_test.Histo1D(("hist_InvMl_good","",100,0,600),"InvMl_good")
hist_InvMl_bad = df_test.Histo1D(("hist_InvMl_bad","",100,0,600),"InvMl_bad")
hist_chi2_good = df_test.Histo1D(("hist_chi2_good","",100,0,10),"chi2_good")
hist_chi2_bad = df_test.Histo1D(("hist_chi2_bad","",100,0,10),"chi2_bad")

#############################
####     DATA SAVING     ####
#############################

path_hist = "/nfs/cms/vazqueze/higgssearch/fromJF/chi_invMl_test.root"
myfile = TFile( path_hist, 'RECREATE' )

hist_InvMl_good.Write()
hist_InvMl_bad.Write()
hist_chi2_good.Write()
hist_chi2_bad.Write()

myfile.Close()

print('Ended succesfully')
