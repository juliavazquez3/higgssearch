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

path1 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/data_further_analysis/btagMM/folder2018/wcs/"

archives_list = [path1+"dataset_wqq_btagMM_fromJF_ttbar_sl2018_charm.root",path1+"dataset_wqq_btagMM_fromJF_ttbar_sl2018_nocharm.root"]

df_test = ROOT.RDataFrame("Events",set(archives_list))

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

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto vector2mu(const int vec0, const int vec1) {
            vector<int> vb;
            vb.push_back(vec0);
            vb.push_back(vec1);
            return vb;
      };
      auto muonEffSelection(Vint qind, Vint good, Vint bgood, Vfloat pt, Vfloat eta, Vfloat phi, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vint mu_good, Vbool mu_id){
            vector<int> vb;
            bool cond1 = false;
            bool cond2 = false;
            bool condb1 = false;
            bool condb2 = false;
            bool condmu = true;
            int indM1 = 0;
            int indM2 = 0;
            int indMb1 = 0;
            int indMb2 = 0;
            float ptM{-10.};
            for (unsigned int i=0; i<mu_pt.size(); ++i){
                cond1 = ROOT::VecOps::DeltaR(mu_eta[i],eta[qind[0]],mu_phi[i],phi[qind[0]]) < 0.4;
                cond2 = ROOT::VecOps::DeltaR(mu_eta[i],eta[qind[1]],mu_phi[i],phi[qind[1]]) < 0.4;
                condb1 = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[bgood[0]]],mu_phi[i],phi[good[bgood[0]]]) < 0.4;
                condb2 = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[bgood[1]]],mu_phi[i],phi[good[bgood[1]]]) < 0.4;
                if (mu_good.size() > 0) condmu = mu_good[0] != i;
                if(condmu && mu_id[i] && mu_pt[i]<25. && cond1 && mu_pt[i]>ptM && mu_pt[i]>5. && fabs(mu_eta[i])<2.4){
                     indM1 = 1;
                } else if (condmu && mu_id[i] && mu_pt[i]<25. && cond2 && mu_pt[i]>ptM && mu_pt[i]>5. && fabs(mu_eta[i])<2.4){
                     indM2 = 1;
                } else if (condmu && mu_id[i] && mu_pt[i]<25. && condb1 && mu_pt[i]>ptM && mu_pt[i]>5. && fabs(mu_eta[i])<2.4){
                     indMb1 = 1;
                } else if (condmu && mu_id[i] && mu_pt[i]<25. && condb2 && mu_pt[i]>ptM && mu_pt[i]>5. && fabs(mu_eta[i])<2.4){
                     indMb2 = 1;
                }
            }            
            vb.push_back(indM1);
            vb.push_back(indM2);
            vb.push_back(indMb1);
            vb.push_back(indMb2);
            return vb;
      };
""")

#############################################################
######################## Definitions ########################
#############################################################

gInterpreter.Declare("""
   Double_t xbins[11] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.};
""")


df_test = df_test.Define('vec_muon_aux','muonEffSelection(JetQInd, JetGoodInd, JetBotInd, Jet_pt, Jet_eta, Jet_phi, Muon_pt, Muon_eta, Muon_phi, MuonGoodInd, Muon_tightId)')
df_test = df_test.Define('jet_1_flav','fabs(Jet_partonFlavour[JetQInd[0]])')
df_test = df_test.Define('jet_2_flav','fabs(Jet_partonFlavour[JetQInd[1]])')
df_test = df_test.Define('jet_bot1_flav','fabs(Jet_partonFlavour[JetGoodInd[JetBotInd[0]]])')
df_test = df_test.Define('jet_bot2_flav','fabs(Jet_partonFlavour[JetGoodInd[JetBotInd[1]]])')

list_jets = ["jet_1","jet_2","jet_bot1","jet_bot2"]
cond_jet = {};cond_jet["jet_1"] = 'vec_muon_aux[0]>0';cond_jet["jet_2"] = 'vec_muon_aux[1]>0';
cond_jet["jet_bot1"] = 'vec_muon_aux[2]>0';cond_jet["jet_bot2"] = 'vec_muon_aux[3]>0';

df_jet_mu = {}
df_jet_charm = {}
df_jet_charm_mu = {}
df_jet_bot = {}
df_jet_bot_mu = {}
df_jet_light = {}
df_jet_light_mu = {}
for l in list_jets:
    df_jet_mu[l] = df_test.Filter(cond_jet[l])
    df_jet_charm[l] = df_test.Filter(str(l)+'_flav == 4')
    df_jet_bot[l] = df_test.Filter(str(l)+'_flav == 5')
    df_jet_light[l] = df_test.Filter(str(l)+'_flav == 1 || '+str(l)+'_flav == 2 || '+str(l)+'_flav == 3')
    df_jet_charm_mu[l] = df_jet_charm[l].Filter(cond_jet[l])
    df_jet_bot_mu[l] = df_jet_bot[l].Filter(cond_jet[l])
    df_jet_light_mu[l] = df_jet_light[l].Filter(cond_jet[l])

histpT = {}
histbtag = {}
histpT_mu = {}
histbtag_mu = {}
histpT_charm = {}
histbtag_charm = {}
histpT_charm_mu = {}
histbtag_charm_mu = {}
histpT_bot = {}
histbtag_bot = {}
histpT_bot_mu = {}
histbtag_bot_mu = {}
histpT_light = {}
histbtag_light = {}
histpT_light_mu = {}
histbtag_light_mu = {}
for l in list_jets:
    histpT[l] = df_test.Histo1D(("histpT_"+str(l),"",50,15,115),str(l)+"_pt")
    histbtag[l] = df_test.Histo1D(("histbtag_"+str(l),"",10,xbins),str(l)+"_btag")
    histpT_mu[l] = df_jet_mu[l].Histo1D(("histpT_"+str(l)+"_mu","",50,15,115),str(l)+"_pt")
    histbtag_mu[l] = df_jet_mu[l].Histo1D(("histbtag_"+str(l)+"_mu","",10,xbins),str(l)+"_btag")
    histpT_charm[l] = df_jet_charm[l].Histo1D(("histpT_charm_"+str(l),"",50,15,115),str(l)+"_pt")
    histbtag_charm[l] = df_jet_charm[l].Histo1D(("histbtag_charm_"+str(l),"",10,xbins),str(l)+"_btag")
    histpT_charm_mu[l] = df_jet_charm_mu[l].Histo1D(("histpT_charm_"+str(l)+"_mu","",50,15,115),str(l)+"_pt")
    histbtag_charm_mu[l] = df_jet_charm_mu[l].Histo1D(("histbtag_charm_"+str(l)+"_mu","",10,xbins),str(l)+"_btag")
    histpT_bot[l] = df_jet_bot[l].Histo1D(("histpT_bot_"+str(l),"",50,15,115),str(l)+"_pt")
    histbtag_bot[l] = df_jet_bot[l].Histo1D(("histbtag_bot_"+str(l),"",10,xbins),str(l)+"_btag")
    histpT_bot_mu[l] = df_jet_bot_mu[l].Histo1D(("histpT_bot_"+str(l)+"_mu","",50,15,115),str(l)+"_pt")
    histbtag_bot_mu[l] = df_jet_bot_mu[l].Histo1D(("histbtag_bot_"+str(l)+"_mu","",10,xbins),str(l)+"_btag")
    histpT_light[l] = df_jet_light[l].Histo1D(("histpT_light_"+str(l),"",50,15,115),str(l)+"_pt")
    histbtag_light[l] = df_jet_light[l].Histo1D(("histbtag_light_"+str(l),"",10,xbins),str(l)+"_btag")
    histpT_light_mu[l] = df_jet_light_mu[l].Histo1D(("histpT_light_"+str(l)+"_mu","",50,15,115),str(l)+"_pt")
    histbtag_light_mu[l] = df_jet_light_mu[l].Histo1D(("histbtag_light_"+str(l)+"_mu","",10,xbins),str(l)+"_btag")

#############################
####     DATA SAVING     ####
#############################

path_hist = "/nfs/cms/vazqueze/higgssearch/fromJF/muon_eff_jetflavour.root"
myfile = TFile( path_hist, 'RECREATE' )

for l in list_jets:
    histpT[l].Write()
    histbtag[l].Write()
    histpT_mu[l].Write()
    histbtag_mu[l].Write()
    histpT_charm[l].Write()
    histbtag_charm[l].Write()
    histpT_charm_mu[l].Write()
    histbtag_charm_mu[l].Write()
    histpT_bot[l].Write()
    histbtag_bot[l].Write()
    histpT_bot_mu[l].Write()
    histbtag_bot_mu[l].Write()
    histpT_light[l].Write()
    histbtag_light[l].Write()
    histpT_light_mu[l].Write()
    histbtag_light_mu[l].Write()

myfile.Close()

print('Ended succesfully')
