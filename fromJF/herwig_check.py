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

#file = uproot.open("/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL18NanoAODv9/TT_TuneCH3_13TeV-powheg-herwig7/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/54FCF727-5A53-244F-93A7-A7A519DD150D.root")
#file = uproot.open("/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL18NanoAODv9/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/280000/02003FED-4DED-2849-BAFB-494212F691EE.root")

archives_list = ["/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL18NanoAODv9/TT_TuneCH3_13TeV-powheg-herwig7/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/54FCF727-5A53-244F-93A7-A7A519DD150D.root"]
#archives_list = ["/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL18NanoAODv9/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/280000/02003FED-4DED-2849-BAFB-494212F691EE.root"]

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
      struct ASCII
      {
                std::string toBinary(int n)
                {
                        std::string r;
                        while(n!=0) {r=(n%2==0 ?"0":"1")+r; n/=2;}
                        return r;
                }

      };
      auto testSF(UInt_t part_ind, UInt_t status, UInt_t n) {
            char ind;
            bool hardP = false;
            std::bitset<15> b(status);
            auto statusflags_string = b.to_string();
            for(unsigned int j=statusflags_string.length()-1; j>=statusflags_string.length()-(n+1); j--)
            {
                   //std::cout << "statusflags bit " << j << " " << statusflags_string[j] <<std::endl;
                   ind = statusflags_string.at(j);
            }
            if(ind=='1') hardP = true;
            return hardP;
      };
      auto vectorHP(UInt_t nPart, Vint status, Vint pdg, UInt_t n) {
            vector<bool> vb;
            for (unsigned int i=0; i<nPart; ++i) {
                vb.push_back(testSF(i,status[i],n));
            }
            return vb;
      };
      auto ttbartop(UInt_t nPart, Vint status, Vint pdg, Vint mother) {
            bool typetop = false;
            int arr1[] = {1,2,3};
            bool cond = false;
            for (unsigned int i=0; i<nPart; ++i) {
                //cond = fabs(pdg[mother[i]]) == 24;
                //if (cond && ((fabs(pdg[i])>400 && fabs(pdg[i])<499) || (fabs(pdg[i])>4000 && fabs(pdg[i])<4999))) {
                cond = (fabs(pdg[mother[i]])>400 && fabs(pdg[mother[i]])<499) || (fabs(pdg[mother[i]])>4000 && fabs(pdg[mother[i]])<4999);
                if (cond && (fabs(pdg[i])==13)) {
                          typetop = true;
                }
            }
            return typetop;
      };
      auto ttbarcharm(UInt_t nPart, Vint status, Vint pdg, Vint mother) {
            vector<bool> vb;
            bool typeC = false;
            bool typeL = false;
            int arr1[] = {1,2,3};
            bool cond = false;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==4 && fabs(pdg[mother[i]])==24) {
                          typeC = true;
                }
                cond = std::find(std::begin(arr1), std::end(arr1), fabs(pdg[i])) != std::end(arr1);
                if (cond && fabs(pdg[mother[i]])==24) {
                          typeL = true;
                }
            }
            vb.push_back(typeC);
            vb.push_back(typeL);
            return vb;
      };
      auto muonmother(Vint pdg, Vint mother, const int muon_id){
            int mother_id_v1 = fabs(pdg[mother[muon_id]]);
            int grandmother_id_v1 = fabs(pdg[mother[mother[muon_id]]]);
            int greatmother_id_v1 = fabs(pdg[mother[mother[mother[muon_id]]]]);
            int greatmother_id_v2 = fabs(pdg[mother[mother[mother[mother[muon_id]]]]]);
            int greatmother_id_v3 = fabs(pdg[mother[mother[mother[mother[mother[muon_id]]]]]]);
            int mother_id_v2 = -10;
            if(mother_id_v1 == 13){
                 mother_id_v1 = fabs(pdg[mother[mother[muon_id]]]);
                 grandmother_id_v1 = fabs(pdg[mother[mother[mother[muon_id]]]]);
                 greatmother_id_v1 = fabs(pdg[mother[mother[mother[mother[muon_id]]]]]);
                 greatmother_id_v2 = fabs(pdg[mother[mother[mother[mother[mother[muon_id]]]]]]);
                 greatmother_id_v3 = fabs(pdg[mother[mother[mother[mother[mother[mother[muon_id]]]]]]]);
            }
            int arr[] = {511,521,531,5122,5,10411,10421,413,423,10413,10423,20413,20423,415,425,10431,433,10433,20433,435,24};
            int arr1[] = {511,521,531,5122,5};
            int arr2[] = {411,421,431,4122};
            bool cond = std::find(std::begin(arr), std::end(arr), mother_id_v1) != std::end(arr);
            bool cond1 = std::find(std::begin(arr1), std::end(arr1), mother_id_v1) != std::end(arr1);
            bool cond2 = std::find(std::begin(arr2), std::end(arr2), mother_id_v1) != std::end(arr2);
            bool cond_aux = std::find(std::begin(arr2), std::end(arr2), grandmother_id_v1) != std::end(arr2);
            bool cond3 = (greatmother_id_v1 == 24) || (greatmother_id_v2 == 24) || (grandmother_id_v1 == 24) || (greatmother_id_v3 == 24);
            vector<bool> vb;
            vb.push_back(cond);
            vb.push_back(cond1);
            vb.push_back(cond2 || cond_aux);
            vb.push_back(cond3);
            return vb;
      };
      auto muonmother_aux(Vint pdg, Vint mother, const int muon_id){
            int mother_id_v1 = mother[muon_id];
            int grandmother_id_v1 = mother[mother[muon_id]];
            int greatmother_id_v1 = mother[mother[mother[muon_id]]];
            int greatmother_id_v2 = mother[mother[mother[mother[muon_id]]]];
            int greatmother_id_v3 = mother[mother[mother[mother[mother[muon_id]]]]];
            int mother_id_v2 = -10;
            if(fabs(pdg[mother_id_v1]) == 13){
                 mother_id_v1 = mother[mother[muon_id]];
                 grandmother_id_v1 = mother[mother[mother[muon_id]]];
                 greatmother_id_v1 = mother[mother[mother[mother[muon_id]]]];
                 greatmother_id_v2 = mother[mother[mother[mother[mother[muon_id]]]]];
                 greatmother_id_v3 = mother[mother[mother[mother[mother[mother[muon_id]]]]]];
            }
            int arr[] = {511,521,531,5122,5,10411,10421,413,423,10413,10423,20413,20423,415,425,10431,433,10433,20433,435,24};
            int arr1[] = {511,521,531,5122,5};
            int arr2[] = {411,421,431,4122};
            bool cond = std::find(std::begin(arr), std::end(arr), fabs(pdg[mother_id_v1])) != std::end(arr);
            bool cond1 = std::find(std::begin(arr1), std::end(arr1), fabs(pdg[mother_id_v1])) != std::end(arr1);
            bool cond2 = std::find(std::begin(arr2), std::end(arr2), fabs(pdg[mother_id_v1])) != std::end(arr2);
            bool cond_aux = std::find(std::begin(arr2), std::end(arr2), fabs(pdg[grandmother_id_v1])) != std::end(arr2);
            bool cond3 = (fabs(pdg[greatmother_id_v1]) == 24) || (fabs(pdg[greatmother_id_v2]) == 24) || (fabs(pdg[grandmother_id_v1]) == 24) || (fabs(pdg[greatmother_id_v3]) == 24);
            vector<int> vb;
            if (cond2) {
               vb.push_back(mother_id_v1);
            } else if (cond_aux) {
               vb.push_back(grandmother_id_v1);
            } else {
               vb.push_back(-1);
            }
            if (fabs(pdg[greatmother_id_v1]) == 24) {
               vb.push_back(greatmother_id_v1);
            } else if (fabs(pdg[greatmother_id_v2]) == 24) {
               vb.push_back(greatmother_id_v2);
            } else if (fabs(pdg[grandmother_id_v1]) == 24) {
               vb.push_back(grandmother_id_v1);
            } else if (fabs(pdg[greatmother_id_v3]) == 24) {
               vb.push_back(greatmother_id_v3);
            } else {
               vb.push_back(-1);
            }
            return vb;
      };
      auto part_id(UInt_t nPart, Vint status, Vint pdg, Vint mother) {
            vector<int> vb;
            int had_ind = -1;
            int muon_ind = -1;
            vector<bool> vb_aux;
            vector<int> ind_aux;
            for (unsigned int i=0; i<nPart; ++i) {
                vb_aux = muonmother(pdg, mother, i);
                if (fabs(pdg[i])==13 &&  vb_aux[2] && vb_aux[3]) {
                          muon_ind = i;
                          ind_aux = muonmother_aux(pdg, mother, i);
                }
            }
            if (muon_ind > -1) {
               had_ind = ind_aux[0];
            }
            vb.push_back(muon_ind);
            vb.push_back(had_ind);
            return vb;
      };
      auto ttbarcharm_particles(UInt_t nPart, Vint status, Vint pdg, Vint mother) {
            vector<int> vb;
            int had_ind = -1;
            int muon_ind = -1;
            int quark_ind = -1;
            bool typeC = false;
            bool typeL = false;
            int arr1[] = {1,2,3};
            bool cond = false;
            vector<int> ind_aux;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==4 && fabs(pdg[mother[i]])==24) {
                          typeC = true;
                          quark_ind = i;
                          ind_aux = part_id( nPart, status, pdg, mother);
                }
                cond = std::find(std::begin(arr1), std::end(arr1), fabs(pdg[i])) != std::end(arr1);
                if (cond && fabs(pdg[mother[i]])==24) {
                          typeL = true;
                }
            }
            if (quark_ind > -1) {
               muon_ind = ind_aux[0];
               if (muon_ind > -1) {
                   had_ind = ind_aux[1];
               }
            }
            vb.push_back(quark_ind);
            vb.push_back(muon_ind);
            vb.push_back(had_ind);
            return vb;
      };
      auto ttsllepton(UInt_t nPart, Vint status, Vint pdg, Vint mother) {
            int typelep = 0;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==13 && fabs(pdg[mother[i]])==24) {
                          typelep = 1;
                } else if (fabs(pdg[i])==11 && fabs(pdg[mother[i]])==24) {
                          typelep = 2;
                } else if (fabs(pdg[i])==15 && fabs(pdg[mother[i]])==24) {
                          typelep = 3;
                }
            }
            return typelep;
      };
""")

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


df_test = df_test.Define('charm_light_aux','ttbarcharm(nGenPart,GenPart_statusFlags,GenPart_pdgId,GenPart_genPartIdxMother)')
df_test = df_test.Define('top_aux','ttbartop(nGenPart,GenPart_statusFlags,GenPart_pdgId,GenPart_genPartIdxMother)')
df_test = df_test.Define('particles_aux','ttbarcharm_particles(nGenPart,GenPart_statusFlags,GenPart_pdgId,GenPart_genPartIdxMother)')
df_test = df_test.Define('quark_ind','particles_aux[0]')
df_test = df_test.Define('muon_ind','particles_aux[1]')
df_test = df_test.Define('had_ind','particles_aux[2]')
df_test = df_test.Define('quark_pt','quark_ind > -1 ? GenPart_pt[quark_ind] : -10.')
df_test = df_test.Define('muon_pt_aux','muon_ind > -1 ? GenPart_pt[muon_ind] : -10.')
df_test = df_test.Define('had_pt','had_ind > -1 ? GenPart_pt[had_ind] : -10.')
df_charm = df_test.Filter('charm_light_aux[0]')
df_nocharm = df_test.Filter('!charm_light_aux[0]')
df_top = df_test.Filter('top_aux')


hist_quark_pt_charm = df_charm.Histo1D(("hist_quark_pt_charm","",50,15,115),"quark_pt")
hist_quark_pt_nocharm = df_nocharm.Histo1D(("hist_quark_pt_nocharm","",50,15,115),"quark_pt")
hist_muon_pt_charm = df_charm.Histo1D(("hist_muon_pt_charm","",50,15,115),"muon_pt_aux")
hist_muon_pt_nocharm = df_nocharm.Histo1D(("hist_muon_pt_nocharm","",50,15,115),"muon_pt_aux")
hist_had_pt_charm = df_charm.Histo1D(("hist_had_pt_charm","",50,15,115),"had_pt")
hist_had_pt_nocharm = df_nocharm.Histo1D(("hist_had_pt_nocharm","",50,15,115),"had_pt")

print("Total events "+str(df_test.Count().GetValue()))
print("Top filtered events "+str(df_top.Count().GetValue()))

#############################
####     DATA SAVING     ####
#############################

path_hist = "/nfs/cms/vazqueze/higgssearch/fromJF/pythia_genpart_pt.root"
myfile = TFile( path_hist, 'RECREATE' )

hist_quark_pt_charm.Write()
hist_quark_pt_nocharm.Write()
hist_muon_pt_charm.Write()
hist_muon_pt_nocharm.Write()
hist_had_pt_charm.Write()
hist_had_pt_nocharm.Write()

myfile.Close()

print('Ended succesfully')
