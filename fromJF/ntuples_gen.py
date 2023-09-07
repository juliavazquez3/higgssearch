######################################                                        
######################################                                        
#######      GOOD  VERSION     #######
######################################                                     
######################################                                     

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join, isdir

import json
import argparse

sys.path.append('/nfs/cms/vazqueze/ttbaranalisis/last_corrections/')

df = ROOT.RDataFrame("Events","/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux/myconfig2018/wjets_3/cat_base/prod_test/data_30.root")

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
            vector<int> vb;
            for (unsigned int i=0; i<nPart; ++i) {
                vb.push_back(testSF(i,status[i],n));
            }
            return vb;
      };
""")

# Attempt to classify with GenPart
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using Vbool = const ROOT::RVec<bool>&;
      using namespace std;
      auto counting_dhadrons(UInt_t nGenPart, Vint pdg, Vint mother, Vfloat gen_eta, Vfloat gen_phi, Vfloat gen_pt, Vfloat jet_pt, Vint jetind, Vbool lastC) {
            vector<int> vb;
            int ind1 = 0; // D+
            int ind2 = 0; // D0
            int ind3 = 0; // D+s
            int ind4 = 0; // lambda
            int ind5 = 0; // rest
            int ind6 = 0; // charm quarks
            vector<int> vdh = {10411,10421,413,423,10413,10423,20413,20423,415,425,10431,433,10433,20433,435};
            for (unsigned int i=0; i<nGenPart; ++i) {
              if (lastC[i]) {
                 if (fabs(pdg[i]) == 411) {
                    ind1 = ind1 +1;
                 } else if (fabs(pdg[i]) == 421) {
                    ind2 = ind2 +1;
                 } else if (fabs(pdg[i]) == 431) {
                    ind3 = ind3 +1;
                 } else if (fabs(pdg[i]) == 4122) {
                    ind4 = ind4 +1;
                 } else if (std::find(vdh.begin(), vdh.end(), fabs(pdg[i])) != vdh.end()) {
                    ind5 = ind5 +1;
                 } else if (fabs(pdg[i]) == 4) {
                    ind6 = ind6+1;
                 }
              }
            }
            vb.push_back(ind1);
            vb.push_back(ind2);
            vb.push_back(ind3);
            vb.push_back(ind4);
            vb.push_back(ind5);
            vb.push_back(ind6);
            return vb;
      };
""")

#############################################################
######################## Definitions ########################
#############################################################


df = df.Define('ishard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,7)')
df = df.Define('fromhard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,8)')
df = df.Define('first_copy','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,12)')
df = df.Define('last_copy','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,13)')

brlist = ["nGenJet", "GenJet_eta", "GenJet_hadronFlavour", "GenJet_mass", "GenJet_partonFlavour", "GenJet_phi", "GenJet_pt",
       "GenMET_phi", "GenMET_pt", "nGenPart","GenPart_eta", "GenPart_genPartIdxMother", "GenPart_mass", "GenPart_pdgId", "GenPart_phi",
       "GenPart_pt", "GenPart_status", "GenPart_statusFlags","ishard","first_copy","last_copy","fromhard"]

path = '/pnfs/ciemat.es/data/cms/store/user/juvazque/data_test/wjets_3'
df.Snapshot("Events", path+".root", brlist)

