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

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

ROOT.EnableImplicitMT()


df = ROOT.RDataFrame("Events","/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_ver2/myconfig2016/ttbar_sl/cat_base/prod_test/data_0.root")

###################################################
################   DEFINITIONS   ##################
###################################################

#### pT top reweight ##########

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
      auto wpluscbool(UInt_t nPart, Vint status, Vint pdg, Vbool hardP, Vbool firstC) {
            int typeC = 0;
            int indC = 0;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==5 && hardP[i]) {
                          typeC = 3;
                } else {
                     if (fabs(pdg[i])==4 && hardP[i] && firstC[i]) indC++;
                }
            }
            if (typeC!=3 && (indC % 2 != 0)) typeC = 2;
            if (typeC!=3 && (indC % 2 == 0) && (indC > 0)) typeC = 1;
            return typeC;
      };
""")

#############################################################
######################## Definitions ########################
#############################################################

df = df.Define('status7','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,7)')
df = df.Define('status8','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,8)')
df = df.Define('status9','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,9)')
df = df.Define('status11','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,11)')
df = df.Define('status12','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,12)')
df = df.Define('status13','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,13)')

#############################
####     DATA SAVING     ####
#############################

brlist = [
'nGenPart','GenPart_pdgId','status7','status8','status9','status11','status12','status13','GenPart_statusFlags'
]

path = '/tmp/hsperfdata_vazqueze/dataset_test_genpart_ttbarsl'
df.Snapshot("Events", path+".root", brlist)

print('Ended succesfully')

