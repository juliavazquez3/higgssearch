import ROOT, os, sys
from ROOT import *

import numpy as np
import awkward as ak
import uproot
import pandas as pd
import pyarrow as pa
import urllib.request


df_init = ROOT.RDataFrame("Events","/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux/myconfig2018/ttbar_dl/cat_base/prod_test/data_38.root")

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
      auto topreweight(UInt_t nPart, Vint status, Vint pdg, Vbool lastC, Vfloat gen_pt) {
            float wei = 1.;
            float pt1 = 0.;
            float pt2 = 0.;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==6 && lastC[i]) {
                          if (pdg[i] == 6) pt1 = gen_pt[i];
                          if (pdg[i] == -6) pt2 = gen_pt[i];
                }
            }
            if (pt1 != 0. && pt2 != 0.) {
                   wei = std::sqrt(exp(0.0615-(0.0005*pt1))*exp(0.0615-(0.0005*pt2)));
            }
            return wei;
      };
""")

df_init = df_init.Define('last_copy','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,13)')
df_init = df_init.Define('ishard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,7)')
df_init = df_init.Define('fromhard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,8)')

brlist = ["nGenPart","GenPart_pdgId","GenPart_genPartIdxMother","ishard","fromhard","last_copy"]

df_init.Snapshot("Events", "file_ttbar_dl.root", brlist)

#file = uproot.open("/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/260000/062085CD-8DF4-1D40-8042-63998B6A3A95.root")
file = uproot.open("file_ttbar_dl.root")
print(file.keys())
print(file.classnames())

events = file['Events']

pdgID = events["GenPart_pdgId"].array()
motherID = events["GenPart_genPartIdxMother"].array()
isHardP = events["ishard"].array()
fromHardP = events["fromhard"].array()

ak.to_pandas(pdgID)
ak.to_pandas(motherID)
ak.to_pandas(isHardP)
ak.to_pandas(fromHardP)

print('--------------------------------')
print(pdgID[1271][0:15])
print(motherID[1271][0:15])
print(isHardP[1271][0:15])
print(fromHardP[1271][0:15])

print('--------------------------------')
print(pdgID[1271][15:27])
print(motherID[1271][15:27])
print(isHardP[1271][15:27])
print(fromHardP[1271][15:27])

print('--------------------------------')
print(pdgID[1581][15:27])
print(motherID[1581][15:27])
print(isHardP[1581][15:27])
print(fromHardP[1581][15:27])

print('--------------------------------')
print(pdgID[11201][15:27])
print(motherID[11201][15:27])
print(isHardP[11201][15:27])
print(fromHardP[11201][15:27])

print('--------------------------------')
print(pdgID[21803][15:27])
print(motherID[21803][15:27])
print(isHardP[21803][15:27])
print(fromHardP[21803][15:27])

print('--------------------------------')
print(pdgID[2104][0:22])
print(motherID[2104][0:22])
print(isHardP[2104][0:22])
print(fromHardP[2104][0:22])

print('--------------------------------')
print(pdgID[108011][0:18])
print(motherID[108011][0:18])
print(isHardP[108011][0:18])
print(fromHardP[108011][0:18])

print('--------------------------------')
print(pdgID[15811][0:22])
print(motherID[15811][0:22])
print(isHardP[15811][0:22])
print(fromHardP[15811][0:22])

print('--------------------------------')
print(pdgID[2702][0:22])
print(motherID[2702][0:22])
print(isHardP[2702][0:22])
print(fromHardP[2702][0:22])

print('--------------------------------')
print(pdgID[2184][0:22])
print(motherID[2184][0:22])
print(isHardP[2184][0:22])
print(fromHardP[2184][0:22])

print('--------------------------------')
print(pdgID[2107][0:22])
print(motherID[2107][0:22])
print(isHardP[2107][0:22])
print(fromHardP[2107][0:22])

print('--------------------------------')
print(pdgID[1082][0:18])
print(motherID[1082][0:18])
print(isHardP[1082][0:18])
print(fromHardP[1082][0:18])

