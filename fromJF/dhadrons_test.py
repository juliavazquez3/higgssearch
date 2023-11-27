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
   #include <vector>
   TFile *file1 = new TFile("higgssearch/fromJF/ratios_muons_botjet/hists_confidenceintervals/confidence_intervals_etabinfive.root");
   TH1F * h1 = (TH1F*)file1->Get("hist_test_bin_1_new");
   std::vector<float> arr;
   int limi = h1->GetNbinsX();
   auto testC(UInt_t nPart) {
            for (unsigned int i=0; i<limi; ++i) {
                arr.push_back(h1->GetBinError(i+1));
            }
            float vb;
            for (unsigned int i=0; i<nPart; ++i) {
                vb = arr[8];
            }
            return vb;
   };  
""")

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

df_init = df_init.Define('testC','testC(nGenPart)')
histd0 = df_init.Histo1D(("d0","",10,0,1),"testC")

print("%s" %histd0.Integral())
