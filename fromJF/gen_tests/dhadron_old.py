import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join, isdir

import json
import argparse

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


gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto dhad_count(UInt_t nPart, Vint mother, Vint pdg, Vbool lastC) {
            vector<int> vb;
            int ind_d0 = 0;
            int ind_dplus = 0;
            int ind_ds = 0;
            int ind_lambda = 0;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==411 && lastC[i]) {
                   ind_dplus = ind_dplus+1;
                } else if (fabs(pdg[i])==421 && lastC[i]) {
                   ind_d0 = ind_d0+1;
                } else if (fabs(pdg[i])==431 && lastC[i]) {
                   ind_ds = ind_ds+1;
                } else if (fabs(pdg[i])==4122 && lastC[i]) {
                   ind_lambda = ind_lambda+1;
                }
            }
            vb.push_back(ind_d0);
            vb.push_back(ind_dplus);
            vb.push_back(ind_ds);
            vb.push_back(ind_lambda);
            return vb;
      };
""")

df_init = df_init.Define('last_copy','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,13)')
df_init = df_init.Define('dhadrons_quants','dhad_count(nGenPart, GenPart_genPartIdxMother, GenPart_pdgId, last_copy)')
df_init = df_init.Define('quant_d0','dhadrons_quants[0]')
df_init = df_init.Define('quant_dplus','dhadrons_quants[1]')
df_init = df_init.Define('quant_ds','dhadrons_quants[2]')
df_init = df_init.Define('quant_lambda','dhadrons_quants[3]')
df_init = df_init.Define('quant_aux','1.')

entries1 = df_init.Count()

hqd0 = df_init.Histo1D(("hqd0","",2,0,2),"quant_aux","quant_d0")
hqdplus = df_init.Histo1D(("hqdplus","",2,0,2),"quant_aux","quant_dplus")
hqds = df_init.Histo1D(("hqds","",2,0,2),"quant_aux","quant_ds")
hqlambda = df_init.Histo1D(("hqlambda","",2,0,2),"quant_aux","quant_lambda")

print("%s total entries" %entries1.GetValue())
print("%s number of D0 hadrons" %hqd0.Integral())
print("%s number of Dplus hadrons" %hqdplus.Integral())
print("%s number of Ds hadrons" %hqds.Integral())
print("%s number of lambda D hadrons" %hqlambda.Integral())
