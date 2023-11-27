import ROOT, os, sys
from ROOT import *

import numpy as np
import awkward as ak
import uproot
import pandas as pd
import pyarrow as pa
import urllib.request


#df_init = ROOT.RDataFrame("Events","/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux/myconfig2018/ttbar_dl/cat_base/prod_test/data_38.root")
df_init = ROOT.RDataFrame("Events","/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux/myconfig2018/st_3/cat_base/prod_test/data_0.root")
#df_init = ROOT.RDataFrame("Events","/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux/myconfig2018/ttbar_sl/cat_base/prod_test/data_38.root")

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
      auto dhad_count(UInt_t nPart, Vint mother, Vint pdg, Vbool lastC) {
            int ind_d0 = 0;
            int ind_dplus = 0;
            int ind_ds = 0;
            int ind_lambda = 0;
            vector<int> vb;
            bool cond = false;
            bool cond1 = false;
            int moth;
            int moth1;
            int arr[] = {511,521,531,5122,5,10411,10421,413,423,10413,10423,20413,20423,415,425,10431,433,10433,20433,435,24};
            int arr1[] = {511,521,531,5122,5};
            for (unsigned int i=0; i<nPart; ++i) {
                moth = fabs(pdg[mother[i]]);
                moth1 = fabs(pdg[mother[mother[i]]]);
                cond = std::find(std::begin(arr), std::end(arr), moth) != std::end(arr);
                cond1 = (std::find(std::begin(arr1), std::end(arr1), moth) != std::end(arr1)) ? true : std::find(std::begin(arr1), std::end(arr1), moth1) != std::end(arr1);
                if (fabs(pdg[i])==421 && lastC[i] && cond&& cond1) {
                     ind_d0 = ind_d0 + 1;
                } else if (fabs(pdg[i])==411 && lastC[i] && cond && cond1) {
                     ind_dplus = ind_dplus+1;
                } else if (fabs(pdg[i])==431 && lastC[i] && cond && cond1) {
                     ind_ds = ind_ds+1;
                } else if (fabs(pdg[i])==4122 && lastC[i] && cond && cond1) {
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

df_init = df_init.Define('first_copy','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,12)')
df_init = df_init.Define('ishard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,7)')
df_init = df_init.Define('fromhard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,8)')
df_init = df_init.Define('dhadrons_quant','dhad_count(nGenPart, GenPart_genPartIdxMother, GenPart_pdgId, first_copy)')

df_init = df_init.Define('d0_quant','dhadrons_quant[0]')
df_init = df_init.Define('dplus_quant','dhadrons_quant[1]')
df_init = df_init.Define('ds_quant','dhadrons_quant[2]')
df_init = df_init.Define('lambda_quant','dhadrons_quant[3]')
df_init = df_init.Define('dhadrons_aux','1.')

histd0 = df_init.Histo1D(("d0","",2,0,2),"dhadrons_aux","d0_quant")
histdplus = df_init.Histo1D(("dplus","",2,0,2),"dhadrons_aux","dplus_quant")
histds = df_init.Histo1D(("ds","",2,0,2),"dhadrons_aux","ds_quant")
histlambda = df_init.Histo1D(("lambda","",2,0,2),"dhadrons_aux","lambda_quant")


ent = df_init.Count().GetValue()
d0q = histd0.Integral()
dplusq = histdplus.Integral()
dsq = histds.Integral()
lambdaq = histlambda.Integral()

print("Number of entries is %s" %ent)
print("Number of D0 hadrons is %s" %d0q)
print("Number of D+ hadrons is %s" %dplusq)
print("Number of Ds hadrons is %s" %dsq)
print("Number of D lambda hadrons is %s" %lambdaq)
print("Total number of D hadrons is %s" %(d0q+dplusq+dsq+lambdaq))
print("Proportion of D hadrons with respect to entries %s" %((d0q+dplusq+dsq+lambdaq)/ent))

