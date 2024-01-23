import ROOT, os, sys
from ROOT import *

import numpy as np
import awkward as ak
import uproot
import pandas as pd
import pyarrow as pa
import urllib.request


#df_init = ROOT.RDataFrame("Events","/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux/myconfig2018/ttbar_dl/cat_base/prod_test/data_38.root")
df_init = ROOT.RDataFrame("Events","/pnfs/ciemat.es/data/cms/store/user/juvazque/data_further_analysis/btagMM/chitest/sl/folder2016/wcs/dataset_wqq_btagMM_fromJF_ttbar_sl2016_charm.root")

histfile = TFile.Open("higgssearch/fromJF/ratios_muons_botjet/hists_confidenceintervals/confidence_intervals_etabins_all.root","READ")

h1 = histfile.Get("etabinone")
limi = h1.GetNbinsX()

isoabs_edges = [2.5,5.,7.5,10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,32.5,35.,37.5,40.,45.,50.,60.,70.,80.,100.]
etabins = [0,1,2,3,4]

dict_to_cpp = "{"
for k,keta in enumerate(etabins):
  dict_to_cpp += "{%s, {" % etabins[k]
  for i,isoed in enumerate(isoabs_edges):
     dict_to_cpp += "{%s, %s, %s} %s" % (h1.GetBinContent(i+1), (h1.GetBinContent(i+1))+(h1.GetBinError(i+1)), (h1.GetBinContent(i+1))-(h1.GetBinError(i+1)),
      ("," if i< len(isoabs_edges)-1 else "}"))
  dict_to_cpp +="}%s" % (", " if k < len(etabins) - 1 else "}")

gInterpreter.Declare("""
   #include <vector>
   //TFile *file1 = new TFile("higgssearch/fromJF/ratios_muons_botjet/hists_confidenceintervals/confidence_intervals_etabins_all.root");
   //TH1F * h1 = (TH1F*)file1->Get("etabinone");
   std::map<int, std::vector<std::vector<float>>> m = %s;
   using Vfloat = const ROOT::RVec<float>&;
   std::vector<float> arr = {2.5,5.,7.5,10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,32.5,35.,37.5,40.,45.,50.,60.,70.,80.,100.};
   auto testC(const float muon_pt, const float muon_eta) {
            int binx = 0;
            float res;
            for (unsigned int i=0; i<arr.size(); ++i) {
               if (muon_pt >= arr[i]) {
                  binx = i+1;
               }
               if (binx == 22) binx = 21;
            }
            res = m[0][binx][0];
            return res;
   };  
""" %dict_to_cpp)

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
""")

df_init = df_init.Define('last_copy','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,13)')
df_init = df_init.Define('ishard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,7)')
df_init = df_init.Define('fromhard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,8)')

df_mu = df_init.Define('muon_jet_mother_aux','muonmother(GenPart_pdgId, GenPart_genPartIdxMother, muon_jet_genindx)')
df_muB = df_mu.Filter('muon_jet_mother_aux[1]')
df_muD = df_mu.Filter('muon_jet_mother_aux[2]')
df_muW = df_mu.Filter('muon_jet_mother_aux[3]')

df_init = df_init.Define('testC','nMuon>0 ? testC(Muon_pt[0], Muon_eta[0]) : 0.')
#df_init = df_init.Define('testC','testC_aux.size() > 0 ? testC_aux[0] : 0.')
histd0 = df_init.Histo1D(("d0","",30,0,1),"testC")

print("%s" %histd0.Integral())

print("Entries are %s" %df_mu.Count().GetValue())
print("Entries with muon from B hadron are %s" %df_muB.Count().GetValue())
print("Entries with muon from D hadron are %s" %df_muD.Count().GetValue())
print("Entries with muon from W boson are %s" %df_muW.Count().GetValue())
