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
            vector<bool> vb;
            for (unsigned int i=0; i<nPart; ++i) {
                vb.push_back(testSF(i,status[i],n));
            }
            return vb;
      };
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
      auto muongen(UInt_t nPart, Vint pdg, Vint mother, Vfloat eta, Vfloat phi, Vfloat pt, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vbool mu_id) {
            vector<int> vb_gen;
            vector<int> muind1;
            vector<int> muind2;
            vector<int> muind3;
            vector<int> muind4;
            vector<int> muind5;
            vector<int> muind6;
            vector<int> muind7;
            vector<int> muind8;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==13 && muonmotherverprev(pdg,mother,i)>0) {
                    if (pt[i] < 3.) {
                       muind1.push_back(i);
                    } else if (pt[i] < 4.) {
                       muind2.push_back(i);
                    } else if (pt[i] < 5.) {
                       muind3.push_back(i);
                    } else if (pt[i] < 7.) {
                       muind4.push_back(i);
                    } else if (pt[i] < 10.) {
                       muind5.push_back(i);
                    } else if (pt[i] < 15.) {
                       muind6.push_back(i);
                    } else if (pt[i] < 20.) {
                       muind7.push_back(i);
                    } else if (pt[i] < 25.) {
                       muind8.push_back(i);
                    }
                }
            }
            vb_gen.push_back(muind1.size());
            vb_gen.push_back(muind2.size());
            vb_gen.push_back(muind3.size());
            vb_gen.push_back(muind4.size());
            vb_gen.push_back(muind5.size());
            vb_gen.push_back(muind6.size());
            vb_gen.push_back(muind7.size());
            vb_gen.push_back(muind8.size());
            return vb_gen;
      };
      auto muonreco(UInt_t nPart, Vint pdg, Vint mother, Vfloat eta, Vfloat phi, Vfloat pt, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vbool mu_id) {
            vector<int> vb_gen;
            vector<int> vb_reco;
            vector<int> muind1;
            vector<int> muind2;
            vector<int> muind3;
            vector<int> muind4;
            vector<int> muind5;
            vector<int> muind6;
            vector<int> muind7;
            vector<int> muind8;
            vector<int> murecind1;
            vector<int> murecind2;
            vector<int> murecind3;
            vector<int> murecind4;
            vector<int> murecind5;
            vector<int> murecind6;
            vector<int> murecind7;
            vector<int> murecind8;
            bool cond1 = false;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==13 && muonmotherverprev(pdg,mother,i)>0) {
                    if (pt[i] < 3.) {
                       muind1.push_back(i);
                    } else if (pt[i] < 4.) {
                       muind2.push_back(i);
       	       	    } else if (pt[i] < 5.) {
                       muind3.push_back(i);
       	       	    } else if (pt[i] < 7.) {
                       muind4.push_back(i);
       	       	    } else if (pt[i] < 10.) {
                       muind5.push_back(i);
       	       	    } else if (pt[i] < 15.) {
                       muind6.push_back(i);
       	       	    } else if (pt[i] < 20.) {
                       muind7.push_back(i);
       	       	    } else if (pt[i] < 25.) {
                       muind8.push_back(i);
                    }
                }
            }
            for (unsigned int k=0; k<mu_pt.size(); ++k) {
                for (unsigned int j=0; j<muind1.size(); ++j) {
                    cond1 = ROOT::VecOps::DeltaR(mu_eta[k],eta[muind1[j]],mu_phi[k],phi[muind1[j]]) < 0.4;
                    if (mu_id[k] && cond1) {
                       murecind1.push_back(k);
                    }
                }
                for (unsigned int j=0; j<muind2.size(); ++j) {
                    cond1 = ROOT::VecOps::DeltaR(mu_eta[k],eta[muind2[j]],mu_phi[k],phi[muind2[j]]) < 0.4;
                    if (mu_id[k] && cond1) {
                       murecind2.push_back(k);
                    }
                }
                for (unsigned int j=0; j<muind3.size(); ++j) {
                    cond1 = ROOT::VecOps::DeltaR(mu_eta[k],eta[muind3[j]],mu_phi[k],phi[muind3[j]]) < 0.4;
                    if (mu_id[k] && cond1) {
                       murecind3.push_back(k);
                    }
                }
                for (unsigned int j=0; j<muind4.size(); ++j) {
                    cond1 = ROOT::VecOps::DeltaR(mu_eta[k],eta[muind4[j]],mu_phi[k],phi[muind4[j]]) < 0.4;
                    if (mu_id[k] && cond1) {
                       murecind4.push_back(k);
                    }
                }
                for (unsigned int j=0; j<muind5.size(); ++j) {
                    cond1 = ROOT::VecOps::DeltaR(mu_eta[k],eta[muind5[j]],mu_phi[k],phi[muind5[j]]) < 0.4;
                    if (mu_id[k] && cond1) {
                       murecind5.push_back(k);
                    }
                }
                for (unsigned int j=0; j<muind6.size(); ++j) {
                    cond1 = ROOT::VecOps::DeltaR(mu_eta[k],eta[muind6[j]],mu_phi[k],phi[muind6[j]]) < 0.4;
                    if (mu_id[k] && cond1) {
                       murecind6.push_back(k);
                    }
                }
                for (unsigned int j=0; j<muind7.size(); ++j) {
                    cond1 = ROOT::VecOps::DeltaR(mu_eta[k],eta[muind7[j]],mu_phi[k],phi[muind7[j]]) < 0.4;
                    if (mu_id[k] && cond1) {
                       murecind7.push_back(k);
                    }
                }
                for (unsigned int j=0; j<muind8.size(); ++j) {
                    cond1 = ROOT::VecOps::DeltaR(mu_eta[k],eta[muind8[j]],mu_phi[k],phi[muind8[j]]) < 0.4;
                    if (mu_id[k] && cond1) {
                       murecind8.push_back(k);
                    }
                }
            }
            vb_reco.push_back(murecind1.size());
            vb_reco.push_back(murecind2.size());
            vb_reco.push_back(murecind3.size());
            vb_reco.push_back(murecind4.size());
            vb_reco.push_back(murecind5.size());
            vb_reco.push_back(murecind6.size());
            vb_reco.push_back(murecind7.size());
            vb_reco.push_back(murecind8.size());
            return vb_reco;
      };
""")

### SL selection
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
      auto SLselection(Vint qind, Vfloat pt, Vfloat eta, Vfloat phi, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vint mu_good, Vbool mu_id){
            vector<int> vb;
            bool cond1 = false;
            bool cond2 = false;
            bool condmu = true;
            int indJ=-1;
            int indM=-1;
            float ptM{-10.};
            for (unsigned int i=0; i<mu_pt.size(); ++i){
                cond1 = ROOT::VecOps::DeltaR(mu_eta[i],eta[qind[0]],mu_phi[i],phi[qind[0]]) < 0.4;
                cond2 = ROOT::VecOps::DeltaR(mu_eta[i],eta[qind[1]],mu_phi[i],phi[qind[1]]) < 0.4;
                if (mu_good.size() > 0) condmu = mu_good[0] != i;
                if(condmu && mu_id[i] && mu_pt[i]<25. && (cond1 || cond2) && mu_pt[i]>ptM && mu_pt[i]>3.){
                     indM = i;
                     ptM = mu_pt[i];
                     if (cond1) {
                        indJ = qind[0];
                     } else if (cond2) {
                        indJ = qind[1];
                     }
                }
            }
            vb.push_back(indM);
            vb.push_back(indJ);
            return vb;
      };
""")

######### number of muons in bottom jets

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto bottomMuons(const int blep, const int bhad, Vfloat pt, Vfloat eta, Vfloat phi, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vint mu_good, Vbool mu_id){
            vector<int> vb;
            bool cond1 = false;
            bool cond2 = false;
            bool condmu = true;
            int muon_mul_lep = 0;
            int indMlep=-1;
            float ptMlep{-10.};
            int muon_mul_had = 0;
            int indMhad=-1;
            float ptMhad{-10.};
            for (unsigned int i=0; i<mu_pt.size(); ++i){
                if (mu_good.size() > 0) condmu = mu_good[0] != i;
                cond1 = ROOT::VecOps::DeltaR(mu_eta[i],eta[blep],mu_phi[i],phi[blep]) < 0.4;
                cond2 = ROOT::VecOps::DeltaR(mu_eta[i],eta[bhad],mu_phi[i],phi[bhad]) < 0.4;
                if (condmu && mu_id[i] && mu_pt[i]<25. && cond1 && mu_pt[i]>3.){
                   muon_mul_lep = muon_mul_lep+1;
                   if (mu_pt[i]>ptMlep) {
                     indMlep = i;
                     ptMlep = mu_pt[i];
                   }
                } else if (condmu && mu_id[i] && mu_pt[i]<25. && cond2 && mu_pt[i]>3.) {
                   muon_mul_had = muon_mul_had+1;
                   if (mu_pt[i]>ptMhad) {
                     indMhad = i;
                     ptMhad = mu_pt[i];
                   }
                }
            }
            vb.push_back(muon_mul_lep);
            vb.push_back(indMlep);
            vb.push_back(muon_mul_had);
            vb.push_back(indMhad);
            return vb;
      };
""")

### mother of muon

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonmother(Vint pdg, Vint mother, Vint muon_id_gen, const int muon_id){
            int mother_id_v1 = fabs(pdg[mother[muon_id_gen[muon_id]]]);
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
""")


#############################################################
######################## Definitions ########################
#############################################################

gInterpreter.Declare("""
   Double_t xbins[10] = {0.,3.,4.,5.,7.,10.,15.,20.,25.,100.};
""")


df_test = df_test.Define('vec_gen','muongen(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, GenPart_pt, Muon_pt, Muon_eta, Muon_phi, Muon_tightId)')
df_test = df_test.Define('vec_reco','muonreco(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_eta, GenPart_phi, GenPart_pt, Muon_pt, Muon_eta, Muon_phi, Muon_tightId)')
df_test = df_test.Define('aux_gen1','vec_gen[0]')
df_test = df_test.Define('aux_gen2','vec_gen[1]')
df_test = df_test.Define('aux_gen3','vec_gen[2]')
df_test = df_test.Define('aux_gen4','vec_gen[3]')
df_test = df_test.Define('aux_gen5','vec_gen[4]')
df_test = df_test.Define('aux_gen6','vec_gen[5]')
df_test = df_test.Define('aux_gen7','vec_gen[6]')
df_test = df_test.Define('aux_gen8','vec_gen[7]')
df_test = df_test.Define('aux_gen9','vec_gen[8]')
df_test = df_test.Define('aux_reco1','vec_reco[0]')
df_test = df_test.Define('aux_reco2','vec_reco[1]')
df_test = df_test.Define('aux_reco3','vec_reco[2]')
df_test = df_test.Define('aux_reco4','vec_reco[3]')
df_test = df_test.Define('aux_reco5','vec_reco[4]')
df_test = df_test.Define('aux_reco6','vec_reco[5]')
df_test = df_test.Define('aux_reco7','vec_reco[6]')
df_test = df_test.Define('aux_reco8','vec_reco[7]')
df_test = df_test.Define('aux_reco9','vec_reco[8]')
df_test = df_test.Define('aux_pt1','1.')
df_test = df_test.Define('aux_pt2','3.5')
df_test = df_test.Define('aux_pt3','4.5')
df_test = df_test.Define('aux_pt4','6.')
df_test = df_test.Define('aux_pt5','8.')
df_test = df_test.Define('aux_pt6','12.')
df_test = df_test.Define('aux_pt7','17.')
df_test = df_test.Define('aux_pt8','22.')
df_test = df_test.Define('aux_pt9','30.')

#histgen1 = df_test.Histo1D(("histgen1","",9,xbins),"aux_pt1","aux_gen1")
#histgen2 = df_test.Histo1D(("histgen2","",9,xbins),"aux_pt2","aux_gen2")
histgen1 = df_test.Histo1D(("histgen1","",10,0,10),"aux_pt1","aux_gen1")
histgen2 = df_test.Histo1D(("histgen2","",10,0,10),"aux_pt2","aux_gen2")


#############################
####     DATA SAVING     ####
#############################

myratiofile = TFile('/nfs/cms/vazqueze/higgssearch/fromJF/ratios_muons_botjet/muon_eff_files.root', 'RECREATE')

myratiofile.WriteObject(histgen1,"first_hist_gen")
myratiofile.WriteObject(histgen2,"second_hist_gen")

myratiofile.Close()



print('Ended succesfully')
