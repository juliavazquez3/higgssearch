######################################
######################################
#######      NAIVE VERSION     #######
######################################
######################################

print('SV CHANNEL')

## Modified version, meant to be run locally, not sent to condor 

## Selection for noth MC and data samples, created to distinguish between years
## Different histogram files are produced for each situation
## It includes an option for SSOS substraction

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join

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

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--process", type=string, default="WW",
                    help="Select type of process to run")
parser.add_argument("--year", type=string, default="2016",
                    help="Select year of process to run")
parser.add_argument("--type", type=string, default="data",
                    help="Selec type of data to run")
parser.add_argument("--notfull", action="store_true", default=False,
                    help="Use just a range of the sample")
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Perform ssos substraction")
parser.add_argument('-l','--list', nargs='+', help='range of sample to use')
# Use like:
# python arg.py -l 1234 2345 3456 4567

args = parser.parse_args()

if args.process == "allMCH": proc = ["higgs_tocs_m80","higgs_tocs_m90","higgs_tocs_m100","higgs_tocs_m120","higgs_tocs_m140","higgs_tocs_m150",
           "higgs_tocs_m155","higgs_tocs_m160"]
elif args.process == "allMC": proc = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6","zjets_7","zjets_8",
        "ttbar_sl","ttbar_dl","ttbar_dh","zz","wz","st_1","st_2","st_3","st_4"]
elif (args.process == "higgs_tocs_m80" or args.process == "higgs_tocs_m90" or args.process == "higgs_tocs_m100" or args.process == "higgs_tocs_m120"
        or args.process == "higgs_tocs_m140" or args.process == "higgs_tocs_m150"
        or args.process == "higgs_tocs_m155" or args.process == "higgs_tocs_m160"
        or args.process == "ww" or args.process == "wjets_1" or args.process == "wjets_2" or args.process == "wjets_3" or args.process == "wjets_4"
        or args.process == "ttbar_sl" or args.process == "wjets_5" or args.process == "wjets_6" or args.process == "wjets_7" or args.process == "wjets_8"
        or args.process == "zjets_1" or args.process == "zjets_2" or args.process == "zjets_3" or args.process == "zjets_4"
        or args.process == "zjets_5" or args.process == "zjets_6" or args.process == "zjets_7" or args.process == "zjets_8"
        or args.process == "DY" or args.process == "WZ"or args.process == "ZZ" or args.process == "ST1"
        or args.process == "ST2" or args.process == "ST3" or args.process == "ST4"  or args.process == "DYjets"  or args.process == "ttbar_dl"
        or args.process == "ttbar_dh" or args.process == "M" or args.process == "E"): proc = [str(args.process)]
else: raise NameError('Incorrect process name')

if args.year == "all": years = ["2016","2016B","2017","2018"]
elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): years = [str(args.year)]
else: raise NameError('Incorrect year')

if args.type == "data": mode = "data"
elif args.type == "mc": mode = "mc"
else: raise NameError('Incorrect type')

special_weights = ["Wjets0J2016","Wjets1J2016","Wjets2J2016","DY0J2016","DY1J2016","DY2J2016", 
        "Wjets0J2016B","Wjets1J2016B","Wjets2J2016B","DY0J2016B","DY1J2016B","DY2J2016B", 
        "Wjets0J2017","Wjets1J2017","Wjets2J2017","DY0J2017","DY1J2017","DY2J2017", 
        "Wjets0J2018","Wjets1J2018","Wjets2J2018","DY0J2018","DY1J2018","DY2J2018"]

samples = []

for p in proc:
  for y in years:
    if mode == "mc": samples.append(p+y)
    else: samples.append(y+p)

print("the samples treated are",samples)

samples_test = samples[:]

if args.notfull:
	if len(args.list) != 2: raise NameError('List has to have 2 elements')
	print("the range of files is",args.list)
 

# Create a ROOT dataframe for each dataset
# Note that we load the filenames from the external json file placed in the same folder than this script.
# Example "python analisisWW/selection_v2.py --process="WW" --notfull -l 0 50"

if mode == "mc": files = json.load(open("/nfs/cms/vazqueze/higgssearch/mcinfo"+years[0]+".json"))
else: files = json.load(open("/nfs/cms/vazqueze/higgssearch/datainfo.json"))

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

if mode == "mc":
  for p in proc:
    # Construct the dataframes
    folder = files[p]["folder_name"] # Folder name
    list1 = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(folder)) for f in fn]
    #print("list1 is %s" %(list1))
    if mode == "mc": num_events = files[p]["events"] # Number of events
    else: num_events = files[p]["events_total"] # Number of events
    num_files = files[p]["files"] # Number of files
    list_files = [f for f in list1 if (isfile(f) and f[-5:] == ".root")] # Lista de archivos
    #print("list_files is %s" %(list_files))
    print(len(list_files))
    #if (num_files == len(list_files)):
    for f in list_files:
          archives[p+years[0]].append(f)
else:
  for p in names_proc[str(args.process)][years[0]]:
    # Construct the dataframes
    folder = files[p]["folder_name"] # Folder name
    list1 = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(folder)) for f in fn]
    #print("list1 is %s" %(list1))
    if mode == "mc": num_events = files[p]["events"] # Number of events
    else: num_events = files[p]["events_total"] # Number of events
    num_files = files[p]["files"] # Number of files
    list_files = [f for f in list1 if (isfile(f) and f[-5:] == ".root")] # Lista de archivos
    #print("list_files is %s" %(list_files))
    print(len(list_files))
    #if (num_files == len(list_files)):
    for f in list_files:
          archives[p[:-1]].append(f)

for s in samples:
  if args.notfull: archives[s]=archives[s][int(args.list[0]):int(args.list[1])]
  df[s] = ROOT.RDataFrame("Events",set(archives[s]))
  print("Number of files for",s,len(archives[s]))

## Cuts per year

cuts_btag = {}

cuts_btag["2016"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2016B"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2017"]=[0.0532, 0.3040, 0.7476]
cuts_btag["2018"]=[0.0490,0.2783,0.7100]

## Muon and Electron pT cuts per year

muon_pt = {}

muon_pt["2016"]=26
muon_pt["2016B"]=26
muon_pt["2017"]=29
muon_pt["2018"]=29

el_pt = {}

el_pt["2016"]=30
el_pt["2016B"]=30
el_pt["2017"]=35
el_pt["2018"]=35

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

#######################################
########      GEN_PART      ###########
#######################################

# Attempt to classify with GenPart
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;
 
// Function to print the
// index of an element
// 1: Full hadronic 2: One W hadronic, the other leptonic (SEMI) 3: Full leptonic
      auto getIndex(UInt_t nGenPart,Vint v, int K){
            int ind = 2;
            while (v[ind]!=K) {
                ++ind;
            }
            return ind;
      }
      auto typeWW(UInt_t nGenPart, Vint partId, Vint motherId) {
	    int type = -1;
	    auto c = (partId == 24)||(partId == -24);
	    std::vector<int> v1(size(partId));
	    std::iota(std::begin(v1),std::end(v1),0);
	    ROOT::VecOps::RVec<int> myRVec(v1.data(), v1.size());
	    int v2 = -1;
	    auto Windexes = ROOT::VecOps::Where(c,myRVec,v2);
	    int idWlast1 = ROOT::VecOps::Max(Windexes);
	    int idWlast2 = idWlast1-1;
	    int idpart1 = getIndex(nGenPart,motherId, idWlast1);
            int idpart2 = getIndex(nGenPart,motherId, idWlast2);
	    
	    if (motherId[idpart1]==motherId[idpart1+1] && motherId[idpart2]==motherId[idpart2+1]){
		if (fabs(partId[idpart1])< 9 && fabs(partId[idpart2])< 9) {
  			type = 1;
		} else if (fabs(partId[idpart1])< 9 || fabs(partId[idpart2])< 9) {
  			type = 2;
		} else {
  			type = 3;
		}
	    }
    	    return type;
      };
""")


# Attempt to classify with GenPart
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;

// Function to print the
// index of an element
// 1: Charm -1: No Charm
      auto typeC(UInt_t nGenPart,Vint partId, Vint motherId){
	    int type = -1;
	    auto c = (partId == 24)||(partId == -24);
            std::vector<int> v1(size(partId));
            std::iota(std::begin(v1),std::end(v1),0);
            ROOT::VecOps::RVec<int> myRVec(v1.data(), v1.size());
            int v2 = -1;
            auto Windexes = ROOT::VecOps::Where(c,myRVec,v2);
            int idWlast1 = ROOT::VecOps::Max(Windexes);
            int idWlast2 = idWlast1-1;
            int idpart1 = getIndex(nGenPart,motherId, idWlast1);
            int idpart2 = getIndex(nGenPart,motherId, idWlast2);

            if (motherId[idpart1]==motherId[idpart1+1] && motherId[idpart2]==motherId[idpart2+1]){
                if ((fabs(partId[idpart1])==4 || fabs(partId[idpart2])==4) || (fabs(partId[idpart1+1])==4 || fabs(partId[idpart2+1])==4)) {
                        type = 1;
                }
            }
            return type;
      };
""")

######## Gen identification for W plus jets

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
      auto ttbarcharm(UInt_t nPart, Vint status, Vint pdg, Vint mother) {
            bool typeC = false;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==4 && fabs(pdg[mother[i]])==24) {
                          typeC = true;
                }
            }
            return typeC;
      };
      auto ttbarcharmbottom(UInt_t nPart, Vint status, Vint pdg, Vint mother) {
            bool typeC = false;
            for (unsigned int i=0; i<nPart; ++i) {
                for (unsigned int j=0; j<nPart; ++j) {
                          //if (fabs(pdg[i])==4 && fabs(pdg[mother[i]])==24 && fabs(pdg[j])==5 && fabs(pdg[mother[j]])==24) {
                          if (fabs(pdg[j])==5 && fabs(pdg[mother[j]])==24) {
                                      typeC = true;
                          }
                }
            }
            return typeC;
      };
""")

#######################    main selection    #########################

###################################################
################   DEFINITIONS   ##################
###################################################

### LEPTONS

## Funciones para seleccionar muones y electrones
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vbool tID, Float_t cutpt) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
		if (pt[i]>cutpt && fabs(eta[i])<2.4 && iso[i]<0.15 && tID[i]){
                	vb.push_back(i);
		}
            }
            return vb;
      };
      auto elInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vint cutB, Vbool mva80 , Vbool mva90, Float_t cutpt) {
            vector<int> vb;
            bool cond_eta = false;
            for (unsigned int i=0; i<nmu; ++i) {
                cond_eta = !(fabs(eta[i])>1.442 && fabs(eta[i])<1.556);
                //if (pt[i]>cutpt && fabs(eta[i])<2.5 && iso[i]<0.15 && mva80[i]){
                if (pt[i]>cutpt && fabs(eta[i])<2.4 && mva80[i] && cond_eta){
                	vb.push_back(i);
                }
            }
            return vb;
      };
""")

## Funciones para seleccionar muones y electrones secundarios
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonIndSec(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vbool tID, Float_t cutpt, Vint mu_good) {
            vector<int> vb;
            bool cond_lep = true;
            for (unsigned int i=0; i<nmu; ++i) {
                if (mu_good.size()>0) cond_lep = i != mu_good[0];
                if (pt[i]>15. && fabs(eta[i])<2.4 && iso[i]<0.2 && tID[i] && cond_lep){
                        vb.push_back(i);
                }
            }
            return vb;
      };
      auto elIndSec(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vint cutB, Vbool mva80 , Vbool mva90, Float_t cutpt, Vint cutbased, Vint el_good) {
            vector<int> vb;
            bool cond_eta = false;
            bool cond_lep = true;
            for (unsigned int i=0; i<nmu; ++i) {
                if (el_good.size()>0) cond_lep = i != el_good[0];
                cond_eta = !(fabs(eta[i])>1.442 && fabs(eta[i])<1.556);
                //if (pt[i]>cutpt && fabs(eta[i])<2.5 && iso[i]<0.15 && mva80[i]){
                if (pt[i]>15. && fabs(eta[i])<2.4 && mva80[i] && cond_eta){
                        vb.push_back(i);
                }
            }
            return vb;
      };
""")

## Numero de SV dentro de un jet
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muoninjet(UInt_t nmu, Vint mu_jetid, Vint mu_good) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
                if (mu_good.size()>0){
                	if (i!=mu_good[0] && mu_jetid[i]!=-1){
                                vb.push_back(i);
                	}
                } else {
                        if (mu_jetid[i]!=-1){
                                vb.push_back(i);
                        }
                }

            }
            return vb;
      };
""")

#######   JETS   #######

## Funciones para seleccionar JETS
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetInd(UInt_t njet, Vfloat pt, Vfloat eta, Vfloat phi, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vint el_good, Vfloat el_eta, Vfloat el_phi, Vint puID, Vint jetID) {
            vector<int> vb;
	    bool cond = false;
	    bool cond1 = false; 
            for (unsigned int i=0; i<njet; ++i) {
		if(mu_good.size()>0){
			pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
			cond = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],eta[i],mu_phi[mu_good[0]],phi[i]) > 0.4;
                        if (pt[i]>25. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                                vb.push_back(i);
                        }
		}
                if(el_good.size()>0){
			pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
                        cond = ROOT::VecOps::DeltaR(el_eta[el_good[0]],eta[i],el_phi[el_good[0]],phi[i]) > 0.4;
                        if (pt[i]>25. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                                vb.push_back(i);
                        }
                }
            }
            return vb;
      };

      auto InvariantM(const float pt, const float eta, const float phi, const float mass, const float pt1, const float eta1, const float phi1, const float mass1) {
            auto x = pt*std::cos(phi);
	    auto x1 = pt1*std::cos(phi1);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto z = pt*std::sinh(eta);
	    auto z1 = pt1*std::sinh(eta1);
            auto e = std::sqrt(x*x+y*y+z*z+mass*mass);
	    auto e1 = std::sqrt(x1*x1+y1*y1+z1*z1+mass1*mass1);

            auto mJet = std::sqrt((e+e1)*(e+e1)-(x+x1)*(x+x1)-(y+y1)*(y+y1)-(z+z1)*(z+z1));
            return mJet;
      };
      auto InvariantM3(const float pt, const float eta, const float phi, const float mass, const float pt1, const float eta1, const float phi1, const float mass1, const float pt2, const float eta2, const float phi2, const float mass2) {
            auto x = pt*std::cos(phi);
            auto x1 = pt1*std::cos(phi1);
            auto x2 = pt2*std::cos(phi2);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto y2 = pt2*std::sin(phi2);
            auto z = pt*std::sinh(eta);
            auto z1 = pt1*std::sinh(eta1);
            auto z2 = pt2*std::sinh(eta2);
            auto e = std::sqrt(x*x+y*y+z*z+mass*mass);
            auto e1 = std::sqrt(x1*x1+y1*y1+z1*z1+mass1*mass1);
            auto e2 = std::sqrt(x2*x2+y2*y2+z2*z2+mass2*mass2);

            auto mJet = std::sqrt((e+e1+e2)*(e+e1+e2)-(x+x1+x2)*(x+x1+x2)-(y+y1+y2)*(y+y1+y2)-(z+z1+z2)*(z+z1+z2));
            return mJet;
      };
""")

######### bottom jets ordering

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      template <typename T>
      vector<size_t> sort_indexes(const vector<T> &v) {

          // initialize original index locations
          vector<size_t> idx(v.size());
          iota(idx.begin(), idx.end(), 0);

          // sort indexes based on comparing values in v
          // using std::stable_sort instead of std::sort
          // to avoid unnecessary index re-orderings
          // when v contains elements of equal values
          stable_sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

          return idx;
      };
      auto bottomjets(UInt_t njetgood, Vint jetgood, UInt_t njet, Vfloat jet_btag) {
          vector<float> vb;
          vector<int> fin;
          if (njetgood>0) {
                for (unsigned int i=0; i<njetgood; ++i){
                    vb.push_back(jet_btag[jetgood[i]]);
                }
          }
          for (auto i: sort_indexes(vb)){
                fin.push_back(i);
          }
          return fin;      
      };
""")

## Muon dentro del jet

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetMuonIndJet(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_pt, Vfloat mu_iso, Vint el_good, Vint mu_jetid, Vint jetbotind) {
            vector<int> vb;
            bool cond = false;
            bool cond1 = true;
            bool condb = true;
            int ind=-1;
            float ptM{-10.};
            float ptJ{-10.};
            if (good.size() == 1){
                for (unsigned int i=0; i<nmu; ++i){
                        cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[0]],mu_phi[i],phi[good[0]]) < 0.4;
                        if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                        if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0] && mu_iso[i] > 0.2){
                        //if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0]){
                                ind = good[0];
                                ptM = mu_pt[i];
                        }
                }
            }
            if (good.size() > 1){
                for (unsigned int j=0; j<good.size(); ++j){
                        for (unsigned int i=0; i<nmu; ++i){
                                cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[j]],mu_phi[i],phi[good[j]]) < 0.4;
                                if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                                if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                                if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && mu_iso[i] > 0.2 && condb){
                                //if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && condb){
                                        ind = good[j];
                                        ptM = mu_pt[i];
                                        ptJ = pt[good[j]];
                                }
                        }
                }
            }
            if (ind>-1) {
		vb.push_back(ind);
	    }
            return vb;
      };
      auto JetMuonInd(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_pt, Vfloat mu_iso, Vint el_good, Vint mu_jetid, Vint jetbotind) {
            vector<int> vb;
            bool cond = false;
	    bool cond1 = true;
            bool condb = true;
            int ind=-1;
            float ptM{-10.};
            float ptJ{-10.};
	    if (good.size() == 1){
		for (unsigned int i=0; i<nmu; ++i){
			cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[0]],mu_phi[i],phi[good[0]]) < 0.4;
			if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                        if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0] && mu_iso[i] > 0.2){
			//if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0]){
                                ind = i;
				ptM = mu_pt[i];
			}
		}
	    }
	    if (good.size() > 1){
		for (unsigned int j=0; j<good.size(); ++j){
                	for (unsigned int i=0; i<nmu; ++i){
                        	cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[j]],mu_phi[i],phi[good[j]]) < 0.4;
				if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                                if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                                if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && mu_iso[i] > 0.2 && condb){
                        	//if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && condb){
                                        ind = i;
                                	ptM = mu_pt[i];
                                        ptJ = pt[good[j]];
                        	}
                	}
		}
	    }
            if (ind>-1) {
		vb.push_back(ind);
            }
	    return vb;
      };
""")

######### pedimos algunas condiciones al muon seleccionado

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonCond( Vint mu_jet,Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_iso, Vbool mu_softid) {
	    bool cond = false;
	    if (mu_jet.size()>0){ 
	    	if (mu_pt[mu_jet[0]]<25. && fabs(mu_eta[mu_jet[0]])<2.4 && mu_softid[mu_jet[0]]) {
			cond = true;
	    	}
	    }
            return cond;
      };
""")

## Secondary vertex

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetSVIndJet(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nSV, Vfloat sv_pt, Vfloat sv_phi, Vfloat sv_eta, Vint jetbotind) {
            vector<int> vb;
            bool cond = false;
            bool condb = true;
            float ptM{-10.};
            float ptJ{-10.};
            int ind = -1;
            for (unsigned int j=0; j<good.size(); ++j){
                        for (unsigned int i=0; i<nSV; ++i){
                                cond = ROOT::VecOps::DeltaR(sv_eta[i],eta[good[j]],sv_phi[i],phi[good[j]]) < 0.4;
                                if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                                if(cond && pt[good[j]] > ptJ && condb){
                                        ind = good[j];
                                        ptJ = pt[good[j]];
                                }
                        }
            }
            if (ind>-1) {
		vb.push_back(ind);
            }
            return vb;
      };
      auto JetSVIndSV(UInt_t njet, Vint good, Vint jet_sv, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nSV, Vfloat sv_pt, Vfloat sv_phi, Vfloat sv_eta) {
            vector<int> vb;
            bool cond = false;
            bool condb = true;
            float ptM{-10.};
            float ptJ{-10.};
            int ind = -1;
            for (unsigned int i=0; i<nSV; ++i){
                  if (jet_sv.size() > 0) cond = ROOT::VecOps::DeltaR(sv_eta[i],eta[jet_sv[0]],sv_phi[i],phi[jet_sv[0]]) < 0.4;
                  if(cond && sv_pt[i] > ptM){
                               ind = i;
                               ptM = sv_pt[i];
                  }
            }
            if (ind>-1) {
                vb.push_back(ind);
            }
            return vb;
      };
""")

######### pedimos algunas condiciones al SV seleccionado

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto svCond( Vint sv_jet,Vfloat sv_pt, Vfloat sv_eta) {
	    bool cond = false;
	    if (sv_jet.size()>0){ 
			cond = true;
	    }
            return cond;
      };
""")

######### Segundo jet

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto secondjet( UInt_t njetgood,Vint good, Vint jetsv, Vint jetbotind, Vfloat jet_pt) {
            int ind = -1;
            float pT = -10.;
            for (unsigned int i=0; i<njetgood; ++i){
              if (jet_pt[good[i]] > pT && good[i]!=jetsv[0] && good[i]!=good[jetbotind[0]] && good[i]!=good[jetbotind[1]]) {
                    ind = good[i];
                    pT = jet_pt[good[i]];
              }
            }
            return ind;
      };
""")

######### Masa invariante con los bottom

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto InvMassBot(Vint jetgood, Vint jetbot, Vint jetmuon, Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi, UInt_t jet_notmuon, Vfloat Jet_mass) {
            int ind = -1;
            vector<float> vb;
            float dR1 = ROOT::VecOps::DeltaR(Jet_eta[jetmuon[0]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetmuon[0]],Jet_phi[jetgood[jetbot[0]]]);
            float dR2 = ROOT::VecOps::DeltaR(Jet_eta[jetmuon[0]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetmuon[0]],Jet_phi[jetgood[jetbot[1]]]);
            if (dR1 > dR2){
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_mass[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_mass[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
            }else{
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_mass[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_mass[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_mass[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon],Jet_mass[jet_notmuon]));
            }
            return vb;
      };
""")

####   SSOS   ####

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto SSOS(Vint mu_good, Vint el_good, Vint sv_jet,Vint el_charge, Vint mu_charge, Vint sv_charge) {
            int ind = 0;
            if(sv_jet.size()>0){
                if(mu_good.size()>0){
                        ind = mu_charge[mu_good[0]]*sv_charge[sv_jet[0]];
                }
                if(el_good.size()>0){
                        ind = el_charge[el_good[0]]*sv_charge[sv_jet[0]];
                }
            }
            return ind;
      };
""")

## pT component calculations

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto pTsum(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            auto ptot = plep1+plep2+phad1+phad2;
            auto plepX = plep.Px(); auto plepY = plep.Py(); auto phadX = phad.Px(); auto phadY = phad.Py(); auto ptotX = ptot.Px(); auto ptotY = ptot.Py();
            float plepmod; float phadmod; float ptotmod;
            plepmod = std::sqrt(plepX*plepX + plepY*plepY);
            phadmod = std::sqrt(phadX*phadX + phadY*phadY);
            ptotmod = std::sqrt(ptotX*ptotX + ptotY*ptotY);
            return ptotmod/(plepmod+phadmod);
      };
      auto pTprod(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            auto plepX = plep.Px(); auto plepY = plep.Py(); auto phadX = phad.Px(); auto phadY = phad.Py();
            float plepmod; float phadmod; float ptotmod;
            plepmod = std::sqrt(plepX*plepX + plepY*plepY);
            phadmod = std::sqrt(phadX*phadX + phadY*phadY);
            ptotmod = plepX*phadX + plepY*phadY;
            return ptotmod/(plepmod*phadmod);
      };
      auto variousSUM(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            vector<float> vb;
            ROOT::Math::PtEtaPhiMVector plep1;
            float deltaRlepjet1;
            float deltaRlepjet2;
            float deltaPhilepjet1;
            float deltaPhilepjet2;
            float deltaEtalepjet1;
            float deltaEtalepjet2;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
                  deltaRlepjet1 = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_muon[0]],mu_phi[mu_good[0]],jet_phi[jet_muon[0]]);
                  deltaRlepjet2 = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_notmuon],mu_phi[mu_good[0]],jet_phi[jet_notmuon]);
                  deltaPhilepjet1 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_notmuon]);
                  deltaEtalepjet1 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_notmuon]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
                  deltaRlepjet1 = ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_muon[0]],el_phi[el_good[0]],jet_phi[jet_muon[0]]);
                  deltaRlepjet2 = ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_notmuon],el_phi[el_good[0]],jet_phi[jet_notmuon]);
                  deltaPhilepjet1 = fabs(el_phi[el_good[0]]-jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = fabs(el_phi[el_good[0]]-jet_phi[jet_notmuon]);
                  deltaEtalepjet1 = fabs(el_eta[el_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(el_eta[el_good[0]]-jet_eta[jet_notmuon]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            vb.push_back(ROOT::VecOps::DeltaR(plep1.Eta(),phad.Eta(),plep1.Phi(),phad.Phi()));
            vb.push_back(fabs(plep2.Phi()-phad.Phi()));
            vb.push_back(phad.Eta());
            vb.push_back(phad.Pt());
            vb.push_back(fabs(plep.Phi()-phad.Phi()));
            vb.push_back(ROOT::VecOps::DeltaR(plep.Eta(),phad.Eta(),plep.Phi(),phad.Phi()));
            vb.push_back(fabs(plep1.Phi()-phad.Phi()));
            vb.push_back(fabs(plep.Eta()-phad.Eta()));
            vb.push_back(fabs(plep1.Eta()-phad.Eta()));
            vb.push_back(fabs(plep.Pt()-phad.Pt()));
            vb.push_back(fabs(plep1.Pt()-phad.Pt()));
            vb.push_back(phad2.Pt()*std::sin(phad1.Phi()-phad2.Phi()));
            vb.push_back(phad.Pt()/(phad1.Pt()+phad2.Pt()));
            vb.push_back(plep.Pt());
            vb.push_back(deltaRlepjet1);
            vb.push_back(deltaRlepjet2);
            vb.push_back(deltaPhilepjet1);
            vb.push_back(deltaPhilepjet2);
            vb.push_back(deltaEtalepjet1);
            vb.push_back(deltaEtalepjet2);
            return vb;
      };
""")

## Trigger function for 2017 electron data

gInterpreter.Declare("""
      #include <iomanip>
      #include <math.h>
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto calDeltaR(float eta1, float phi1, float eta2, float phi2) {
          float dphi = phi1-phi2;
          float deta = eta1-eta2;
          if (dphi <= - M_PI) {
             dphi = dphi + M_PI;
          }
          if (dphi >= M_PI) {
             dphi = dphi - M_PI;
          }
          float prod = deta*deta + dphi*dphi;
          return prod;
      };
      auto triggeremulator(UInt_t nel, UInt_t ntrig, Vfloat el_eta, Vfloat el_deltaeta, Vfloat el_phi, Vfloat trig_eta, Vfloat trig_phi, Vint trig_bits, Vint trig_id) {
          bool cond = false;
          bool cond1 = false;
          bool cond2 = false;
          for (unsigned int j=0; j<nel; ++j) {
            for (unsigned int i=0; i<ntrig; ++i) {
                cond1 = calDeltaR(el_eta[j]+el_deltaeta[j],el_phi[j], trig_eta[i], trig_phi[i]) < 0.01;
                cond2 = trig_bits[i] & (0x1 << 10);
                if( trig_id[i]==11 && cond1 && cond2) {
                    cond = true;
                }
            }
          }
          return cond;
      };
""")

## Funciones para seleccionar muones y electrones
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonSecInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vbool tID, Float_t cutpt) {
            vector<int> vb;
            float ptest = -10.;
            int ind= -1;
            for (unsigned int i=0; i<nmu; ++i) {
                if (fabs(eta[i])<2.4 && iso[i]<0.15 && tID[i] && pt[i] < cutpt && pt[i]>ptest){
                        ptest = pt[i];
                        ind = i;
                }
            }
            if (ind > -1) vb.push_back(ind);
            return vb;
      };
      auto elSecInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vint cutB, Vbool mva80 , Vbool mva90, Float_t cutpt) {
            vector<int> vb;
            float ptest = -10.;
            int ind= -1;
            for (unsigned int i=0; i<nmu; ++i) {
                //if (pt[i]>cutpt && fabs(eta[i])<2.5 && iso[i]<0.15 && mva80[i]){
                if (fabs(eta[i])<2.5 && mva80[i] && pt[i] < cutpt && pt[i]>ptest){
                        ptest = pt[i];
                        ind = i;
                }
            }
            if (ind > -1) vb.push_back(ind);
            return vb;
      };
""")

#############################################################
######################## Definitions ########################
#############################################################

df_muon = {}
df_electron = {}
df_test = {}

for s in samples:
      if s in special_weights:
             df[s] = df[s].Define('weight_aux','fabs(genWeight) > 0 ? genWeight/fabs(genWeight) : 0')
      else:
             df[s] = df[s].Define('weight_aux','1')

for s in samples:
	df[s] = df[s].Define('MuonGoodInd','muonInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all, Muon_tightId,'+str(muon_pt[years[0]])+')')
	df[s] = df[s].Define('ElectronGoodInd','elInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_cutBased, Electron_mvaFall17V2Iso_WP80, Electron_mvaFall17V2Iso_WP90,'+str(el_pt[years[0]])+')')
	df[s] = df[s].Define('nMuonGood','MuonGoodInd.size()')
	df[s] = df[s].Define('nElectronGood','ElectronGoodInd.size()')
	df[s] = df[s].Define('JetGoodInd','JetInd(nJet, Jet_pt, Jet_eta, Jet_phi, MuonGoodInd, Muon_eta, Muon_phi, ElectronGoodInd, Electron_eta, Electron_phi, Jet_puId, Jet_jetId)')
	df[s] = df[s].Define('nJetGood','JetGoodInd.size()')
	df[s] = df[s].Define('JetBotInd','bottomjets(nJetGood, JetGoodInd, nJet, Jet_btagDeepFlavB)')
	df[s] = df[s].Define('JetSVInd','JetSVIndJet(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nSV, SV_pt, SV_phi, SV_eta, JetBotInd)')
	df[s] = df[s].Define('SVJetInd','JetSVIndSV(nJet, JetGoodInd, JetSVInd, Jet_pt, Jet_eta, Jet_phi, nSV, SV_pt, SV_phi, SV_eta)')
	df[s] = df[s].Define('SVJetGood','svCond(SVJetInd, SV_pt, SV_eta)')
	df[s] = df[s].Define('MuonJetInd','JetMuonInd(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx, JetBotInd)')
	df[s] = df[s].Define('JetMuonInd','JetMuonIndJet(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx, JetBotInd)')
	df[s] = df[s].Define('MuonJetGood','muonCond( MuonJetInd, Muon_pt, Muon_eta, Muon_pfRelIso04_all,Muon_softId)')
	df[s] = df[s].Define('MuonGoodSec','muonIndSec(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all, Muon_mediumId,'+str(muon_pt[years[0]])+', MuonGoodInd)')
	df[s] = df[s].Define('ElectronGoodSec','elIndSec(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_cutBased, Electron_mvaFall17V2Iso_WP80, Electron_mvaFall17V2Iso_WP90,'+str(el_pt[years[0]])+', Electron_cutBased, ElectronGoodInd)')
	if args.ssos: df[s] = df[s].Define('weightSSOS','weight_aux*(-1)*SVLepSign/std::abs(SVLepSign)')
	else: df[s] = df[s].Define('weightSSOS','weight_aux')
	df_test[s] = df[s]
	if s[0:3] == "hig":
		df[s] = df[s].Define('SVLepSign','-1')
	else:
		df[s] = df[s].Define('SVLepSign','SSOS(MuonGoodInd, ElectronGoodInd, SVJetInd, Electron_charge, Muon_charge, SV_charge)')



#################     FILTERS     #######################

##### New cuts, compared to verison 0 and 1
##### Exactly 2 jets, not more or less
##### ETA restrictions: 2.4 for jets and muons and 2.5 for electrons

###########    Extra definitions

for s in samples:
	df[s] = df[s].Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood >= 4').Filter('SVJetInd.size() >= 1').Filter('SVJetGood')
	df[s] = df[s].Filter('!(MuonJetInd.size() >= 1 && MuonJetGood)')
	#df[s] = df[s].Filter('MuonGoodSec.size() < 1 && ElectronGoodSec.size() < 1')
	df[s] = df[s].Define('JetnotSVInd','secondjet(nJetGood, JetGoodInd, JetSVInd, JetBotInd, Jet_pt)')
	df[s] = df[s].Define('SVmaxpTInd','std::distance(SV_pt.begin(),std::max_element(SV_pt.begin(),SV_pt.end()))')
	### hists definitions
	df[s] = df[s].Define('jet_muon_pt','Jet_pt[JetSVInd[0]]')
	df[s] = df[s].Define('jet_muon_nmu','Jet_nMuons[JetSVInd[0]]')
	df[s] = df[s].Define('jet_muon_mass','Jet_mass[JetSVInd[0]]')
	df[s] = df[s].Define('jet_notmuon_pt','nJetGood>1 ? Jet_pt[JetnotSVInd] : 0')
	df[s] = df[s].Define('jet_muon_eta','Jet_eta[JetSVInd[0]]')
	df[s] = df[s].Define('jet_notmuon_eta','nJetGood>1 ? Jet_eta[JetnotSVInd] : 0')
	df[s] = df[s].Define('jet_notmuon_mass','nJetGood>1 ? Jet_mass[JetnotSVInd] : 0')
	df[s] = df[s].Define('jet_notmuon_qgl','nJetGood>1 ? Jet_qgl[JetnotSVInd] : 0')
	df[s] = df[s].Define('jet_notmuon_nmu','nJetGood>1 ? Jet_nMuons[JetnotSVInd] : 0')
	df[s] = df[s].Define('sv_jet_pt','SV_pt[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_eta','SV_eta[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_mass','SV_mass[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_pangle','SV_pAngle[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_chi','SV_chi2[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_ndof','SV_ndof[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_relpt','SV_pt[SVJetInd[0]]/Jet_pt[JetSVInd[0]]')
	df[s] = df[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt[JetSVInd[0]],Jet_eta[JetSVInd[0]],Jet_phi[JetSVInd[0]],0.,Jet_pt[JetnotSVInd],Jet_eta[JetnotSVInd],Jet_phi[JetnotSVInd],0.) : 0')
	df[s] = df[s].Define('deltaR_jetM_jetNM','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotSVInd], Jet_eta[JetSVInd[0]] , Jet_phi[JetnotSVInd], Jet_phi[JetSVInd[0]])  : 10')
	df[s] = df[s].Define('deltaphi_jetM_jetNM','fabs(Jet_phi[JetSVInd[0]]-Jet_phi[JetnotSVInd])')
	df[s] = df[s].Define('deltaeta_jetM_jetNM','fabs(Jet_eta[JetSVInd[0]]-Jet_eta[JetnotSVInd])')
	df[s] = df[s].Define('deltapt_jetM_jetNM','fabs(Jet_pt[JetSVInd[0]]-Jet_pt[JetnotSVInd])')
	df[s] = df[s].Define('tracks_jetM','Jet_nConstituents[JetSVInd[0]]')
	df[s] = df[s].Define('tracks_jetNM','nJetGood>1 ? Jet_nConstituents[JetnotSVInd] : 0')
	df[s] = df[s].Define('EMN_jetM','Jet_neEmEF[JetSVInd[0]]')
	df[s] = df[s].Define('EMC_jetM','Jet_chEmEF[JetSVInd[0]]')
	df[s] = df[s].Define('EMtotal_jetM','Jet_chEmEF[JetSVInd[0]]+Jet_neEmEF[JetSVInd[0]]')
	df[s] = df[s].Define('MuoninJetAux','muoninjet(nMuon, Muon_jetIdx, MuonGoodInd)')
	df[s] = df[s].Define('nMuoninJet','MuoninJetAux.size()')
	df[s] = df[s].Define('sv_jet_xy','SV_dxy[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_sigxy','SV_dxySig[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_r','SV_dlen[SVJetInd[0]]')
	df[s] = df[s].Define('sv_jet_sigr','SV_dlenSig[SVJetInd[0]]')
	df[s] = df[s].Define('pT_sum','pTsum(MuonGoodInd, ElectronGoodInd, JetSVInd, JetnotSVInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
	df[s] = df[s].Define('pT_product','pTprod(MuonGoodInd, ElectronGoodInd, JetSVInd, JetnotSVInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
	df[s] = df[s].Define('aux_various','variousSUM(MuonGoodInd, ElectronGoodInd, JetSVInd, JetnotSVInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
	df[s] = df[s].Define('deltaR_lep_2jets','aux_various[0]')
	df[s] = df[s].Define('deltaphi_MET_2jets','aux_various[1]')
	df[s] = df[s].Define('eta_2jets','aux_various[2]')
	df[s] = df[s].Define('pt_2jets','aux_various[3]')
	df[s] = df[s].Define('deltaphi_lephad','aux_various[4]')
	df[s] = df[s].Define('deltaR_lephad','aux_various[5]')
	df[s] = df[s].Define('deltaphi_lep_2jets','aux_various[6]')
	df[s] = df[s].Define('deltaeta_lephad','aux_various[7]')
	df[s] = df[s].Define('deltaeta_lep_2jets','aux_various[8]')
	df[s] = df[s].Define('deltapt_lephad','aux_various[9]')
	df[s] = df[s].Define('deltapt_lep_2jets','aux_various[10]')
	df[s] = df[s].Define('pT_proy','aux_various[11]')
	df[s] = df[s].Define('pT_sum_2J','aux_various[12]')
	df[s] = df[s].Define('pT_Wlep','aux_various[13]')
	df[s] = df[s].Define('deltaR_lep_jet1','aux_various[14]')
	df[s] = df[s].Define('deltaR_lep_jet2','aux_various[15]')
	df[s] = df[s].Define('deltaPhi_lep_jet1','aux_various[16]')
	df[s] = df[s].Define('deltaPhi_lep_jet2','aux_various[17]')
	df[s] = df[s].Define('deltaEta_lep_jet1','aux_various[18]')
	df[s] = df[s].Define('deltaEta_lep_jet2','aux_various[19]')
	df[s] = df[s].Define('jet_bot1_btag','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[0]]]')
	df[s] = df[s].Define('jet_bot2_btag','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[1]]]')
	df[s] = df[s].Define('jet_muon_btag','Jet_btagDeepFlavB[JetSVInd[0]]')
	df[s] = df[s].Define('jet_notmuon_btag','Jet_btagDeepFlavB[JetnotSVInd]')
	df[s] = df[s].Define('nLooseLepton','nMuon+nElectron-1')
	df[s] = df[s].Define('InvM3_aux','InvMassBot(JetGoodInd, JetBotInd, JetSVInd, Jet_pt, Jet_eta, Jet_phi, JetnotSVInd, Jet_mass)')
	df[s] = df[s].Define('InvM_bot_closer','InvM3_aux[0]')
	df[s] = df[s].Define('InvM_bot_farther','InvM3_aux[1]')
	df[s] = df[s].Define('second_muon_aux','muonSecInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all, Muon_tightId,'+str(muon_pt[years[0]])+')')
	df[s] = df[s].Define('second_electron_aux','elSecInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_cutBased, Electron_mvaFall17V2Iso_WP80, Electron_mvaFall17V2Iso_WP90,'+str(el_pt[years[0]])+')')
	df[s] = df[s].Define('second_muon_pt','second_muon_aux.size() > 0 ? Muon_pt[second_muon_aux[0]] : -1')
	df[s] = df[s].Define('second_electron_pt','second_electron_aux.size() > 0 ? Electron_pt[second_electron_aux[0]] : -1')
	if s[0:3] == "hig":
		df[s] = df[s].Define('sv_jet_ntracks','3')
	else:
		df[s] = df[s].Define('sv_jet_ntracks','SV_ntracks[SVJetInd[0]]')


############ Gen level definitions

if mode == "mc":
	for s in samples:
		df[s] = df[s].Define('jet_muon_flavourH','Jet_hadronFlavour[JetSVInd[0]]')
		df[s] = df[s].Define('jet_notmuon_flavourH','Jet_hadronFlavour[JetnotSVInd]')
		df[s] = df[s].Define('jet_muon_flavourP','Jet_partonFlavour[JetSVInd[0]]')
		df[s] = df[s].Define('jet_notmuon_flavourP','Jet_partonFlavour[JetnotSVInd]')


## Differentiated definitions
for s in samples:
        df[s] = df[s].Filter(met_filter[years[0]])
        df_muon[s] = df[s].Filter('nMuonGood>0')
        df_electron[s] = df[s].Filter('nElectronGood >0')
        df[s] = df[s].Define('lepton_pt','nMuonGood>0 ? Muon_pt[MuonGoodInd[0]] : Electron_pt[ElectronGoodInd[0]]')
        df[s] = df[s].Define('lepton_eta','nMuonGood>0 ? Muon_eta[MuonGoodInd[0]] : Electron_eta[ElectronGoodInd[0]]')
        df[s] = df[s].Define('transverse_mass', 'nMuonGood>0 ? std::sqrt(2*Muon_pt[MuonGoodInd[0]]*MET_pt*(1-std::cos(Muon_phi[MuonGoodInd[0]]-MET_phi))): std::sqrt(2*Electron_pt[ElectronGoodInd[0]]*MET_pt*(1-std::cos(Electron_phi[ElectronGoodInd[0]]-MET_phi)))')
        df[s] = df[s].Define('lepton_iso', 'nMuonGood>0 ? Muon_pfRelIso04_all[MuonGoodInd[0]] : Electron_pfRelIso03_all[ElectronGoodInd[0]]')
        ## New cuts
        df[s] = df[s].Filter('transverse_mass > 50')
        df[s] = df[s].Filter('sv_jet_ntracks>2')
        #df[s] = df[s].Filter('sv_jet_sigr > 5')
        #df[s] = df[s].Filter('EMtotal_jetM<0.4')
        ## BWP: cutting on btagging working points
        df[s] = df[s].Filter('jet_bot1_btag >'+str(cuts_btag[years[0]][1]))
        df[s] = df[s].Filter('jet_bot2_btag >'+str(cuts_btag[years[0]][0]))

df_M_ss = {}
df_M_os = {}
df_E_ss = {}
df_E_os = {}

for s in samples:
        df_M_ss[s] = df[s].Filter('SVLepSign > 0 && nMuonGood>0')
        df_M_os[s] = df[s].Filter('SVLepSign < 0 && nMuonGood>0')
        df_E_ss[s] = df[s].Filter('SVLepSign > 0 && nElectronGood>0')
        df_E_os[s] = df[s].Filter('SVLepSign < 0 && nElectronGood>0')
        if args.type == "mc":
                df_M_ss[s] = df_M_ss[s].Filter(muon_trig[years[0]])
                df_M_os[s] = df_M_os[s].Filter(muon_trig[years[0]])
                df_E_ss[s] = df_E_ss[s].Filter(el_trig[years[0]])
                df_E_os[s] = df_E_os[s].Filter(el_trig[years[0]])
        else:
                if args.process == "M":
                        df_M_ss[s] = df_M_ss[s].Filter(muon_trig[years[0]])
                        df_M_os[s] = df_M_os[s].Filter(muon_trig[years[0]])
                        df_E_ss[s] = df_E_ss[s].Filter(muon_trig[years[0]])
                        df_E_os[s] = df_E_os[s].Filter(muon_trig[years[0]])
                if args.process == "E":
                        if years[0] == "2017":
                                df_M_ss[s] = df_M_ss[s].Filter('triggeremulator(nElectron, nTrigObj, Electron_eta, Electron_deltaEtaSC, Electron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)')
                                df_M_os[s] = df_M_os[s].Filter('triggeremulator(nElectron, nTrigObj, Electron_eta, Electron_deltaEtaSC, Electron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)')
                                df_E_ss[s] = df_E_ss[s].Filter('triggeremulator(nElectron, nTrigObj, Electron_eta, Electron_deltaEtaSC, Electron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)')
                                df_E_os[s] = df_E_os[s].Filter('triggeremulator(nElectron, nTrigObj, Electron_eta, Electron_deltaEtaSC, Electron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)')
                        else:
                                df_M_ss[s] = df_M_ss[s].Filter(el_trig[years[0]])
                                df_M_os[s] = df_M_os[s].Filter(el_trig[years[0]])
                                df_E_ss[s] = df_E_ss[s].Filter(el_trig[years[0]])
                                df_E_os[s] = df_E_os[s].Filter(el_trig[years[0]])


############################################################
####################     HISTS    ##########################
############################################################

hist_InvM_2jets_M_ss = {}
hist_InvM_2jets_M_os = {}
hist_nJetGood_M_ss = {}
hist_nJetGood_M_os = {}
hist_jet_muon_pt_M_ss = {}
hist_jet_muon_pt_M_os = {}
hist_InvM_2jets_E_ss = {}
hist_InvM_2jets_E_os = {}
hist_nJetGood_E_ss = {}
hist_nJetGood_E_os = {}
hist_jet_muon_pt_E_ss = {}
hist_jet_muon_pt_E_os = {}

for s in samples:
        hist_nJetGood_M_ss[s] = df_M_ss[s].Histo1D(("nJetGood_M_ss","",10,0,10),"nJetGood")
        hist_nJetGood_M_os[s] = df_M_os[s].Histo1D(("nJetGood_M_os","",10,0,10),"nJetGood")

        hist_jet_muon_pt_M_ss[s] = df_M_ss[s].Histo1D(("jet_muon_pt_M_ss","",50,20,120),"jet_muon_pt")
        hist_jet_muon_pt_M_os[s] = df_M_os[s].Histo1D(("jet_muon_pt_M_os","",50,20,120),"jet_muon_pt")

        hist_InvM_2jets_M_ss[s] = df_M_ss[s].Histo1D(("InvM_2jets_M_ss","",108,30,300),"InvM_2jets")
        hist_InvM_2jets_M_os[s] = df_M_os[s].Histo1D(("InvM_2jets_M_os","",108,30,300),"InvM_2jets")

        hist_nJetGood_E_ss[s] = df_E_ss[s].Histo1D(("nJetGood_E_ss","",10,0,10),"nJetGood")
        hist_nJetGood_E_os[s] = df_E_os[s].Histo1D(("nJetGood_E_os","",10,0,10),"nJetGood")

        hist_jet_muon_pt_E_ss[s] = df_E_ss[s].Histo1D(("jet_muon_pt_E_ss","",50,20,120),"jet_muon_pt")
        hist_jet_muon_pt_E_os[s] = df_E_os[s].Histo1D(("jet_muon_pt_E_os","",50,20,120),"jet_muon_pt")

        hist_InvM_2jets_E_ss[s] = df_E_ss[s].Histo1D(("InvM_2jets_E_ss","",108,30,300),"InvM_2jets")
        hist_InvM_2jets_E_os[s] = df_E_os[s].Histo1D(("InvM_2jets_E_os","",108,30,300),"InvM_2jets")

#############################
####     DATA SAVING     ####
#############################

for s in samples:
                if args.ssos: path_hist = '/nfs/cms/vazqueze/hists_ttbar/hists/higgs/ssos/histstt_SV_v1v2vBWPSSOS_'+s+'.root'
                else: path_hist = '/nfs/cms/vazqueze/hists_ttbar/hists/higgs/sv/histstt_SV_v1v2vBWP_'+s+'.root'
                myfile = TFile( path_hist, 'RECREATE' )

                hist_nJetGood_M_ss[s].Write()
                hist_nJetGood_M_os[s].Write()
                hist_jet_muon_pt_M_ss[s].Write()
                hist_jet_muon_pt_M_os[s].Write()
                hist_InvM_2jets_M_ss[s].Write()
                hist_InvM_2jets_M_os[s].Write()
                hist_nJetGood_E_ss[s].Write()
                hist_nJetGood_E_os[s].Write()
                hist_jet_muon_pt_E_ss[s].Write()
                hist_jet_muon_pt_E_os[s].Write()
                hist_InvM_2jets_E_ss[s].Write()
                hist_InvM_2jets_E_os[s].Write()

                myfile.Close()

if args.ssos: print('SSOS version')

print('Ended succesfully')
