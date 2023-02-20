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

if args.process == "allMC": proc = ["higgs_tocs_m80","higgs_tocs_m90","higgs_tocs_m100","higgs_tocs_m120","higgs_tocs_m140","higgs_tocs_m150",
           "higgs_tocs_m155","higgs_tocs_m160"] 
elif (args.process == "higgs_tocs_m80" or args.process == "higgs_tocs_m90" or args.process == "higgs_tocs_m100" or args.process == "higgs_tocs_m120"
        or args.process == "higgs_tocs_m140" or args.process == "higgs_tocs_m150" 
        or args.process == "higgs_tocs_m155" or args.process == "higgs_tocs_m160"): proc = [str(args.process)]
else: raise NameError('Incorrect process name')

if args.year == "all": years = ["2016","2016B","2017","2018"]
elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): years = [str(args.year)]
else: raise NameError('Incorrect year')

if args.type == "data": mode = "data"
elif args.type == "mc": mode = "mc"
else: raise NameError('Incorrect type')

samples = []

for p in proc:
  for y in years:
    if mode == "mc": samples.append(p+y)
    else: samples.append(y+p)

samples_test = samples[:]

print("the samples treated are",samples)

if args.notfull:
	if len(args.list) != 2: raise NameError('List has to have 2 elements')
	print("the range of files is",args.list)
 

# Create a ROOT dataframe for each dataset
# Note that we load the filenames from the external json file placed in the same folder than this script.
# Example "python analisisWW/selection_v2.py --process="WW" --notfull -l 0 50"

files = json.load(open("/nfs/cms/vazqueze/higgssearch/mcinfo.json"))

df = {}
xsecs = {}
sumws = {}
archives = {}

for s in samples:
  archives[s]=[]

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
    if (num_files == len(list_files)):
       for f in list_files:
          archives[p+years[0]].append(f)


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

muon_pt["2016"]=30
muon_pt["2016B"]=30
muon_pt["2017"]=30
muon_pt["2018"]=30

el_pt = {}

el_pt["2016"]=35
el_pt["2016B"]=35
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
                //if (pt[i]>cutpt && fabs(eta[i])<2.4 && iso[i]<0.15 && tID[i]){
                if (iso[i]<0.15 && tID[i]){
                	vb.push_back(i);
		}
            }
            return vb;
      };
      auto elInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vint cutB, Vbool mva80 , Vbool mva90, Float_t cutpt) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
                //if (pt[i]>cutpt && fabs(eta[i])<2.5 && mva80[i]){
                if (mva80[i]){
                	vb.push_back(i);
                }
            }
            return vb;
      };
""")

## Numero de muones dentro de un jet
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
      auto JetInd(UInt_t njet, Vfloat pt, Vfloat eta, Vfloat phi, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vint el_good, Vfloat el_eta, Vfloat el_phi) {
            vector<int> vb;
	    bool cond = false;
	    bool cond1 = false;
            for (unsigned int i=0; i<njet; ++i) {
		if(mu_good.size()>0){
			//pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
			cond = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],eta[i],mu_phi[mu_good[0]],phi[i]) > 0.4;
                        //if (pt[i]>30. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                        if (cond){
                                vb.push_back(i);
                        }
		}
                if(el_good.size()>0){
			//pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
                        cond = ROOT::VecOps::DeltaR(el_eta[el_good[0]],eta[i],el_phi[el_good[0]],phi[i]) > 0.4;
                        //if (pt[i]>30. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                        if (cond){
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
      auto bottomjets(UInt_t njetgood, Vint jetgood, Vfloat partonflav, Vint mu_good, Vint el_good, Vint muon_charge, Vint el_charge) {
            vector<int> vb;
            int indB1 = -1;
            int indB2 = -1;
            for (unsigned int i=0; i<jetgood.size(); ++i) {
                if (fabs(partonflav[i]) == 5){
                   if (mu_good.size()>0){
                      if (partonflav[i]*muon_charge[mu_good[0]] > 0){
                         indB1 = i;
                      } else  if (partonflav[i]*muon_charge[mu_good[0]] < 0) {
       	       	       	 indB2 = i;
                      }
                   }
       	       	   if (el_good.size()>0){
                      if (partonflav[i]*el_charge[el_good[0]] > 0){
       	       	       	 indB1 = i;
                      } else  if (partonflav[i]*el_charge[el_good[0]] < 0) {
                         indB2 = i;
                      }
       	       	   }
                }
            }
            if (indB1>-1) vb.push_back(indB1);
            if (indB2>-1) vb.push_back(indB2);
            return vb;
      };
      auto charmjets(UInt_t njetgood, Vint jetgood, Vfloat partonflav, Vint mu_good, Vint el_good, Vint muon_charge, Vint el_charge) {
            vector<int> vb;
            int indC = -1;
            for (unsigned int i=0; i<jetgood.size(); ++i) {
                if (fabs(partonflav[i]) == 4){
                   if (mu_good.size()>0){
                      if (partonflav[i]*muon_charge[mu_good[0]] < 0){
                         indC = jetgood[i];
                      }
                   }
                   if (el_good.size()>0){
                      if (partonflav[i]*el_charge[el_good[0]] < 0){
                         indC = jetgood[i];
                      }
                   }
                }
            }
            if (indC>-1) vb.push_back(indC);
            return vb;
      };
      auto strjets(UInt_t njetgood, Vint jetgood, Vfloat partonflav, Vfloat pt, Vint mu_good, Vint el_good, Vint muon_charge, Vint el_charge) {
            vector<int> vb;
            int indC = -1;
            float ptJ = 0.;
            for (unsigned int i=0; i<jetgood.size(); ++i) {
                if (fabs(partonflav[i]) == 3 && pt[jetgood[i]] > ptJ){
                         indC = jetgood[i];
                         ptJ = pt[jetgood[i]];
                }
            }
            if (indC>-1) vb.push_back(indC);
            return vb;
      };
""")

####### Cantidad de muones total

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto nMUtotal(Vint jetbot, Vint jetgood, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nMuonG, Vint muongood, Vfloat muon_phi, Vfloat muon_eta, Vint mu_jetid) {
            vector<int> vb;
            bool condb = false;
            bool cond = false;
            bool condm = false;
            int ind = -1;
            for (unsigned int i=0; i<muon_phi.size(); ++i){
                  ind = -1;
                  if (nMuonG>0) condm = i != muongood[0];
                  for (unsigned int j=0; j<jetgood.size(); ++j){
                             if (jetbot.size() > 1) condb = (jetgood[j] != jetgood[jetbot[0]] && jetgood[j] != jetgood[jetbot[1]]);
                             cond = ROOT::VecOps::DeltaR(muon_eta[i],eta[jetgood[j]],muon_phi[i],phi[jetgood[j]]) < 0.4;
                             if(cond && condb && condm && mu_jetid[i]==jetgood[j]){
                                     ind = i;
                             }

                  }
                  if(ind>-1){
                             vb.push_back(ind);
                  }

            }
            return vb;
      };
""")

######### Segundo jet

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto secondjet( UInt_t njetgood,Vint good, Vint jetmuon, Vint jetbotind, Vfloat jet_pt) {
            int ind = -1;
            float pT = -10.;
            for (unsigned int i=0; i<njetgood; ++i){
              if (jet_pt[good[i]] > pT && good[i]!=jetmuon[0] && good[i]!=good[jetbotind[0]] && good[i]!=good[jetbotind[1]]) {
                    ind = good[i];
                    pT = jet_pt[good[i]];
              }
            }
            return ind;
      };
""")


####   SSOS   ####

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto SSOS(Vint mu_good, Vint el_good, Vint mu_jet,Vint el_charge, Vint mu_charge) {
	    int ind = 0;
	    if(mu_jet.size()>0){
            	if(mu_good.size()>0){
			ind = mu_charge[mu_good[0]]*mu_charge[mu_jet[0]];
            	}
            	if(el_good.size()>0){
                	ind = el_charge[el_good[0]]*mu_charge[mu_jet[0]];
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
       	       	  deltaRlepjet2	= ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_notmuon],mu_phi[mu_good[0]],jet_phi[jet_notmuon]);
                  deltaPhilepjet1 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_muon[0]]);
       	       	  deltaPhilepjet2 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_notmuon]);
       	       	  deltaEtalepjet1 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_notmuon]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
       	       	  deltaRlepjet1	= ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_muon[0]],el_phi[el_good[0]],jet_phi[jet_muon[0]]);
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
      auto muonSecInd(UInt_t nmu, Vint mu_good, Vfloat pt, Vfloat eta, Vfloat iso, Vbool tID, Float_t cutpt) {
            vector<int> vb;
            float ptest = -10.;
            int ind= -1;
            for (unsigned int i=0; i<nmu; ++i) {
                   if (mu_good.size()>0){
                      if (fabs(eta[i])<2.4 && iso[i]<0.15 && tID[i] && i != mu_good[0] && pt[i]>ptest){
                        ptest = pt[i];
                        ind = i;
                      }
                   }
            }
            if (ind > -1) vb.push_back(ind);
            return vb;
      };
      auto elSecInd(UInt_t nmu, Vint mu_good, Vfloat pt, Vfloat eta, Vfloat iso, Vint cutB, Vbool mva80 , Vbool mva90, Float_t cutpt) {
            vector<int> vb;
            float ptest = -10.;
            int ind= -1;
            for (unsigned int i=0; i<nmu; ++i) {
                   if (mu_good.size()>0){
                      if (fabs(eta[i])<2.5 && mva80[i] && i != mu_good[0] && pt[i]>ptest){
                        ptest = pt[i];
                        ind = i;
                      }
                   }
            }
            if (ind > -1) vb.push_back(ind);
            return vb;
      };
""")

## Funciones para seleccionar jets segun parton flavour
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetFlavAux( Vint jetgood , Vint partonflav) {
            vector<int> vb;
            int indB= 0;
            int indC= 0;
            int indS= 0;
            for (unsigned int i=0; i<jetgood.size(); ++i) {
                if (fabs(partonflav[i]) == 5){
                        ++indB;
                } else if (fabs(partonflav[i]) == 4) {
                        ++indC;
                } else if (fabs(partonflav[i]) == 3) {
                        ++indS;
                }
            }
            vb.push_back(indB);
            vb.push_back(indC);
            vb.push_back(indS);
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
	df[s] = df[s].Define('MuonGoodInd','muonInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all, Muon_tightId,'+str(muon_pt[years[0]])+')')
	df[s] = df[s].Define('ElectronGoodInd','elInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_cutBased, Electron_mvaFall17V2Iso_WP80, Electron_mvaFall17V2Iso_WP90,'+str(el_pt[years[0]])+')')
	df[s] = df[s].Define('nMuonGood','MuonGoodInd.size()')
	df[s] = df[s].Define('nElectronGood','ElectronGoodInd.size()')
	df[s] = df[s].Define('GenJetGoodInd','JetInd(nGenJet, GenJet_pt, GenJet_eta, GenJet_phi, MuonGoodInd, Muon_eta, Muon_phi, ElectronGoodInd, Electron_eta, Electron_phi)')
	df[s] = df[s].Define('nGenJetGood','GenJetGoodInd.size()')
	df[s] = df[s].Define('GenJetFlav_aux','JetFlavAux(GenJetGoodInd,GenJet_partonFlavour)')
	df[s] = df[s].Define('GenJetBotInd','bottomjets(nGenJetGood, GenJetGoodInd, GenJet_partonFlavour, MuonGoodInd, ElectronGoodInd, Muon_charge, Electron_charge)')
	df[s] = df[s].Define('GenJetCharmInd','charmjets(nGenJetGood, GenJetGoodInd, GenJet_partonFlavour, MuonGoodInd, ElectronGoodInd, Muon_charge, Electron_charge)')
	df[s] = df[s].Define('GenJetStrInd','strjets(nGenJetGood, GenJetGoodInd, GenJet_partonFlavour, GenJet_pt, MuonGoodInd, ElectronGoodInd, Muon_charge, Electron_charge)')
	df[s] = df[s].Define('weightSSOS','1')
	df_test[s] = df[s]

#################     FILTERS     #######################

##### New cuts, compared to verison 0 and 1
##### Exactly 2 jets, not more or less
##### ETA restrictions: 2.4 for jets and muons and 2.5 for electrons

###########    Extra definitions

for s in samples:
	df[s] = df[s].Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)')
	df[s] = df[s].Filter('GenJetFlav_aux[0] == 2 && GenJetFlav_aux[1] == 1').Filter('GenJetBotInd.size() == 2').Filter('GenJetCharmInd.size() == 1')
	df[s] = df[s].Define('JetnotMuonInd','GenJetStrInd.size() > 0 ? GenJetStrInd[0] : secondjet(nGenJetGood, GenJetGoodInd, GenJetCharmInd, GenJetBotInd, GenJet_pt)')
	### hists definitions
	df[s] = df[s].Define('jet_muon_pt','GenJet_pt[GenJetCharmInd[0]]')
	df[s] = df[s].Define('jet_muon_mass','GenJet_mass[GenJetCharmInd[0]]')
	df[s] = df[s].Define('jet_notmuon_pt','JetnotMuonInd>-1 ? GenJet_pt[JetnotMuonInd] : 0')
	df[s] = df[s].Define('jet_muon_eta','GenJet_eta[GenJetCharmInd[0]]')
	df[s] = df[s].Define('jet_notmuon_eta','JetnotMuonInd>-1 ? GenJet_eta[JetnotMuonInd] : 0')
	df[s] = df[s].Define('jet_notmuon_mass','JetnotMuonInd>-1 ? GenJet_mass[JetnotMuonInd] : 0')
	df[s] = df[s].Define('jet_bot1_pt','GenJet_pt[GenJetGoodInd[GenJetBotInd[0]]]')
	df[s] = df[s].Define('jet_bot1_mass','GenJet_mass[GenJetGoodInd[GenJetBotInd[0]]]')
	df[s] = df[s].Define('jet_bot2_pt','GenJet_pt[GenJetGoodInd[GenJetBotInd[1]]]')
	df[s] = df[s].Define('jet_bot1_eta','GenJet_eta[GenJetGoodInd[GenJetBotInd[0]]]')
	df[s] = df[s].Define('jet_bot2_eta','GenJet_eta[GenJetGoodInd[GenJetBotInd[1]]]')
	df[s] = df[s].Define('jet_bot2_mass','GenJet_mass[GenJetGoodInd[GenJetBotInd[1]]]')
	df[s] = df[s].Define('InvM_2jets','JetnotMuonInd>-1 ? InvariantM(GenJet_pt[GenJetCharmInd[0]],GenJet_eta[GenJetCharmInd[0]],GenJet_phi[GenJetCharmInd[0]],GenJet_mass[GenJetCharmInd[0]],GenJet_pt[JetnotMuonInd],GenJet_eta[JetnotMuonInd],GenJet_phi[JetnotMuonInd],GenJet_mass[JetnotMuonInd]) : 0')
	df[s] = df[s].Define('nLooseLepton','nMuon+nElectron-1')
	df[s] = df[s].Define('second_muon_aux','muonSecInd(nMuon, MuonGoodInd, Muon_pt,Muon_eta,Muon_pfRelIso04_all, Muon_tightId,'+str(muon_pt[years[0]])+')')
	df[s] = df[s].Define('second_electron_aux','elSecInd(nElectron, ElectronGoodInd, Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_cutBased, Electron_mvaFall17V2Iso_WP80, Electron_mvaFall17V2Iso_WP90,'+str(el_pt[years[0]])+')')
	df[s] = df[s].Define('second_muon_pt','second_muon_aux.size() > 0 ? Muon_pt[second_muon_aux[0]] : -1')
	df[s] = df[s].Define('second_electron_pt','second_electron_aux.size() > 0 ? Electron_pt[second_electron_aux[0]] : -1')
	df[s] = df[s].Define('second_muon_eta','second_muon_aux.size() > 0 ? Muon_eta[second_muon_aux[0]] : -1')
	df[s] = df[s].Define('second_electron_eta','second_electron_aux.size() > 0 ? Electron_eta[second_electron_aux[0]] : -1')

#################################
#######    GEN FILTER    ########
#################################

######### WW sectioning in events of interest (semi charm) and not 

if mode == "mc":
	for s in samples:
		if (s[0]+s[1] == "WW" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
			df[s] = df[s].Define('typeWW','typeWW(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother)')
			df[s] = df[s].Define('typeC','typeC(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother)')
			df[s+"_hadronic"] = df[s].Filter('typeWW == 1')
			df[s+"_leptonic"] = df[s].Filter('typeWW == 3')
			df[s+"_semi_charm"] = df[s].Filter('typeWW == 2 && typeC == 1')
			df[s+"_semi_nocharm"] = df[s].Filter('typeWW == 2 && typeC != 1')
			## Samples correction
			samples.append(s+"_hadronic")
			samples.append(s+"_leptonic")
			samples.append(s+"_semi_charm")
			samples.append(s+"_semi_nocharm")

samples = [s for s in samples if not (s[0]+s[1] == "WW" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

########## Wjets sectioning for W plus c discrimination

if mode	== "mc":
	for s in samples:
        	if (s[0]+s[1]+s[2]+s[3]+s[4] == "Wjets" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                	df[s] = df[s].Define('ishard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,7)')
                	df[s] = df[s].Define('first_copy','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,12)')
                	df[s] = df[s].Define('isWplusc','wpluscbool(nGenPart,GenPart_statusFlags,GenPart_pdgId,ishard,first_copy)')
                	df[s+"_charm"] = df[s].Filter('isWplusc == 2')
                	df[s+"_bottom"] = df[s].Filter('isWplusc == 3')
                	df[s+"_doublecharm"] = df[s].Filter('isWplusc == 1')
                	df[s+"_light"] = df[s].Filter('isWplusc == 0')
                	## Samples correction
                	samples.append(s+"_charm")
                	samples.append(s+"_bottom")
                	samples.append(s+"_doublecharm")
                	samples.append(s+"_light")

samples = [s for s in samples if not (s[0]+s[1]+s[2]+s[3]+s[4] == "Wjets" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

########## ttbar sectioning for charm discrimination

if mode	== "mc":
	for s in samples:
        	if (s[0]+s[1]+s[2]+s[3]+s[4] == "ttbar" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                	df[s] = df[s].Define('isttbarC','ttbarcharm(nGenPart,GenPart_statusFlags,GenPart_pdgId,GenPart_genPartIdxMother)')
                	df[s] = df[s].Define('isttbarCB','ttbarcharmbottom(nGenPart,GenPart_statusFlags,GenPart_pdgId,GenPart_genPartIdxMother)')
                	df[s+"_charmbottom"] = df[s].Filter('isttbarCB')
                	df[s+"_charm"] = df[s].Filter('isttbarC && !isttbarCB')
                	df[s+"_nocharm"] = df[s].Filter('!isttbarC && !isttbarCB')
                	## Samples correction
                	samples.append(s+"_charmbottom")
                	samples.append(s+"_charm")
                	samples.append(s+"_nocharm")

samples = [s for s in samples if not (s[0]+s[1]+s[2]+s[3]+s[4] == "ttbar"  and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

########## ST sectioning for charm discrimination

if mode == "mc":
        for s in samples:
                if (s[0]+s[1] == "ST" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                        df[s] = df[s].Define('isSTC','ttbarcharm(nGenPart,GenPart_statusFlags,GenPart_pdgId,GenPart_genPartIdxMother)')
                        df[s+"_charm"] = df[s].Filter('isSTC')
                        df[s+"_nocharm"] = df[s].Filter('!isSTC')
                        ## Samples correction
                        samples.append(s+"_charm")
                        samples.append(s+"_nocharm")

samples = [s for s in samples if not (s[0]+s[1] == "ST"  and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

print(samples)

## Differentiated definitions
for s in samples:
	df[s] = df[s].Filter(met_filter[years[0]])
	df_muon[s] = df[s].Filter('nMuonGood>0')
	df_electron[s] = df[s].Filter('nElectronGood >0')
	if args.type == "mc":
		df_muon[s] = df_muon[s].Filter(muon_trig[years[0]])
		df_electron[s] = df_electron[s].Filter(el_trig[years[0]])
	else:
		if args.process == "M": 
			df_muon[s] = df_muon[s].Filter(muon_trig[years[0]])
			df_electron[s] = df_electron[s].Filter(muon_trig[years[0]])
		if args.process == "E":
			if years[0] == "2017":
				df_muon[s] = df_muon[s].Filter('triggeremulator(nElectron, nTrigObj, Electron_eta, Electron_deltaEtaSC, Electron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)')
				df_electron[s] = df_electron[s].Filter('triggeremulator(nElectron, nTrigObj, Electron_eta, Electron_deltaEtaSC, Electron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)')
			else:
				df_muon[s] = df_muon[s].Filter(el_trig[years[0]])
				df_electron[s] = df_electron[s].Filter(el_trig[years[0]])
	df_muon[s] = df_muon[s].Define('lepton_pt','Muon_pt[MuonGoodInd[0]]')
	df_muon[s] = df_muon[s].Define('lepton_eta','Muon_eta[MuonGoodInd[0]]')
	df_electron[s] = df_electron[s].Define('lepton_pt','Electron_pt[ElectronGoodInd[0]]')
	df_electron[s] = df_electron[s].Define('lepton_eta','Electron_eta[ElectronGoodInd[0]]')
	df_muon[s] = df_muon[s].Define('transverse_mass','std::sqrt(2*Muon_pt[MuonGoodInd[0]]*MET_pt*(1-std::cos(Muon_phi[MuonGoodInd[0]]-MET_phi)))')
	df_electron[s] = df_electron[s].Define('transverse_mass','std::sqrt(2*Electron_pt[ElectronGoodInd[0]]*MET_pt*(1-std::cos(Electron_phi[ElectronGoodInd[0]]-MET_phi)))')

############################################################
####################     HISTS    ##########################
############################################################

hist_nJetGood_M = {}
hist_nJetGood_E = {}
hist_nLooseLepton_M = {}
hist_nLooseLepton_E = {}
hist_jet_muon_pt_M = {}
hist_jet_muon_mass_M = {}
hist_jet_notmuon_pt_M = {}
hist_jet_muon_eta_M = {}
hist_jet_notmuon_eta_M = {}
hist_jet_notmuon_mass_M = {}
hist_jet_muon_pt_E = {}
hist_jet_muon_mass_E = {}
hist_jet_notmuon_pt_E = {}
hist_jet_muon_eta_E = {}
hist_jet_notmuon_eta_E = {}
hist_jet_notmuon_mass_E = {}
hist_jet_bot1_pt_M = {}
hist_jet_bot1_mass_M = {}
hist_jet_bot2_pt_M = {}
hist_jet_bot1_eta_M = {}
hist_jet_bot2_eta_M = {}
hist_jet_bot2_mass_M = {}
hist_jet_bot1_pt_E = {}
hist_jet_bot1_mass_E = {}
hist_jet_bot2_pt_E = {}
hist_jet_bot1_eta_E = {}
hist_jet_bot2_eta_E = {}
hist_jet_bot2_mass_E = {}
hist_lepton_pt_M = {}
hist_lepton_eta_M = {}
hist_lepton_pt_E = {}
hist_lepton_eta_E = {}
hist_second_muon_pt_M = {}
hist_second_muon_pt_E = {}
hist_second_electron_pt_M = {}
hist_second_electron_pt_E = {}
hist_second_muon_eta_M = {}
hist_second_muon_eta_E = {}
hist_second_electron_eta_M = {}
hist_second_electron_eta_E = {}
hist_InvM_2jets_M = {}
hist_InvM_2jets_E = {}
hist_transverse_mass_M = {}
hist_transverse_mass_E = {}

for s in samples:
        hist_nJetGood_M[s] = df_muon[s].Histo1D(("nGenJetGood_M","",10,0,10),"nGenJetGood","weightSSOS")
        hist_nJetGood_E[s] = df_electron[s].Histo1D(("nGenJetGood_E","",10,0,10),"nGenJetGood","weightSSOS")

        hist_nLooseLepton_M[s] = df_muon[s].Histo1D(("nLooseLepton_M","",10,0,10),"nLooseLepton","weightSSOS")
        hist_nLooseLepton_E[s] = df_electron[s].Histo1D(("nLooseLepton_E","",10,0,10),"nLooseLepton","weightSSOS")

        hist_second_muon_pt_M[s] = df_muon[s].Histo1D(("second_muon_pt_M","",101,-1,100),"second_muon_pt","weightSSOS")
        hist_second_muon_pt_E[s] = df_electron[s].Histo1D(("second_muon_pt_E","",101,-1,100),"second_muon_pt","weightSSOS")
        hist_second_electron_pt_M[s] = df_muon[s].Histo1D(("second_electron_pt_M","",101,-1,100),"second_electron_pt","weightSSOS")
        hist_second_electron_pt_E[s] = df_electron[s].Histo1D(("second_electron_pt_E","",101,-1,100),"second_electron_pt","weightSSOS")
        hist_second_muon_eta_M[s] = df_muon[s].Histo1D(("second_muon_eta_M","",40,-4,4),"second_muon_eta","weightSSOS")
        hist_second_muon_eta_E[s] = df_electron[s].Histo1D(("second_muon_eta_E","",40,-4,4),"second_muon_eta","weightSSOS")
        hist_second_electron_eta_M[s] = df_muon[s].Histo1D(("second_electron_eta_M","",40,-4,4),"second_electron_eta","weightSSOS")
        hist_second_electron_eta_E[s] = df_electron[s].Histo1D(("second_electron_eta_E","",40,-4,4),"second_electron_eta","weightSSOS")

        hist_jet_muon_pt_M[s] = df_muon[s].Histo1D(("jet_muon_pt_M","",60,0,120),"jet_muon_pt","weightSSOS")
        hist_jet_muon_mass_M[s] = df_muon[s].Histo1D(("jet_muon_mass_M","",40,0,40),"jet_muon_mass","weightSSOS")
        hist_jet_notmuon_pt_M[s] = df_muon[s].Histo1D(("jet_not_muon_pt_M","",60,0,120),"jet_notmuon_pt","weightSSOS")
        hist_jet_muon_eta_M[s] = df_muon[s].Histo1D(("jet_muon_eta_M","",80,-4,4),"jet_muon_eta","weightSSOS")
        hist_jet_notmuon_eta_M[s] = df_muon[s].Histo1D(("jet_not_muon_eta_M","",80,-4,4),"jet_notmuon_eta","weightSSOS")
        hist_jet_notmuon_mass_M[s] = df_muon[s].Histo1D(("jet_notmuon_mass_M","",40,0,40),"jet_notmuon_mass","weightSSOS")

        hist_jet_muon_pt_E[s] = df_electron[s].Histo1D(("jet_muon_pt_E","",60,0,120),"jet_muon_pt","weightSSOS")
        hist_jet_muon_mass_E[s] = df_electron[s].Histo1D(("jet_muon_mass_E","",40,0,40),"jet_muon_mass","weightSSOS")
        hist_jet_notmuon_pt_E[s] = df_electron[s].Histo1D(("jet_not_muon_pt_E","",60,0,120),"jet_notmuon_pt","weightSSOS")
        hist_jet_muon_eta_E[s] = df_electron[s].Histo1D(("jet_muon_eta_E","",80,-4,4),"jet_muon_eta","weightSSOS")
        hist_jet_notmuon_eta_E[s] = df_electron[s].Histo1D(("jet_not_muon_eta_E","",80,-4,4),"jet_notmuon_eta","weightSSOS")
        hist_jet_notmuon_mass_E[s] = df_electron[s].Histo1D(("jet_notmuon_mass_E","",40,0,40),"jet_notmuon_mass","weightSSOS")

        hist_jet_bot1_pt_M[s] = df_muon[s].Histo1D(("jet_bot1_pt_M","",60,0,120),"jet_bot1_pt","weightSSOS")
        hist_jet_bot1_mass_M[s] = df_muon[s].Histo1D(("jet_bot1_mass_M","",40,0,40),"jet_bot1_mass","weightSSOS")
        hist_jet_bot2_pt_M[s] = df_muon[s].Histo1D(("jet_bot2_pt_M","",60,0,120),"jet_bot2_pt","weightSSOS")
        hist_jet_bot1_eta_M[s] = df_muon[s].Histo1D(("jet_bot1_eta_M","",80,-4,4),"jet_bot1_eta","weightSSOS")
        hist_jet_bot2_eta_M[s] = df_muon[s].Histo1D(("jet_bot2_eta_M","",80,-4,4),"jet_bot2_eta","weightSSOS")
        hist_jet_bot2_mass_M[s] = df_muon[s].Histo1D(("jet_bot2_mass_M","",40,0,40),"jet_bot2_mass","weightSSOS")

        hist_jet_bot1_pt_E[s] = df_electron[s].Histo1D(("jet_bot1_pt_E","",60,0,120),"jet_bot1_pt","weightSSOS")
        hist_jet_bot1_mass_E[s] = df_electron[s].Histo1D(("jet_bot1_mass_E","",40,0,40),"jet_bot1_mass","weightSSOS")
        hist_jet_bot2_pt_E[s] = df_electron[s].Histo1D(("jet_bot2_pt_E","",60,0,120),"jet_bot2_pt","weightSSOS")
        hist_jet_bot1_eta_E[s] = df_electron[s].Histo1D(("jet_bot1_eta_E","",80,-4,4),"jet_bot1_eta","weightSSOS")
        hist_jet_bot2_eta_E[s] = df_electron[s].Histo1D(("jet_bot2_eta_E","",80,-4,4),"jet_bot2_eta","weightSSOS")
        hist_jet_bot2_mass_E[s] = df_electron[s].Histo1D(("jet_bot2_mass_E","",40,0,40),"jet_bot2_mass","weightSSOS")

        hist_lepton_pt_M[s] = df_muon[s].Histo1D(("lepton_pt_M","",60,0,120),"lepton_pt","weightSSOS")
        hist_lepton_eta_M[s] = df_muon[s].Histo1D(("lepton_eta_M","",80,-4,4),"lepton_eta","weightSSOS")
        hist_lepton_pt_E[s] = df_electron[s].Histo1D(("lepton_pt_E","",60,0,120),"lepton_pt","weightSSOS")
        hist_lepton_eta_E[s] = df_electron[s].Histo1D(("lepton_eta_E","",80,-4,4),"lepton_eta","weightSSOS")

        hist_InvM_2jets_M[s] = df_muon[s].Histo1D(("InvM_2jets_M","",100,0,300),"InvM_2jets","weightSSOS")
        hist_InvM_2jets_E[s] = df_electron[s].Histo1D(("InvM_2jets_E","",100,0,300),"InvM_2jets","weightSSOS")

        hist_transverse_mass_M[s] = df_muon[s].Histo1D(("transverse_massM","",50,0,150),"transverse_mass","weightSSOS")
        hist_transverse_mass_E[s] = df_electron[s].Histo1D(("transverse_massE","",50,0,150),"transverse_mass","weightSSOS")

#############################
####     DATA SAVING     ####
#############################

if args.notfull:
	for s in samples:
		if args.ssos : path_hist = '/nfs/cms/vazqueze/higgssearch/hists/ssos/histstt_v1v2SSOS_'+s+'_range_'+str(args.list[0])+'_'+str(args.list[1])+'.root'
		else: path_hist = '/nfs/cms/vazqueze/higgssearch/hists/histstt_v1v2_'+s+'_range_'+str(args.list[0])+'_'+str(args.list[1])+'.root'
		myfile = TFile( path_hist, 'RECREATE' )

		hist_nJetGood_M[s].Write()
		hist_nJetGood_E[s].Write()
		hist_nLooseLepton_M[s].Write()
		hist_nLooseLepton_E[s].Write()
		hist_second_muon_pt_M[s].Write()
		hist_second_muon_pt_E[s].Write()
		hist_second_electron_pt_M[s].Write()
		hist_second_electron_pt_E[s].Write()
		hist_second_muon_eta_M[s].Write()
		hist_second_muon_eta_E[s].Write()
		hist_second_electron_eta_M[s].Write()
		hist_second_electron_eta_E[s].Write()
		hist_jet_muon_pt_M[s].Write()
		hist_jet_muon_mass_M[s].Write()
		hist_jet_muon_eta_M[s].Write()
		hist_jet_muon_pt_E[s].Write()
		hist_jet_muon_eta_E[s].Write()
		hist_jet_muon_mass_E[s].Write()
		hist_jet_notmuon_pt_M[s].Write()
		hist_jet_notmuon_eta_M[s].Write()
		hist_jet_notmuon_mass_M[s].Write()
		hist_jet_notmuon_pt_E[s].Write()
		hist_jet_notmuon_eta_E[s].Write()
		hist_jet_notmuon_mass_E[s].Write()
		hist_jet_bot1_pt_M[s].Write()
		hist_jet_bot1_mass_M[s].Write()
		hist_jet_bot1_eta_M[s].Write()
		hist_jet_bot1_pt_E[s].Write()
		hist_jet_bot1_eta_E[s].Write()
		hist_jet_bot1_mass_E[s].Write()
		hist_jet_bot2_pt_M[s].Write()
		hist_jet_bot2_eta_M[s].Write()
		hist_jet_bot2_mass_M[s].Write()
		hist_jet_bot2_pt_E[s].Write()
		hist_jet_bot2_eta_E[s].Write()
		hist_jet_bot2_mass_E[s].Write()
		hist_lepton_pt_M[s].Write()
		hist_lepton_eta_M[s].Write()
		hist_lepton_pt_E[s].Write()
		hist_lepton_eta_E[s].Write()
		hist_InvM_2jets_M[s].Write()
		hist_InvM_2jets_E[s].Write()
		hist_transverse_mass_M[s].Write()
		hist_transverse_mass_E[s].Write()

		myfile.Close()

else:
	for s in samples:
                if args.ssos: path_hist = '/nfs/cms/vazqueze/higgssearch/hists/ssos/hists_higgsSSOS_'+s+'.root'
                else: path_hist = '/nfs/cms/vazqueze/higgssearch/hists/hists_higgs_'+s+'.root'
                myfile = TFile( path_hist, 'RECREATE' )

                hist_nJetGood_M[s].Write()
                hist_nJetGood_E[s].Write()
                hist_nLooseLepton_M[s].Write()
                hist_nLooseLepton_E[s].Write()
                hist_second_muon_pt_M[s].Write()
                hist_second_muon_pt_E[s].Write()
                hist_second_electron_pt_M[s].Write()
                hist_second_electron_pt_E[s].Write()
                hist_second_muon_eta_M[s].Write()
                hist_second_muon_eta_E[s].Write()
                hist_second_electron_eta_M[s].Write()
                hist_second_electron_eta_E[s].Write()
                hist_jet_muon_pt_M[s].Write()
                hist_jet_muon_mass_M[s].Write()
                hist_jet_muon_eta_M[s].Write()
                hist_jet_muon_pt_E[s].Write()
                hist_jet_muon_eta_E[s].Write()
                hist_jet_muon_mass_E[s].Write()
                hist_jet_notmuon_pt_M[s].Write()
                hist_jet_notmuon_eta_M[s].Write()
                hist_jet_notmuon_mass_M[s].Write()
                hist_jet_notmuon_pt_E[s].Write()
                hist_jet_notmuon_eta_E[s].Write()
                hist_jet_notmuon_mass_E[s].Write()
                hist_jet_bot1_pt_M[s].Write()
                hist_jet_bot1_mass_M[s].Write()
                hist_jet_bot1_eta_M[s].Write()
                hist_jet_bot1_pt_E[s].Write()
                hist_jet_bot1_eta_E[s].Write()
                hist_jet_bot1_mass_E[s].Write()
                hist_jet_bot2_pt_M[s].Write()
                hist_jet_bot2_eta_M[s].Write()
                hist_jet_bot2_mass_M[s].Write()
                hist_jet_bot2_pt_E[s].Write()
                hist_jet_bot2_eta_E[s].Write()
                hist_jet_bot2_mass_E[s].Write()
                hist_lepton_pt_M[s].Write()
                hist_lepton_eta_M[s].Write()
                hist_lepton_pt_E[s].Write()
                hist_lepton_eta_E[s].Write()
                hist_InvM_2jets_M[s].Write()
                hist_InvM_2jets_E[s].Write()
                hist_transverse_mass_M[s].Write()
                hist_transverse_mass_E[s].Write()

                myfile.Close()

for s in samples_test:
	print("Number of events analyzed for %s sample are %s" %(s,df_test[s].Count().GetValue()))
	print("Final events for %s sample are %s" %(s,df[s].Count().GetValue()))


if args.ssos: print('SSOS version')

print('Ended succesfully')

