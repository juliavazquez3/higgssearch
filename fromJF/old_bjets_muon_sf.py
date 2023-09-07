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
parser.add_argument("--syst", type=string, default="nom",
                    help="Selec type of systematic to run")
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
        or args.process == "DY" or args.process == "wz"or args.process == "zz" or args.process == "st_1"
        or args.process == "st_2" or args.process == "st_3" or args.process == "st_4"  or args.process == "DYjets"  or args.process == "ttbar_dl"
        or args.process == "ttbar_dh" or args.process == "M" or args.process == "E"): proc = [str(args.process)]
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

if mode == "mc": 
   if years[0]=="2016B": 
      files = json.load(open("/nfs/cms/vazqueze/higgssearch/mcinfo2016.json"))
   else:
      files = json.load(open("/nfs/cms/vazqueze/higgssearch/mcinfo"+years[0]+".json"))
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

files_data = {"M":{"2016":230,"2016B":157,"2017":436,"2018":393},
   "E":{"2016":335,"2016B":156,"2017":245,"2018":738}}

term1 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_ver2/"
term2 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_down/"
term3 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux/"
#term3 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_nosmearing/"
for p in proc:
    # Construct the dataframes
    if args.syst == "nom": folder = term1+"myconfig"+years[0]+"/{sample}/cat_base/prod_test/" # Folder name
    if args.syst == "down": folder = term2+"myconfig"+years[0]+"/{sample}/cat_base/prod_test/" # Folder name
    folder = term3 + "myconfig"+years[0]+"/{sample}/cat_base/prod_test/" # Folder name
    if args.type == "mc":
        if years[0] == "2016B": folder = term3 + "myconfig2016/{sample}/cat_base/prod_test/"
        folder = folder.format(sample = p) # Sample name
        num_files = files[p]["files"] # Number of files
        list_files = [f for f in listdir(folder) if isfile(join(folder, f))] # Lista de archivos
        if (num_files == len(list_files)):
           for f in list_files:
               file_r = TFile(join(folder, f))
               if file_r.GetListOfKeys().Contains("Events"):
                  archives[p+years[0]].append(join(folder,f))
    else:
        if p == "M":
          #folder1 = term1+"myconfig"+years[0]+"/"
          folder1 = term3+"myconfig"+years[0]+"/"
          for f in [f1 for f1 in listdir(folder1) if (isdir(join(folder1, f1)) and f1[0:7] == "data_mu")]:
             folder = folder1 + f + "/cat_base/prod_test/"
             num_files = files_data["M"][years[0]] # Number of files
             list_files = [f1 for f1 in listdir(folder) if isfile(join(folder, f1))] # Lista de archivos
             #if (num_files == len(list_files)):
             for fil in list_files:
                   file_r = TFile(join(folder, fil))
                   if file_r.GetListOfKeys().Contains("Events"):
                     archives[years[0]+"M"].append(join(folder,fil))
        if p == "E":
          #folder1 = term1+"myconfig"+years[0]+"/"
          folder1 = term3+"myconfig"+years[0]+"/"
          for f in [f1 for f1 in listdir(folder1) if (isdir(join(folder1, f1)) and f1[0:7] == "data_el")]:
             folder = folder1 + f + "/cat_base/prod_test/"
             num_files = files_data["E"][years[0]] # Number of files
             list_files = [f1 for f1 in listdir(folder) if isfile(join(folder, f1))] # Lista de archivos
             #if (num_files == len(list_files)):
             for fil in list_files:
                   file_r = TFile(join(folder, fil))
                   if file_r.GetListOfKeys().Contains("Events"):
                     archives[years[0]+"E"].append(join(folder,fil))

for s in samples:
  if args.notfull: archives[s]=archives[s][int(args.list[0]):int(args.list[1])]
  df[s] = ROOT.RDataFrame("Events",set(archives[s]))
  print("Number of files for",s,len(archives[s]))

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

#######   JETS   #######

## Funciones para seleccionar JETS
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
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
      auto InvariantM3(const float pt, const float eta, const float phi, const float pt1, const float eta1, const float phi1, const float pt2, const float eta2, const float phi2) {
            auto x = pt*std::cos(phi);
            auto x1 = pt1*std::cos(phi1);
            auto x2 = pt2*std::cos(phi2);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto y2 = pt2*std::sin(phi2);
            auto z = pt*std::sinh(eta);
            auto z1 = pt1*std::sinh(eta1);
            auto z2 = pt2*std::sinh(eta2);
            auto e = pt*std::cosh(eta);
            auto e1 =  pt1*std::cosh(eta1);
            auto e2 = pt2*std::cosh(eta2);

            auto mJet = std::sqrt((e+e1+e2)*(e+e1+e2)-(x+x1+x2)*(x+x1+x2)-(y+y1+y2)*(y+y1+y2)-(z+z1+z2)*(z+z1+z2));
            return mJet;
      };
""")

## Eleccion de los jets que no sean bottom

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetInds(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, Vint jetbotind, Vfloat qgl) {
            vector<int> vb;
            bool condb = true;
            int ind1 = -1;
            int ind2 = -1;
            float ptJ{-10.};
            for (unsigned int j=0; j<4; ++j){
                        if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                        //if(pt[good[j]]>ptJ && condb && qgl[good[j]]>0.25){
                        if(pt[good[j]]>ptJ && condb){
                                ind1 = good[j];
                                ptJ = pt[good[j]];
                        }
            }
            if (ind1 > -1) {
                vb.push_back(ind1);
            }
            ptJ = -10.;
            for (unsigned int j=0; j<4; ++j){
                        if (good.size() > 2) condb = (good[j] != good[jetbotind[0]] && good[j] != good[jetbotind[1]]);
                        //if(pt[good[j]]>ptJ && condb && good[j] != ind1 && qgl[good[j]]>0.25){
                        if(pt[good[j]]>ptJ && condb && good[j] != ind1){
                                ind2 = good[j];
                                ptJ = pt[good[j]];
                        }
            }
            if (ind2 > -1) {
                vb.push_back(ind2);
            }
            return vb;
      };
""")

######### Masa invariante con los bottom
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto InvMassBot(Vint jetbot, Vint jetmuon, Vint jetgood, Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi) {
            int ind = -1;
            vector<float> vb;
            ROOT::Math::PtEtaPhiMVector jet1(Jet_pt[jetmuon[0]], Jet_eta[jetmuon[0]], Jet_phi[jetmuon[0]], 0.);
            ROOT::Math::PtEtaPhiMVector jet2(Jet_pt[jetmuon[1]], Jet_eta[jetmuon[1]], Jet_phi[jetmuon[1]], 0.);
            auto had = jet1+jet2;
            float had_eta = had.Eta();
       	    float had_phi = had.Phi();
            float dR1 = ROOT::VecOps::DeltaR(had_eta,Jet_eta[jetgood[jetbot[0]]],had_phi,Jet_phi[jetgood[jetbot[0]]]);
            float dR2 = ROOT::VecOps::DeltaR(had_eta,Jet_eta[jetgood[jetbot[1]]],had_phi,Jet_phi[jetgood[jetbot[1]]]);
            if (dR1 > dR2){
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_pt[jetmuon[1]],Jet_eta[jetmuon[1]],Jet_phi[jetmuon[1]]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_pt[jetmuon[1]],Jet_eta[jetmuon[1]],Jet_phi[jetmuon[1]]));           
            }else{
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_pt[jetmuon[1]],Jet_eta[jetmuon[1]],Jet_phi[jetmuon[1]]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_pt[jetmuon[1]],Jet_eta[jetmuon[1]],Jet_phi[jetmuon[1]]));
            }
            return vb;
      };
      auto InvMassVarious(Vint jetbot, Vint jetmuon, Vint jetgood, Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi, Vint muongood, Vint elgood, Vfloat muon_pt, Vfloat muon_eta, Vfloat muon_phi, Vfloat muon_mass, Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass) {
            int ind = -1;
            vector<float> vb;
            ROOT::Math::PtEtaPhiMVector jet1(Jet_pt[jetmuon[0]], Jet_eta[jetmuon[0]], Jet_phi[jetmuon[0]], 0.);
            ROOT::Math::PtEtaPhiMVector jet2(Jet_pt[jetmuon[1]], Jet_eta[jetmuon[1]], Jet_phi[jetmuon[1]], 0.);
            auto had = jet1+jet2;
            float had_eta = had.Eta();
       	    float had_phi = had.Phi();
            vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_pt[jetmuon[1]],Jet_eta[jetmuon[1]],Jet_phi[jetmuon[1]]));
            vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_pt[jetmuon[1]],Jet_eta[jetmuon[1]],Jet_phi[jetmuon[1]]));
            if (muongood.size() > 0) {
               vb.push_back(InvariantM(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],0., muon_pt[muongood[0]], muon_eta[muongood[0]], muon_phi[muongood[0]], muon_mass[muongood[0]]));
               vb.push_back(InvariantM(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],0., muon_pt[muongood[0]], muon_eta[muongood[0]], muon_phi[muongood[0]], muon_mass[muongood[0]]));
            } else {
               vb.push_back(InvariantM(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],0., el_pt[elgood[0]], el_eta[elgood[0]], el_phi[elgood[0]], el_mass[elgood[0]]));
               vb.push_back(InvariantM(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],0., el_pt[elgood[0]], el_eta[elgood[0]], el_phi[elgood[0]], el_mass[elgood[0]]));
            }
            return vb;
      };
""")


######### CHI 2 TEST ############
gInterpreter.Declare("""
      auto chi2calc(const float mw, const float mt){
        Double_t xx2 = 0.;

        xx2  = (mw-77.5)*(mw-77.5)/11.5/11.5;
        xx2 += (mt-161.)*(mt-161.)/20./20.;
        xx2 += -0.6*2.*(mw-77.5)*(mt-161.)/11.5/20.;
        xx2 = xx2/(1-0.6*0.6);

        if (xx2<0) return -999.;

        return xx2;
      };
      auto chibottomindex(const float imdijet, const float im30, const float im31, const float iml0, const float iml1) {
            vector<int> vb;
            bool cond_final = true;
            bool cond1 = true;
            Double_t xi20 = chi2calc(imdijet, im30);
            Double_t xi21 = chi2calc(imdijet, im31);
            if ( iml0 > 150. && iml1 > 150. ){
                  cond_final = false;
            } else if (iml0 > 150. && iml1 < 150.) {
                  if (xi20 > 3.2) {
                     cond_final = false;
                  } else {
                     vb.push_back(1);
                     vb.push_back(0);
                  }
            } else if (iml1 > 150. && iml0 < 150.) {
                  if (xi21 > 3.2) {
                     cond_final = false;
                  } else {
                     vb.push_back(0);
                     vb.push_back(1);
                  }
            } else {
                  if (xi20 > xi21) {
                     if (xi21 >3.2) {
                        cond_final = false;
                     } else {
                        vb.push_back(0);
                        vb.push_back(1);
                     }
                  } else {
                     if (xi20 >3.2) {
                        cond_final = false;
                     } else {
                        vb.push_back(1);
                        vb.push_back(0);
                     }
                  }
            }
            return vb;
      };
      auto kinfittest(const float imdijet, const float im30, const float im31, const float iml0, const float iml1) {
            bool cond_final = true;
            bool cond1 = true;
            Double_t xi20 = chi2calc(imdijet, im30);
            Double_t xi21 = chi2calc(imdijet, im31);
            if ( iml0 > 150. && iml1 > 150. ){
                  cond_final = false;
            } else if (iml0 > 150. && iml1 < 150.) {
                  if (xi20 > 3.2) {
                     cond_final = false;
                  }
            } else if (iml1 > 150. && iml0 < 150.) {
                  if (xi21 > 3.2) {
                     cond_final = false;
                  }
            } else {
                  if (xi20 > xi21) {
                     if (xi21 >3.2) {
                        cond_final = false;
                     }
                  } else { 
                     if (xi20 >3.2) {
                        cond_final = false;
                     }
                  }
            }
            return cond_final;
      };
      auto kinfitmasses(const float imdijet, const float im30, const float im31, const float iml0, const float iml1) {
            bool cond_final = true;
            bool cond1 = true;
            vector<float> vb;
            Double_t xi20 = chi2calc(imdijet, im30);
            Double_t xi21 = chi2calc(imdijet, im31);
            if ( iml0 > 150. && iml1 > 150. ){
                  cond_final = false;
                  vb.push_back(-999.);
                  vb.push_back(-999.);
                  vb.push_back(-999.);
                  vb.push_back(-999.);
                  vb.push_back(-999.);
                  vb.push_back(-999.);
            } else if (iml0 > 150. && iml1 < 150.) {
                  if (xi20 > 3.2) {
                     cond_final = false;
                     vb.push_back(-999.);
                     vb.push_back(-999.);
                     vb.push_back(-999.);
                     vb.push_back(-999.);
                     vb.push_back(-999.);
                     vb.push_back(-999.);
                  } else {
                     cond_final = true;
                     vb.push_back(iml1);
                     vb.push_back(iml0);
                     vb.push_back(im30);
                     vb.push_back(im31);
                     vb.push_back(xi20);
                     vb.push_back(xi21);
                  }
            } else if (iml1 > 150. && iml0 < 150.) {
                  if (xi21 > 3.2) {
                     cond_final = false;
                     vb.push_back(-999.);
                     vb.push_back(-999.);
                     vb.push_back(-999.);
                     vb.push_back(-999.);
                     vb.push_back(-999.);
                     vb.push_back(-999.);
                  } else {
                     cond_final = true;
                     vb.push_back(iml0);
                     vb.push_back(iml1);
                     vb.push_back(im31);
                     vb.push_back(im30);
                     vb.push_back(xi21);
                     vb.push_back(xi20);
                  }
            } else {
                  if (xi20 > xi21) {
       	       	     if	(xi21 >3.2) {
                        cond_final = false;
                        vb.push_back(-999.);
                        vb.push_back(-999.);
                        vb.push_back(-999.);
                        vb.push_back(-999.);
                        vb.push_back(-999.);
                        vb.push_back(-999.);
                     } else {
                        cond_final = true;
                        vb.push_back(iml0);
                        vb.push_back(iml1);
                        vb.push_back(im31);
                        vb.push_back(im30);
                        vb.push_back(xi21);
                        vb.push_back(xi20);
                     }
                  } else {
                     if (xi20 >3.2) {
                        cond_final = false;
                        vb.push_back(-999.);
                        vb.push_back(-999.);
                        vb.push_back(-999.);
                        vb.push_back(-999.);
                        vb.push_back(-999.);
                        vb.push_back(-999.);
       	       	     } else {
                        cond_final = true;
                        vb.push_back(iml1);
                        vb.push_back(iml0);
                        vb.push_back(im30);
                        vb.push_back(im31);
                        vb.push_back(xi20);
                        vb.push_back(xi21);
                     }
                  }
            }
            return vb;
      };
""")

######### SL function

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
      auto SLselection(const int blep, const int bhad, Vfloat pt, Vfloat eta, Vfloat phi, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vint mu_good, Vbool mu_id){
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

### tau selector

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto tau_selection(Vint xind, Vfloat pt, Vfloat eta, Vfloat phi, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi){
            bool cond1 = false;
            int indM=-1;
            float ptM{-10.};
            for (unsigned int i=0; i<mu_pt.size(); ++i){
                cond1 = ROOT::VecOps::DeltaR(mu_eta[i],eta[xind[0]],mu_phi[i],phi[xind[0]]) < 0.4;
                if(cond1 && mu_pt[i]>ptM){
                     indM = i;
                     ptM = mu_pt[i];
                }
            }
            return indM;
      };
""")

## pT component calculations

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto pTsum(Vint mu_good, Vint el_good, Vint jet_muon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_muon[1]], jet_eta[jet_muon[1]], jet_phi[jet_muon[1]], jet_mass[jet_muon[1]]);
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
      auto pTprod(Vint mu_good, Vint el_good, Vint jet_muon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_muon[1]], jet_eta[jet_muon[1]], jet_phi[jet_muon[1]], jet_mass[jet_muon[1]]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            auto plepX = plep.Px(); auto plepY = plep.Py(); auto phadX = phad.Px(); auto phadY = phad.Py();
            float plepmod; float phadmod; float ptotmod;
            plepmod = std::sqrt(plepX*plepX + plepY*plepY);
            phadmod = std::sqrt(phadX*phadX + phadY*phadY);
            ptotmod = plepX*phadX + plepY*phadY;
            return ptotmod/(plepmod*phadmod);
      };
      auto mydeltaphi(const float phi1, const float phi2 ) {
         float diff = 0.;
         float PI = 3.14159265358979323846;
         if (fabs(phi1 - phi2) > PI) {
            diff = 2*PI - fabs(phi1 - phi2);
         } else {
            diff = fabs(phi1 - phi2);
         }
         return diff;
      };
      auto variousSUM(Vint mu_good, Vint el_good, Vint jet_muon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
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
       	       	  deltaRlepjet2	= ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_muon[1]],mu_phi[mu_good[0]],jet_phi[jet_muon[1]]);
                  deltaPhilepjet1 = mydeltaphi(mu_phi[mu_good[0]],jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = mydeltaphi(mu_phi[mu_good[0]],jet_phi[jet_muon[1]]);
       	       	  deltaEtalepjet1 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_muon[1]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
       	       	  deltaRlepjet1	= ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_muon[0]],el_phi[el_good[0]],jet_phi[jet_muon[0]]);
                  deltaRlepjet2 = ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_muon[1]],el_phi[el_good[0]],jet_phi[jet_muon[1]]);
                  deltaPhilepjet1 = mydeltaphi(el_phi[el_good[0]],jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = mydeltaphi(el_phi[el_good[0]],jet_phi[jet_muon[1]]);
                  deltaEtalepjet1 = fabs(el_eta[el_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(el_eta[el_good[0]]-jet_eta[jet_muon[1]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_muon[1]], jet_eta[jet_muon[1]], jet_phi[jet_muon[1]], jet_mass[jet_muon[1]]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            vb.push_back(ROOT::VecOps::DeltaR(plep1.Eta(),phad.Eta(),plep1.Phi(),phad.Phi()));
            vb.push_back(mydeltaphi(plep2.Phi(),phad.Phi()));
            vb.push_back(phad.Eta());
            vb.push_back(phad.Pt());
            vb.push_back(mydeltaphi(plep.Phi(),phad.Phi()));
            vb.push_back(ROOT::VecOps::DeltaR(plep.Eta(),phad.Eta(),plep.Phi(),phad.Phi()));
            vb.push_back(mydeltaphi(plep1.Phi(),phad.Phi()));
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
            vb.push_back(mydeltaphi(plep2.Phi(),jet_phi[jet_muon[0]]));
            vb.push_back(mydeltaphi(plep2.Phi(),jet_phi[jet_muon[1]]));
            return vb;
      };
""")

#############################################################
######################## Definitions ########################
#############################################################

if mode == "data":
        for s in samples:
           df[s] = df[s].Define('Jet_pt_nom','Jet_pt')
           df[s] = df[s].Define('Jet_mass_nom','Jet_mass')
           df[s] = df[s].Define('MET_smeared_phi','MET_phi')
           df[s] = df[s].Define('MET_smeared_pt','MET_pt')

df_muon = {}
df_electron = {}
df_test = {}

##### Systematic run

if mode == "mc":
   if args.syst == "nom":
        for s in samples:
           df[s] = df[s].Define('Jet_pt_aux','0.99*Jet_pt_nom')
           #df[s] = df[s].Define('Jet_pt_aux','Jet_pt')
           df[s] = df[s].Define('MET_pt_aux','MET_smeared_pt')
           #df[s] = df[s].Define('MET_pt_aux','MET_pt')
           df[s] = df[s].Define('MET_phi_aux','MET_smeared_phi')
           #df[s] = df[s].Define('MET_phi_aux','MET_phi')
   elif args.syst == "down":
        for s in samples:
           df[s] = df[s].Define('Jet_pt_aux','Jet_pt_down')
           df[s] = df[s].Define('MET_pt_aux','MET_smeared_pt')
           df[s] = df[s].Define('MET_phi_aux','MET_smeared_phi')
   elif args.syst == "up":
        for s in samples:
           df[s] = df[s].Define('Jet_pt_aux','Jet_pt_up')
           df[s] = df[s].Define('MET_pt_aux','MET_smeared_pt')
           df[s] = df[s].Define('MET_phi_aux','MET_smeared_phi')
   else: raise NameError('Incorrect systematic option')
elif mode == "data":
   df[s] = df[s].Define('Jet_pt_aux','Jet_pt_nom')
   #df[s] = df[s].Define('Jet_pt_aux','Jet_pt')
   df[s] = df[s].Define('MET_pt_aux','MET_smeared_pt')
   #df[s] = df[s].Define('MET_pt_aux','MET_pt')
   df[s] = df[s].Define('MET_phi_aux','MET_smeared_phi')
   #df[s] = df[s].Define('MET_phi_aux','MET_phi')

###########    Extra definitions

for s in samples:
        df[s] = df[s].Define('lepton_charge','nMuonGood>0 ? Muon_charge[MuonGoodInd[0]] : Electron_charge[ElectronGoodInd[0]]')
        df[s] = df[s].Define('JetQInd','JetInds(nJet, JetGoodInd, Jet_pt_aux, Jet_eta, Jet_phi, JetBotInd, Jet_qgl)')
        ########### Filtering and further definitions
        df[s] = df[s].Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>=4')
        df[s] = df[s].Filter('JetQInd.size() > 1')
        df[s] = df[s].Define('JetXInd','JetQInd')
        ### hists definitions
        df[s] = df[s].Define('jet_1_pt','Jet_pt_aux[JetXInd[0]]')
        df[s] = df[s].Define('jet_1_nmu','Jet_nMuons[JetXInd[0]]')
        df[s] = df[s].Define('jet_1_qgl','Jet_qgl[JetXInd[0]]')
        df[s] = df[s].Define('jet_2_pt','nJetGood>1 ? Jet_pt_aux[JetXInd[1]] : 0')
        df[s] = df[s].Define('jet_1_eta','Jet_eta[JetXInd[0]]')
        df[s] = df[s].Define('jet_2_eta','nJetGood>1 ? Jet_eta[JetXInd[1]] : 0')
        df[s] = df[s].Define('jet_1_phi','Jet_phi[JetXInd[0]]')
        df[s] = df[s].Define('jet_2_phi','nJetGood>1 ? Jet_phi[JetXInd[1]] : 0')
        df[s] = df[s].Define('jet_1_mass','Jet_mass_nom[JetXInd[0]]')
        df[s] = df[s].Define('jet_2_mass','nJetGood>1 ? Jet_mass_nom[JetXInd[1]] : 0')
        df[s] = df[s].Define('jet_2_qgl','nJetGood>1 ? Jet_qgl[JetXInd[1]] : 0')
        df[s] = df[s].Define('jet_2_nmu','nJetGood>1 ? Jet_nMuons[JetXInd[1]] : 0')
        df[s] = df[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt_aux[JetXInd[0]],Jet_eta[JetXInd[0]],Jet_phi[JetXInd[0]],0.,Jet_pt_aux[JetXInd[1]],Jet_eta[JetXInd[1]],Jet_phi[JetXInd[1]],0.) : 0')
        df[s] = df[s].Define('deltaR_jet1_jet2','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetXInd[1]], Jet_eta[JetXInd[0]] , Jet_phi[JetXInd[1]], Jet_phi[JetXInd[0]])  : 10')
        df[s] = df[s].Define('deltaphi_jet1_jet2','mydeltaphi(Jet_phi[JetXInd[0]],Jet_phi[JetXInd[1]])')
        df[s] = df[s].Define('deltaeta_jet1_jet2','fabs(Jet_eta[JetXInd[0]]-Jet_eta[JetXInd[1]])')
        df[s] = df[s].Define('deltapt_jet1_jet2','fabs(Jet_pt_aux[JetXInd[0]]-Jet_pt_aux[JetXInd[1]])')
        df[s] = df[s].Define('tracks_jet1','Jet_nConstituents[JetXInd[0]]')
        df[s] = df[s].Define('tracks_jet2','nJetGood>1 ? Jet_nConstituents[JetXInd[1]] : 0')
        df[s] = df[s].Define('EMN_jet1','Jet_neEmEF[JetXInd[0]]')
        df[s] = df[s].Define('EMC_jet1','Jet_chEmEF[JetXInd[0]]')
        df[s] = df[s].Define('EMtotal_jet1','Jet_chEmEF[JetXInd[0]]+Jet_neEmEF[JetXInd[0]]')
        df[s] = df[s].Define('pT_sum','pTsum(MuonGoodInd, ElectronGoodInd, JetXInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt_aux, MET_phi_aux, Jet_pt_aux, Jet_eta, Jet_phi, Jet_mass_nom)')
        df[s] = df[s].Define('pT_product','pTprod(MuonGoodInd, ElectronGoodInd, JetXInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt_aux, MET_phi_aux, Jet_pt_aux, Jet_eta, Jet_phi, Jet_mass_nom)')
        df[s] = df[s].Define('aux_various','variousSUM(MuonGoodInd, ElectronGoodInd, JetXInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt_aux, MET_phi_aux, Jet_pt_aux, Jet_eta, Jet_phi, Jet_mass_nom)')
        df[s] = df[s].Define('deltaR_lep_2jets','aux_various[0]')
        df[s] = df[s].Define('deltaphi_MET_2jets','aux_various[1]')
        df[s] = df[s].Define('deltaphi_MET_jets_1','aux_various[20]')
        df[s] = df[s].Define('deltaphi_MET_jets_2','aux_various[21]')
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
        df[s] = df[s].Define('jet_bot1_pt','Jet_pt_aux[JetGoodInd[JetBotInd[0]]]')
        df[s] = df[s].Define('jet_bot2_pt','Jet_pt_aux[JetGoodInd[JetBotInd[1]]]')
        df[s] = df[s].Define('jet_bot1_eta','Jet_eta[JetGoodInd[JetBotInd[0]]]')
        df[s] = df[s].Define('jet_bot2_eta','Jet_eta[JetGoodInd[JetBotInd[1]]]')
        df[s] = df[s].Define('jet_bot1_phi','Jet_phi[JetGoodInd[JetBotInd[0]]]')
        df[s] = df[s].Define('jet_bot2_phi','Jet_phi[JetGoodInd[JetBotInd[1]]]')
        df[s] = df[s].Define('jet_1_btag','Jet_btagDeepFlavB[JetXInd[0]]')
        df[s] = df[s].Define('jet_2_btag','Jet_btagDeepFlavB[JetXInd[1]]')
        df[s] = df[s].Define('nLooseLepton','nMuon+nElectron-1')
        df[s] = df[s].Define('InvM3_aux','InvMassBot(JetBotInd, JetQInd, JetGoodInd, Jet_pt_aux, Jet_eta, Jet_phi)')
        df[s] = df[s].Define('InvM_various_aux','InvMassVarious(JetBotInd, JetQInd, JetGoodInd, Jet_pt_aux, Jet_eta, Jet_phi, MuonGoodInd, ElectronGoodInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass)')
        df[s] = df[s].Define('InvM_bot_closer','InvM3_aux[0]')
        df[s] = df[s].Define('InvM_bot_farther','InvM3_aux[1]')
        df[s] = df[s].Define('InvM30','InvM_various_aux[0]')
        df[s] = df[s].Define('InvM31','InvM_various_aux[1]')
        df[s] = df[s].Define('InvMl0','InvM_various_aux[2]')
        df[s] = df[s].Define('InvMl1','InvM_various_aux[3]')
        df[s] = df[s].Define('jet_bot1_btagnumber','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[0]]] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetGoodInd[JetBotInd[0]]] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')
        df[s] = df[s].Define('jet_bot2_btagnumber','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[1]]] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetGoodInd[JetBotInd[1]]] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')
        df[s] = df[s].Define('jet_1_btagnumber','Jet_btagDeepFlavB[JetXInd[0]] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetXInd[0]] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')
        df[s] = df[s].Define('jet_2_btagnumber','Jet_btagDeepFlavB[JetXInd[1]] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetXInd[1]] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')
        df[s] = df[s].Define('jet_1_cvltag','Jet_btagDeepFlavCvL[JetXInd[0]]')
        df[s] = df[s].Define('jet_2_cvltag','Jet_btagDeepFlavCvL[JetXInd[1]]')
        df[s] = df[s].Define('jet_1_cvltag_csv','Jet_btagDeepCvL[JetXInd[0]]')
        df[s] = df[s].Define('jet_2_cvltag_csv','Jet_btagDeepCvL[JetXInd[1]]')
        df[s] = df[s].Define('jet_1_cvbtag','Jet_btagDeepFlavCvB[JetXInd[0]]')
        df[s] = df[s].Define('jet_2_cvbtag','Jet_btagDeepFlavCvB[JetXInd[1]]')
        df[s] = df[s].Define('jet_1_cvbtag_csv','Jet_btagDeepCvB[JetXInd[0]]')
        df[s] = df[s].Define('jet_2_cvbtag_csv','Jet_btagDeepCvB[JetXInd[1]]')
        df[s] = df[s].Define('jet_max_cvltag','jet_1_cvltag > jet_2_cvltag ? jet_1_cvltag : jet_2_cvltag')
        df[s] = df[s].Define('jet_min_cvltag','jet_1_cvltag > jet_2_cvltag ? jet_2_cvltag : jet_1_cvltag')
        df[s] = df[s].Define('jet_max_cvbtag','jet_1_cvbtag > jet_2_cvbtag ? jet_1_cvbtag : jet_2_cvbtag')
        df[s] = df[s].Define('jet_min_cvbtag','jet_1_cvbtag > jet_2_cvbtag ? jet_2_cvbtag : jet_1_cvbtag')
        df[s] = df[s].Define('chi2_test0','chi2calc(InvM_2jets, InvM30)')
        df[s] = df[s].Define('chi2_test1','chi2calc(InvM_2jets, InvM31)')
        df[s] = df[s].Define('nJetMuonInd','JetMuonInd.size()')
        df[s] = df[s].Define('nJetSVInd','JetSVInd.size()')
        df[s] = df[s].Define('nJetBotInd','JetBotInd.size()')
        df[s] = df[s].Define('nMuonJetInd','MuonJetInd.size()')
        df[s] = df[s].Define('nSVJetInd','SVJetInd.size()')
        df[s] = df[s].Define('nJetQInd','JetQInd.size()')
        df[s] = df[s].Define('kfmasses_aux','kinfitmasses(InvM_2jets, InvM30, InvM31, InvMl0, InvMl1)')
        df[s] = df[s].Define('InvMl_good','kfmasses_aux[0]')
        df[s] = df[s].Define('InvMl_bad','kfmasses_aux[1]')
        df[s] = df[s].Define('InvM3_good','kfmasses_aux[2]')
        df[s] = df[s].Define('InvM3_bad','kfmasses_aux[3]')
        df[s] = df[s].Define('chi2_test_good','kfmasses_aux[4]')
        df[s] = df[s].Define('chi2_test_bad','kfmasses_aux[5]')

############################################################
#### Distinguishing between same sign and opposite sign ####
############################################################

#################### Final definitions and filters ###############################

for s in samples:
        df[s] = df[s].Define('lepton_phi','nMuonGood>0 ? Muon_phi[MuonGoodInd[0]] : Electron_phi[ElectronGoodInd[0]]')
        df[s] = df[s].Define('transverse_mass_old','std::sqrt(2*lepton_pt*MET_pt_aux*(1-std::cos(lepton_phi-MET_phi_aux)))')
        df[s] = df[s].Define('weight_aux','1.')
        df[s] = df[s].Define('deltaR_jetM_lep','nMuonGood>0 ? ROOT::VecOps::DeltaR(Muon_eta[MuonGoodInd[0]],Jet_eta[JetQInd[0]] , Muon_phi[MuonGoodInd[0]], Jet_phi[JetQInd[0]]) : ROOT::VecOps::DeltaR(Electron_eta[ElectronGoodInd[0]],Jet_eta[JetQInd[0]] , Electron_phi[ElectronGoodInd[0]], Jet_phi[JetQInd[0]])')
        df[s] = df[s].Define('InvM_jetM_lep', 'nMuonGood>0 ? InvariantM(Jet_pt_aux[JetXInd[0]],Jet_eta[JetXInd[0]],Jet_phi[JetXInd[0]],0.,Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]]) : InvariantM(Jet_pt_aux[JetXInd[0]],Jet_eta[JetXInd[0]],Jet_phi[JetXInd[0]],0.,Electron_pt[ElectronGoodInd[0]],Electron_eta[ElectronGoodInd[0]],Electron_phi[ElectronGoodInd[0]], Electron_mass[ElectronGoodInd[0]])')
        df[s] = df[s].Define('InvM_muon_jet','nMuonGood>0 ? InvariantM(Muon_pt[MuonJetInd[0]],Muon_eta[MuonJetInd[0]],Muon_phi[MuonJetInd[0]],Muon_mass[MuonJetInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]]) : 50.')
        df[s] = df[s].Define('deltaphi_MET_lep','mydeltaphi(lepton_phi,MET_phi_aux)')
        df[s] = df[s].Define('MET_my_significance','(MET_pt_aux*MET_pt_aux)/((0.62*0.62)*MET_sumEt)')
        ### Filters
        df[s] = df[s].Filter('transverse_mass > 20.')
        df[s] = df[s].Filter('MET_pt_aux > 20.')
        df[s] = df[s].Filter('nMuonGood>0 ? lepton_pt > 30. : lepton_pt > 35.')
        ### Jet pT filters
        #df[s] = df[s].Filter('jet_bot1_pt > 25.')
        #df[s] = df[s].Filter('jet_bot2_pt > 25.')
        #df[s] = df[s].Filter('jet_1_pt > 25.')
        #df[s] = df[s].Filter('jet_2_pt > 25.')
        ## Invariant mass filter
        df[s] = df[s].Filter('InvM_2jets > 30. && InvM_2jets < 300.')
        ### Cur based for electron, EL ID tests
        #df[s] = df[s].Filter('nElectronGood > 0 ? Electron_cutBased[ElectronGoodInd[0]] > 3 : 1')
        ### Chi2 test
        df[s] = df[s].Define('chi_bool','kinfittest(InvM_2jets, InvM30, InvM31, InvMl0, InvMl1)')
        df[s] = df[s].Filter('jet_bot1_btag >'+str(cuts_btag[years[0]][1]))
        df[s] = df[s].Filter('jet_bot2_btag >'+str(cuts_btag[years[0]][1]))
        df[s] = df[s].Filter('kinfittest(InvM_2jets, InvM30, InvM31, InvMl0, InvMl1)')
        ####### SL stuff
        df[s] = df[s].Define('bottom_chi_aux','chibottomindex(InvM_2jets, InvM30, InvM31, InvMl0, InvMl1)')
        df[s] = df[s].Define('bottom_chi_lep','JetGoodInd[JetBotInd[bottom_chi_aux[0]]]')
        df[s] = df[s].Define('bottom_chi_had','JetGoodInd[JetBotInd[bottom_chi_aux[1]]]')
        df[s] = df[s].Define('sl_selection_aux','SLselection(JetGoodInd[JetBotInd[0]], JetGoodInd[JetBotInd[1]] ,Jet_pt, Jet_eta, Jet_phi, Muon_pt, Muon_eta, Muon_phi, MuonGoodInd, Muon_tightId)')
        df[s] = df[s].Define('bot1_muons','sl_selection_aux[0]')
        df[s] = df[s].Define('bot2_muons','sl_selection_aux[2]')
        df[s] = df[s].Define('bot1_muon_id','sl_selection_aux[1]')
        df[s] = df[s].Define('bot2_muon_id','sl_selection_aux[3]')
        df[s] = df[s].Filter('bot1_muons > 0 && bot2_muons > 0')
        df[s] = df[s].Define('muon_bot1_eta','Muon_eta[bot1_muon_id]')
        df[s] = df[s].Define('muon_bot2_eta','Muon_eta[bot2_muon_id]')
        df[s] = df[s].Define('muon_bot1_pt','Muon_pt[bot1_muon_id]')
        df[s] = df[s].Define('muon_bot2_pt','Muon_pt[bot2_muon_id]')
        df[s] = df[s].Define('muon_bot1_relpt','jet_bot1_pt > 0. ? muon_bot1_pt/jet_bot1_pt : -10.')
        df[s] = df[s].Define('muon_bot2_relpt','jet_bot2_pt > 0. ? muon_bot2_pt/jet_bot2_pt : -10.')

############ Trigger scale factors ##############

from trigger_sf import *

if mode == "mc":
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= trigger_mu_sfRDF(**kwargs)
                df[s] = b().run(df[s])

from trigger_sf_el import *

if mode == "mc":
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= trigger_el_sfRDF(**kwargs)
                df[s] = b().run(df[s])

####### Other scale factors

from low_pt_muonid import *

if mode == "mc":
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= displaced_mu_idRDF(**kwargs)
                df[s] = b().run(df[s])

from complete_btag_weight import *

if mode == "mc":
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= btag_weights_totRDF(**kwargs)
                df[s] = b().run(df[s])

############ BR corrections ##############

from branchingfractions_corrections import *

if mode == "mc":
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= branchingfractions_SL_corRDF(**kwargs)
                df[s] = b().run(df[s])

############ Gen level definitions

if mode == "mc":
        for s in samples:
                df[s] = df[s].Define('jet_1_flavourH','Jet_hadronFlavour[JetXInd[0]]')
                df[s] = df[s].Define('jet_2_flavourH','Jet_hadronFlavour[JetXInd[1]]')
                df[s] = df[s].Define('jet_1_flavourP','Jet_partonFlavour[JetXInd[0]]')
                df[s] = df[s].Define('jet_2_flavourP','Jet_partonFlavour[JetXInd[1]]')
                df[s] = df[s].Define('jet_bot1_flavourP','Jet_partonFlavour[JetGoodInd[JetBotInd[0]]]')
                df[s] = df[s].Define('jet_bot2_flavourP','Jet_partonFlavour[JetGoodInd[JetBotInd[1]]]')
                df[s] = df[s].Define('last_Copy','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,13)')
                df[s] = df[s].Define('top_weight','topreweight(nGenPart, GenPart_statusFlags, GenPart_pdgId, last_Copy, GenPart_pt)')
                df[s] = df[s].Define('btag_sf','Aux_btag_weight[1]/Aux_btag_weight[0]')
                df[s] = df[s].Define('btag_sf_light_up','Aux_btag_weight_light_up[1]/Aux_btag_weight_light_up[0]')
                df[s] = df[s].Define('btag_sf_heavy_up','Aux_btag_weight_heavy_up[1]/Aux_btag_weight_heavy_up[0]')
                df[s] = df[s].Define('btag_sf_light_down','Aux_btag_weight_light_down[1]/Aux_btag_weight_light_down[0]')
                df[s] = df[s].Define('btag_sf_heavy_down','Aux_btag_weight_heavy_down[1]/Aux_btag_weight_heavy_down[0]')
                df[s] = df[s].Define('muon_bot1_mother','fabs(GenPart_pdgId[GenPart_genPartIdxMother[Muon_genPartIdx[bot1_muon_id]]])')
                df[s] = df[s].Define('muon_bot2_mother','fabs(GenPart_pdgId[GenPart_genPartIdxMother[Muon_genPartIdx[bot2_muon_id]]])')

########## ttbar sectioning for charm discrimination

cond1 = '(fabs(jet_1_flavourP) == 4 || fabs(jet_2_flavourP) == 4)'
cond2 = '(fabs(jet_1_flavourP) == 5 || fabs(jet_2_flavourP) == 5)'
cond3 = '(fabs(jet_1_flavourP) == 3 || fabs(jet_1_flavourP) == 2 || fabs(jet_1_flavourP) == 1)'
cond4 =	'(fabs(jet_2_flavourP) == 3 || fabs(jet_2_flavourP) == 2 || fabs(jet_2_flavourP) == 1)'
cond5 = '(fabs(jet_1_flavourP)>0 && fabs(jet_1_flavourP)<6 && fabs(jet_2_flavourP)>0 && fabs(jet_2_flavourP)<6)'
cond6 = '((fabs(jet_1_flavourP)==0 || fabs(jet_1_flavourP)>6) && (fabs(jet_2_flavourP)==0 || fabs(jet_2_flavourP)>6))'

if mode == "mc":
        for s in samples:
                if (s[0:8] == "ttbar_sl" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                        df[s+"_charm"] = df[s].Filter(cond1+' && '+cond5)
                        df[s+"_bottom"] = df[s].Filter(cond2+' && !'+cond1+' && '+cond5)
                        df[s+"_light"] = df[s].Filter('('+cond3+') && ('+cond4+') && !('+cond1+') && !('+cond2+')')
                        df[s+"_charmgluon"] = df[s].Filter(cond1+' && !'+cond5)
                        df[s+"_bottomgluon"] = df[s].Filter(cond2+' && !'+cond1+' && !'+cond5)
                        #df[s+"_gluongluon"] = df[s].Filter(cond6)
                        #df[s+"_else"] = df[s].Filter('!('+cond1+') && !('+cond2+') && !(('+cond3+') && ('+cond4+') && !('+cond6+')')
                        df[s+"_else"] = df[s].Filter('!('+cond1+') && !('+cond2+') && !(('+cond3+') && ('+cond4+'))')
                        ## Samples correction
                        samples.append(s+"_charm")
                        samples.append(s+"_bottom")
                        samples.append(s+"_light")
                        samples.append(s+"_charmgluon")
                        samples.append(s+"_bottomgluon")
                        #samples.append(s+"_gluongluon")
                        samples.append(s+"_else")

samples = [s for s in samples if not (s[0:8]=="ttbar_sl"  and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

#if mode == "mc":
#    for s in samples:
#        df[s] = df[s].Filter('!(fabs(jet_bot1_flavourP)==5 && fabs(jet_bot2_flavourP)==4) && !(fabs(jet_bot1_flavourP)==4 && fabs(jet_bot2_flavourP)==5) && !(fabs(jet_bot1_flavourP)==5 && fabs(jet_bot2_flavourP)==5)')

#if mode == "mc":
#    for s in samples:
#        df[s] = df[s].Filter('(fabs(jet_bot1_flavourP)==5 && fabs(jet_bot2_flavourP)==5)')

################# channel partitions

df_M = {}
df_E = {}


for s in samples:
        if args.type == "mc":
               if years[0] == "2017":
                     df[s] = df[s].Define('l1_prefw','L1PreFiringWeight_Nom')
               else:
                     df[s] = df[s].Define('l1_prefw','1.')
               df[s] = df[s].Define('lep_id_sf','nMuonGood>0 ? musf_tight_id[MuonGoodInd[0]] : elesf_wp80iso[ElectronGoodInd[0]]')
               df[s] = df[s].Define('lep_iso_sf','nMuonGood>0 ? musf_tight_reliso[MuonGoodInd[0]] : 1.')
               #df_E[s] = df_E[s].Define('lep_id_sf_2','elesf_Tight[ElectronGoodInd[0]]')
               df[s] = df[s].Define('lep_trig_sf','nMuonGood>0 ? trigger_sf_mu_aux[MuonGoodInd[0]] : trigger_sf_el_aux[ElectronGoodInd[0]]')
               #df_M[s] = df_M[s].Define('weightSSOS_final','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*puWeight*PUjetID_SF*lep_trig_sf*top_weight')
               #df_E[s] = df_E[s].Define('weightSSOS_final','weight_aux*lep_id_sf*btag_sf*puWeight*PUjetID_SF*lep_trig_sf*top_weight')
               df[s] = df[s].Define('weightSSOS_final','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw')
               df[s] = df[s].Define('weightSSOS_final_btaglightup','weight_aux*btag_sf_light_up*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw')
               df[s] = df[s].Define('weightSSOS_final_btaglightdown','weight_aux*btag_sf_light_down*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw')
               df[s] = df[s].Define('weightSSOS_final_btagheavyup','weight_aux*btag_sf_heavy_up*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw')
               df[s] = df[s].Define('weightSSOS_final_btagheavydown','weight_aux*btag_sf_heavy_down*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw')
        else:
               df[s] = df[s].Define('weightSSOS_final','weight_aux')

for s in samples:
        df_M[s] = df[s].Filter('nMuonGood>0')
        df_E[s] = df[s].Filter('nElectronGood>0')

############################################################
####################     HISTS    ##########################
############################################################

hist_nom_M = {}
hist_nom_E = {}
hist_sfbtag_light_up_M = {}
hist_sfbtag_light_up_E = {}
hist_sfbtag_light_down_M = {}
hist_sfbtag_light_down_E = {}
hist_sfbtag_heavy_up_M = {}
hist_sfbtag_heavy_up_E = {}
hist_sfbtag_heavy_down_M = {}
hist_sfbtag_heavy_down_E = {}

observable_names = ["nJetGood", "jet_1_pt", "jet_1_nmu", "jet_1_eta", "jet_2_pt", "jet_2_eta", "jet_2_mass", "jet_2_qgl", "jet_2_nmu", "jet_1_qgl",
   "lepton_pt", "lepton_eta", "lepton_pt_detail", "lepton_eta_thick", "InvM_2jets", "InvM_bot_closer", "InvM_bot_farther", 
   "deltaR_jet1_jet2", "deltaphi_jet1_jet2", "deltaeta_jet1_jet2", "MET_pt", "MET_sig", "MET_my_sig",
   "transverse_mass", "tracks_jet1", "tracks_jet2", "EMN_jet1", "EMtotal_jet1", "EMC_jet1", "pT_sum", "pT_product", "deltaR_lep_2jets",
   "deltaphi_MET_2jets", "deltaphi_MET_jets_1", "deltaphi_MET_jets_2", "deltaphi_lephad", "eta_2jets", "pt_2jets", "pT_Wlep", 
   "jet_1_btag", "jet_2_btag", "deltaR_lep_jet1", "deltaPhi_lep_jet1", "deltaEta_lep_jet1", "deltaR_lep_jet2", "deltaPhi_lep_jet2", "deltaEta_lep_jet2",
   "deltaphi_MET_lep", "pT_proy", "pT_sum_2J", "jet_bot1_btag", "jet_bot2_btag", "jet_bot1_pt", "jet_bot1_eta", "jet_bot2_pt", "jet_bot2_eta",
   "jet_bot1_btag_thick", "jet_bot2_btag_thick", "jet_1_btag_thick", "jet_2_btag_thick",
   "jet_bot1_btagnumber", "jet_bot2_btagnumber", "jet_1_btagnumber", "jet_2_btagnumber",
   "jet_1_cvltag_csv", "jet_2_cvltag_csv", "jet_1_cvltag", "jet_2_cvltag", "jet_max_cvltag", "jet_min_cvltag",
   "jet_1_cvbtag_csv", "jet_2_cvbtag_csv", "jet_1_cvbtag", "jet_2_cvbtag", "jet_max_cvbtag", "jet_min_cvbtag"]

observable_names = ["nJetGood", "jet_1_pt", "jet_1_nmu", "jet_1_eta", "jet_2_pt", "jet_2_eta", "jet_2_mass", "jet_2_qgl", "jet_2_nmu", "jet_1_qgl",
   "lepton_pt", "lepton_eta", "lepton_pt_detail", "lepton_eta_thick", "InvM_2jets", "InvM_bot_closer", "InvM_bot_farther",
   "deltaR_jet1_jet2", "deltaphi_jet1_jet2", "deltaeta_jet1_jet2", "MET_pt_aux", "MET_sig", "MET_my_sig",
   "transverse_mass", "tracks_jet1", "tracks_jet2", "deltaphi_MET_2jets", "deltaphi_MET_jets_1", "deltaphi_MET_jets_2", "pT_Wlep",
   "jet_1_btag", "jet_2_btag", "deltaphi_MET_lep", "jet_bot1_btag", "jet_bot2_btag", "jet_bot1_pt", "jet_bot1_eta", "jet_bot2_pt", "jet_bot2_eta",
   "jet_bot1_btag_thick", "jet_bot2_btag_thick", "jet_1_btag_thick", "jet_2_btag_thick",
   "jet_bot1_btagnumber", "jet_bot2_btagnumber", "jet_1_btagnumber", "jet_2_btagnumber",
   "jet_1_cvltag_csv", "jet_2_cvltag_csv", "jet_1_cvltag", "jet_2_cvltag", "InvM30","InvM31","InvMl0","InvMl1", "chi2_test0", "chi2_test1",
   "InvMl_good", "InvMl_bad", "InvM3_good", "InvM3_bad", "chi2_test_good", "chi2_test_bad", "jet_max_cvltag", "jet_min_cvltag",
   "jet_1_cvbtag_csv", "jet_2_cvbtag_csv", "jet_1_cvbtag", "jet_2_cvbtag", "jet_max_cvbtag", "jet_min_cvbtag",
   "jet_1_eta_thick","jet_2_eta_thick","jet_bot1_eta_thick","jet_bot2_eta_thick",
   "InvM_2jets_thick","InvM_2jets_short","bot1_muons","bot2_muons"]

observable_names = observable_names + ["muon_bot1_eta","muon_bot2_eta","muon_bot1_pt","muon_bot2_pt",
   "muon_bot1_relpt","muon_bot2_relpt"]

column_names = {}

for name in observable_names:
   column_names[name] = name

column_names["lepton_pt_detail"] = "lepton_pt"; column_names["jet_1_btag_thick"] = "jet_1_btag"; column_names["jet_2_btag_thick"] = "jet_2_btag";
column_names["jet_bot1_btag_thick"] = "jet_bot1_btag"; column_names["jet_bot2_btag_thick"] = "jet_bot2_btag"; column_names["lepton_eta_thick"] = "lepton_eta";
column_names["MET_sig"] = "MET_significance"; column_names["MET_my_sig"] = "MET_my_significance";
column_names["jet_1_eta_thick"] = "jet_1_eta";column_names["jet_2_eta_thick"] = "jet_2_eta";column_names["jet_bot1_eta_thick"] = "jet_bot1_eta";column_names["jet_bot2_eta_thick"] = "jet_bot2_eta";
column_names["InvM_2jets_thick"] = "InvM_2jets";column_names["InvM_2jets_short"] = "InvM_2jets";

dict_binlim = {}
dict_binlim["nJetGood"] = [7,4,11]; dict_binlim["jet_2_mass"] = [40,0,40];
dict_binlim["jet_1_pt"] = [50,15,115]; dict_binlim["jet_1_eta"] = [60,-3,3]; dict_binlim["jet_1_nmu"] = [10,0,10]; dict_binlim["jet_1_qgl"] = [50,0,1];
dict_binlim["jet_2_pt"] = [50,15,115]; dict_binlim["jet_2_eta"] = [60,-3,3]; dict_binlim["jet_2_nmu"] = [10,0,10]; dict_binlim["jet_2_qgl"] = [50,0,1];
dict_binlim["jet_1_eta_thick"] = [18,-2.7,2.7];dict_binlim["jet_2_eta_thick"] = [18,-2.7,2.7];dict_binlim["jet_bot1_eta_thick"] = [18,-2.7,2.7];dict_binlim["jet_bot2_eta_thick"] = [18,-2.7,2.7];
dict_binlim["lepton_pt"] = [50,20,120]; dict_binlim["lepton_eta"] = [80,-4,4]; dict_binlim["lepton_eta_thick"] = [18,-2.7,2.7]; 
dict_binlim["lepton_pt_detail"] = [80,40,60]; dict_binlim["InvM_2jets"] = [108,30,300]; dict_binlim["InvM_bot_closer"] = [100,0,300]; 
dict_binlim["InvM_bot_farther"] = [100,0,300]; dict_binlim["MET_pt_aux"] = [60,10,130]; dict_binlim["MET_sig"] = [50,0,18]; dict_binlim["MET_my_sig"] = [50,0,18]; 
dict_binlim["deltaR_jet1_jet2"] = [40,0,4]; dict_binlim["deltaphi_jet1_jet2"] = [100,0,5]; dict_binlim["deltaeta_jet1_jet2"] = [100,0,5]; 
dict_binlim["transverse_mass"] = [35,0,140]; dict_binlim["tracks_jet1"] = [60,0,60]; dict_binlim["tracks_jet2"] = [60,0,60]; 
dict_binlim["EMN_jet1"] = [60,0,1]; dict_binlim["EMC_jet1"] = [60,0,1]; dict_binlim["EMtotal_jet1"] = [60,0,1]; 
dict_binlim["pT_sum"] = [100,0,1]; dict_binlim["pT_product"] = [100,-1,1]; dict_binlim["deltaR_lep_2jets"] = [100,0,5]; dict_binlim["deltaphi_MET_2jets"] = [100,0,5];
dict_binlim["deltaphi_MET_lep"] = [100,0,5]; dict_binlim["deltaphi_MET_jets_1"] = [100,0,5]; dict_binlim["deltaphi_MET_jets_2"] = [100,0,5];
dict_binlim["deltaphi_lephad"] = [100,0,5]; dict_binlim["eta_2jets"] = [50,-5,5]; dict_binlim["pt_2jets"] = [100,0,200]; dict_binlim["pT_Wlep"] = [50,0,200]; 
dict_binlim["deltaR_lep_jet1"] = [100,0,5]; dict_binlim["deltaR_lep_jet2"] = [100,0,5]; dict_binlim["deltaPhi_lep_jet1"] = [100,0,5]; 
dict_binlim["deltaPhi_lep_jet2"] = [100,0,5]; dict_binlim["deltaEta_lep_jet1"] = [100,0,5]; dict_binlim["deltaEta_lep_jet2"] = [100,0,5]; 
dict_binlim["jet_bot1_btag"] = [50,0,1]; dict_binlim["jet_bot1_pt"] = [50,15,115]; dict_binlim["jet_bot1_eta"] = [60,-3,3]; dict_binlim["jet_bot1_btag_thick"] = [10,0,1];
dict_binlim["jet_bot2_btag"] = [50,0,1]; dict_binlim["jet_bot2_pt"] = [50,15,115]; dict_binlim["jet_bot2_eta"] = [60,-3,3]; dict_binlim["jet_bot2_btag_thick"] = [10,0,1]; 
dict_binlim["jet_1_btag"] = [50,0,1]; dict_binlim["jet_1_btag_thick"] = [10,0,1]; dict_binlim["jet_2_btag"] = [50,0,1]; dict_binlim["jet_2_btag_thick"] = [10,0,1]; 
dict_binlim["jet_bot1_btagnumber"] = [3,0,3]; dict_binlim["jet_bot2_btagnumber"] = [3,0,3]; dict_binlim["jet_1_btagnumber"] = [3,0,3]; 
dict_binlim["jet_2_btagnumber"] = [3,0,3]; dict_binlim["pT_proy"] = [100,-100,100]; dict_binlim["pT_sum_2J"] = [100,0,1]; 
dict_binlim["jet_1_flavourP"] = [28,-6,22]; dict_binlim["jet_bot1_flavourP"] = [28,-6,22];
dict_binlim["jet_2_flavourP"] = [28,-6,22]; dict_binlim["jet_bot2_flavourP"] = [28,-6,22];
dict_binlim["btag_sf"] = [100,0.5,1.5]; dict_binlim["lep_id_sf"] = [100,0.5,1.5]; dict_binlim["lep_trig_sf"] = [100,0.5,1.5]; 
dict_binlim["lep_iso_sf"] = [100,0.5,1.5]; dict_binlim["puWeight"] = [100,0.5,1.5]; dict_binlim["PUjetID_SF"] = [100,0.5,1.5]; 
dict_binlim["top_weight"] = [100,0.5,1.5]; dict_binlim["jet_1_cvltag_csv"]=[50,0,1]; dict_binlim["jet_2_cvltag_csv"]=[50,0,1];
dict_binlim["jet_1_cvltag"]=[50,0,1]; dict_binlim["jet_2_cvltag"]=[50,0,1];
dict_binlim["InvM30"]=[50,50,300];dict_binlim["InvM31"]=[50,50,300];dict_binlim["InvMl0"]=[50,50,300];dict_binlim["InvMl1"]=[50,50,300];
dict_binlim["chi2_test0"]=[50,0,10];dict_binlim["chi2_test1"]=[50,0,10];dict_binlim["chi2_test_good"]=[20,0,4];dict_binlim["chi2_test_bad"]=[50,0,10];
dict_binlim["InvM3_good"]=[36,120,210];dict_binlim["InvM3_bad"]=[50,50,300];dict_binlim["InvMl_good"]=[40,40,160];dict_binlim["InvMl_bad"]=[50,50,300];
dict_binlim["jet_max_cvltag"]=[50,0,1]; dict_binlim["jet_min_cvltag"]=[50,0,1];dict_binlim["jet_max_cvbtag"]=[50,0,1]; dict_binlim["jet_min_cvbtag"]=[50,0,1];
dict_binlim["jet_1_cvbtag_csv"]=[50,0,1]; dict_binlim["jet_2_cvbtag_csv"]=[50,0,1];dict_binlim["jet_1_cvbtag"]=[50,0,1]; dict_binlim["jet_2_cvbtag"]=[50,0,1];
dict_binlim["deltaR_jet1_muon"]=[50,0,5];dict_binlim["deltaR_jet1_tau"]=[50,0,5];dict_binlim["deltaR_muon_tau"]=[50,0,5];
dict_binlim["Frag_weight_sl"]=[40,0.4,2.4];dict_binlim["Br_weight_sl"]=[40,0.4,2.4];
dict_binlim["InvM_2jets_thick"] = [54,30,300];dict_binlim["InvM_2jets_short"] = [30,50,110];
dict_binlim["bot1_muons"]=[10,0,10];dict_binlim["bot2_muons"]=[10,0,10];dict_binlim["muon_bot1_eta"]=[18,-2.7,2.7];dict_binlim["muon_bot2_eta"]=[18,-2.7,2.7];
dict_binlim["muon_bot1_pt"]=[30,0,30];dict_binlim["muon_bot2_pt"]=[30,0,30];dict_binlim["muon_bot1_relpt"]=[20,0,1];dict_binlim["muon_bot2_relpt"]=[20,0,1];
dict_binlim["muon_bot1_mother"]=[500,0,500];dict_binlim["muon_bot2_mother"]=[500,0,500];

dict_binlim_M = dict(dict_binlim); dict_binlim_E = dict(dict_binlim);
dict_binlim_M["lepton_pt"] = [50,20,120];dict_binlim_E["lepton_pt"] = [50,25,125];

for name in observable_names:
        hist_nom_M[name] = {}
        hist_nom_E[name] = {}
        hist_sfbtag_light_up_M[name] = {}
        hist_sfbtag_light_up_E[name] = {}
        hist_sfbtag_light_down_M[name] = {}
        hist_sfbtag_light_down_E[name] = {}
        hist_sfbtag_heavy_up_M[name] = {}
        hist_sfbtag_heavy_up_E[name] = {}
        hist_sfbtag_heavy_down_M[name] = {}
        hist_sfbtag_heavy_down_E[name] = {}
        if mode == "mc":
           if args.syst == "down":
              for s in samples:
                 hist_nom_M[name][s] = df_M[s].Histo1D((s+"_"+name+"_M_ptscaledown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final")
       	         hist_nom_E[name][s] = df_E[s].Histo1D((s+"_"+name+"_E_ptscaledown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final")
           if args.syst == "nom":           
              for s in samples:
                 hist_nom_M[name][s] = df_M[s].Histo1D((s+"_"+name+"_M","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final")
       	         hist_nom_E[name][s] = df_E[s].Histo1D((s+"_"+name+"_E","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final")
                 #if name == "InvM_2jets" and (chan == "btagMM_chitest" or chan == "lepton50_chitest"):
                 #   hist_nom_M[name][s].Sumw2()
       	         #   hist_nom_E[name][s].Sumw2()
                 hist_sfbtag_light_up_M[name][s] = df_M[s].Histo1D((s+"_"+name+"_M_btaglightup","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_btaglightup")
       	         hist_sfbtag_light_up_E[name][s] = df_E[s].Histo1D((s+"_"+name+"_E_btaglightup","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_btaglightup")
                 hist_sfbtag_light_down_M[name][s] = df_M[s].Histo1D((s+"_"+name+"_M_btaglightdown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_btaglightdown")
       	         hist_sfbtag_light_down_E[name][s] = df_E[s].Histo1D((s+"_"+name+"_E_btaglightdown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_btaglightdown")
                 hist_sfbtag_heavy_up_M[name][s] = df_M[s].Histo1D((s+"_"+name+"_M_btagheavyup","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_btagheavyup")
       	         hist_sfbtag_heavy_up_E[name][s] = df_E[s].Histo1D((s+"_"+name+"_E_btagheavyup","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_btagheavyup")
                 hist_sfbtag_heavy_down_M[name][s] = df_M[s].Histo1D((s+"_"+name+"_M_btagheavydown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_btagheavydown")
       	         hist_sfbtag_heavy_down_E[name][s] = df_E[s].Histo1D((s+"_"+name+"_E_btagheavydown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_btagheavydown")
        else:
           for s in samples:
              hist_nom_M[name][s] = df_M[s].Histo1D(("data"+s+"_"+name+"_M","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final")
       	      hist_nom_E[name][s] = df_E[s].Histo1D(("data"+s+"_"+name+"_E","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final")

########## gen hists

hist_mc_M = {}
hist_mc_E = {}

observable_mc_names = ["jet_1_flavourP", "jet_2_flavourP", "jet_bot1_flavourP", "jet_bot2_flavourP", "btag_sf", "lep_id_sf", "lep_trig_sf", "lep_iso_sf",
     "puWeight", "PUjetID_SF", "top_weight","Frag_weight_sl","Br_weight_sl"]

if mode == "mc":
    for name in observable_mc_names:
        hist_mc_M[name] = {}
        hist_mc_E[name] = {}
        for s in samples:
              hist_mc_M[name][s] = df_M[s].Histo1D((s+"_"+name+"_M","",dict_binlim[name][0],dict_binlim[name][1],dict_binlim[name][2]),name,"weightSSOS_final")
       	      hist_mc_E[name][s] = df_E[s].Histo1D((s+"_"+name+"_E","",dict_binlim[name][0],dict_binlim[name][1],dict_binlim[name][2]),name,"weightSSOS_final")

#############################
####     DATA SAVING     ####
#############################

term1 = "botjets_muons/"

#observable_names = ["transverse_mass","MET_pt_aux"]

if mode == "mc":
   for name in observable_names:
        path_hist = '/nfs/cms/vazqueze/new_hists/fromJF/wqq/'+term1+'hist_wqqfromJF_MC_'+years[0]+'_'+name+'.root'

        if args.syst == "down":
           myfile = TFile( path_hist, 'UPDATE' )
           for s in samples:
              hist_nom_M[name][s].Write()
              hist_nom_E[name][s].Write()

           myfile.Close()
        elif args.syst == "nom":
           myfile = TFile( path_hist, 'RECREATE' )
           for s in samples:
              hist_nom_M[name][s].Write()
              hist_nom_E[name][s].Write()
              hist_sfbtag_light_up_M[name][s].Write()
              hist_sfbtag_light_up_E[name][s].Write()
              hist_sfbtag_light_down_M[name][s].Write()
              hist_sfbtag_light_down_E[name][s].Write()
              hist_sfbtag_heavy_up_M[name][s].Write()
              hist_sfbtag_heavy_up_E[name][s].Write()
              hist_sfbtag_heavy_down_M[name][s].Write()
              hist_sfbtag_heavy_down_E[name][s].Write()

           myfile.Close()

   for name in observable_mc_names:
        path_hist = '/nfs/cms/vazqueze/new_hists/fromJF/wqq/'+term1+'hist_wqqfromJF_MC_'+years[0]+'_'+name+'.root'
        myfile = TFile( path_hist, 'RECREATE' )
        for s in samples:
          hist_mc_M[name][s].Write()
          hist_mc_E[name][s].Write()
        myfile.Close()

if (mode == "data" and samples[0][-1] == "M"):
   for name in observable_names:
        path_hist = '/nfs/cms/vazqueze/new_hists/fromJF/wqq/'+term1+'hist_wqqfromJF_dataM_'+years[0]+'_'+name+'.root'
        myfile = TFile( path_hist, 'RECREATE' )

        hist_nom_M[name][samples[0]].Write()

        myfile.Close()

if (mode == "data" and samples[0][-1] == "E"):
   for name in observable_names:
        path_hist = '/nfs/cms/vazqueze/new_hists/fromJF/wqq/'+term1+'hist_wqqfromJF_dataE_'+years[0]+'_'+name+'.root'
        myfile = TFile( path_hist, 'RECREATE' )

        hist_nom_E[name][samples[0]].Write()

        myfile.Close()

print('Ended succesfully')
