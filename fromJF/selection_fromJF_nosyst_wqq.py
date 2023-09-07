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
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Perform ssos substraction")
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

files_data = {"M":{"2016":230,"2016B":157,"2017":436,"2018":393},
   "E":{"2016":335,"2016B":156,"2017":245,"2018":738}}

term1 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_ver2/"
term2 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_down/"
term3 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID/"
#term3 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_nosmearing/"
for p in proc:
    # Construct the dataframes
    if args.syst == "nom": folder = term1+"myconfig"+years[0]+"/{sample}/cat_base/prod_test/" # Folder name
    if args.syst == "down": folder = term2+"myconfig"+years[0]+"/{sample}/cat_base/prod_test/" # Folder name
    folder = term3 + "myconfig"+years[0]+"/{sample}/cat_base/prod_test/" # Folder name
    if args.type == "mc":
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

## Cuts per year

cuts_btag = {}

cuts_btag["2016"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2016B"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2017"]=[0.0532, 0.3040, 0.7476]
cuts_btag["2018"]=[0.0490,0.2783,0.7100]

cuts_btag["2016"]=[0.0614, 0.3093, 0.7221]
cuts_btag["2016B"]=[0.0614, 0.3093, 0.7221]

## Muon and Electron pT cuts per year

muon_pt = {}

muon_pt["2016"]=26
muon_pt["2016B"]=26
muon_pt["2017"]=29
muon_pt["2018"]=26

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
      auto InvMassBot(Vint jetbot, Vint jetmuon, Vint jetgood, Vfloat Jet_pt, Vfloat Jet_eta, Vfloat Jet_phi, UInt_t jet_notmuon) {
            int ind = -1;
            vector<float> vb;
            ROOT::Math::PtEtaPhiMVector jet1(Jet_pt[jetmuon[0]], Jet_eta[jetmuon[0]], Jet_phi[jetmuon[0]], 0.);
            ROOT::Math::PtEtaPhiMVector jet2(Jet_pt[jet_notmuon], Jet_eta[jet_notmuon], Jet_phi[jet_notmuon], 0.);
            auto had = jet1+jet2;
            float had_eta = had.Eta();
       	    float had_phi = had.Phi();
            float dR1 = ROOT::VecOps::DeltaR(had_eta,Jet_eta[jetgood[jetbot[0]]],had_phi,Jet_phi[jetgood[jetbot[0]]]);
            float dR2 = ROOT::VecOps::DeltaR(had_eta,Jet_eta[jetgood[jetbot[1]]],had_phi,Jet_phi[jetgood[jetbot[1]]]);
            if (dR1 > dR2){
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon]));           
            }else{
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[0]]],Jet_eta[jetgood[jetbot[0]]],Jet_phi[jetgood[jetbot[0]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon]));
                   vb.push_back(InvariantM3(Jet_pt[jetgood[jetbot[1]]],Jet_eta[jetgood[jetbot[1]]],Jet_phi[jetgood[jetbot[1]]],Jet_pt[jetmuon[0]],Jet_eta[jetmuon[0]],Jet_phi[jetmuon[0]],Jet_pt[jet_notmuon],Jet_eta[jet_notmuon],Jet_phi[jet_notmuon]));
            }
            return vb;
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
                  deltaPhilepjet1 = mydeltaphi(mu_phi[mu_good[0]],jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = mydeltaphi(mu_phi[mu_good[0]],jet_phi[jet_notmuon]);
       	       	  deltaEtalepjet1 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_notmuon]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
       	       	  deltaRlepjet1	= ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_muon[0]],el_phi[el_good[0]],jet_phi[jet_muon[0]]);
                  deltaRlepjet2 = ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_notmuon],el_phi[el_good[0]],jet_phi[jet_notmuon]);
                  deltaPhilepjet1 = mydeltaphi(el_phi[el_good[0]],jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = mydeltaphi(el_phi[el_good[0]],jet_phi[jet_notmuon]);
                  deltaEtalepjet1 = fabs(el_eta[el_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(el_eta[el_good[0]]-jet_eta[jet_notmuon]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
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
            vb.push_back(mydeltaphi(plep2.Phi(),jet_phi[jet_notmuon]));
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
           df[s] = df[s].Define('Jet_pt_aux','Jet_pt_nom')
           #df[s] = df[s].Define('Jet_pt_aux','Jet_pt')
           df[s] = df[s].Define('MET_pt_aux','MET_smeared_pt')
           #df[s] = df[s].Define('MET_pt_aux','MET_pt')
           df[s] = df[s].Define('MET_phi_aux','MET_smeared_phi')
           #df[s] = df[s].Define('MET_phi_aux','MET_phi')
   elif args.syst == "down":
        for s in samples:
           df[s] = df[s].Define('Jet_pt_aux','Jet_pt_down')
   elif args.syst == "up":
        for s in samples:
           df[s] = df[s].Define('Jet_pt_aux','Jet_pt_up')
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
        df[s] = df[s].Define('JetQInd','JetInds(nJet, JetGoodInd, Jet_pt_aux, Jet_eta, Jet_phi, JetBotInd, Jet_qgl)')
        ########### Filtering and further definitions
        df[s] = df[s].Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>=4')
        df[s] = df[s].Filter('JetQInd.size() > 1')
        df[s] = df[s].Define('JetnotMuonInd','JetQInd[1]')
        ### hists definitions
        df[s] = df[s].Define('jet_1_pt','Jet_pt_aux[JetQInd[0]]')
        df[s] = df[s].Define('jet_1_nmu','Jet_nMuons[JetQInd[0]]')
        df[s] = df[s].Define('jet_1_qgl','Jet_qgl[JetQInd[0]]')
        df[s] = df[s].Define('jet_2_pt','nJetGood>1 ? Jet_pt_aux[JetnotMuonInd] : 0')
        df[s] = df[s].Define('jet_1_eta','Jet_eta[JetQInd[0]]')
        df[s] = df[s].Define('jet_2_eta','nJetGood>1 ? Jet_eta[JetnotMuonInd] : 0')
        df[s] = df[s].Define('jet_1_phi','Jet_phi[JetQInd[0]]')
        df[s] = df[s].Define('jet_2_phi','nJetGood>1 ? Jet_phi[JetnotMuonInd] : 0')
        df[s] = df[s].Define('jet_1_mass','Jet_mass_nom[JetQInd[0]]')
        df[s] = df[s].Define('jet_2_mass','nJetGood>1 ? Jet_mass_nom[JetnotMuonInd] : 0')
        df[s] = df[s].Define('jet_2_qgl','nJetGood>1 ? Jet_qgl[JetnotMuonInd] : 0')
        df[s] = df[s].Define('jet_2_nmu','nJetGood>1 ? Jet_nMuons[JetnotMuonInd] : 0')
        df[s] = df[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt_aux[JetQInd[0]],Jet_eta[JetQInd[0]],Jet_phi[JetQInd[0]],0.,Jet_pt_aux[JetnotMuonInd],Jet_eta[JetnotMuonInd],Jet_phi[JetnotMuonInd],0.) : 0')
        df[s] = df[s].Define('deltaR_jet1_jet2','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotMuonInd], Jet_eta[JetQInd[0]] , Jet_phi[JetnotMuonInd], Jet_phi[JetQInd[0]])  : 10')
        df[s] = df[s].Define('deltaphi_jet1_jet2','mydeltaphi(Jet_phi[JetQInd[0]],Jet_phi[JetnotMuonInd])')
        df[s] = df[s].Define('deltaeta_jet1_jet2','fabs(Jet_eta[JetQInd[0]]-Jet_eta[JetnotMuonInd])')
        df[s] = df[s].Define('deltapt_jet1_jet2','fabs(Jet_pt_aux[JetQInd[0]]-Jet_pt_aux[JetnotMuonInd])')
        df[s] = df[s].Define('tracks_jet1','Jet_nConstituents[JetQInd[0]]')
        df[s] = df[s].Define('tracks_jet2','nJetGood>1 ? Jet_nConstituents[JetnotMuonInd] : 0')
        df[s] = df[s].Define('EMN_jet1','Jet_neEmEF[JetQInd[0]]')
        df[s] = df[s].Define('EMC_jet1','Jet_chEmEF[JetQInd[0]]')
        df[s] = df[s].Define('EMtotal_jet1','Jet_chEmEF[JetQInd[0]]+Jet_neEmEF[JetQInd[0]]')
        df[s] = df[s].Define('pT_sum','pTsum(MuonGoodInd, ElectronGoodInd, JetQInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt_aux, MET_phi_aux, Jet_pt_aux, Jet_eta, Jet_phi, Jet_mass_nom)')
        df[s] = df[s].Define('pT_product','pTprod(MuonGoodInd, ElectronGoodInd, JetQInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt_aux, MET_phi_aux, Jet_pt_aux, Jet_eta, Jet_phi, Jet_mass_nom)')
        df[s] = df[s].Define('aux_various','variousSUM(MuonGoodInd, ElectronGoodInd, JetQInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt_aux, MET_phi_aux, Jet_pt_aux, Jet_eta, Jet_phi, Jet_mass_nom)')
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
        df[s] = df[s].Define('jet_1_btag','Jet_btagDeepFlavB[JetQInd[0]]')
        df[s] = df[s].Define('jet_2_btag','Jet_btagDeepFlavB[JetnotMuonInd]')
        df[s] = df[s].Define('nLooseLepton','nMuon+nElectron-1')
        df[s] = df[s].Define('InvM3_aux','InvMassBot(JetBotInd, JetQInd, JetGoodInd, Jet_pt_aux, Jet_eta, Jet_phi, JetnotMuonInd)')
        df[s] = df[s].Define('InvM_bot_closer','InvM3_aux[0]')
        df[s] = df[s].Define('InvM_bot_farther','InvM3_aux[1]')
        df[s] = df[s].Define('jet_bot1_btagnumber','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[0]]] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetGoodInd[JetBotInd[0]]] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')
        df[s] = df[s].Define('jet_bot2_btagnumber','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[1]]] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetGoodInd[JetBotInd[1]]] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')
        df[s] = df[s].Define('jet_1_btagnumber','Jet_btagDeepFlavB[JetQInd[0]] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetQInd[0]] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')
        df[s] = df[s].Define('jet_2_btagnumber','Jet_btagDeepFlavB[JetnotMuonInd] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetnotMuonInd] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')

############################################################
#### Distinguishing between same sign and opposite sign ####
############################################################

#################### Final definitions and filters ###############################

for s in samples:
        df[s] = df[s].Define('lepton_phi','nMuonGood>0 ? Muon_phi[MuonGoodInd[0]] : Electron_phi[ElectronGoodInd[0]]')
        df[s] = df[s].Define('transverse_mass_old','std::sqrt(2*lepton_pt*MET_pt_aux*(1-std::cos(lepton_phi-MET_phi_aux)))')
        df[s] = df[s].Define('weight_aux','1.')
        df[s] = df[s].Define('deltaR_jetM_lep','nMuonGood>0 ? ROOT::VecOps::DeltaR(Muon_eta[MuonGoodInd[0]],Jet_eta[JetQInd[0]] , Muon_phi[MuonGoodInd[0]], Jet_phi[JetQInd[0]]) : ROOT::VecOps::DeltaR(Electron_eta[ElectronGoodInd[0]],Jet_eta[JetQInd[0]] , Electron_phi[ElectronGoodInd[0]], Jet_phi[JetQInd[0]])')
        df[s] = df[s].Define('InvM_jetM_lep', 'nMuonGood>0 ? InvariantM(Jet_pt_aux[JetQInd[0]],Jet_eta[JetQInd[0]],Jet_phi[JetQInd[0]],0.,Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]]) : InvariantM(Jet_pt_aux[JetQInd[0]],Jet_eta[JetQInd[0]],Jet_phi[JetQInd[0]],0.,Electron_pt[ElectronGoodInd[0]],Electron_eta[ElectronGoodInd[0]],Electron_phi[ElectronGoodInd[0]], Electron_mass[ElectronGoodInd[0]])')
        df[s] = df[s].Define('InvM_muon_jet','nMuonGood>0 ? InvariantM(Muon_pt[MuonJetInd[0]],Muon_eta[MuonJetInd[0]],Muon_phi[MuonJetInd[0]],Muon_mass[MuonJetInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]]) : 50.')
        df[s] = df[s].Define('deltaphi_MET_lep','mydeltaphi(lepton_phi,MET_phi_aux)')
        df[s] = df[s].Define('MET_my_significance','(MET_pt_aux*MET_pt_aux)/((0.62*0.62)*MET_sumEt)')
        ### Filters
        df[s] = df[s].Filter('jet_bot1_btag >'+str(cuts_btag[years[0]][1]))
        df[s] = df[s].Filter('jet_bot2_btag >'+str(cuts_btag[years[0]][0]))
        #df[s] = df[s].Filter('InvM_muon_jet >12').Filter('InvM_muon_jet > 110 || InvM_muon_jet < 70')
        df[s] = df[s].Filter('transverse_mass > 50.')
        df[s] = df[s].Filter('lepton_pt > 50.')
        ### Jet pT filters
        df[s] = df[s].Filter('jet_bot1_pt > 30.')
        df[s] = df[s].Filter('jet_bot2_pt > 30.')
        df[s] = df[s].Filter('jet_1_pt > 30.')
        df[s] = df[s].Filter('jet_2_pt > 30.')
        ## Invariant mass filter
        df[s] = df[s].Filter('InvM_2jets > 30. && InvM_2jets < 300.')
        ### Cur based for electron, EL ID tests
        #df[s] = df[s].Filter('nElectronGood > 0 ? Electron_cutBased[ElectronGoodInd[0]] > 3 : 1')

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

############ Gen level definitions

if mode == "mc":
        for s in samples:
                df[s] = df[s].Define('jet_1_flavourH','Jet_hadronFlavour[JetQInd[0]]')
                df[s] = df[s].Define('jet_2_flavourH','Jet_hadronFlavour[JetnotMuonInd]')
                df[s] = df[s].Define('jet_1_flavourP','Jet_partonFlavour[JetQInd[0]]')
                df[s] = df[s].Define('jet_2_flavourP','Jet_partonFlavour[JetnotMuonInd]')
                df[s] = df[s].Define('jet_bot1_flavourP','Jet_partonFlavour[JetGoodInd[JetBotInd[0]]]')
                df[s] = df[s].Define('jet_bot2_flavourP','Jet_partonFlavour[JetGoodInd[JetBotInd[1]]]')
                df[s] = df[s].Define('last_Copy','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,13)')
                df[s] = df[s].Define('top_weight','topreweight(nGenPart, GenPart_statusFlags, GenPart_pdgId, last_Copy, GenPart_pt)')
                df[s] = df[s].Define('btag_sf','Aux_btag_weight[1]/Aux_btag_weight[0]')

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
        df_M[s] = df[s].Filter('nMuonGood>0')
        df_E[s] = df[s].Filter('nElectronGood>0')

for s in samples:
        if args.type == "mc":
               if years[0] == "2017":
                     df_M[s] = df_M[s].Define('l1_prefw','L1PreFiringWeight_Nom')
                     df_E[s] = df_E[s].Define('l1_prefw','L1PreFiringWeight_Nom')
               else:
                     df_M[s] = df_M[s].Define('l1_prefw','1.')
                     df_E[s] = df_E[s].Define('l1_prefw','1.')
               df_M[s] = df_M[s].Define('lep_id_sf','musf_tight_id[MuonGoodInd[0]]')
               df_M[s] = df_M[s].Define('lep_iso_sf','musf_tight_reliso[MuonGoodInd[0]]')
               df_E[s] = df_E[s].Define('lep_id_sf','elesf_wp80iso[ElectronGoodInd[0]]')
               df_E[s] = df_E[s].Define('lep_iso_sf','1.')
               #df_E[s] = df_E[s].Define('lep_id_sf_2','elesf_Tight[ElectronGoodInd[0]]')
               df_M[s] = df_M[s].Define('lep_trig_sf','trigger_sf_mu_aux[MuonGoodInd[0]]')
               df_E[s] = df_E[s].Define('lep_trig_sf','trigger_sf_el_aux[ElectronGoodInd[0]]')
               df_M[s] = df_M[s].Define('lep_id_lowpt_sf','musf_tight_id_lowpt[MuonJetInd[0]]')
               df_E[s] = df_E[s].Define('lep_id_lowpt_sf','musf_tight_id_lowpt[MuonJetInd[0]]')
               #df_M[s] = df_M[s].Define('weightSSOS_final','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*puWeight*PUjetID_SF*lep_trig_sf*top_weight')
               #df_E[s] = df_E[s].Define('weightSSOS_final','weight_aux*lep_id_sf*btag_sf*puWeight*PUjetID_SF*lep_trig_sf*top_weight')
               df_M[s] = df_M[s].Define('weightSSOS_final','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw')
               df_E[s] = df_E[s].Define('weightSSOS_final','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw')
        else:
               df_M[s] = df_M[s].Define('weightSSOS_final','weight_aux')
               df_E[s] = df_E[s].Define('weightSSOS_final','weight_aux')

############################################################
####################     HISTS    ##########################
############################################################

hist_nJetGood_M = {}
hist_nJetGood_E = {}
hist_jet_1_pt_M = {}
hist_jet_2_pt_M = {}
hist_jet_1_eta_M = {}
hist_jet_2_eta_M = {}
hist_jet_2_qgl_M = {}
hist_jet_1_qgl_M = {}
hist_jet_2_mass_M = {}
hist_jet_1_nmu_M = {}
hist_jet_2_nmu_M = {}
hist_jet_1_pt_E = {}
hist_jet_2_pt_E = {}
hist_jet_1_eta_E = {}
hist_jet_2_eta_E = {}
hist_jet_2_qgl_E = {}
hist_jet_1_qgl_E = {}
hist_jet_2_mass_E = {}
hist_jet_1_nmu_E = {}
hist_jet_2_nmu_E = {}
hist_lepton_pt_M = {}
hist_lepton_pt_detail_M = {}
hist_lepton_eta_M = {}
hist_lepton_eta_thick_M = {}
hist_lepton_pt_E = {}
hist_lepton_pt_detail_E = {}
hist_lepton_eta_E = {}
hist_lepton_eta_thick_E = {}
hist_InvM_2jets_M = {}
hist_InvM_2jets_E = {}
hist_InvM_bot_closer_M = {}
hist_InvM_bot_closer_E = {}
hist_InvM_bot_farther_M = {}
hist_InvM_bot_farther_E = {}
hist_deltaR_jet1_jet2_M = {}
hist_deltaR_jet1_jet2_E = {}
hist_deltaphi_jet1_jet2_M = {}
hist_deltaphi_jet1_jet2_E = {}
hist_deltaeta_jet1_jet2_M = {}
hist_deltaeta_jet1_jet2_E = {}
hist_deltapt_jet1_jet2_M = {}
hist_deltapt_jet1_jet2_E = {}
hist_MET_M = {}
hist_MET_E = {}
hist_MET_sig_M = {}
hist_MET_sig_E = {}
hist_MET_my_sig_M = {}
hist_MET_my_sig_E = {}
hist_tranverse_mass_M = {}
hist_tranverse_mass_E = {}
hist_tracks_jet1_M = {}
hist_tracks_jet2_M = {}
hist_tracks_jet1_E = {}
hist_tracks_jet2_E = {}
hist_EMN_jet1_M = {}
hist_EMC_jet1_M = {}
hist_EMN_jet1_E = {}
hist_EMC_jet1_E = {}
hist_EMtotal_jet1_M = {}
hist_EMtotal_jet1_E = {}
hist_pT_sum_M = {}
hist_pT_sum_E = {}
hist_pT_product_M = {}
hist_pT_product_E = {}
hist_deltaR_lep_2jets_M = {}
hist_deltaR_lep_2jets_E = {}
hist_deltaphi_MET_2jets_M = {}
hist_deltaphi_MET_2jets_E = {}
hist_deltaphi_MET_jets_1_M = {}
hist_deltaphi_MET_jets_1_E = {}
hist_deltaphi_MET_jets_2_M = {}
hist_deltaphi_MET_jets_2_E = {}
hist_deltaphi_lephad_M = {}
hist_deltaphi_lephad_E = {}
hist_eta_2jets_M = {}
hist_eta_2jets_E = {}
hist_pt_2jets_M = {}
hist_pt_2jets_E = {}
hist_pt_Wlep_M = {}
hist_pt_Wlep_E = {}
hist_jet_1_btag_M = {}
hist_jet_1_btag_E = {}
hist_jet_2_btag_M = {}
hist_jet_2_btag_E = {}
hist_jet_1_flavourP_M = {}
hist_jet_1_flavourP_E = {}
hist_jet_2_flavourP_M = {}
hist_jet_2_flavourP_E = {}
hist_jet_bot1_flavourP_M = {}
hist_jet_bot1_flavourP_E = {}
hist_jet_bot2_flavourP_M = {}
hist_jet_bot2_flavourP_E = {}
hist_deltaR_lep_jet1_M = {}
hist_deltaR_lep_jet2_M = {}
hist_deltaPhi_lep_jet1_M = {}
hist_deltaPhi_lep_jet2_M = {}
hist_deltaEta_lep_jet1_M = {}
hist_deltaEta_lep_jet2_M = {}
hist_deltaR_lep_jet1_E = {}
hist_deltaR_lep_jet2_E = {}
hist_deltaPhi_lep_jet1_E = {}
hist_deltaPhi_lep_jet2_E = {}
hist_deltaEta_lep_jet1_E = {}
hist_deltaEta_lep_jet2_E = {}
hist_deltaPhi_MET_lep_M = {}
hist_deltaPhi_MET_lep_E = {}
hist_pT_proy_M = {}
hist_pT_proy_E = {}
hist_pT_sum_2J_M = {}
hist_pT_sum_2J_E = {}
hist_jet_bot1_btag_M = {}
hist_jet_bot1_btag_E = {}
hist_jet_bot2_btag_M = {}
hist_jet_bot2_btag_E = {}
hist_jet_bot1_pt_M = {}
hist_jet_bot1_pt_E = {}
hist_jet_bot2_pt_M = {}
hist_jet_bot2_pt_E = {}
hist_jet_bot1_eta_M = {}
hist_jet_bot1_eta_E = {}
hist_jet_bot2_eta_M = {}
hist_jet_bot2_eta_E = {}
hist_jet_bot1_btag_thick_M = {}
hist_jet_bot1_btag_thick_E = {}
hist_jet_bot2_btag_thick_M = {}
hist_jet_bot2_btag_thick_E = {}
hist_jet_1_btag_thick_M = {}
hist_jet_1_btag_thick_E = {}
hist_jet_2_btag_thick_M = {}
hist_jet_2_btag_thick_E = {}
hist_jet_bot1_btagnumber_M = {}
hist_jet_bot1_btagnumber_E = {}
hist_jet_bot2_btagnumber_M = {}
hist_jet_bot2_btagnumber_E = {}
hist_jet_1_btagnumber_M = {}
hist_jet_1_btagnumber_E = {}
hist_jet_2_btagnumber_M = {}
hist_jet_2_btagnumber_E = {}
## just for mc
hist_btag_sf_M = {}
hist_btag_sf_E = {}
hist_lep_id_sf_M = {}
hist_lep_id_sf_E = {}
hist_lep_iso_sf_M = {}
hist_puWeight_M = {}
hist_puWeight_E = {}
hist_PUjetID_SF_M = {}
hist_PUjetID_SF_E = {}
hist_lep_trig_sf_M = {}
hist_lep_trig_sf_E = {}
hist_topweight_M = {}
hist_topweight_E = {}
### 2D hists
hist_lepton_pt_eta_M = {}
hist_lepton_pt_eta_E = {}

for s in samples:

        hist_nJetGood_M[s] = df_M[s].Histo1D(("nJetGood_M","",10,0,10),"nJetGood","weightSSOS_final")
        hist_nJetGood_E[s] = df_E[s].Histo1D(("nJetGood_E","",10,0,10),"nJetGood","weightSSOS_final")

        hist_jet_1_pt_M[s] = df_M[s].Histo1D(("jet_1_pt_M","",50,20,120),"jet_1_pt","weightSSOS_final")
        hist_jet_1_nmu_M[s] = df_M[s].Histo1D(("jet_1_nmu_M","",10,0,10),"jet_1_nmu","weightSSOS_final")
        hist_jet_2_pt_M[s] = df_M[s].Histo1D(("jet_2_pt_M","",50,20,120),"jet_2_pt","weightSSOS_final")
        hist_jet_1_eta_M[s] = df_M[s].Histo1D(("jet_1_eta_M","",80,-4,4),"jet_1_eta","weightSSOS_final")
        hist_jet_2_eta_M[s] = df_M[s].Histo1D(("jet_2_eta_M","",80,-4,4),"jet_2_eta","weightSSOS_final")
        hist_jet_2_nmu_M[s] = df_M[s].Histo1D(("jet_2_nmu_M","",10,0,10),"jet_2_nmu","weightSSOS_final")
        hist_jet_2_mass_M[s] = df_M[s].Histo1D(("jet_2_mass_M","",40,0,40),"jet_2_mass","weightSSOS_final")
        hist_jet_2_qgl_M[s] = df_M[s].Histo1D(("jet_2_qgl_M","",100,0,1),"jet_2_qgl","weightSSOS_final")
        hist_jet_1_qgl_M[s] = df_M[s].Histo1D(("jet_1_qgl_M","",100,0,1),"jet_1_qgl","weightSSOS_final")

        hist_jet_1_pt_E[s] = df_E[s].Histo1D(("jet_1_pt_E","",50,20,120),"jet_1_pt","weightSSOS_final")
        hist_jet_1_nmu_E[s] = df_E[s].Histo1D(("jet_1_nmu_E","",10,0,10),"jet_1_nmu","weightSSOS_final")
        hist_jet_2_pt_E[s] = df_E[s].Histo1D(("jet_2_pt_E","",50,20,120),"jet_2_pt","weightSSOS_final")
        hist_jet_1_eta_E[s] = df_E[s].Histo1D(("jet_1_eta_E","",80,-4,4),"jet_1_eta","weightSSOS_final")
        hist_jet_2_eta_E[s] = df_E[s].Histo1D(("jet_2_eta_E","",80,-4,4),"jet_2_eta","weightSSOS_final")
        hist_jet_2_nmu_E[s] = df_E[s].Histo1D(("jet_2_nmu_E","",10,0,10),"jet_2_nmu","weightSSOS_final")
        hist_jet_2_mass_E[s] = df_E[s].Histo1D(("jet_2_mass_E","",40,0,40),"jet_2_mass","weightSSOS_final")
        hist_jet_2_qgl_E[s] = df_E[s].Histo1D(("jet_2_qgl_E","",100,0,1),"jet_2_qgl","weightSSOS_final")
        hist_jet_1_qgl_E[s] = df_E[s].Histo1D(("jet_1_qgl_E","",100,0,1),"jet_1_qgl","weightSSOS_final")

        hist_lepton_pt_M[s] = df_M[s].Histo1D(("lepton_pt_M","",50,20,120),"lepton_pt","weightSSOS_final")
        hist_lepton_eta_M[s] = df_M[s].Histo1D(("lepton_eta_M","",80,-4,4),"lepton_eta","weightSSOS_final")
        hist_lepton_eta_thick_M[s] = df_M[s].Histo1D(("lepton_eta_thick_M","",18,-2.7,2.7),"lepton_eta","weightSSOS_final")
        hist_lepton_pt_detail_M[s] = df_M[s].Histo1D(("lepton_pt_detail_M","",80,40,60),"lepton_pt","weightSSOS_final")

        hist_lepton_pt_E[s] = df_E[s].Histo1D(("lepton_pt_E","",50,20,120),"lepton_pt","weightSSOS_final")
        hist_lepton_eta_E[s] = df_E[s].Histo1D(("lepton_eta_E","",80,-4,4),"lepton_eta","weightSSOS_final")
        hist_lepton_eta_thick_E[s] = df_E[s].Histo1D(("lepton_eta_thick_E","",18,-2.7,2.7),"lepton_eta","weightSSOS_final")
        hist_lepton_pt_detail_E[s] = df_E[s].Histo1D(("lepton_pt_detail_E","",80,40,60),"lepton_pt","weightSSOS_final")

        hist_InvM_2jets_M[s] = df_M[s].Histo1D(("InvM_2jets_M","",108,30,300),"InvM_2jets","weightSSOS_final")
        hist_InvM_2jets_E[s] = df_E[s].Histo1D(("InvM_2jets_E","",108,30,300),"InvM_2jets","weightSSOS_final")

        hist_InvM_bot_closer_M[s] = df_M[s].Histo1D(("InvM_bot_closer_M","",100,0,300),"InvM_bot_closer","weightSSOS_final")
        hist_InvM_bot_closer_E[s] = df_E[s].Histo1D(("InvM_bot_closer_E","",100,0,300),"InvM_bot_closer","weightSSOS_final")
        hist_InvM_bot_farther_M[s] = df_M[s].Histo1D(("InvM_bot_farther_M","",100,0,300),"InvM_bot_farther","weightSSOS_final")
        hist_InvM_bot_farther_E[s] = df_E[s].Histo1D(("InvM_bot_farther_E","",100,0,300),"InvM_bot_farther","weightSSOS_final")

        hist_MET_M[s] = df_M[s].Histo1D(("MET_pt_M","",100,0,150),"MET_pt","weightSSOS_final")
        hist_MET_E[s] = df_E[s].Histo1D(("MET_pt_E","",100,0,150),"MET_pt","weightSSOS_final")
        hist_MET_sig_M[s] = df_M[s].Histo1D(("MET_sig_M","",50,0,18),"MET_significance","weightSSOS_final")
        hist_MET_sig_E[s] = df_E[s].Histo1D(("MET_sig_E","",50,0,18),"MET_significance","weightSSOS_final")
        hist_MET_my_sig_M[s] = df_M[s].Histo1D(("MET_my_sig_M","",50,0,18),"MET_my_significance","weightSSOS_final")
        hist_MET_my_sig_E[s] = df_E[s].Histo1D(("MET_my_sig_E","",50,0,18),"MET_my_significance","weightSSOS_final")

        hist_deltaR_jet1_jet2_M[s] = df_M[s].Histo1D(("deltaR_jet1_jet2_M","",100,0,5),"deltaR_jet1_jet2","weightSSOS_final")
        hist_deltaR_jet1_jet2_E[s] = df_E[s].Histo1D(("deltaR_jet1_jet2_E","",100,0,5),"deltaR_jet1_jet2","weightSSOS_final")
        hist_deltaphi_jet1_jet2_M[s] = df_M[s].Histo1D(("deltaphi_jet1_jet2_M","",100,0,5),"deltaphi_jet1_jet2","weightSSOS_final")
        hist_deltaphi_jet1_jet2_E[s] = df_E[s].Histo1D(("deltaphi_jet1_jet2_E","",100,0,5),"deltaphi_jet1_jet2","weightSSOS_final")
        hist_deltaeta_jet1_jet2_M[s] = df_M[s].Histo1D(("deltaeta_jet1_jet2_M","",100,0,5),"deltaeta_jet1_jet2","weightSSOS_final")
        hist_deltaeta_jet1_jet2_E[s] = df_E[s].Histo1D(("deltaeta_jet1_jet2_E","",100,0,5),"deltaeta_jet1_jet2","weightSSOS_final")

        hist_tranverse_mass_M[s] = df_M[s].Histo1D(("transverse_mass_M","",50,0,150),"transverse_mass_old","weightSSOS_final")
        hist_tranverse_mass_E[s] = df_E[s].Histo1D(("transverse_mass_E","",50,0,150),"transverse_mass_old","weightSSOS_final")

        hist_tracks_jet1_M[s] = df_M[s].Histo1D(("tracks_jet1_M","",60,0,60),"tracks_jet1","weightSSOS_final")
        hist_tracks_jet2_M[s] = df_M[s].Histo1D(("tracks_jet2_M","",60,0,60),"tracks_jet2","weightSSOS_final")
        hist_tracks_jet1_E[s] = df_E[s].Histo1D(("tracks_jet1_E","",60,0,60),"tracks_jet1","weightSSOS_final")
        hist_tracks_jet2_E[s] = df_E[s].Histo1D(("tracks_jet2_E","",60,0,60),"tracks_jet2","weightSSOS_final")

        hist_EMN_jet1_M[s] = df_M[s].Histo1D(("EMN_jet1_M","",60,0,1),"EMN_jet1","weightSSOS_final")
        hist_EMC_jet1_M[s] = df_M[s].Histo1D(("EMC_jet1_M","",60,0,1),"EMC_jet1","weightSSOS_final")
        hist_EMN_jet1_E[s] = df_E[s].Histo1D(("EMN_jet1_E","",60,0,1),"EMN_jet1","weightSSOS_final")
        hist_EMC_jet1_E[s] = df_E[s].Histo1D(("EMC_jet1_E","",60,0,1),"EMC_jet1","weightSSOS_final")
        hist_EMtotal_jet1_M[s] = df_M[s].Histo1D(("EMtotal_jet1_M","",60,0,1),"EMtotal_jet1","weightSSOS_final")
        hist_EMtotal_jet1_E[s] = df_E[s].Histo1D(("EMtotal_jet1_E","",60,0,1),"EMtotal_jet1","weightSSOS_final")

        hist_pT_sum_M[s] = df_M[s].Histo1D(("pT_sum_M","",100,0,1),"pT_sum","weightSSOS_final")
        hist_pT_sum_E[s] = df_E[s].Histo1D(("pT_sum_E","",100,0,1),"pT_sum","weightSSOS_final")

        hist_pT_product_M[s] = df_M[s].Histo1D(("pT_product_M","",100,-1,1),"pT_product","weightSSOS_final")
        hist_pT_product_E[s] = df_E[s].Histo1D(("pT_product_E","",100,-1,1),"pT_product","weightSSOS_final")

        hist_deltaR_lep_2jets_M[s] = df_M[s].Histo1D(("deltaR_lep_2jets_M","",100,0,5),"deltaR_lep_2jets","weightSSOS_final")
        hist_deltaR_lep_2jets_E[s] = df_E[s].Histo1D(("deltaR_lep_2jets_E","",100,0,5),"deltaR_lep_2jets","weightSSOS_final")
        hist_deltaphi_MET_2jets_M[s] = df_M[s].Histo1D(("deltaphi_MET_2jets_M","",100,0,5),"deltaphi_MET_2jets","weightSSOS_final")
        hist_deltaphi_MET_2jets_E[s] = df_E[s].Histo1D(("deltaphi_MET_2jets_E","",100,0,5),"deltaphi_MET_2jets","weightSSOS_final")

        hist_deltaPhi_MET_lep_M[s] = df_M[s].Histo1D(("deltaphi_MET_lep_M","",100,0,5),"deltaphi_MET_lep","weightSSOS_final")
       	hist_deltaPhi_MET_lep_E[s] = df_E[s].Histo1D(("deltaphi_MET_lep_E","",100,0,5),"deltaphi_MET_lep","weightSSOS_final")

        hist_deltaphi_MET_jets_1_M[s] = df_M[s].Histo1D(("deltaphi_MET_jets_1_M","",100,0,5),"deltaphi_MET_jets_1","weightSSOS_final")
        hist_deltaphi_MET_jets_1_E[s] = df_E[s].Histo1D(("deltaphi_MET_jets_1_E","",100,0,5),"deltaphi_MET_jets_1","weightSSOS_final")
        hist_deltaphi_MET_jets_2_M[s] = df_M[s].Histo1D(("deltaphi_MET_jets_2_M","",100,0,5),"deltaphi_MET_jets_2","weightSSOS_final")
        hist_deltaphi_MET_jets_2_E[s] = df_E[s].Histo1D(("deltaphi_MET_jets_2_E","",100,0,5),"deltaphi_MET_jets_2","weightSSOS_final")
        hist_deltaphi_lephad_M[s] = df_M[s].Histo1D(("deltaphi_lephad_M","",100,0,5),"deltaphi_lephad","weightSSOS_final")
        hist_deltaphi_lephad_E[s] = df_E[s].Histo1D(("deltaphi_lephad_E","",100,0,5),"deltaphi_lephad","weightSSOS_final")

        hist_eta_2jets_M[s] = df_M[s].Histo1D(("eta_2jets_M","",50,-5,5),"eta_2jets","weightSSOS_final")
        hist_eta_2jets_E[s] = df_E[s].Histo1D(("eta_2jets_E","",50,-5,5),"eta_2jets","weightSSOS_final")
        hist_pt_2jets_M[s] = df_M[s].Histo1D(("pt_2jets_M","",100,0,200),"pt_2jets","weightSSOS_final")
        hist_pt_2jets_E[s] = df_E[s].Histo1D(("pt_2jets_E","",100,0,200),"pt_2jets","weightSSOS_final")

        hist_pt_Wlep_M[s] = df_M[s].Histo1D(("pt_Wlep_M","",100,0,200),"pT_Wlep","weightSSOS_final")
        hist_pt_Wlep_E[s] = df_E[s].Histo1D(("pt_Wlep_E","",100,0,200),"pT_Wlep","weightSSOS_final")

        hist_deltaR_lep_jet1_M[s] = df_M[s].Histo1D(("deltaR_lep_jet1_M","",100,0,5),"deltaR_lep_jet1","weightSSOS_final")
        hist_deltaR_lep_jet2_M[s] = df_M[s].Histo1D(("deltaR_lep_jet2_M","",100,0,5),"deltaR_lep_jet2","weightSSOS_final")
        hist_deltaPhi_lep_jet1_M[s] = df_M[s].Histo1D(("deltaPhi_lep_jet1_M","",100,0,5),"deltaPhi_lep_jet1","weightSSOS_final")
        hist_deltaPhi_lep_jet2_M[s] = df_M[s].Histo1D(("deltaPhi_lep_jet2_M","",100,0,5),"deltaPhi_lep_jet2","weightSSOS_final")
        hist_deltaEta_lep_jet1_M[s] = df_M[s].Histo1D(("deltaEta_lep_jet1_M","",100,0,5),"deltaEta_lep_jet1","weightSSOS_final")
        hist_deltaEta_lep_jet2_M[s] = df_M[s].Histo1D(("deltaEta_lep_jet2_M","",100,0,5),"deltaEta_lep_jet2","weightSSOS_final")

        hist_deltaR_lep_jet1_E[s] = df_E[s].Histo1D(("deltaR_lep_jet1_E","",100,0,5),"deltaR_lep_jet1","weightSSOS_final")
        hist_deltaR_lep_jet2_E[s] = df_E[s].Histo1D(("deltaR_lep_jet2_E","",100,0,5),"deltaR_lep_jet2","weightSSOS_final")
        hist_deltaPhi_lep_jet1_E[s] = df_E[s].Histo1D(("deltaPhi_lep_jet1_E","",100,0,5),"deltaPhi_lep_jet1","weightSSOS_final")
        hist_deltaPhi_lep_jet2_E[s] = df_E[s].Histo1D(("deltaPhi_lep_jet2_E","",100,0,5),"deltaPhi_lep_jet2","weightSSOS_final")
        hist_deltaEta_lep_jet1_E[s] = df_E[s].Histo1D(("deltaEta_lep_jet1_E","",100,0,5),"deltaEta_lep_jet1","weightSSOS_final")
        hist_deltaEta_lep_jet2_E[s] = df_E[s].Histo1D(("deltaEta_lep_jet2_E","",100,0,5),"deltaEta_lep_jet2","weightSSOS_final")

        hist_jet_bot1_btag_M[s] = df_M[s].Histo1D(("jet_bot1_btag_M","",50,0,1),"jet_bot1_btag","weightSSOS_final")
        hist_jet_bot1_btag_E[s] = df_E[s].Histo1D(("jet_bot1_btag_E","",50,0,1),"jet_bot1_btag","weightSSOS_final")
        hist_jet_bot2_btag_M[s] = df_M[s].Histo1D(("jet_bot2_btag_M","",50,0,1),"jet_bot2_btag","weightSSOS_final")
        hist_jet_bot2_btag_E[s] = df_E[s].Histo1D(("jet_bot2_btag_E","",50,0,1),"jet_bot2_btag","weightSSOS_final")

        hist_jet_bot1_pt_M[s] = df_M[s].Histo1D(("jet_bot1_pt_M","",50,20,120),"jet_bot1_pt","weightSSOS_final")
        hist_jet_bot1_eta_M[s] = df_M[s].Histo1D(("jet_bot1_eta_M","",80,-4,4),"jet_bot1_eta","weightSSOS_final")
        hist_jet_bot1_pt_E[s] = df_E[s].Histo1D(("jet_bot1_pt_E","",50,20,120),"jet_bot1_pt","weightSSOS_final")
        hist_jet_bot1_eta_E[s] = df_E[s].Histo1D(("jet_bot1_eta_E","",80,-4,4),"jet_bot1_eta","weightSSOS_final")
        hist_jet_bot2_pt_M[s] = df_M[s].Histo1D(("jet_bot2_pt_M","",50,20,120),"jet_bot2_pt","weightSSOS_final")
        hist_jet_bot2_eta_M[s] = df_M[s].Histo1D(("jet_bot2_eta_M","",80,-4,4),"jet_bot2_eta","weightSSOS_final")
        hist_jet_bot2_pt_E[s] = df_E[s].Histo1D(("jet_bot2_pt_E","",50,20,120),"jet_bot2_pt","weightSSOS_final")
        hist_jet_bot2_eta_E[s] = df_E[s].Histo1D(("jet_bot2_eta_E","",80,-4,4),"jet_bot2_eta","weightSSOS_final")

        hist_jet_bot1_btag_thick_M[s] = df_M[s].Histo1D(("jet_bot1_btag_thick_M","",10,0,1),"jet_bot1_btag","weightSSOS_final")
        hist_jet_bot1_btag_thick_E[s] = df_E[s].Histo1D(("jet_bot1_btag_thick_E","",10,0,1),"jet_bot1_btag","weightSSOS_final")
        hist_jet_bot2_btag_thick_M[s] = df_M[s].Histo1D(("jet_bot2_btag_thick_M","",10,0,1),"jet_bot2_btag","weightSSOS_final")
        hist_jet_bot2_btag_thick_E[s] = df_E[s].Histo1D(("jet_bot2_btag_thick_E","",10,0,1),"jet_bot2_btag","weightSSOS_final")

        hist_jet_1_btag_M[s] = df_M[s].Histo1D(("jet_1_btag_M","",50,0,1),"jet_1_btag","weightSSOS_final")
        hist_jet_1_btag_E[s] = df_E[s].Histo1D(("jet_1_btag_E","",50,0,1),"jet_1_btag","weightSSOS_final")
        hist_jet_2_btag_M[s] = df_M[s].Histo1D(("jet_2_btag_M","",50,0,1),"jet_2_btag","weightSSOS_final")
        hist_jet_2_btag_E[s] = df_E[s].Histo1D(("jet_2_btag_E","",50,0,1),"jet_2_btag","weightSSOS_final")
        hist_jet_1_btag_thick_M[s] = df_M[s].Histo1D(("jet_1_btag_thick_M","",10,0,1),"jet_1_btag","weightSSOS_final")
        hist_jet_1_btag_thick_E[s] = df_E[s].Histo1D(("jet_1_btag_thick_E","",10,0,1),"jet_1_btag","weightSSOS_final")
        hist_jet_2_btag_thick_M[s] = df_M[s].Histo1D(("jet_2_btag_thick_M","",10,0,1),"jet_2_btag","weightSSOS_final")
        hist_jet_2_btag_thick_E[s] = df_E[s].Histo1D(("jet_2_btag_thick_E","",10,0,1),"jet_2_btag","weightSSOS_final")

        hist_jet_bot1_btagnumber_M[s] = df_M[s].Histo1D(("jet_bot1_btagnumber_M","",3,0,3),"jet_bot1_btagnumber","weightSSOS_final")
        hist_jet_bot1_btagnumber_E[s] = df_E[s].Histo1D(("jet_bot1_btagnumber_E","",3,0,3),"jet_bot1_btagnumber","weightSSOS_final")
        hist_jet_bot2_btagnumber_M[s] = df_M[s].Histo1D(("jet_bot2_btagnumber_M","",3,0,3),"jet_bot2_btagnumber","weightSSOS_final")
        hist_jet_bot2_btagnumber_E[s] = df_E[s].Histo1D(("jet_bot2_btagnumber_E","",3,0,3),"jet_bot2_btagnumber","weightSSOS_final")
        hist_jet_1_btagnumber_M[s] = df_M[s].Histo1D(("jet_1_btagnumber_M","",3,0,3),"jet_1_btagnumber","weightSSOS_final")
        hist_jet_1_btagnumber_E[s] = df_E[s].Histo1D(("jet_1_btagnumber_E","",3,0,3),"jet_1_btagnumber","weightSSOS_final")
        hist_jet_2_btagnumber_M[s] = df_M[s].Histo1D(("jet_2_btagnumber_M","",3,0,3),"jet_2_btagnumber","weightSSOS_final")
        hist_jet_2_btagnumber_E[s] = df_E[s].Histo1D(("jet_2_btagnumber_E","",3,0,3),"jet_2_btagnumber","weightSSOS_final")

        hist_pT_proy_M[s] = df_M[s].Histo1D(("pT_proy_M","",100,-100,100),"pT_proy","weightSSOS_final")
        hist_pT_proy_E[s] = df_E[s].Histo1D(("pT_proy_E","",100,-100,100),"pT_proy","weightSSOS_final")

        hist_pT_sum_2J_M[s] = df_M[s].Histo1D(("pT_sum_2J_M","",100,0,1),"pT_sum_2J","weightSSOS_final")
        hist_pT_sum_2J_E[s] = df_E[s].Histo1D(("pT_sum_2J_E","",100,0,1),"pT_sum_2J","weightSSOS_final")

        hist_lepton_pt_eta_M[s] = df_M[s].Histo2D(("lepton_pt_eta_M","",80,-4,4,50,40,60),"lepton_eta","lepton_pt","weightSSOS_final")
        hist_lepton_pt_eta_E[s] = df_E[s].Histo2D(("lepton_pt_eta_E","",80,-4,4,50,40,60),"lepton_eta","lepton_pt","weightSSOS_final")

########## gen hists

if mode == "mc":
        for s in samples:
               hist_jet_1_flavourP_M[s] = df_M[s].Histo1D(("jet_1_flavourP_M","",28,-6,22),"jet_1_flavourP","weightSSOS_final")
               hist_jet_1_flavourP_E[s] = df_E[s].Histo1D(("jet_1_flavourP_E","",28,-6,22),"jet_1_flavourP","weightSSOS_final")
               hist_jet_2_flavourP_M[s] = df_M[s].Histo1D(("jet_2_flavourP_M","",28,-6,22),"jet_2_flavourP","weightSSOS_final")
               hist_jet_2_flavourP_E[s] = df_E[s].Histo1D(("jet_2_flavourP_E","",28,-6,22),"jet_2_flavourP","weightSSOS_final")
               hist_jet_bot1_flavourP_M[s] = df_M[s].Histo1D(("jet_bot1_flavourP_M","",28,-6,22),"jet_bot1_flavourP","weightSSOS_final")
               hist_jet_bot1_flavourP_E[s] = df_E[s].Histo1D(("jet_bot1_flavourP_E","",28,-6,22),"jet_bot1_flavourP","weightSSOS_final")
               hist_jet_bot2_flavourP_M[s] = df_M[s].Histo1D(("jet_bot2_flavourP_M","",28,-6,22),"jet_bot2_flavourP","weightSSOS_final")
               hist_jet_bot2_flavourP_E[s] = df_E[s].Histo1D(("jet_bot2_flavourP_E","",28,-6,22),"jet_bot2_flavourP","weightSSOS_final")

               hist_btag_sf_M[s] = df_M[s].Histo1D(("btag_sf_M","",100,0.5,1.5),"btag_sf","weightSSOS_final")
               hist_btag_sf_E[s] = df_E[s].Histo1D(("btag_sf_E","",100,0.5,1.5),"btag_sf","weightSSOS_final")
               hist_lep_id_sf_M[s] = df_M[s].Histo1D(("lep_id_sf_M","",100,0.5,1.5),"lep_id_sf","weightSSOS_final")
               hist_lep_id_sf_E[s] = df_E[s].Histo1D(("lep_id_sf_E","",100,0.5,1.5),"lep_id_sf","weightSSOS_final")
               hist_lep_iso_sf_M[s] = df_M[s].Histo1D(("lep_iso_sf_M","",100,0.5,1.5),"lep_iso_sf","weightSSOS_final")
               hist_lep_trig_sf_M[s] = df_M[s].Histo1D(("lep_trig_sf_M","",100,0.5,1.5),"lep_trig_sf","weightSSOS_final")
               hist_lep_trig_sf_E[s] = df_E[s].Histo1D(("lep_trig_sf_E","",100,0.5,1.5),"lep_trig_sf","weightSSOS_final")
               hist_puWeight_M[s] = df_M[s].Histo1D(("puWeight_M","",100,0.5,1.5),"puWeight","weightSSOS_final")
               hist_puWeight_E[s] = df_E[s].Histo1D(("puWeight_E","",100,0.5,1.5),"puWeight","weightSSOS_final")
               hist_PUjetID_SF_M[s] = df_M[s].Histo1D(("PUjetID_SF_M","",100,0.5,1.5),"PUjetID_SF","weightSSOS_final")
               hist_PUjetID_SF_E[s] = df_E[s].Histo1D(("PUjetID_SF_E","",100,0.5,1.5),"PUjetID_SF","weightSSOS_final")
               hist_topweight_M[s] = df_M[s].Histo1D(("topweight_M","",100,0.5,1.5),"top_weight","weightSSOS_final")
               hist_topweight_E[s] = df_E[s].Histo1D(("topweight_E","",100,0.5,1.5),"top_weight","weightSSOS_final")

#############################
####     DATA SAVING     ####
#############################

for s in samples:
                path_hist = '/nfs/cms/vazqueze/hists_ttbar/hists/higgs/fromJF/wqq/lepton50/ptjets30/histstt_wqq_fromJF_'+s+'.root'
                myfile = TFile( path_hist, 'RECREATE' )

                hist_nJetGood_M[s].Write()
                hist_nJetGood_E[s].Write()
                hist_jet_1_pt_M[s].Write()
                hist_jet_1_nmu_M[s].Write()
                hist_jet_1_eta_M[s].Write()
                hist_jet_1_pt_E[s].Write()
                hist_jet_1_eta_E[s].Write()
                hist_jet_1_nmu_E[s].Write()
                hist_jet_2_pt_M[s].Write()
                hist_jet_2_eta_M[s].Write()
                hist_jet_2_nmu_M[s].Write()
                hist_jet_2_mass_M[s].Write()
                hist_jet_2_qgl_M[s].Write()
                hist_jet_1_qgl_M[s].Write()
                hist_jet_2_pt_E[s].Write()
                hist_jet_2_eta_E[s].Write()
                hist_jet_2_nmu_E[s].Write()
                hist_jet_2_nmu_E[s].Write()
                hist_jet_2_mass_E[s].Write()
                hist_jet_2_qgl_E[s].Write()
                hist_jet_1_qgl_E[s].Write()
                hist_lepton_pt_M[s].Write()
                hist_lepton_pt_detail_M[s].Write()
                hist_lepton_eta_M[s].Write()
                hist_lepton_eta_thick_M[s].Write()
                hist_lepton_pt_E[s].Write()
                hist_lepton_pt_detail_E[s].Write()
                hist_lepton_eta_E[s].Write()
                hist_lepton_eta_thick_E[s].Write()
                hist_InvM_2jets_M[s].Write()
                hist_InvM_2jets_E[s].Write()
                hist_InvM_bot_closer_M[s].Write()
                hist_InvM_bot_closer_E[s].Write()
                hist_InvM_bot_farther_M[s].Write()
                hist_InvM_bot_farther_E[s].Write()
                hist_deltaR_jet1_jet2_M[s].Write()
                hist_deltaR_jet1_jet2_E[s].Write()
                hist_deltaphi_jet1_jet2_M[s].Write()
                hist_deltaphi_jet1_jet2_E[s].Write()
                hist_deltaeta_jet1_jet2_M[s].Write()
                hist_deltaeta_jet1_jet2_E[s].Write()
                hist_MET_M[s].Write()
                hist_MET_E[s].Write()
                hist_MET_sig_M[s].Write()
                hist_MET_sig_E[s].Write()
                hist_MET_my_sig_M[s].Write()
                hist_MET_my_sig_E[s].Write()
                hist_tranverse_mass_M[s].Write()
                hist_tranverse_mass_E[s].Write()
                hist_tracks_jet1_M[s].Write()
                hist_tracks_jet2_M[s].Write()
                hist_tracks_jet1_E[s].Write()
                hist_tracks_jet2_E[s].Write()
                hist_EMN_jet1_M[s].Write()
                hist_EMC_jet1_M[s].Write()
                hist_EMN_jet1_E[s].Write()
                hist_EMC_jet1_E[s].Write()
                hist_EMtotal_jet1_M[s].Write()
                hist_EMtotal_jet1_E[s].Write()
                hist_pT_sum_M[s].Write()
                hist_pT_sum_E[s].Write()
                hist_pT_product_M[s].Write()
                hist_pT_product_E[s].Write()
                hist_deltaR_lep_2jets_M[s].Write()
                hist_deltaR_lep_2jets_E[s].Write()
                hist_deltaphi_MET_2jets_M[s].Write()
                hist_deltaphi_MET_2jets_E[s].Write()
                hist_deltaphi_MET_jets_1_M[s].Write()
                hist_deltaphi_MET_jets_1_E[s].Write()
                hist_deltaphi_MET_jets_2_M[s].Write()
                hist_deltaphi_MET_jets_2_E[s].Write()
                hist_deltaphi_lephad_M[s].Write()
                hist_deltaphi_lephad_E[s].Write()
                hist_deltaPhi_MET_lep_M[s].Write()
       	       	hist_deltaPhi_MET_lep_E[s].Write()
                hist_eta_2jets_M[s].Write()
                hist_eta_2jets_E[s].Write()
                hist_pt_2jets_M[s].Write()
                hist_pt_2jets_E[s].Write()
                hist_pt_Wlep_M[s].Write()
                hist_pt_Wlep_E[s].Write()
                hist_deltaR_lep_jet1_M[s].Write()
                hist_deltaR_lep_jet2_M[s].Write()
                hist_deltaPhi_lep_jet1_M[s].Write()
                hist_deltaPhi_lep_jet2_M[s].Write()
                hist_deltaEta_lep_jet1_M[s].Write()
                hist_deltaEta_lep_jet2_M[s].Write()
                hist_deltaR_lep_jet1_E[s].Write()
                hist_deltaR_lep_jet2_E[s].Write()
                hist_deltaPhi_lep_jet1_E[s].Write()
                hist_deltaPhi_lep_jet2_E[s].Write()
                hist_deltaEta_lep_jet1_E[s].Write()
                hist_deltaEta_lep_jet2_E[s].Write()
                hist_jet_1_btag_M[s].Write()
                hist_jet_1_btag_E[s].Write()
                hist_jet_2_btag_M[s].Write()
                hist_jet_2_btag_E[s].Write()
                hist_jet_bot1_btag_M[s].Write()
                hist_jet_bot1_btag_E[s].Write()
                hist_jet_bot2_btag_M[s].Write()
                hist_jet_bot2_btag_E[s].Write()
                hist_jet_bot1_pt_M[s].Write()
                hist_jet_bot1_pt_E[s].Write()
                hist_jet_bot2_pt_M[s].Write()
                hist_jet_bot2_pt_E[s].Write()
                hist_jet_bot1_eta_M[s].Write()
                hist_jet_bot1_eta_E[s].Write()
                hist_jet_bot2_eta_M[s].Write()
                hist_jet_bot2_eta_E[s].Write()
                hist_jet_bot1_btag_thick_M[s].Write()
                hist_jet_bot1_btag_thick_E[s].Write()
                hist_jet_bot2_btag_thick_M[s].Write()
                hist_jet_bot2_btag_thick_E[s].Write()
                hist_jet_1_btag_thick_M[s].Write()
                hist_jet_1_btag_thick_E[s].Write()
                hist_jet_2_btag_thick_M[s].Write()
                hist_jet_2_btag_thick_E[s].Write()
                hist_jet_bot1_btagnumber_M[s].Write()
       	       	hist_jet_bot1_btagnumber_E[s].Write()
       	       	hist_jet_bot2_btagnumber_M[s].Write()
       	       	hist_jet_bot2_btagnumber_E[s].Write()
                hist_jet_1_btagnumber_M[s].Write()
                hist_jet_1_btagnumber_E[s].Write()
                hist_jet_2_btagnumber_M[s].Write()
       	       	hist_jet_2_btagnumber_E[s].Write()
                hist_pT_proy_M[s].Write()
                hist_pT_proy_E[s].Write()
                hist_pT_sum_2J_M[s].Write()
                hist_pT_sum_2J_E[s].Write()

                if mode=="mc":
                        hist_jet_bot1_flavourP_M[s].Write()
                        hist_jet_bot1_flavourP_E[s].Write()
                        hist_jet_bot2_flavourP_M[s].Write()
                        hist_jet_bot2_flavourP_E[s].Write()
                        hist_jet_1_flavourP_M[s].Write()
                        hist_jet_1_flavourP_E[s].Write()
                        hist_jet_2_flavourP_M[s].Write()
                        hist_jet_2_flavourP_E[s].Write()
                        hist_btag_sf_M[s].Write()
                        hist_btag_sf_E[s].Write()
                        hist_lep_id_sf_M[s].Write()
                        hist_lep_id_sf_E[s].Write()
                        hist_lep_iso_sf_M[s].Write()
                        hist_lep_trig_sf_M[s].Write()
                        hist_lep_trig_sf_E[s].Write()
                        hist_puWeight_M[s].Write()
                        hist_puWeight_E[s].Write()
                        hist_PUjetID_SF_M[s].Write()
                        hist_PUjetID_SF_E[s].Write()
                        hist_topweight_M[s].Write()
                        hist_topweight_E[s].Write()

                myfile.Close()

print('Ended succesfully')

