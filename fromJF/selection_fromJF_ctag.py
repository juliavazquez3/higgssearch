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
parser.add_argument("--wcs", action="store_true", default=False,
                    help="Wcs clasiffication or not")
parser.add_argument("--charmtag", type=string, default="no", 
                    help="Perform ssos substraction")
parser.add_argument("--presel", type=string, default="btagMM_chitest",
                    help="preselection")
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

if args.charmtag == "ctag": channel = "ctag"
elif args.charmtag == "noctag": channel = "noctag"
else: raise NameError('Incorrect channel')
print(channel)

if args.presel == "nobtag": chan = "nobtag"
elif args.presel == "btagMM": chan = "btagMM"
elif args.presel == "lepton50": chan = "lepton50"
elif args.presel == "btagMM_chitest": chan = "btagMM_chitest"
elif args.presel == "lepton50_chitest": chan = "lepton50_chitest"
else: raise NameError('Incorrect channel')
print(chan)

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

#list_jetsyst = ["nom","smearup","smeardown","jecup","jecdown"]
list_jetsyst = ["nom"]
#list_jetsyst = ["nom","smearup","smeardown"]
#list_jetsyst = ["nom","jecup","jecdown"]

df = {}
xsecs = {}
sumws = {}
archives = {}
event_test = {}

for syst in list_jetsyst:
  df[syst] = {}
  archives[syst] = {}
  for s in samples:
    archives[syst][s]=[]

names_proc = {"M":{"2016":["2016M1","2016M2","2016M3","2016M4","2016M5","2016M6"],"2016B":["2016BM1","2016BM2","2016BM3"],
   "2017":["2017M1","2017M2","2017M3","2017M4","2017M5"],"2018":["2018M1","2018M2","2018M3","2018M4"]},
   "E":{"2016":["2016E1","2016E2","2016E3","2016E4","2016E5","2016E6"],"2016B":["2016BE1","2016BE2","2016BE3"],
   "2017":["2017E1","2017E2","2017E3","2017E4","2017E5"],"2018":["2018E1","2018E2","2018E3","2018E4"]}}

files_data = {"M":{"2016":230,"2016B":157,"2017":436,"2018":393},
   "E":{"2016":335,"2016B":156,"2017":245,"2018":738}}

term3 = {}

term1 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_ver2/"
term2 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_down/"
term3["nom"] = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux/"
term3["smearup"] = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux_smearup/"
term3["smeardown"] = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux_smeardown/"
term3["jecup"] = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux_jecup/"
term3["jecdown"] = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_aux_jecdown/"
#term3 = "/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID_nosmearing/"
for p in proc:
    # Construct the dataframes
  if args.type == "mc":
    for syst in list_jetsyst:  
        folder = term3[syst] + "myconfig"+years[0]+"/{sample}/cat_base/prod_test/" # Folder name
        if years[0] == "2016B": folder = term3[syst] + "myconfig2016/{sample}/cat_base/prod_test/"
        folder = folder.format(sample = p) # Sample name
        num_files = files[p]["files"] # Number of files
        list_files = [f for f in listdir(folder) if isfile(join(folder, f))] # Lista de archivos
        if (num_files == len(list_files)):
           for f in list_files:
               file_r = TFile(join(folder, f))
               if file_r.GetListOfKeys().Contains("Events"):
                  archives[syst][p+years[0]].append(join(folder,f))
  else:
        if p == "M":
          folder1 = term3["nom"]+"myconfig"+years[0]+"/"
          for f in [f1 for f1 in listdir(folder1) if (isdir(join(folder1, f1)) and f1[0:7] == "data_mu")]:
             folder = folder1 + f + "/cat_base/prod_test/"
             num_files = files_data["M"][years[0]] # Number of files
             list_files = [f1 for f1 in listdir(folder) if isfile(join(folder, f1))] # Lista de archivos
             #if (num_files == len(list_files)):
             for fil in list_files:
                   file_r = TFile(join(folder, fil))
                   if file_r.GetListOfKeys().Contains("Events"):
                     for syst in list_jetsyst:
                       archives[syst][years[0]+"M"].append(join(folder,fil))
        if p == "E":
          folder1 = term3["nom"]+"myconfig"+years[0]+"/"
          for f in [f1 for f1 in listdir(folder1) if (isdir(join(folder1, f1)) and f1[0:7] == "data_el")]:
             folder = folder1 + f + "/cat_base/prod_test/"
             num_files = files_data["E"][years[0]] # Number of files
             list_files = [f1 for f1 in listdir(folder) if isfile(join(folder, f1))] # Lista de archivos
             #if (num_files == len(list_files)):
             for fil in list_files:
                   file_r = TFile(join(folder, fil))
                   if file_r.GetListOfKeys().Contains("Events"):
                     for syst in list_jetsyst:
                       archives[syst][years[0]+"E"].append(join(folder,fil))

for s in samples:
  for syst in list_jetsyst:
    if args.notfull: archives[syst][s]=archives[syst][s][int(args.list[0]):int(args.list[1])]
    df[syst][s] = ROOT.RDataFrame("Events",set(archives[syst][s]))
  print("Number of files for",s,len(archives["nom"][s]))

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

cuts_cvl_deepj = {}
cuts_cvl_deepj["2016"]=[0.039,0.098,0.270];cuts_cvl_deepj["2016B"]=[0.039,0.099,0.269];
cuts_cvl_deepj["2017"]=[0.030,0.085,0.520];cuts_cvl_deepj["2018"]=[0.038,0.099,0.282];
cuts_cvb_deepj = {}
cuts_cvb_deepj["2016"]=[0.327,0.370,0.256];cuts_cvb_deepj["2016B"]=[0.305,0.353,0.247];
cuts_cvb_deepj["2017"]=[0.400,0.340,0.050];cuts_cvb_deepj["2018"]=[0.246,0.325,0.267];

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

######## LUMI Info #############

lumi = {}
xsecs = {}
nevents = {}

if mode == "mc":
   for data_op in years:
        if data_op=="2016B":
           files = json.load(open("/nfs/cms/vazqueze/higgssearch/mcinfo2016.json"))
        else:
           files = json.load(open("/nfs/cms/vazqueze/higgssearch/mcinfo"+data_op+".json"))
        lumi[data_op] = {}
        xsecs[data_op] = {}
        for p in samples:
                if years[0]=="2016B":
                  p = p[:-5]
                else:
                  p = p[:-4]
                xsecval = files[p]["xsec"]
                num_files = files[p]["files"] # Number of files
                luminosity = files[p]["lumi"] # Luminosity
                #print(files[p]["type"])
                lumi[data_op][p] = luminosity
                xsecs[data_op][p] = xsecval
        if (args.process == "allMC" or args.process == "ttbar_sl"):
          lumi[data_op]["ttbar_sl_nomuon"] = lumi[data_op]["ttbar_sl"]
          lumi[data_op]["ttbar_sl_muon_charm"] = lumi[data_op]["ttbar_sl"]
          lumi[data_op]["ttbar_sl_muon_bottom"] = lumi[data_op]["ttbar_sl"]
          lumi[data_op]["ttbar_sl_muon_else"] = lumi[data_op]["ttbar_sl"]
          xsecs[data_op]["ttbar_sl_nomuon"] = xsecs[data_op]["ttbar_sl"]
          xsecs[data_op]["ttbar_sl_muon_charm"] = xsecs[data_op]["ttbar_sl"]
          xsecs[data_op]["ttbar_sl_muon_bottom"] = xsecs[data_op]["ttbar_sl"]
          xsecs[data_op]["ttbar_sl_muon_else"] = xsecs[data_op]["ttbar_sl"]
        if (args.process == "allMC" or args.process == "ttbar_dl"):
          lumi[data_op]["ttbar_dl_nomuon"] = lumi[data_op]["ttbar_dl"]
          lumi[data_op]["ttbar_dl_muon_charm"] = lumi[data_op]["ttbar_dl"]
          lumi[data_op]["ttbar_dl_muon_bottom"] = lumi[data_op]["ttbar_dl"]
          lumi[data_op]["ttbar_dl_muon_else"] = lumi[data_op]["ttbar_dl"]
          xsecs[data_op]["ttbar_dl_nomuon"] = xsecs[data_op]["ttbar_dl"]
          xsecs[data_op]["ttbar_dl_muon_charm"] = xsecs[data_op]["ttbar_dl"]
          xsecs[data_op]["ttbar_dl_muon_bottom"] = xsecs[data_op]["ttbar_dl"]
          xsecs[data_op]["ttbar_dl_muon_else"] = xsecs[data_op]["ttbar_dl"]

listsampl = ["ww","wjets_1","wjets_2","wjets_3","wjets_4","wjets_5","wjets_6","wjets_7","wjets_8",
        "ttbar_sl","ttbar_dl","ttbar_dh","zjets_1","zjets_2","zjets_3","zjets_4","zjets_5","zjets_6",
        "zjets_7","zjets_8","st_1","st_2","st_3","st_4","zz","wz", "ttbar_sl_nomuon", "ttbar_sl_muon_charm",
        "ttbar_sl_muon_bottom", "ttbar_sl_muon_else","ttbar_dl_nomuon", "ttbar_dl_muon_charm",
        "ttbar_dl_muon_bottom", "ttbar_dl_muon_else","ttbar_sl","ttbar_dl"]

lumi_d = {}
lumi_d["2016"] = 19.5
lumi_d["2016B"] = 16.8
lumi_d["2017"] = 41.5
lumi_d["2018"] = 59.8

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
      auto ttbarcharm(UInt_t nPart, Vint status, Vint pdg, Vint mother) {
            vector<bool> vb;
            bool typeC = false;
            bool typeL = false;
            int arr1[] = {1,2,3};
            bool cond = false;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==4 && fabs(pdg[mother[i]])==24) {
                          typeC = true;
                }
                cond = std::find(std::begin(arr1), std::end(arr1), fabs(pdg[i])) != std::end(arr1);
                if (cond && fabs(pdg[mother[i]])==24) {
                          typeL = true;
                }
            }
            vb.push_back(typeC);
            vb.push_back(typeL);
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
                if(condmu && mu_id[i] && mu_pt[i]<25. && (cond1 || cond2) && mu_pt[i]>ptM && mu_pt[i]>3. && fabs(mu_eta[i])<2.4){
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

### tau selector

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto tau_selection(Vint goodind, Vint botind, Vint xind, Vfloat pt, Vfloat eta, Vfloat phi, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi){
            vector<int> vb;
            bool cond1 = false;
            bool cond2 = false;
            bool cond3 = false;
            bool cond4 = false;
            int indM1 =-1;
            int indM2 =-1;
            int indM3 =-1;
            int indM4 =-1;
            float ptM1{-10.};
            float ptM2{-10.};
            float ptM3{-10.};
            float ptM4{-10.};
            for (unsigned int i=0; i<mu_pt.size(); ++i){
                cond1 = ROOT::VecOps::DeltaR(mu_eta[i],eta[xind[0]],mu_phi[i],phi[xind[0]]) < 0.4;
                cond2 = ROOT::VecOps::DeltaR(mu_eta[i],eta[xind[1]],mu_phi[i],phi[xind[1]]) < 0.4;
                cond3 = ROOT::VecOps::DeltaR(mu_eta[i],eta[goodind[botind[0]]],mu_phi[i],phi[goodind[botind[0]]]) < 0.4;
                cond4 = ROOT::VecOps::DeltaR(mu_eta[i],eta[goodind[botind[1]]],mu_phi[i],phi[goodind[botind[1]]]) < 0.4;
                if(cond1 && mu_pt[i]>ptM1){
                     indM1 = i;
                     ptM1 = mu_pt[i];
                } else if (cond2 && mu_pt[i]>ptM2) {
                     indM2 = i;
                     ptM2 = mu_pt[i];
                } else if (cond3 && mu_pt[i]>ptM3) {
                     indM3 = i;
                     ptM3 = mu_pt[i];
                } else if (cond4 && mu_pt[i]>ptM4) {
                     indM4 = i;
                     ptM4 = mu_pt[i];
                }
            }
            vb.push_back(indM1);
            vb.push_back(indM2);
            vb.push_back(indM3);
            vb.push_back(indM4);
            return vb;
      };
""")

### mother of muon

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonmotherdhad(Vint pdg, Vint mother, Vint muon_gen, const int muon_id){
            int mother_id_v1 = fabs(pdg[mother[muon_gen[muon_id]]]);
            int gmother_id_v1 = fabs(pdg[mother[mother[muon_gen[muon_id]]]]);
            int gmother_id_v2 = fabs(pdg[mother[mother[mother[muon_gen[muon_id]]]]]);
            int gmother_id_v3 = fabs(pdg[mother[mother[mother[mother[muon_gen[muon_id]]]]]]);
            int gmother_id_v4 = fabs(pdg[mother[mother[mother[mother[mother[muon_gen[muon_id]]]]]]]);
            int mother_id_v2 = -10;
            if(mother_id_v1 == 13){
                 mother_id_v1 = fabs(pdg[mother[mother[muon_gen[muon_id]]]]);
                 gmother_id_v1 = fabs(pdg[mother[mother[mother[muon_gen[muon_id]]]]]);
                 gmother_id_v2 = fabs(pdg[mother[mother[mother[mother[muon_gen[muon_id]]]]]]);
                 gmother_id_v3 = fabs(pdg[mother[mother[mother[mother[mother[muon_gen[muon_id]]]]]]]);
                 gmother_id_v4 = fabs(pdg[mother[mother[mother[mother[mother[mother[muon_gen[muon_id]]]]]]]]);
            }
            int arr[] = {511,521,531,5122,5,10411,10421,413,423,10413,10423,20413,20423,415,425,10431,433,10433,20433,435,24};
            int arr1[] = {511,521,531,5122,5};
            int arr2[] = {411,421,431,4122};
            bool cond = std::find(std::begin(arr), std::end(arr), mother_id_v1) != std::end(arr);
            bool cond1 = std::find(std::begin(arr1), std::end(arr1), mother_id_v1) != std::end(arr1);
            bool cond2 = std::find(std::begin(arr2), std::end(arr2), mother_id_v1) != std::end(arr2);
            bool cond3 = ((gmother_id_v1==24) || (gmother_id_v2==24)) || ((gmother_id_v3==24) || (gmother_id_v4==24));
            bool cond4 = cond2 && cond3;
            vector<bool> vb;
            vb.push_back(cond);
            vb.push_back(cond1);
            vb.push_back(cond2);
            return cond4;
      };
""")

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
      auto muonmother_ttdilep(Vint pdg, Vint mother, Vint muon_id_gen, const int muon_id){
            int mother_id_v1 = fabs(pdg[mother[muon_id_gen[muon_id]]]);
            int mother_id_v2 = -10;
            if(mother_id_v1 == 13){
                 mother_id_v1 = fabs(pdg[mother[mother[muon_id_gen[muon_id]]]]);
            }
            int arr1[] = {111,211,130,310,311,321,313,323};
            if(mother_id_v1 == 15){
                 mother_id_v2 = 1;
            } else if (mother_id_v1 == 21) {
                 mother_id_v2 = 2;
            } else if (std::find(std::begin(arr1), std::end(arr1), mother_id_v1) != std::end(arr1)) {
                 mother_id_v2 = 3;
            } else {
                 mother_id_v2 = 0;
            }
            return mother_id_v2;
      };
""")

#### muon calculations

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto scalprod(const float x1, const float y1, const float z1, const float x2, const float y2, const float z2) {
          float pr;
          pr = x1*x2+y1*y2+z1*z2;
          return pr;
      };
      auto muons_zvals(Vfloat mu_iso, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass, const int muind, const int jetind) {
            vector<float> vb;
            ROOT::Math::PtEtaPhiMVector muvec;
            ROOT::Math::PtEtaPhiMVector jetvec;
            muvec.SetPt(mu_pt[muind]);muvec.SetEta(mu_eta[muind]);muvec.SetPhi(mu_phi[muind]);muvec.SetM(mu_mass[muind]);
            jetvec.SetPt(jet_pt[jetind]);jetvec.SetEta(jet_eta[jetind]);jetvec.SetPhi(jet_phi[jetind]);jetvec.SetM(jet_mass[jetind]);
            float z2 = (muvec.Dot(jetvec))/(jetvec.Dot(jetvec)); 
            float z3 = (1+mu_iso[muind])*(muvec.Pt()/jetvec.Pt());
            float z4 = (scalprod(muvec.X(),muvec.Y(),muvec.Z(),jetvec.X(),jetvec.Y(),jetvec.Z()))/(scalprod(jetvec.X(),jetvec.Y(),jetvec.Z(),jetvec.X(),jetvec.Y(),jetvec.Z()));
            vb.push_back(z2);
            vb.push_back(z3);
            vb.push_back(z4);
            return vb;
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
      auto metrecalculation(Float_t met_pt, Float_t met_phi, Vfloat jetpt, Vfloat jetptnom, Vfloat jetphi) {
           vector<float> vb;
           auto x = std::cos(jetphi[0])*(jetpt[0]-jetptnom[0]);       
           auto y = std::sin(jetphi[0])*(jetpt[0]-jetptnom[0]);       
           for (unsigned int i=1; i<jetpt.size(); ++i){
                x = x+std::cos(jetphi[i])*(jetpt[i]-jetptnom[i]);
                y = y+std::sin(jetphi[i])*(jetpt[i]-jetptnom[i]);
           }
           auto xmet = met_pt*std::cos(met_phi)-x;
           auto ymet = met_pt*std::sin(met_phi)-y;
           auto new_pt = std::sqrt(xmet*xmet+ymet*ymet);
           auto new_phi = (fabs(ymet)>0. ? std::atan(xmet/ymet) : 0.);
           vb.push_back(new_pt);
           vb.push_back(new_phi);
           return vb;
      };
""")

#############################################################
######################## Definitions ########################
#############################################################

if mode == "data":
   for syst in list_jetsyst:
        for s in samples:
           df[syst][s] = df[syst][s].Define('Jet_pt_nom','Jet_pt')
           df[syst][s] = df[syst][s].Define('Jet_mass_nom','Jet_mass')
           df[syst][s] = df[syst][s].Define('MET_smeared_phi','MET_phi')
           df[syst][s] = df[syst][s].Define('MET_smeared_pt','MET_pt')

df_muon = {}
df_electron = {}
df_test = {}
for syst in list_jetsyst:
    df_muon[syst] = {}
    df_electron[syst] = {}
    df_test[syst] = {}

##### Systematic run

if mode == "mc":
   for s in samples:
     for syst in list_jetsyst:
       if syst == "nom":
           df[syst][s] = df[syst][s].Define('MET_pt_aux','MET_smeared_pt')
           df[syst][s] = df[syst][s].Define('MET_phi_aux','MET_smeared_phi')
           df[syst][s] = df[syst][s].Define('Jet_pt_aux','0.99*Jet_pt_nom')
       else:
           df[syst][s] = df[syst][s].Define('Jet_pt_aux','Jet_pt_FJ')
           df[syst][s] = df[syst][s].Define('aux4met','metrecalculation(MET_smeared_pt, MET_smeared_phi, Jet_pt_FJ, Jet_pt_nom, Jet_phi)')
           df[syst][s] = df[syst][s].Define('MET_pt_aux','aux4met[0]')
           df[syst][s] = df[syst][s].Define('MET_phi_aux','aux4met[1]')

elif mode == "data":
 for syst in list_jetsyst:
   for s in samples:
     df[syst][s] = df[syst][s].Define('Jet_pt_aux','Jet_pt_nom')
     df[syst][s] = df[syst][s].Define('MET_pt_aux','MET_smeared_pt')
     df[syst][s] = df[syst][s].Define('MET_phi_aux','MET_smeared_phi')

###########    Extra definitions

for s in samples:
    for syst in list_jetsyst:
        df[syst][s] = df[syst][s].Define('lepton_charge','nMuonGood>0 ? Muon_charge[MuonGoodInd[0]] : Electron_charge[ElectronGoodInd[0]]')
        df[syst][s] = df[syst][s].Define('JetQInd','JetInds(nJet, JetGoodInd, Jet_pt_aux, Jet_eta, Jet_phi, JetBotInd, Jet_qgl)')
        ########### Filtering and further definitions
        df[syst][s] = df[syst][s].Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>=4')
        df[syst][s] = df[syst][s].Filter('JetQInd.size() > 1')
        df[syst][s] = df[syst][s].Define('sl_selection_aux','SLselection(JetQInd, Jet_pt, Jet_eta, Jet_phi, Muon_pt, Muon_eta, Muon_phi, MuonGoodInd, Muon_tightId)')
        df[syst][s] = df[syst][s].Define('sl_bool','sl_selection_aux[0] > -1')
        df[syst][s] = df[syst][s].Define('MuonSLInd','sl_selection_aux[0]')
        df[syst][s] = df[syst][s].Define('Jetauxmuon','sl_selection_aux[1]')
        df[syst][s] = df[syst][s].Define('JetnotMuonInd','Jetauxmuon == JetQInd[0] ? JetQInd[1] : JetQInd[0]')
        df[syst][s] = df[syst][s].Define('JetSLInd','vector2mu(Jetauxmuon, JetnotMuonInd)')
        df[syst][s] = df[syst][s].Define('nJetSLInd','JetSLInd.size()')
        df[syst][s] = df[syst][s].Define('JetXInd','JetQInd')
        ### hists definitions
        df[syst][s] = df[syst][s].Define('jet_1_pt','Jet_pt_aux[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('jet_1_nmu','Jet_nMuons[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('jet_1_qgl','Jet_qgl[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('jet_2_pt','nJetGood>1 ? Jet_pt_aux[JetXInd[1]] : 0')
        df[syst][s] = df[syst][s].Define('jet_1_eta','Jet_eta[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('jet_2_eta','nJetGood>1 ? Jet_eta[JetXInd[1]] : 0')
        df[syst][s] = df[syst][s].Define('jet_1_phi','Jet_phi[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('jet_2_phi','nJetGood>1 ? Jet_phi[JetXInd[1]] : 0')
        df[syst][s] = df[syst][s].Define('jet_1_mass','Jet_mass_nom[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('jet_2_mass','nJetGood>1 ? Jet_mass_nom[JetXInd[1]] : 0')
        df[syst][s] = df[syst][s].Define('jet_2_qgl','nJetGood>1 ? Jet_qgl[JetXInd[1]] : 0')
        df[syst][s] = df[syst][s].Define('jet_2_nmu','nJetGood>1 ? Jet_nMuons[JetXInd[1]] : 0')
        df[syst][s] = df[syst][s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt_aux[JetXInd[0]],Jet_eta[JetXInd[0]],Jet_phi[JetXInd[0]],0.,Jet_pt_aux[JetXInd[1]],Jet_eta[JetXInd[1]],Jet_phi[JetXInd[1]],0.) : 0')
        df[syst][s] = df[syst][s].Define('deltaR_jet1_jet2','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetXInd[1]], Jet_eta[JetXInd[0]] , Jet_phi[JetXInd[1]], Jet_phi[JetXInd[0]])  : 10')
        df[syst][s] = df[syst][s].Define('deltaphi_jet1_jet2','mydeltaphi(Jet_phi[JetXInd[0]],Jet_phi[JetXInd[1]])')
        df[syst][s] = df[syst][s].Define('deltaeta_jet1_jet2','fabs(Jet_eta[JetXInd[0]]-Jet_eta[JetXInd[1]])')
        df[syst][s] = df[syst][s].Define('deltapt_jet1_jet2','fabs(Jet_pt_aux[JetXInd[0]]-Jet_pt_aux[JetXInd[1]])')
        df[syst][s] = df[syst][s].Define('tracks_jet1','Jet_nConstituents[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('tracks_jet2','nJetGood>1 ? Jet_nConstituents[JetXInd[1]] : 0')
        df[syst][s] = df[syst][s].Define('EMN_jet1','Jet_neEmEF[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('EMC_jet1','Jet_chEmEF[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('EMtotal_jet1','Jet_chEmEF[JetXInd[0]]+Jet_neEmEF[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('pT_sum','pTsum(MuonGoodInd, ElectronGoodInd, JetXInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt_aux, MET_phi_aux, Jet_pt_aux, Jet_eta, Jet_phi, Jet_mass_nom)')
        df[syst][s] = df[syst][s].Define('pT_product','pTprod(MuonGoodInd, ElectronGoodInd, JetXInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt_aux, MET_phi_aux, Jet_pt_aux, Jet_eta, Jet_phi, Jet_mass_nom)')
        df[syst][s] = df[syst][s].Define('aux_various','variousSUM(MuonGoodInd, ElectronGoodInd, JetXInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt_aux, MET_phi_aux, Jet_pt_aux, Jet_eta, Jet_phi, Jet_mass_nom)')
        df[syst][s] = df[syst][s].Define('deltaR_lep_2jets','aux_various[0]')
        df[syst][s] = df[syst][s].Define('deltaphi_MET_2jets','aux_various[1]')
        df[syst][s] = df[syst][s].Define('deltaphi_MET_jets_1','aux_various[20]')
        df[syst][s] = df[syst][s].Define('deltaphi_MET_jets_2','aux_various[21]')
        df[syst][s] = df[syst][s].Define('eta_2jets','aux_various[2]')
        df[syst][s] = df[syst][s].Define('pt_2jets','aux_various[3]')
        df[syst][s] = df[syst][s].Define('deltaphi_lephad','aux_various[4]')
        df[syst][s] = df[syst][s].Define('deltaR_lephad','aux_various[5]')
        df[syst][s] = df[syst][s].Define('deltaphi_lep_2jets','aux_various[6]')
        df[syst][s] = df[syst][s].Define('deltaeta_lephad','aux_various[7]')
        df[syst][s] = df[syst][s].Define('deltaeta_lep_2jets','aux_various[8]')
        df[syst][s] = df[syst][s].Define('deltapt_lephad','aux_various[9]')
        df[syst][s] = df[syst][s].Define('deltapt_lep_2jets','aux_various[10]')
        df[syst][s] = df[syst][s].Define('pT_proy','aux_various[11]')
        df[syst][s] = df[syst][s].Define('pT_sum_2J','aux_various[12]')
        df[syst][s] = df[syst][s].Define('pT_Wlep','aux_various[13]')
        df[syst][s] = df[syst][s].Define('deltaR_lep_jet1','aux_various[14]')
        df[syst][s] = df[syst][s].Define('deltaR_lep_jet2','aux_various[15]')
        df[syst][s] = df[syst][s].Define('deltaPhi_lep_jet1','aux_various[16]')
        df[syst][s] = df[syst][s].Define('deltaPhi_lep_jet2','aux_various[17]')
        df[syst][s] = df[syst][s].Define('deltaEta_lep_jet1','aux_various[18]')
        df[syst][s] = df[syst][s].Define('deltaEta_lep_jet2','aux_various[19]')
        df[syst][s] = df[syst][s].Define('jet_bot1_btag','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[0]]]')
        df[syst][s] = df[syst][s].Define('jet_bot2_btag','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[1]]]')
        df[syst][s] = df[syst][s].Define('jet_bot1_bcorr','Jet_bRegCorr[JetGoodInd[JetBotInd[0]]]')
        df[syst][s] = df[syst][s].Define('jet_bot2_bcorr','Jet_bRegCorr[JetGoodInd[JetBotInd[1]]]')
        df[syst][s] = df[syst][s].Define('jet_bot1_pt','Jet_pt_aux[JetGoodInd[JetBotInd[0]]]')
        df[syst][s] = df[syst][s].Define('jet_bot2_pt','Jet_pt_aux[JetGoodInd[JetBotInd[1]]]')
        df[syst][s] = df[syst][s].Define('jet_bot1_eta','Jet_eta[JetGoodInd[JetBotInd[0]]]')
        df[syst][s] = df[syst][s].Define('jet_bot2_eta','Jet_eta[JetGoodInd[JetBotInd[1]]]')
        df[syst][s] = df[syst][s].Define('jet_bot1_phi','Jet_phi[JetGoodInd[JetBotInd[0]]]')
        df[syst][s] = df[syst][s].Define('jet_bot2_phi','Jet_phi[JetGoodInd[JetBotInd[1]]]')
        df[syst][s] = df[syst][s].Define('jet_1_btag','Jet_btagDeepFlavB[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('jet_2_btag','Jet_btagDeepFlavB[JetXInd[1]]')
        df[syst][s] = df[syst][s].Define('jet_bot1_tracks','Jet_nConstituents[JetGoodInd[JetBotInd[0]]]')
        df[syst][s] = df[syst][s].Define('jet_bot2_tracks','Jet_nConstituents[JetGoodInd[JetBotInd[1]]]')
        df[syst][s] = df[syst][s].Define('nLooseLepton','nMuon+nElectron-1')
        df[syst][s] = df[syst][s].Define('InvM3_aux','InvMassBot(JetBotInd, JetQInd, JetGoodInd, Jet_pt_aux, Jet_eta, Jet_phi)')
        df[syst][s] = df[syst][s].Define('InvM_various_aux','InvMassVarious(JetBotInd, JetQInd, JetGoodInd, Jet_pt_aux, Jet_eta, Jet_phi, MuonGoodInd, ElectronGoodInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass)')
        df[syst][s] = df[syst][s].Define('InvM_bot_closer','InvM3_aux[0]')
        df[syst][s] = df[syst][s].Define('InvM_bot_farther','InvM3_aux[1]')
        df[syst][s] = df[syst][s].Define('InvM30','InvM_various_aux[0]')
        df[syst][s] = df[syst][s].Define('InvM31','InvM_various_aux[1]')
        df[syst][s] = df[syst][s].Define('InvMl0','InvM_various_aux[2]')
        df[syst][s] = df[syst][s].Define('InvMl1','InvM_various_aux[3]')
        df[syst][s] = df[syst][s].Define('jet_bot1_btagnumber','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[0]]] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetGoodInd[JetBotInd[0]]] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')
        df[syst][s] = df[syst][s].Define('jet_bot2_btagnumber','Jet_btagDeepFlavB[JetGoodInd[JetBotInd[1]]] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetGoodInd[JetBotInd[1]]] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')
        df[syst][s] = df[syst][s].Define('jet_1_btagnumber','Jet_btagDeepFlavB[JetXInd[0]] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetXInd[0]] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')
        df[syst][s] = df[syst][s].Define('jet_2_btagnumber','Jet_btagDeepFlavB[JetXInd[1]] < '+str(cuts_btag[years[0]][0])+' ? 0 : (Jet_btagDeepFlavB[JetXInd[1]] < '+str(cuts_btag[years[0]][1])+' ? 1 : 2)')
        df[syst][s] = df[syst][s].Define('jet_1_cvltag','Jet_btagDeepFlavCvL[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('jet_2_cvltag','Jet_btagDeepFlavCvL[JetXInd[1]]')
        df[syst][s] = df[syst][s].Define('jet_1_cvltag_csv','Jet_btagDeepCvL[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('jet_2_cvltag_csv','Jet_btagDeepCvL[JetXInd[1]]')
        df[syst][s] = df[syst][s].Define('jet_1_cvbtag','Jet_btagDeepFlavCvB[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('jet_2_cvbtag','Jet_btagDeepFlavCvB[JetXInd[1]]')
        df[syst][s] = df[syst][s].Define('jet_1_cvbtag_csv','Jet_btagDeepCvB[JetXInd[0]]')
        df[syst][s] = df[syst][s].Define('jet_2_cvbtag_csv','Jet_btagDeepCvB[JetXInd[1]]')
        df[syst][s] = df[syst][s].Define('jet_max_ctag_indx','jet_1_cvltag > jet_2_cvltag ? JetXInd[0] : JetXInd[1]')
        df[syst][s] = df[syst][s].Define('jet_max_cvltag','jet_1_cvltag > jet_2_cvltag ? jet_1_cvltag : jet_2_cvltag')
        df[syst][s] = df[syst][s].Define('jet_min_cvltag','jet_1_cvltag > jet_2_cvltag ? jet_2_cvltag : jet_1_cvltag')
        df[syst][s] = df[syst][s].Define('jet_max_cvbtag','jet_1_cvbtag > jet_2_cvbtag ? jet_1_cvbtag : jet_2_cvbtag')
        df[syst][s] = df[syst][s].Define('jet_min_cvbtag','jet_1_cvbtag > jet_2_cvbtag ? jet_2_cvbtag : jet_1_cvbtag')
        df[syst][s] = df[syst][s].Define('jet_int_cvltag','Jet_btagDeepFlavCvL[jet_max_ctag_indx]')
        df[syst][s] = df[syst][s].Define('jet_int_cvbtag','Jet_btagDeepFlavCvB[jet_max_ctag_indx]')
        df[syst][s] = df[syst][s].Define('chi2_test0','chi2calc(InvM_2jets, InvM30)')
        df[syst][s] = df[syst][s].Define('chi2_test1','chi2calc(InvM_2jets, InvM31)')
        df[syst][s] = df[syst][s].Define('nJetMuonInd','JetMuonInd.size()')
        df[syst][s] = df[syst][s].Define('nJetSVInd','JetSVInd.size()')
        df[syst][s] = df[syst][s].Define('nJetBotInd','JetBotInd.size()')
        df[syst][s] = df[syst][s].Define('nMuonJetInd','MuonJetInd.size()')
        df[syst][s] = df[syst][s].Define('nSVJetInd','SVJetInd.size()')
        df[syst][s] = df[syst][s].Define('nJetQInd','JetQInd.size()')
        df[syst][s] = df[syst][s].Define('kfmasses_aux','kinfitmasses(InvM_2jets, InvM30, InvM31, InvMl0, InvMl1)')
        df[syst][s] = df[syst][s].Define('InvMl_good','kfmasses_aux[0]')
        df[syst][s] = df[syst][s].Define('InvMl_bad','kfmasses_aux[1]')
        df[syst][s] = df[syst][s].Define('InvM3_good','kfmasses_aux[2]')
        df[syst][s] = df[syst][s].Define('InvM3_bad','kfmasses_aux[3]')
        df[syst][s] = df[syst][s].Define('chi2_test_good','kfmasses_aux[4]')
        df[syst][s] = df[syst][s].Define('chi2_test_bad','kfmasses_aux[5]')
        df[syst][s] = df[syst][s].Define('tau_indx','tau_selection(JetGoodInd,JetBotInd,JetXInd, Jet_pt, Jet_eta, Jet_phi, Tau_pt, Tau_eta, Tau_phi)')
        df[syst][s] = df[syst][s].Define('tau_discr_jet1','tau_indx[0] > -1 ? Tau_rawDeepTau2017v2p1VSjet[tau_indx[0]] : -1')
        df[syst][s] = df[syst][s].Define('tau_discr_jet2','tau_indx[1] > -1 ? Tau_rawDeepTau2017v2p1VSjet[tau_indx[1]] : -1')
        df[syst][s] = df[syst][s].Define('tau_discr_jetbot1','tau_indx[2] > -1 ? Tau_rawDeepTau2017v2p1VSjet[tau_indx[2]] : -1')
        df[syst][s] = df[syst][s].Define('tau_discr_jetbot2','tau_indx[3] > -1 ? Tau_rawDeepTau2017v2p1VSjet[tau_indx[3]] : -1')
        df[syst][s] = df[syst][s].Define('deltaR_jet1_tau','tau_indx[0] > -1 ? ROOT::VecOps::DeltaR(Tau_eta[tau_indx[0]],Jet_eta[JetXInd[0]],Tau_phi[tau_indx[0]],Jet_phi[JetXInd[0]]) : -1.')
        df[syst][s] = df[syst][s].Define('bot_mu_aux','bottomMuons(JetGoodInd[JetBotInd[0]], JetGoodInd[JetBotInd[1]] ,Jet_pt, Jet_eta, Jet_phi, Muon_pt, Muon_eta, Muon_phi, MuonGoodInd, Muon_tightId)')
        df[syst][s] = df[syst][s].Define('bot1_muons','bot_mu_aux[0]')
        df[syst][s] = df[syst][s].Define('bot2_muons','bot_mu_aux[2]')
        df[syst][s] = df[syst][s].Define('bot1_muon_id','bot_mu_aux[1]')
        df[syst][s] = df[syst][s].Define('bot2_muon_id','bot_mu_aux[3]')
        df[syst][s] = df[syst][s].Define('muon_bot1_eta','bot1_muon_id > -1 ? Muon_eta[bot1_muon_id] : -10')
        df[syst][s] = df[syst][s].Define('muon_bot2_eta','bot2_muon_id > -1 ? Muon_eta[bot2_muon_id] : -10')
        df[syst][s] = df[syst][s].Define('muon_bot1_pt','bot1_muon_id > -1 ? Muon_pt[bot1_muon_id] : -10')
        df[syst][s] = df[syst][s].Define('muon_bot2_pt','bot2_muon_id > -1 ? Muon_pt[bot2_muon_id] : -10')
        df[syst][s] = df[syst][s].Define('muon_bot1_relpt','jet_bot1_pt > 0. ? muon_bot1_pt/jet_bot1_pt : -10.')
        df[syst][s] = df[syst][s].Define('muon_bot2_relpt','jet_bot2_pt > 0. ? muon_bot2_pt/jet_bot2_pt : -10.')
        df[syst][s] = df[syst][s].Define('muon_bot1_iso','bot1_muon_id > -1 ? Muon_pfRelIso04_all[bot1_muon_id] : -10')
        df[syst][s] = df[syst][s].Define('muon_bot2_iso','bot2_muon_id > -1 ? Muon_pfRelIso04_all[bot2_muon_id] : -10')

############################################################
#### Distinguishing between same sign and opposite sign ####
############################################################

#################### Final definitions and filters ###############################

for s in samples:
    for syst in list_jetsyst:
        df[syst][s] = df[syst][s].Define('lepton_phi_aux','nMuonGood>0 ? Muon_phi[MuonGoodInd[0]] : Electron_phi[ElectronGoodInd[0]]')
        df[syst][s] = df[syst][s].Define('transverse_mass_old','std::sqrt(2*lepton_pt*MET_pt_aux*(1-std::cos(lepton_phi_aux-MET_phi_aux)))')
        df[syst][s] = df[syst][s].Define('weight_aux','1.')
        df[syst][s] = df[syst][s].Define('deltaR_jetM_lep','nMuonGood>0 ? ROOT::VecOps::DeltaR(Muon_eta[MuonGoodInd[0]],Jet_eta[JetQInd[0]] , Muon_phi[MuonGoodInd[0]], Jet_phi[JetQInd[0]]) : ROOT::VecOps::DeltaR(Electron_eta[ElectronGoodInd[0]],Jet_eta[JetQInd[0]] , Electron_phi[ElectronGoodInd[0]], Jet_phi[JetQInd[0]])')
        df[syst][s] = df[syst][s].Define('InvM_jetM_lep', 'nMuonGood>0 ? InvariantM(Jet_pt_aux[JetXInd[0]],Jet_eta[JetXInd[0]],Jet_phi[JetXInd[0]],0.,Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]]) : InvariantM(Jet_pt_aux[JetXInd[0]],Jet_eta[JetXInd[0]],Jet_phi[JetXInd[0]],0.,Electron_pt[ElectronGoodInd[0]],Electron_eta[ElectronGoodInd[0]],Electron_phi[ElectronGoodInd[0]], Electron_mass[ElectronGoodInd[0]])')
        df[syst][s] = df[syst][s].Define('InvM_muon_jet','nMuonGood>0 ? InvariantM(Muon_pt[MuonJetInd[0]],Muon_eta[MuonJetInd[0]],Muon_phi[MuonJetInd[0]],Muon_mass[MuonJetInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]]) : 50.')
        df[syst][s] = df[syst][s].Define('deltaphi_MET_lep','mydeltaphi(lepton_phi_aux,MET_phi_aux)')
        df[syst][s] = df[syst][s].Define('MET_my_significance','(MET_pt_aux*MET_pt_aux)/((0.62*0.62)*MET_sumEt)')
        ### Filters
        df[syst][s] = df[syst][s].Filter('transverse_mass > 20.')
        df[syst][s] = df[syst][s].Filter('MET_pt_aux > 20.')
        df[syst][s] = df[syst][s].Filter('nMuonGood>0 ? lepton_pt > 30. : lepton_pt > 35.')
        ### Jet pT filters
        #df[syst][s] = df[syst][s].Filter('jet_bot1_pt > 25.')
        #df[syst][s] = df[syst][s].Filter('jet_bot2_pt > 25.')
        #df[syst][s] = df[syst][s].Filter('jet_1_pt > 25.')
        #df[syst][s] = df[syst][s].Filter('jet_2_pt > 25.')
        ## Invariant mass filter
        df[syst][s] = df[syst][s].Filter('InvM_2jets > 30. && InvM_2jets < 300.')
        ### Cur based for electron, EL ID tests
        #df[syst][s] = df[syst][s].Filter('nElectronGood > 0 ? Electron_cutBased[ElectronGoodInd[0]] > 3 : 1')
        ### Chi2 test
        df[syst][s] = df[syst][s].Define('chi_bool','kinfittest(InvM_2jets, InvM30, InvM31, InvMl0, InvMl1)')
        if chan == "btagMM":
           df[syst][s] = df[syst][s].Filter('jet_bot1_btag >'+str(cuts_btag[years[0]][1]))
           df[syst][s] = df[syst][s].Filter('jet_bot2_btag >'+str(cuts_btag[years[0]][1]))
        elif chan == "lepton50":
           df[syst][s] = df[syst][s].Filter('jet_bot1_btag >'+str(cuts_btag[years[0]][1]))
           df[syst][s] = df[syst][s].Filter('jet_bot2_btag >'+str(cuts_btag[years[0]][0]))
           df[syst][s] = df[syst][s].Filter('lepton_pt > 50.')
        elif chan == "btagMM_chitest":
           df[syst][s] = df[syst][s].Filter('jet_bot1_btag >'+str(cuts_btag[years[0]][1]))
           df[syst][s] = df[syst][s].Filter('jet_bot2_btag >'+str(cuts_btag[years[0]][1]))
           df[syst][s] = df[syst][s].Filter('kinfittest(InvM_2jets, InvM30, InvM31, InvMl0, InvMl1)')
        elif chan == "lepton50_chitest":
           df[syst][s] = df[syst][s].Filter('jet_bot1_btag >'+str(cuts_btag[years[0]][1]))
           df[syst][s] = df[syst][s].Filter('jet_bot2_btag >'+str(cuts_btag[years[0]][0]))
           df[syst][s] = df[syst][s].Filter('lepton_pt > 50.')
           df[syst][s] = df[syst][s].Filter('kinfittest(InvM_2jets, InvM30, InvM31, InvMl0, InvMl1)')
        #### Charm cuts
        if channel == "ctag":
           df[syst][s] = df[syst][s].Filter('(jet_int_cvltag >'+str(cuts_cvl_deepj[years[0]][2])+') && (jet_int_cvbtag >'+str(cuts_cvb_deepj[years[0]][2])+')')
        if channel == "noctag":
           df[syst][s] = df[syst][s].Filter('!((jet_int_cvltag >'+str(cuts_cvl_deepj[years[0]][2])+') && (jet_int_cvbtag >'+str(cuts_cvb_deepj[years[0]][2])+'))')

########### XSEC vars ############

if mode == "mc":
   for s in samples:
     for syst in list_jetsyst:
        if years[0]=="2016B":
               df[syst][s] = df[syst][s].Define('var_xsec',str(xsecs[years[0]][s[:-5]]))
               df[syst][s] = df[syst][s].Define('var_lumi',str(lumi[years[0]][s[:-5]]))
               df[syst][s] = df[syst][s].Define('lumi_data',str(lumi_d[years[0]]))
               df[syst][s] = df[syst][s].Define('weight_lumi','lumi_data/var_lumi')
        else:
               df[syst][s] = df[syst][s].Define('var_xsec',str(xsecs[years[0]][s[:-4]]))
               df[syst][s] = df[syst][s].Define('var_lumi',str(lumi[years[0]][s[:-4]]))
               df[syst][s] = df[syst][s].Define('lumi_data',str(lumi_d[years[0]]))
               df[syst][s] = df[syst][s].Define('weight_lumi','lumi_data/var_lumi')
        df[syst][s] = df[syst][s].Define('wcs_var_auxx','ttbarcharm(nGenPart,GenPart_statusFlags,GenPart_pdgId,GenPart_genPartIdxMother)')
        df[syst][s] = df[syst][s].Define('wcs_var','wcs_var_auxx[0]')

############ Trigger scale factors ##############

from trigger_sf import *

if mode == "mc":
   for syst in list_jetsyst:
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= trigger_mu_sfRDF(**kwargs)
                df[syst][s] = b().run(df[syst][s])

from trigger_sf_el import *

if mode == "mc":
   for syst in list_jetsyst:
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= trigger_el_sfRDF(**kwargs)
                df[syst][s] = b().run(df[syst][s])

####### Other scale factors

from low_pt_muonid import *

if mode == "mc":
   for syst in list_jetsyst:
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= displaced_mu_idRDF(**kwargs)
                df[syst][s] = b().run(df[syst][s])

from complete_btag_weight import *

if mode == "mc":
   for syst in list_jetsyst:
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= btag_weights_totRDF(**kwargs)
                df[syst][s] = b().run(df[syst][s])

############ BR corrections ##############

from branchingfractions_corrections import *

if mode == "mc":
   for syst in list_jetsyst:
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= branchingfractions_SL_corRDF(**kwargs)
                df[syst][s] = b().run(df[syst][s])

############ Muon in bot jet SF ##############

from muon_in_bot_sf import *

if (mode=="mc" and channel[0:2]=="sl"):
   for syst in list_jetsyst:
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True, "canal":"selection_sl"}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True, "canal":"selection_sl"}
                       #print(kwargs)
                b= muon_frombot_sfRDF(**kwargs)
                df[syst][s] = b().run(df[syst][s])

########### ctag SF

from ctag_sf import *

if mode == "mc":
   for syst in list_jetsyst:
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= ctag_weights_totRDF(**kwargs)
                df[syst][s] = b().run(df[syst][s])

from vetoctag_sf import *

if mode == "mc":
   for syst in list_jetsyst:
        for s in samples:
                if s[-1]=="B":
                       kwargs = {"year":s[-5:],"isMC":True, "isUL":True}
                else:
                       kwargs = {"year":s[-4:],"isMC":True, "isUL":True}
                       #print(kwargs)
                b= vetoctag_weights_totRDF(**kwargs)
                df[syst][s] = b().run(df[syst][s])

############ Gen level definitions

if mode == "mc":
   for syst in list_jetsyst:
        for s in samples:
                df[syst][s] = df[syst][s].Define('jet_1_flavourH','Jet_hadronFlavour[JetXInd[0]]')
                df[syst][s] = df[syst][s].Define('jet_2_flavourH','Jet_hadronFlavour[JetXInd[1]]')
                df[syst][s] = df[syst][s].Define('jet_1_flavourP','Jet_partonFlavour[JetXInd[0]]')
                df[syst][s] = df[syst][s].Define('jet_2_flavourP','Jet_partonFlavour[JetXInd[1]]')
                df[syst][s] = df[syst][s].Define('jet_bot1_flavourP','Jet_partonFlavour[JetGoodInd[JetBotInd[0]]]')
                df[syst][s] = df[syst][s].Define('jet_bot2_flavourP','Jet_partonFlavour[JetGoodInd[JetBotInd[1]]]')
                df[syst][s] = df[syst][s].Define('last_Copy','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,13)')
                df[syst][s] = df[syst][s].Define('top_weight','topreweight(nGenPart, GenPart_statusFlags, GenPart_pdgId, last_Copy, GenPart_pt)')
                if chan == "nobtag": 
                   df[syst][s] = df[syst][s].Define('btag_sf','1.')
                else:               
                   df[syst][s] = df[syst][s].Define('btag_sf','Aux_btag_weight[1]/Aux_btag_weight[0]')
                df[syst][s] = df[syst][s].Define('btag_sf_light_up','Aux_btag_weight_light_up[1]/Aux_btag_weight_light_up[0]')
                df[syst][s] = df[syst][s].Define('btag_sf_heavy_up','Aux_btag_weight_heavy_up[1]/Aux_btag_weight_heavy_up[0]')
                df[syst][s] = df[syst][s].Define('btag_sf_light_down','Aux_btag_weight_light_down[1]/Aux_btag_weight_light_down[0]')
                df[syst][s] = df[syst][s].Define('btag_sf_heavy_down','Aux_btag_weight_heavy_down[1]/Aux_btag_weight_heavy_down[0]')
                df[syst][s] = df[syst][s].Define('muon_bot1_mother','bot1_muon_id > -1 ? fabs(GenPart_pdgId[GenPart_genPartIdxMother[Muon_genPartIdx[bot1_muon_id]]]) : -10')
                df[syst][s] = df[syst][s].Define('muon_bot2_mother','bot2_muon_id > -1 ? fabs(GenPart_pdgId[GenPart_genPartIdxMother[Muon_genPartIdx[bot2_muon_id]]]) : -10')
                df[syst][s] = df[syst][s].Define('muon_bot1_mother_mine','bot1_muon_id > -1 ? muonmother(GenPart_pdgId, GenPart_genPartIdxMother, Muon_genPartIdx, bot1_muon_id) : -10')
                df[syst][s] = df[syst][s].Define('muon_bot2_mother_mine','bot2_muon_id > -1 ? muonmother(GenPart_pdgId, GenPart_genPartIdxMother, Muon_genPartIdx, bot2_muon_id) : -10')
                if (channel[0:2]=="sl"): df[syst][s] = df[syst][s].Define('muon_jet_genindx','Muon_genPartIdx[MuonSLInd]')
                if (channel[0:2]=="sl"): df[syst][s] = df[syst][s].Define('muon_jet_mother','muonmother_ttdilep(GenPart_pdgId,GenPart_genPartIdxMother,Muon_genPartIdx,MuonSLInd)')

########## ttbar sectioning for charm discrimination

cond1 = '(fabs(jet_1_flavourP) == 4 || fabs(jet_2_flavourP) == 4)'
cond2 = '(fabs(jet_1_flavourP) == 5 || fabs(jet_2_flavourP) == 5)'
cond3 = '(fabs(jet_1_flavourP) == 3 || fabs(jet_1_flavourP) == 2 || fabs(jet_1_flavourP) == 1)'
cond4 =	'(fabs(jet_2_flavourP) == 3 || fabs(jet_2_flavourP) == 2 || fabs(jet_2_flavourP) == 1)'
cond5 = '(fabs(jet_1_flavourP)>0 && fabs(jet_1_flavourP)<6 && fabs(jet_2_flavourP)>0 && fabs(jet_2_flavourP)<6)'
cond6 = '((fabs(jet_1_flavourP)==0 || fabs(jet_1_flavourP)>6) && (fabs(jet_2_flavourP)==0 || fabs(jet_2_flavourP)>6))'

if args.wcs:
   ### alternative distinction
   if mode == "mc":
      for s in samples:
            if (s[0:8]=="ttbar_sl" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                  for syst in list_jetsyst:
                           df[syst][s] = df[syst][s].Define('charm_light_aux','ttbarcharm(nGenPart,GenPart_statusFlags,GenPart_pdgId,GenPart_genPartIdxMother)')
                           df[syst][s+"_charm"] = df[syst][s].Filter('charm_light_aux[0]')
                           df[syst][s+"_nocharm"] = df[syst][s].Filter('!charm_light_aux[0]')
                           ## Samples correction
                           if (channel[0:2]=="sl"):
                              df[syst][s] = df[syst][s].Define('muondhad_bool','muonmotherdhad(GenPart_pdgId,GenPart_genPartIdxMother,Muon_genPartIdx,MuonSLInd)')
                              df[syst][s+"_charm_mudhad"] = df[syst][s].Filter('charm_light_aux[0] && muondhad_bool')
                  samples.append(s+"_charm")
                  samples.append(s+"_nocharm")
                  if (channel[0:2]=="sl"): samples.append(s+"_charm_mudhad")

   samples = [s for s in samples if not (s[0:8]=="ttbar_sl"  and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]
   ######## ST sectioning for charm discrimination
   if mode == "mc":
       for s in samples:
           if (s[0]+s[1] == "st" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                   for syst in list_jetsyst:
                           df[syst][s] = df[syst][s].Define('charm_light_aux','ttbarcharm(nGenPart,GenPart_statusFlags,GenPart_pdgId,GenPart_genPartIdxMother)')
                           df[syst][s+"_charm"] = df[syst][s].Filter('charm_light_aux[0]')
                           df[syst][s+"_nocharm"] = df[syst][s].Filter('!charm_light_aux[0] && charm_light_aux[1]')
                           df[syst][s+"_else"] = df[syst][s].Filter('!charm_light_aux[0] && !charm_light_aux[1]')
                   ## Samples correction
                   samples.append(s+"_charm")
                   samples.append(s+"_nocharm")
                   samples.append(s+"_else")

   samples = [s for s in samples if not (s[0]+s[1] == "st"  and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]
else:
   if mode == "mc":
      for s in samples:
          if (s[0:8] == "ttbar_sl" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                for syst in list_jetsyst:
                           df[syst][s+"_charm"] = df[syst][s].Filter(cond1+' && '+cond5)
                           df[syst][s+"_bottom"] = df[syst][s].Filter(cond2+' && !'+cond1+' && '+cond5)
                           df[syst][s+"_light"] = df[syst][s].Filter('('+cond3+') && ('+cond4+') && !('+cond1+') && !('+cond2+')')
                           df[syst][s+"_charmgluon"] = df[syst][s].Filter(cond1+' && !'+cond5)
                           df[syst][s+"_bottomgluon"] = df[syst][s].Filter(cond2+' && !'+cond1+' && !'+cond5)
                           #df[syst][s+"_gluongluon"] = df[syst][s].Filter(cond6)
                           #df[syst][s+"_else"] = df[syst][s].Filter('!('+cond1+') && !('+cond2+') && !(('+cond3+') && ('+cond4+') && !('+cond6+')')
                           df[syst][s+"_else"] = df[syst][s].Filter('!('+cond1+') && !('+cond2+') && !(('+cond3+') && ('+cond4+'))')
                ## Samples correction
                samples.append(s+"_charm")
                samples.append(s+"_bottom")
                samples.append(s+"_light")
                samples.append(s+"_charmgluon")
                samples.append(s+"_bottomgluon")
                #samples.append(s+"_gluongluon")
                samples.append(s+"_else")

   samples = [s for s in samples if not (s[0:8]=="ttbar_sl"  and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]


################# channel partitions

df_M = {}
df_E = {}
for syst in list_jetsyst:
    df_M[syst] = {}
    df_E[syst] = {}

if args.type == "mc":
    for s in samples:
        for syst in list_jetsyst:
               if years[0] == "2017":
                     df[syst][s] = df[syst][s].Define('l1_prefw','L1PreFiringWeight_Nom')
               else:
                     df[syst][s] = df[syst][s].Define('l1_prefw','1.')
               df[syst][s] = df[syst][s].Define('lep_id_sf','nMuonGood>0 ? musf_tight_id[MuonGoodInd[0]] : elesf_wp80iso[ElectronGoodInd[0]]')
               df[syst][s] = df[syst][s].Define('lep_id_sf_up','nMuonGood>0 ? musf_tight_id_up[MuonGoodInd[0]] : elesf_wp80iso_up[ElectronGoodInd[0]]')
               df[syst][s] = df[syst][s].Define('lep_id_sf_down','nMuonGood>0 ? musf_tight_id_down[MuonGoodInd[0]] : elesf_wp80iso_down[ElectronGoodInd[0]]')
               df[syst][s] = df[syst][s].Define('lep_iso_sf','nMuonGood>0 ? musf_tight_reliso[MuonGoodInd[0]] : 1.')
               df[syst][s] = df[syst][s].Define('lep_iso_sf_up','nMuonGood>0 ? musf_tight_reliso_up[MuonGoodInd[0]] : 1.')
               df[syst][s] = df[syst][s].Define('lep_iso_sf_down','nMuonGood>0 ? musf_tight_reliso_down[MuonGoodInd[0]] : 1.')
               df[syst][s] = df[syst][s].Define('lep_trig_sf','nMuonGood>0 ? trigger_sf_mu_aux[MuonGoodInd[0]] : trigger_sf_el_aux[ElectronGoodInd[0]]')
               df[syst][s] = df[syst][s].Define('lep_trig_sf_aux','nMuonGood>0 ? std::sqrt(trigger_sf_mu_aux_syst[MuonGoodInd[0]]*trigger_sf_mu_aux_syst[MuonGoodInd[0]]+trigger_sf_mu_aux_stat[MuonGoodInd[0]]*trigger_sf_mu_aux_stat[MuonGoodInd[0]]): trigger_sf_el_aux_err[ElectronGoodInd[0]]')
               df[syst][s] = df[syst][s].Define('lep_trig_sf_up','lep_trig_sf+lep_trig_sf_aux')
               df[syst][s] = df[syst][s].Define('lep_trig_sf_down','lep_trig_sf-lep_trig_sf_aux')
               df[syst][s] = df[syst][s].Define('lep_id_lowpt_sf','1.')
               df[syst][s] = df[syst][s].Define('lep_id_lowpt_sf_up','1.')
               df[syst][s] = df[syst][s].Define('lep_id_lowpt_sf_down','1.')
               df[syst][s] = df[syst][s].Define('frag_weight','1.')
               if channel == "ctag":
                  #df[syst][s] = df[syst][s].Define('ctag_cmistag','ctag_weight*cmistag_weight')
                  df[syst][s] = df[syst][s].Define('ctag_cmistag','1.')
               else:
                  #df[syst][s] = df[syst][s].Define('ctag_cmistag','vetoctag_weight')
                  df[syst][s] = df[syst][s].Define('ctag_cmistag','1.')
               df[syst][s] = df[syst][s].Define('ssos_weight','1.')
               df[syst][s] = df[syst][s].Define('weightSSOS_final','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_seclepup','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf_up*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_seclepdown','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf_down*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_btaglightup','weight_aux*btag_sf_light_up*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_btaglightdown','weight_aux*btag_sf_light_down*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_btagheavyup','weight_aux*btag_sf_heavy_up*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_btagheavydown','weight_aux*btag_sf_heavy_down*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_lepidup','weight_aux*btag_sf*lep_id_sf_up*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_lepiddown','weight_aux*btag_sf*lep_id_sf_down*lep_iso_sf*lep_trig_sf*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_lepisoup','weight_aux*btag_sf*lep_id_sf*lep_iso_sf_up*lep_trig_sf*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_lepisodown','weight_aux*btag_sf*lep_id_sf*lep_iso_sf_down*lep_trig_sf*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_leptrigup','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*lep_trig_sf_up*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_leptrigdown','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*lep_trig_sf_down*puWeight*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_notoppt','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_puwup','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeightUp*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')
        df["nom"][s] = df["nom"][s].Define('weightSSOS_final_puwdown','weight_aux*btag_sf*lep_id_sf*lep_iso_sf*lep_trig_sf*puWeightDown*top_weight*l1_prefw*lep_id_lowpt_sf*ssos_weight*frag_weight*ctag_cmistag')

else:
     for syst in list_jetsyst:
            if (channel=="sl_ssos"):
                 df[syst][s] = df[syst][s].Define('ssos_weight','muon_ssos<0 ? 1. : -1.')
            else:
                 df[syst][s] = df[syst][s].Define('ssos_weight','1.')
            df[syst][s] = df[syst][s].Define('weightSSOS_final','weight_aux*ssos_weight')

for s in samples:
    for syst in list_jetsyst:
        df_M[syst][s] = df[syst][s].Filter('nMuonGood>0')
        df_E[syst][s] = df[syst][s].Filter('nElectronGood>0')

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
hist_seclepup_M = {}
hist_seclepup_E = {}
hist_seclepdown_M = {}
hist_seclepdown_E = {}
hist_lepidup_M = {}
hist_lepidup_E = {}
hist_lepiddown_M = {}
hist_lepiddown_E = {}
hist_lepisoup_M = {}
hist_lepisoup_E = {}
hist_lepisodown_M = {}
hist_lepisodown_E = {}
hist_leptrigup_M = {}
hist_leptrigup_E = {}
hist_leptrigdown_M = {}
hist_leptrigdown_E = {}
hist_notoppt_M = {}
hist_notoppt_E = {}
hist_puwup_M = {}
hist_puwup_E = {}
hist_puwdown_M = {}
hist_puwdown_E = {}
hist_smearup_M = {}
hist_smearup_E = {}
hist_smeardown_M = {}
hist_smeardown_E = {}
hist_jecup_M = {}
hist_jecup_E = {}
hist_jecdown_M = {}
hist_jecdown_E = {}

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
   "InvMl_good", "InvMl_bad", "InvM3_good", "InvM3_good_short", "InvM3_bad", "chi2_test_good", "chi2_test_bad", "jet_max_cvltag", "jet_min_cvltag",
   "jet_1_cvbtag_csv", "jet_2_cvbtag_csv", "jet_1_cvbtag", "jet_2_cvbtag", "jet_max_cvbtag", "jet_min_cvbtag",
   "deltaR_jet1_tau","jet_1_eta_thick","jet_2_eta_thick","jet_bot1_eta_thick","jet_bot2_eta_thick",
   "InvM_2jets_thick","InvM_2jets_short","bot1_muons","bot2_muons","muon_bot1_eta","muon_bot2_eta","muon_bot1_pt","muon_bot2_pt",
   "jet_bot1_tracks", "jet_bot2_tracks","tau_discr_jet1","tau_discr_jet2","tau_discr_jetbot1","tau_discr_jetbot2"]

if channel[0:2] == "sl": observable_names = observable_names + ["muon_jet_pt","muon_jet_eta","deltaR_jet1_muon","deltaR_muon_tau","muon_jet_z","muon_jet_iso",
          "muon_jet_pt_rel","muon_jet_iso_log","muon_jet_z_short","muon_jet_sigr","muon_jet_sigxy","muon_jet_sigdz","muon_jet_r","muon_jet_z2","muon_jet_z3","muon_jet_iso_abs",
          "muon_jet_z2_v2"]

column_names = {}

for name in observable_names:
   column_names[name] = name

column_names["lepton_pt_detail"] = "lepton_pt"; column_names["jet_1_btag_thick"] = "jet_1_btag"; column_names["jet_2_btag_thick"] = "jet_2_btag";
column_names["jet_bot1_btag_thick"] = "jet_bot1_btag"; column_names["jet_bot2_btag_thick"] = "jet_bot2_btag"; column_names["lepton_eta_thick"] = "lepton_eta";
column_names["MET_sig"] = "MET_significance"; column_names["MET_my_sig"] = "MET_my_significance";
column_names["jet_1_eta_thick"] = "jet_1_eta";column_names["jet_2_eta_thick"] = "jet_2_eta";column_names["jet_bot1_eta_thick"] = "jet_bot1_eta";column_names["jet_bot2_eta_thick"] = "jet_bot2_eta";
column_names["InvM_2jets_thick"] = "InvM_2jets";column_names["InvM_2jets_short"] = "InvM_2jets";column_names["muon_jet_z_short"] = "muon_jet_z";
column_names["InvM3_good_short"] = "InvM3_good";

dict_binlim = {}
dict_binlim["nJetGood"] = [5,4,9]; dict_binlim["jet_2_mass"] = [40,0,40];
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
dict_binlim["deltaphi_MET_lep"] = [70,0,3.5]; dict_binlim["deltaphi_MET_jets_1"] = [70,0,3.5]; dict_binlim["deltaphi_MET_jets_2"] = [70,0,3.5];
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
dict_binlim["muon_jet_pt"]=[30,0,30];dict_binlim["muon_jet_relpt"]=[20,0,1];dict_binlim["muon_jet_eta"]=[18,-2.7,2.7];
dict_binlim["deltaR_jet1_muon"]=[40,0,0.4];dict_binlim["deltaR_jet1_tau"]=[50,0,5];dict_binlim["deltaR_muon_tau"]=[50,0,5];
dict_binlim["Frag_weight_sl"]=[40,0.4,2.4];dict_binlim["Br_weight_sl"]=[40,0.4,2.4];
dict_binlim["InvM_2jets_thick"] = [54,30,300];dict_binlim["InvM_2jets_short"] = [30,50,110];
dict_binlim["bot1_muons"]=[10,0,10];dict_binlim["bot2_muons"]=[10,0,10];dict_binlim["muon_bot1_eta"]=[18,-2.7,2.7];dict_binlim["muon_bot2_eta"]=[18,-2.7,2.7];
dict_binlim["muon_bot1_pt"]=[30,0,30];dict_binlim["muon_bot2_pt"]=[30,0,30];dict_binlim["muon_bot1_relpt"]=[20,0,1];dict_binlim["muon_bot2_relpt"]=[20,0,1];
dict_binlim["muon_bot1_mother"]=[500,300,800];dict_binlim["muon_bot2_mother"]=[500,300,800];
dict_binlim["muon_bot1_mother_mine"]=[10,0,10];dict_binlim["muon_bot2_mother_mine"]=[10,0,10];
dict_binlim["muon_bot1_iso"]=[20,0,1];dict_binlim["muon_bot2_iso"]=[20,0,1];
dict_binlim["muon_jet_iso_log"]=[30,0,3];dict_binlim["muon_jet_iso"]=[35,0,14];dict_binlim["muon_jet_pt_rel"]=[50,0,3];dict_binlim["muon_jet_z"]=[40,0,1];
dict_binlim["jet_bot1_tracks"] = [60,0,60]; dict_binlim["jet_bot2_tracks"] = [60,0,60];
dict_binlim["muon_from_bot_sf_z"] = [50,0.5,1.1];dict_binlim["muon_from_bot_sf_iso"] = [50,0.7,1.2];
dict_binlim["tau_discr_jet1"]=[25,0,1];dict_binlim["tau_discr_jet2"]=[25,0,1];dict_binlim["tau_discr_jetbot1"]=[25,0,1];dict_binlim["tau_discr_jetbot2"]=[25,0,1];
dict_binlim["muon_jet_z_short"]=[20,0,0.5];dict_binlim["InvM3_good_short"]=[12,120,210];
dict_binlim["muon_jet_sigxy"]=[40,-4,4];dict_binlim["muon_jet_sigdz"]=[40,-4,4];dict_binlim["muon_jet_sigr"]=[40,0,20];dict_binlim["muon_jet_r"]=[30,0,0.05];
dict_binlim["muon_jet_z2"]=[40,0,1];dict_binlim["muon_jet_z3"]=[40,0,1];dict_binlim["muon_jet_iso_abs"]=[40,0,100];
dict_binlim["muon_jet_z2_v2"]=[40,0,1];
dict_binlim["lep_id_lowpt_sf_up"]=[40,0.5,1.5];dict_binlim["lep_id_lowpt_sf_down"]=[40,0.5,1.5];dict_binlim["lep_id_lowpt_sf"]=[40,0.5,1.5];
dict_binlim["muon_jet_mother"]=[4,0,4];

dict_binlim_M = dict(dict_binlim); dict_binlim_E = dict(dict_binlim);

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
        hist_seclepup_M[name] = {}
        hist_seclepup_E[name] = {}
        hist_seclepdown_M[name] = {}
        hist_seclepdown_E[name] = {}
        hist_lepidup_M[name] = {}
        hist_lepidup_E[name] = {}
        hist_lepiddown_M[name] = {}
        hist_lepiddown_E[name] = {}
        hist_lepisoup_M[name] = {}
        hist_lepisoup_E[name] = {}
        hist_lepisodown_M[name] = {}
        hist_lepisodown_E[name] = {}
        hist_leptrigup_M[name] = {}
        hist_leptrigup_E[name] = {}
        hist_leptrigdown_M[name] = {}
        hist_leptrigdown_E[name] = {}
        hist_notoppt_M[name] = {}
        hist_notoppt_E[name] = {}
        hist_puwup_M[name] = {}
        hist_puwup_E[name] = {}
        hist_puwdown_M[name] = {}
        hist_puwdown_E[name] = {}
        hist_smearup_M[name] = {}
        hist_smearup_E[name] = {}
        hist_smeardown_M[name] = {}
        hist_smeardown_E[name] = {}
        hist_jecup_M[name] = {}
        hist_jecup_E[name] = {}
        hist_jecdown_M[name] = {}
        hist_jecdown_E[name] = {}
        if mode == "mc":
              for s in samples:
                 hist_nom_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final")
       	         hist_nom_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final")
                 #if name == "InvM_2jets" and (chan == "btagMM_chitest" or chan == "lepton50_chitest"):
                 #   hist_nom_M[name][s].Sumw2()
       	         #   hist_nom_E[name][s].Sumw2()
                 hist_sfbtag_light_up_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_btaglightup","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_btaglightup")
       	         hist_sfbtag_light_up_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_btaglightup","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_btaglightup")
                 hist_sfbtag_light_down_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_btaglightdown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_btaglightdown")
       	         hist_sfbtag_light_down_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_btaglightdown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_btaglightdown")
                 hist_sfbtag_heavy_up_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_btagheavyup","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_btagheavyup")
       	         hist_sfbtag_heavy_up_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_btagheavyup","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_btagheavyup")
                 hist_sfbtag_heavy_down_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_btagheavydown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_btagheavydown")
       	         hist_sfbtag_heavy_down_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_btagheavydown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_btagheavydown")
                 hist_seclepup_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_seclepup","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_seclepup")
       	         hist_seclepup_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_seclepup","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_seclepup")                 
                 hist_seclepdown_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_seclepdown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_seclepdown")
       	         hist_seclepdown_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_seclepdown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_seclepdown")                 
                 hist_lepidup_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_lepidup","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_lepidup")
       	         hist_lepidup_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_lepidup","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_lepidup")                 
                 hist_lepiddown_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_lepiddown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_lepiddown")
       	         hist_lepiddown_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_lepiddown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_lepiddown")                 
                 hist_lepisoup_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_lepisoup","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_lepisoup")
       	         hist_lepisoup_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_lepisoup","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_lepisoup")                 
                 hist_lepisodown_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_lepisodown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_lepisodown")
       	         hist_lepisodown_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_lepisodown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_lepisodown")                 
                 hist_leptrigup_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_leptrigup","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_leptrigup")
       	         hist_leptrigup_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_leptrigup","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_leptrigup")                 
                 hist_leptrigdown_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_leptrigdown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_leptrigdown")
       	         hist_leptrigdown_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_leptrigdown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_leptrigdown")                 
                 hist_notoppt_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_notoppt","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_notoppt")
       	         hist_notoppt_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_notoppt","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_notoppt")    
                 hist_puwup_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_puwup","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_puwup")
                 hist_puwup_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_puwup","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_puwup")
                 hist_puwdown_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M_puwdown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final_puwdown")
                 hist_puwdown_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E_puwdown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final_puwdown")
                 if set(["smearup","smeardown"]).issubset(list_jetsyst):
                    hist_smearup_M[name][s] = df_M["smearup"][s].Histo1D((s+"_"+name+"_M_smearup","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final")
                    hist_smearup_E[name][s] = df_E["smearup"][s].Histo1D((s+"_"+name+"_E_smearup","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final")
                    hist_smeardown_M[name][s] = df_M["smeardown"][s].Histo1D((s+"_"+name+"_M_smeardown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final")
                    hist_smeardown_E[name][s] = df_E["smeardown"][s].Histo1D((s+"_"+name+"_E_smeardown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final")
                 if set(["jecup","jecdown"]).issubset(list_jetsyst):
                    hist_jecup_M[name][s] = df_M["jecup"][s].Histo1D((s+"_"+name+"_M_jecup","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final")
                    hist_jecup_E[name][s] = df_E["jecup"][s].Histo1D((s+"_"+name+"_E_jecup","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final")
                    hist_jecdown_M[name][s] = df_M["jecdown"][s].Histo1D((s+"_"+name+"_M_jecdown","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final")
                    hist_jecdown_E[name][s] = df_E["jecdown"][s].Histo1D((s+"_"+name+"_E_jecdown","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final")

        else:
           for s in samples:
              hist_nom_M[name][s] = df_M["nom"][s].Histo1D(("data"+s+"_"+name+"_M","",dict_binlim_M[name][0],dict_binlim_M[name][1],dict_binlim_M[name][2]),column_names[name],"weightSSOS_final")
       	      hist_nom_E[name][s] = df_E["nom"][s].Histo1D(("data"+s+"_"+name+"_E","",dict_binlim_E[name][0],dict_binlim_E[name][1],dict_binlim_E[name][2]),column_names[name],"weightSSOS_final")

hist_2D_M = {}
hist_2D_E = {}

obs_2D = ["jet_max_cvltag","jet_max_cvbtag"]

for i in range(int(len(obs_2D)/2)):
      name1 = obs_2D[2*i]
      name2 = obs_2D[2*i+1]
      hist_2D_M[name1] = {}
      hist_2D_E[name1] = {}
      for s in samples:
              hist_2D_M[name1][s] = df_M["nom"][s].Histo2D(("hist2d"+s+"_"+name1+"_M","",dict_binlim[name1][0],dict_binlim[name1][1],dict_binlim[name1][2],dict_binlim[name2][0],dict_binlim[name2][1],dict_binlim[name2][2]),column_names[name1],column_names[name2],"weightSSOS_final")
       	      hist_2D_E[name1][s] = df_E["nom"][s].Histo2D(("hist2d"+s+"_"+name1+"_E","",dict_binlim[name1][0],dict_binlim[name1][1],dict_binlim[name1][2],dict_binlim[name2][0],dict_binlim[name2][1],dict_binlim[name2][2]),column_names[name1],column_names[name2],"weightSSOS_final")


########## gen hists

hist_mc_M = {}
hist_mc_E = {}

observable_mc_names = ["jet_1_flavourP", "jet_2_flavourP", "jet_bot1_flavourP", "jet_bot2_flavourP", "btag_sf", "lep_id_sf", "lep_trig_sf", "lep_iso_sf",
     "puWeight", "PUjetID_SF", "top_weight","Frag_weight_sl","Br_weight_sl","muon_bot1_mother","muon_bot2_mother","muon_bot1_mother_mine","muon_bot2_mother_mine",
     "lep_id_lowpt_sf_up","lep_id_lowpt_sf_down","lep_id_lowpt_sf"]

if channel[0:2] == "sl": observable_mc_names = observable_mc_names + ["muon_from_bot_sf_z","muon_from_bot_sf_iso","muon_jet_mother"]

if mode == "mc":
    for name in observable_mc_names:
        hist_mc_M[name] = {}
        hist_mc_E[name] = {}
        for s in samples:
              hist_mc_M[name][s] = df_M["nom"][s].Histo1D((s+"_"+name+"_M","",dict_binlim[name][0],dict_binlim[name][1],dict_binlim[name][2]),name,"weightSSOS_final")
       	      hist_mc_E[name][s] = df_E["nom"][s].Histo1D((s+"_"+name+"_E","",dict_binlim[name][0],dict_binlim[name][1],dict_binlim[name][2]),name,"weightSSOS_final")

##################################################
################ NTUPLES SAVING ##################
##################################################

first_list_aux = df["nom"][samples[0]].GetColumnNames()
first_list = list(first_list_aux)

brlist = []

for el_st in first_list:
    el = str(el_st)
    if not (el[0:3]=="HLT" or el[0:3]=="HTX" or el[0:2]=="L1" or el[0:3]=="LHE" or el[0:4]=="Soft" or el[0:3]=="Cor"):
       brlist.append(el)       

#brlist = [el for el in first_list if not (el[0:3]=="HLT" or el[0:3]=="HTX" or e[0:2]=="L1" or el[0:3]=="LHE" or el[0:4]=="Soft")]
brlist_wqq = []

brlist = ["nElectronGood","ElectronGoodInd", "nElectron","Electron_charge", "Electron_cleanmask", "Electron_convVeto", "Electron_cutBased", "Electron_cutBased_HEEP", "Electron_dEscaleDown", 
  "Electron_dEscaleUp", "Electron_dEsigmaDown", "Electron_dEsigmaUp", "Electron_deltaEtaSC", "Electron_dr03EcalRecHitSumEt", "Electron_dr03HcalDepth1TowerSumEt", 
  "Electron_dr03TkSumPt", "Electron_dr03TkSumPtHEEP", "Electron_dxy", "Electron_dxyErr", "Electron_dz", "Electron_dzErr", "Electron_eCorr", "Electron_eInvMinusPInv", 
  "Electron_energyErr", "Electron_eta", "Electron_hoe", "Electron_ip3d", "Electron_isPFcand", "Electron_jetIdx", 
  "Electron_jetNDauCharged", "Electron_jetPtRelv2", "Electron_jetRelIso", "Electron_lostHits", "Electron_mass", "Electron_miniPFRelIso_all", 
  "Electron_miniPFRelIso_chg", "Electron_mvaFall17V2Iso", "Electron_mvaFall17V2Iso_WP80", "Electron_mvaFall17V2Iso_WP90", "Electron_mvaFall17V2Iso_WPL", 
  "Electron_mvaFall17V2noIso", "Electron_mvaFall17V2noIso_WP80", "Electron_mvaFall17V2noIso_WP90", "Electron_mvaFall17V2noIso_WPL", "Electron_mvaTTH", 
  "Electron_pdgId", "Electron_pfRelIso03_all", "Electron_pfRelIso03_chg", "Electron_phi", "Electron_photonIdx", "Electron_pt", "Electron_r9", "Electron_scEtOverPt", 
  "Electron_seedGain", "Electron_sieie", "Electron_sip3d", "Electron_tightCharge", "Electron_vidNestedWPBitmap", "Electron_vidNestedWPBitmapHEEP",
  "nJet", "Jet_area", "Jet_bRegCorr", "Jet_bRegRes", "Jet_btagCSVV2", "Jet_btagDeepB", "Jet_btagDeepCvB", "Jet_btagDeepCvL", "Jet_btagDeepFlavB", 
  "Jet_btagDeepFlavCvB", "Jet_btagDeepFlavCvL", "Jet_btagDeepFlavQG", "Jet_cRegCorr", "Jet_cRegRes", "Jet_chEmEF", "Jet_chFPV0EF", "Jet_chHEF", 
  "Jet_cleanmask", "Jet_electronIdx1", "Jet_electronIdx2", "Jet_eta", "Jet_hfadjacentEtaStripsSize", 
  "Jet_hfcentralEtaStripSize", "Jet_hfsigmaEtaEta", "Jet_hfsigmaPhiPhi", "Jet_jetId", "Jet_mass", "Jet_mass_nom", 
  "Jet_muEF", "Jet_muonIdx1", "Jet_muonIdx2", "Jet_muonSubtrFactor", "Jet_nConstituents", "Jet_nElectrons", "Jet_nMuons", "Jet_neEmEF", 
  "Jet_neHEF", "Jet_phi", "Jet_pt", "Jet_pt_nom", "Jet_puId", "Jet_puIdDisc", "Jet_qgl", 
  "Jet_rawFactor", "nJetGood","JetGoodInd", "nJetMuonInd", "JetMuonInd", "nJetSVInd", "JetSVInd", "nJetBotInd", "JetBotInd", "MET_MetUnclustEnUpDeltaX", 
  "MET_MetUnclustEnUpDeltaY", "MET_covXX", "MET_covXY", "MET_covYY", "MET_phi", "MET_pt", "MET_significance", "nJetQInd", "JetQInd",
  "MET_smeared_phi", "MET_smeared_pt", "MET_sumEt", "MET_sumPtUnclustered", "nMuonGood","MuonGoodInd", "MuonJetGood", "nMuonJetInd", "MuonJetInd", "MuonLepSign", 
  "nMuon","Muon_charge", "Muon_cleanmask", "Muon_dxy", "Muon_dxyErr", "Muon_dxybs", "Muon_dz", "Muon_dzErr", "Muon_eta", "Muon_fsrPhotonIdx", 
  "Muon_highPtId", "Muon_highPurity", "Muon_inTimeMuon", "Muon_ip3d", "Muon_isGlobal", "Muon_isPFcand", "Muon_isStandalone", "Muon_isTracker", 
  "Muon_jetIdx", "Muon_jetNDauCharged", "Muon_jetPtRelv2", "Muon_jetRelIso", "Muon_looseId", "Muon_mass", "Muon_mediumId", "Muon_mediumPromptId", "Muon_miniIsoId", 
  "Muon_miniPFRelIso_all", "Muon_miniPFRelIso_chg", "Muon_multiIsoId", "Muon_mvaId", "Muon_mvaLowPt", "Muon_mvaLowPtId", "Muon_mvaTTH", "Muon_nStations", 
  "Muon_nTrackerLayers", "Muon_pdgId", "Muon_pfIsoId", "Muon_pfRelIso03_all", "Muon_pfRelIso03_chg", "Muon_pfRelIso04_all", "Muon_phi", "Muon_pt", "Muon_ptErr", 
  "Muon_puppiIsoId", "Muon_segmentComp", "Muon_sip3d", "Muon_softId", "Muon_softMva", "Muon_softMvaId", "Muon_tightCharge", "Muon_tightId", "Muon_tkIsoId", 
  "Muon_tkRelIso", "Muon_triggerIdLoose", "Muon_tunepRelPt", "SVJetGood", "SVLepSign", "nSVJetInd", "SVJetInd", "nSV","SV_charge", "SV_chi2", "SV_dlen", 
  "SV_dlenSig", "SV_dxy", "SV_dxySig", "SV_eta", "SV_mass", "SV_ndof", "SV_ntracks", "SV_pAngle", "SV_phi", "SV_pt", "SV_x", "SV_y", "SV_z", "nTau","Tau_charge", 
  "Tau_chargedIso", "Tau_cleanmask", "Tau_decayMode", "Tau_dxy", "Tau_dz", "Tau_eta", "Tau_idAntiEleDeadECal", "Tau_idAntiMu", 
  "Tau_idDecayModeOldDMs", "Tau_idDeepTau2017v2p1VSe", "Tau_idDeepTau2017v2p1VSjet", "Tau_idDeepTau2017v2p1VSmu", "Tau_jetIdx", "Tau_leadTkDeltaEta", 
  "Tau_leadTkDeltaPhi", "Tau_leadTkPtOverTauPt", "Tau_mass", "Tau_neutralIso", "Tau_phi", "Tau_photonsOutsideSignalCone", "Tau_pt", "Tau_puCorr", 
  "Tau_rawDeepTau2017v2p1VSe", "Tau_rawDeepTau2017v2p1VSjet", "Tau_rawDeepTau2017v2p1VSmu", "Tau_rawIso", "Tau_rawIsodR03", "jet_1_pt","jet_1_nmu","jet_1_qgl",
  "jet_2_pt","jet_1_eta","jet_2_eta", "jet_1_phi","jet_2_phi","jet_1_mass","jet_2_mass","jet_2_qgl","jet_2_nmu", "InvM_2jets","deltaR_jet1_jet2","deltaphi_jet1_jet2",
  "deltaeta_jet1_jet2", "deltapt_jet1_jet2","tracks_jet1","tracks_jet2","EMN_jet1","EMC_jet1", "EMtotal_jet1","pT_sum","pT_product","deltaR_lep_2jets",
  "deltaphi_MET_2jets", "deltaphi_MET_jets_1","deltaphi_MET_jets_2","eta_2jets","pt_2jets", "deltaphi_lephad","deltaR_lephad","deltaphi_lep_2jets","deltaeta_lephad",
  "deltaeta_lep_2jets","pT_proy","pT_sum_2J","pT_Wlep","deltaR_lep_jet1","deltaR_lep_jet2", "jet_bot1_btag","jet_bot2_btag","jet_bot1_pt","jet_bot2_pt","jet_bot1_eta",
  "jet_bot2_eta","jet_bot1_phi","jet_bot2_phi","jet_1_btag","jet_2_btag", "InvM_bot_closer","InvM_bot_farther","jet_bot1_btagnumber","jet_bot2_btagnumber",
  "jet_1_btagnumber","jet_2_btagnumber","jet_1_cvltag","jet_2_cvltag","jet_1_cvltag_csv","jet_2_cvltag_csv","weightSSOS_final","InvM30","InvM31","InvMl0","InvMl1",
  "chi2_test0", "chi2_test1", "InvMl_good", "InvMl_bad", "InvM3_good", "InvM3_bad", "chi2_test_good", "chi2_test_bad", "jet_max_cvltag", "jet_min_cvltag","chi_bool",
  "jet_1_cvbtag","jet_2_cvbtag","jet_1_cvbtag_csv","jet_2_cvbtag_csv", "jet_max_cvbtag","jet_min_cvbtag","sl_bool","nJetSLInd","JetSLInd","MuonSLInd",
  "jet_bot1_bcorr","jet_bot2_bcorr"]

if (channel[0:2]=="sl"): brlist = brlist + ["muon_ssos"]

if mode == "mc":
   brlist = brlist + ["nGenJet", "GenJet_eta", "GenJet_hadronFlavour", "GenJet_mass", "GenJet_partonFlavour", "GenJet_phi", "GenJet_pt",
       "GenMET_phi", "GenMET_pt", "nGenPart","GenPart_eta", "GenPart_genPartIdxMother", "GenPart_mass", "GenPart_pdgId", "GenPart_phi",
       "GenPart_pt", "GenPart_status", "GenPart_statusFlags","Jet_hadronFlavour", "Jet_mass_smeared_down", "Jet_mass_smeared_up",
       "Jet_partonFlavour", "Jet_pt_smeared_down", "Jet_pt_smeared_up","lep_id_lowpt_sf","Br_weight_sl","Frag_weight_sl",
       "var_xsec","var_lumi","lumi_data","weight_lumi"]
   if (channel[0:2]=="sl"): brlist = brlist + ["muon_jet_genindx"]

#for s in samples:
#   if args.wcs:
#       if (channel[0:2]=="sl"):
#           path = '/pnfs/ciemat.es/data/cms/store/user/juvazque/data_further_analysis/btagMM/chitest/sl/folder'+years[0]+'/wcs/'
#       else:
#           path = '/pnfs/ciemat.es/data/cms/store/user/juvazque/data_further_analysis/btagMM/folder'+years[0]+'/wcs/'
#   else:
#       if (channel[0:2]=="sl"):
#           path = '/pnfs/ciemat.es/data/cms/store/user/juvazque/data_further_analysis/btagMM/chitest/sl/folder'+years[0]+'/'
#       else:
#           path = '/pnfs/ciemat.es/data/cms/store/user/juvazque/data_further_analysis/btagMM/folder'+years[0]+'/'
#   term = 'dataset_wqq_btagMM_fromJF_'+s
#   df["nom"][s].Snapshot("Events", path+term+".root", brlist)

#############################
####     DATA SAVING     ####
#############################

term1 = ""
if chan == "nobtag": term1 = "nobtag/"
elif chan == "btagMM": term1 = "btagMM/"
elif chan == "lepton50": term1 = "lepton50/"
elif chan == "btagMM_chitest": 
  if args.wcs:
     term1 = "btagMM/chi_test/wcs_classes/"
     #term1 = "btagMM/chi_test/wcs_classes_aux/test_aux/"
  else:
     term1 = "btagMM/chi_test/"
elif chan == "lepton50_chitest": term1 = "lepton50/chi_test/"

if (channel=="ctag"): termm = "ctag/"
if (channel=="noctag"): termm = "noctag/"

#observable_names = ["transverse_mass","MET_pt_aux"]

if mode == "mc":
   for name in observable_names:
        path_hist = '/nfs/cms/vazqueze/new_hists/fromJF/wqq/'+term1+termm+'hist_wqqfromJF_MC_'+years[0]+'_'+name+'.root'
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
              hist_seclepup_M[name][s].Write()
              hist_seclepup_E[name][s].Write()
              hist_seclepdown_M[name][s].Write()
              hist_seclepdown_E[name][s].Write()
              hist_lepidup_M[name][s].Write() 
              hist_lepidup_E[name][s].Write() 
              hist_lepiddown_M[name][s].Write() 
              hist_lepiddown_E[name][s].Write()
              hist_lepisoup_M[name][s].Write()
              hist_lepisoup_E[name][s].Write()
              hist_lepisodown_M[name][s].Write()
              hist_lepisodown_E[name][s].Write()
              hist_leptrigup_M[name][s].Write()
              hist_leptrigup_E[name][s].Write()
              hist_leptrigdown_M[name][s].Write()
              hist_leptrigdown_E[name][s].Write()
              hist_notoppt_M[name][s].Write()
              hist_notoppt_E[name][s].Write()
              hist_puwup_M[name][s].Write()
              hist_puwup_E[name][s].Write()
              hist_puwdown_M[name][s].Write()
              hist_puwdown_E[name][s].Write()
              if set(["smearup","smeardown"]).issubset(list_jetsyst):
                 hist_smearup_M[name][s].Write()
                 hist_smearup_E[name][s].Write()
                 hist_smeardown_M[name][s].Write()
                 hist_smeardown_E[name][s].Write()
              if set(["jecup","jecdown"]).issubset(list_jetsyst):
                 hist_jecup_M[name][s].Write()
                 hist_jecup_E[name][s].Write()
                 hist_jecdown_M[name][s].Write()
                 hist_jecdown_E[name][s].Write()
        myfile.Close()

   for name in observable_mc_names:
        path_hist = '/nfs/cms/vazqueze/new_hists/fromJF/wqq/'+term1+termm+'hist_wqqfromJF_MC_'+years[0]+'_'+name+'.root'
        myfile = TFile( path_hist, 'RECREATE' )
        for s in samples:
          hist_mc_M[name][s].Write()
          hist_mc_E[name][s].Write()
        myfile.Close()

if (mode == "data" and samples[0][-1] == "M"):
   for name in observable_names:
        path_hist = '/nfs/cms/vazqueze/new_hists/fromJF/wqq/'+term1+termm+'hist_wqqfromJF_dataM_'+years[0]+'_'+name+'.root'
        myfile = TFile( path_hist, 'RECREATE' )

        hist_nom_M[name][samples[0]].Write()

        myfile.Close()

if (mode == "data" and samples[0][-1] == "E"):
   for name in observable_names:
        path_hist = '/nfs/cms/vazqueze/new_hists/fromJF/wqq/'+term1+termm+'hist_wqqfromJF_dataE_'+years[0]+'_'+name+'.root'
        myfile = TFile( path_hist, 'RECREATE' )

        hist_nom_E[name][samples[0]].Write()

        myfile.Close()

for i in range(int(len(obs_2D)/2)):
      name1 = obs_2D[2*i]
      name2 = obs_2D[2*i+1]
      path_hist = '/nfs/cms/vazqueze/new_hists/fromJF/wqq/'+term1+termm+'hist_wqqfromJF_2D_'+years[0]+'_'+name1+'_'+name2+'.root'
      myfile = TFile( path_hist, 'RECREATE' )
      for s in samples:
           hist_2D_M[name1][s].Write()
           hist_2D_E[name1][s].Write()
      myfile.Close()

print('Ended succesfully')

