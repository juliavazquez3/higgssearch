#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <bitset>
#include <math.h>
#include <string.h>
#include <vector>
#include <map>
#include <iterator>
#include <unistd.h>
#include <stdio.h>
#include <algorithm>

#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TLine.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TMath.h"


void muon_eff_genreco_second()
{
  TChain * chain_ttbar_sl = new TChain("Events","");
  chain_ttbar_sl->Add((const char*)(TString("/pnfs/ciemat.es/data/cms/store/user/juvazque/PreprocessRDF_nojetID/myconfig2016/ttbar_sl/cat_base/prod_test/data_*.root")));
  TCanvas *c1 = new TCanvas("c1","Profile histogram example",200,10,700,500);
  TProfile *hprof  = new TProfile("hprof","Profile of pz versus px",100,-4,4,0,20);
  Float_t px, py, pz;
  for ( Int_t i=0; i<25000; i++) {
     gRandom->Rannor(px,py);
     pz = px*px + py*py;
     hprof->Fill(px,pz,1);
  }
  hprof->Draw();
  c1->SaveAs("test_muoneff.png");
  TProfile *hprof2  = new TProfile("hprof2","Profile of muon eta vs pt",100,0,50,50,-3,3);
  Float_t px, py, pz;
  for ( Int_t i=0; i<25000; i++) {
     gRandom->Rannor(px,py);
     pz = px*px + py*py;
     hprof->Fill(px,pz,1);
  }
  hprof->Draw();
  c1->SaveAs("test_muoneff2.png");
}

