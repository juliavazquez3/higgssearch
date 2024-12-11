import ROOT, os, sys
from ROOT import *

import numpy as np
import awkward as ak
import uproot
import pandas as pd
import pyarrow as pa
import urllib.request


#Events_set=events.arrays(brlist)

#print(Events_set)

#events_pd=ak.to_pandas(Events_set)

#print(events_pd.head())

#events_pd.to_csv('myWW.csv',header=True)

#print(events_pd[events_pd["typeWW"]==2].groupby(level=[0,1])["GenPart_pdgId"])

#file = uproot.open("/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/260000/062085CD-8DF4-1D40-8042-63998B6A3A95.root")
file = uproot.open("/pnfs/ciemat.es/data/cms/store/user/juvazque/data_test/ttbar_sl.root")
file = uproot.open("/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL18NanoAODv9/TT_TuneCH3_13TeV-powheg-herwig7/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/54FCF727-5A53-244F-93A7-A7A519DD150D.root")
file = uproot.open("/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL18NanoAODv9/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/280000/02003FED-4DED-2849-BAFB-494212F691EE.root")
print(file.keys())
print(file.classnames())

events = file['Events']

pdgID = events["GenPart_pdgId"].array()
motherID = events["GenPart_genPartIdxMother"].array()
#isHardP = events["ishard"].array()
#fromHardP = events["fromhard"].array()
#firstC = events["first_copy"].array()
#isWplusc = events["last_copy"].array()

ak.to_pandas(pdgID)
ak.to_pandas(motherID)
#ak.to_pandas(isHardP)
#ak.to_pandas(fromHardP)
#ak.to_pandas(firstC)
#ak.to_pandas(isWplusc)

print('--------------------------------')
print(pdgID[0][0:20])
print(motherID[0][0:20])
print(pdgID[0][20:40])
print(motherID[0][20:40])
#print(isHardP[0][0:22])
#print(fromHardP[0][0:22])
#print(firstC[0][0:22])
#print(isWplusc[0][0:22])

print('--------------------------------')
print(pdgID[1][0:22])
print(motherID[1][0:22])
print(pdgID[1][22:42])
print(motherID[1][22:42])

print('--------------------------------')
print(pdgID[2][0:17])
print(motherID[2][0:17])
print(pdgID[2][17:37])
print(motherID[2][17:37])
