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
file = uproot.open("/nfs/cms/vazqueze/higgssearch/fromJF/dataset_test_genpart_ttbarsl.root")
print(file.keys())
print(file.classnames())

events = file['Events']

pdgID = events["GenPart_pdgId"].array()
statusF = events["GenPart_statusFlags"].array()
status7 = events["status7"].array()
status8 = events["status8"].array()
status9 = events["status9"].array()
status11 = events["status11"].array()
status12 = events["status12"].array()
status13 = events["status13"].array()

ak.to_pandas(pdgID)
ak.to_pandas(statusF)
ak.to_pandas(status7)
ak.to_pandas(status8)
ak.to_pandas(status9)
ak.to_pandas(status11)
ak.to_pandas(status12)
ak.to_pandas(status13)

print('--------------------------------')
print(pdgID[0][0:15])
print(statusF[0][0:15])
print(status7[0][0:15])
print(status8[0][0:15])
print(status9[0][0:15])
print(status11[0][0:15])
print(status12[0][0:15])
print(status13[0][0:15])

print('--------------------------------')
print(pdgID[1][0:22])
print(statusF[1][0:22])
print(status7[1][0:22])
print(status8[1][0:22])
print(status9[1][0:22])
print(status11[1][0:22])
print(status12[1][0:22])
print(status13[1][0:22])

print('--------------------------------')
print(pdgID[2][0:22])
print(statusF[2][0:22])
print(status7[2][0:22])
print(status8[2][0:22])
print(status9[2][0:22])
print(status11[2][0:22])
print(status12[2][0:22])
print(status13[2][0:22])

