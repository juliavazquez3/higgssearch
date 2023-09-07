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
print(file.keys())
print(file.classnames())

events = file['Events']

pdgID = events["GenPart_pdgId"].array()
motherID = events["GenPart_genPartIdxMother"].array()
isHardP = events["ishard"].array()
fromHardP = events["fromhard"].array()
firstC = events["first_copy"].array()
isWplusc = events["last_copy"].array()

ak.to_pandas(pdgID)
ak.to_pandas(motherID)
ak.to_pandas(isHardP)
ak.to_pandas(fromHardP)
ak.to_pandas(firstC)
ak.to_pandas(isWplusc)

d = {'pdgID': pdgID, 'motherID': motherID, 'isHardP': isHardP, 'fromHardP': fromHardP, 'first_copy': firstC, 'last_copy': isWplusc}

#df = pd.DataFrame(data=d)

list_num = []
my_list = [411,421,431,4122]
other_list = [10411,10421,413,423,10413,10423,20413,20423,415,425,10431,433,10433,20433,435]

for count, value in enumerate(pdgID[0:30]):
    my_bool = False
    #print(my_bool)
    for id in pdgID[count]:
        #print(id)
        if id in my_list:
           my_bool = True
    if my_bool:
       list_num.append(count)

print(list_num)

#df1 = df[df]

print('--------------------------------')
print(pdgID[0][0:22])
print(motherID[0][0:22])
print(isHardP[0][0:22])
print(fromHardP[0][0:22])
print(firstC[0][0:22])
print(isWplusc[0][0:22])

