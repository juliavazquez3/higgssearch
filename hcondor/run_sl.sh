#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt/setup.sh

python /nfs/cms/vazqueze/higgssearch/selection_sl.py --process="$1" --year="$2" --type="$3"
