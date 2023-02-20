#!/bin/bash

python higgssearch/fromJF/hist_ttbarlikefromJF_sl.py --stack --ratio --png --linear --type="sl"
python higgssearch/fromJF/hist_ttbarlikefromJF_sl.py --stack --ratio --png --linear --type="sv"
python higgssearch/fromJF/hist_just_one_year.py --stack --ratio --png --linear --type="sl" --year="2016"
python higgssearch/fromJF/hist_just_one_year.py --stack --ratio --png --linear --type="sl" --year="2016B"
python higgssearch/fromJF/hist_just_one_year.py --stack --ratio --png --linear --type="sl" --year="2017"
python higgssearch/fromJF/hist_just_one_year.py --stack --ratio --png --linear --type="sl" --year="2018"
python higgssearch/fromJF/hist_just_one_year.py --stack --ratio --png --linear --type="sv" --year="2016"
python higgssearch/fromJF/hist_just_one_year.py --stack --ratio --png --linear --type="sv" --year="2016B"
python higgssearch/fromJF/hist_just_one_year.py --stack --ratio --png --linear --type="sv" --year="2017"
python higgssearch/fromJF/hist_just_one_year.py --stack --ratio --png --linear --type="sv" --year="2018"

