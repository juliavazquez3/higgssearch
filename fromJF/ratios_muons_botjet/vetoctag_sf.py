print('CTAG scale factors')

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join
import pandas as pd

import json
import argparse

## first is flavour 0, second 4 and third 5
## values correspond to first coeff [0], second [1] and third [3] in equation [0]+[1]*exp([2]*x)

mistag1 = "{{0,{{1.40787e-02,3.68374e-03,6.92609e-03},{3.46265e-01,-5.54932e-01,-4.80450e-02},{0.19,0.0,1.0}}},"
mistag2 = "{1,{{1.84855e-02,1.05117e-03,1.06155e-02},{3.57763e-01,-5.56522e-01,-4.86462e-02},{1.89828e-01,1.24003e-02,-2.46414e-02}}},"
mistag3 = "{2,{{0.0058,0.00008,0.016},{3.44864e-01,-7.02519e-01,-5.98220e-02},{2.78430e-01,7.61335e-02,-1.41138e-02}}},"
mistag4 = "{3,{{1.96041e-02,5.61596e-04,1.28072e-02},{4.12100e-01,-8.39918e-01,-6.40708e-0},{1.69046e-01,5.95589e-02,-1.50555e-02}}}}"

dict_mistag = mistag1+mistag2+mistag3+mistag4

# Funciones para los pesos completos
gInterpreter.Declare("""
   #include <iostream>
   #include <cstring>
   #include <string>
   using Vfloat = const ROOT::RVec<float>&;
   using Vint = const ROOT::RVec<int>&;
   auto myfunc_exp11(const float aa, const float bb, const float cc, const float xx) {
            float res;
            res = aa + bb*std::exp(cc*xx);
            return res;
   };
   auto total_veto_weight(const int jet_ind, Vint flav, Vfloat jetpt, const int yearT, const float ctagsf, const float cmistagsf) {
     std::map <int, std::vector<std::vector<float>>> dict_sf = %s;
     std::map <int, std::vector<std::vector<float>>>::iterator it2;
     float eff_mc = 1.;
     float weight_mc = 1.;
     auto flavranges = dict_sf[yearT];
     if (fabs(flav[jet_ind]) != 4 && fabs(flav[jet_ind]) != 5) {
        if (jetpt[jet_ind] > 20. && jetpt[jet_ind] <= 250.) {
            eff_mc = myfunc_exp11(flavranges[0][0], flavranges[0][1], flavranges[0][2], jetpt[jet_ind]);
        } else if (jetpt[jet_ind] > 250.) {
            eff_mc = myfunc_exp11(flavranges[0][0], flavranges[0][1], flavranges[0][2], 250.);
        } else {
            eff_mc =1.;
        }
        fabs(eff_mc) > 0. ? weight_mc = (1-cmistagsf*eff_mc)/(1-eff_mc) : weight_mc = 1.;
     } else if (fabs(flav[jet_ind]) == 4) {
        if (jetpt[jet_ind] > 20. && jetpt[jet_ind] <= 250.) {
            eff_mc = myfunc_exp11(flavranges[1][0], flavranges[1][1], flavranges[1][2], jetpt[jet_ind]);
        } else if (jetpt[jet_ind] > 250.) {
            eff_mc = myfunc_exp11(flavranges[1][0], flavranges[1][1], flavranges[1][2], 250.);
        } else {
            eff_mc =1.;
        }
        fabs(eff_mc) > 0. ? weight_mc = (1-ctagsf*eff_mc)/(1-eff_mc) : weight_mc = 1.;
     } else if (fabs(flav[jet_ind]) == 5) {
        if (jetpt[jet_ind] > 20. && jetpt[jet_ind] <= 250.) {
            eff_mc = myfunc_exp11(flavranges[2][0], flavranges[2][1], flavranges[2][2], jetpt[jet_ind]);
        } else if (jetpt[jet_ind] > 250.) {
            eff_mc = myfunc_exp11(flavranges[2][0], flavranges[2][1], flavranges[2][2], 250.);
        } else {
            eff_mc =1.;
        }
        fabs(eff_mc) > 0. ? weight_mc = (1-ctagsf*eff_mc)/(1-eff_mc) : weight_mc = 1.;
     }
     return weight_mc;
   }
""" % dict_mistag)

####################################################################################################################################################################################

class vetoctag_weights_tot():
    def __init__(self, *args, **kwargs):
        #self.isUL = kwargs.pop("isUL")
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

    def run(self, df):
        if self.isMC:
           if self.year == "2016":
               df = df.Define('vetoctag_weight','total_veto_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 0, ctag_weight, cmistag_weight)')
           elif self.year == "2016B":
               df = df.Define('vetoctag_weight','total_veto_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 1, ctag_weight, cmistag_weight)')
           elif self.year == "2017":
               df = df.Define('vetoctag_weight','total_veto_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 2, ctag_weight, cmistag_weight)')
           elif self.year == "2018":
               df = df.Define('vetoctag_weight','total_veto_weight(jet_max_ctag_indx, Jet_hadronFlavour, Jet_pt_aux, 3, ctag_weight, cmistag_weight)')

        variables = ['vetoctag_weight']

        branches = variables

        return df

def vetoctag_weights_totRDF(**kwargs):
    """
    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: jetVarRDF
            path: Base.Modules.smearing
            parameters:
                isMC: self.dataset.process.isMC
                year: self.config.year
                isUL: self.dataset.has_tag('ul')
                ¿¿ proc: self.dataset.process ??
    """
    return lambda: vetoctag_weights_tot(**kwargs)

