#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 13:54:58 2022

@author: matin
"""



import tellurium as te
import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
from tools import dirs, common

LPS_model_str = """
import "Zhao_sbml.xml";
model LPS_model()
    PP: pad_mac();
    compartment comp1;
    species aLPS in comp1
    species LPS in comp1
    IKB is PP.IKB;
    IL1b is PP.IL1b;
    IL10 is PP.IL10;
    TNFa is PP.TNFa;
    NFKB_n is PP.NFKB_n;
    // adjustments to zhao's model:
    ### LPS upregulates IRAK recruitment/production ###
    PP.IKB_NFKB + LPS -> PP.IKB + PP.NFKB + LPS; (F_lps-1)*(F_lps>1)*PP.IKB_NFKB
    ## new reaction
    LPS -> aLPS; k_lps_a*LPS;
    aLPS -> ; k_lpsa_d*aLPS;
    LPS -> ; k_lps_d*LPS*(LPS>0);
    
    LPS_0 = 0;
    LPS = 0;
    aLPS = 0
    
    // params
    k_lps_a = 0.5734580614584536;
     k_lps_d = 0.21398196410597586;
     k_lpsa_d = 0.006555185561568644;
     k_lps = 1419.6323301780503;
     kd_lps = 70095.63979887372

    
    
    IKB_0 = IKB;
    IL1b_0 = IL1b;
    IL10_0 = IL10;
    TNFa_0 = TNFa;
    NFKB_n_0 = NFKB_n;
    at (time > 0): IKB_0 = IKB;
    at (time > 0): IL1b_0 = IL1b;
    at (time > 0): IL10_0 = IL10;
    at (time > 0): TNFa_0 = TNFa;
    at (time > 0): NFKB_n_0 = NFKB_n;
    
    nIKB := IKB/IKB_0;
    nIL1b := IL1b/IL1b_0;
    nIL10 := IL10/IL10_0;
    nTNFa := TNFa/TNFa_0;
    nNFKB_n := (NFKB_n)/(NFKB_n_0);
    
    F_lps := 1+k_lps*aLPS/(aLPS+kd_lps)
end

"""

LPS_model = te.loada(LPS_model_str)
Zhao_model = te.loadSBMLModel(dirs.dir_Zhao_model)
species_IDs = Zhao_model.getFloatingSpeciesIds()
LPS_model_m = common.assign_surrogate_names(LPS_model,species_IDs)
LPS_model_m.exportToSBML(dirs.dir_LPS_model)


