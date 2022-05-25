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
from tools import dirs, tools

LPS_model_str = """
import "Zhao_sbml.xml";
model LPS_model()
    PP: pad_mac();
    compartment comp1;
    species LPS in comp1
    IKB is PP.IKB;
    IL1b is PP.IL1b;
    IL10 is PP.IL10;
    TNFa is PP.TNFa;
    // adjustments to zhao's model:
    ### LPS upregulates IRAK recruitment/production ###
     PP.v246: $PP.irak4_prod + PP.m93 + LPS => PP.IRAK4 + PP.m93 + LPS; PP.k246*PP.irak4_prod*(1 - PP.m93/(PP.m93 + PP.ka246))* F_LPS_irak;
    #PP.v248: PP.IL1b_R + PP.TRAF6 + PP.A20 + PP.SOCS1 + PP.m146b + PP.IRAK4 + PP.pIRAK4 => PP.IL1b_R + PP.aTRAF6 + PP.A20 + PP.SOCS1 + PP.m146b + PP.IRAK4 + PP.pIRAK4; PP.k248*PP.IL1b_R*(PP.IRAK4 + PP.pIRAK4)*PP.TRAF6*(1.01 - PP.A20/(PP.A20 + PP.ka248))*(1 - PP.SOCS1/(PP.SOCS1 + PP.kb248))*(1.2 - PP.m146b/(PP.m146b + PP.kc248))*F_LPS_irak;
        ## new reaction
    LPS -> deg; k_lps_d*LPS*(LPS>0);
    
    LPS_0 = 0;
    LPS = 0;
    
    // params
    kd_lps_irak_p  = 270;
     k_lps_d  = 0.199;
     k_lps_irak_p  = 631249
    F_LPS_irak := (1+k_lps_irak_p*(LPS)/(LPS+kd_lps_irak_p))
    
    IKB_0 = IKB;
    IL1b_0 = IL1b;
    IL10_0 = IL10;
    TNFa_0 = TNFa;
    at (time>0): IKB_0 = IKB;
    at (time > 0): IL1b_0 = IL1b;
    at (time > 0): IL10_0 = IL10;
    at (time > 0): TNFa_0 = TNFa;
    
    nIKB := IKB/IKB_0;
    nIL1b := IL1b/IL1b_0;
    nIL10 := IL10/IL10_0;
    nTNFa := TNFa/TNFa_0;
end

"""

LPS_model = te.loada(LPS_model_str)
Zhao_model = te.loadSBMLModel(dirs.dir_Zhao_model)
species_IDs = Zhao_model.getFloatingSpeciesIds()
LPS_model_m = tools.assign_surrogate_names(LPS_model,species_IDs)
LPS_model_m.exportToSBML(dirs.dir_LPS_model)


