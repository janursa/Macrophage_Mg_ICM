#!/usr/bin/env python3
# -*- coding: utf-6 -*-
"""
Created on Fri May 27 20:52:13 2022

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

ILs_model_str = """
import "LPS_sbml.xml";
model IL6_model()
    PP: LPS_model();
    compartment comp1;
    IL6 in comp1;
    IFNGR is PP.IFNGR;
    IL4R is PP.IL4R;
    IRAK4 is PP.IRAK4;
    aTRAF6 is PP.aTRAF6;
    NFKB_n is PP.NFKB_n;
    pSTAT3 is PP.pSTAT3;
    // adjustments to zhao's model:

    // IL6 signaling pathway 
    ## IL6/R/JAK/STAT3
    # IL6_v11: $IL6R_prod => IL6R; k_il6r_p;
    # IL6_v12: IL6R => deg; k_il6r_d*IL6R;
    # IL6_v13: IL6R + PP.JAK => IL6R_JAK; k_il6r_jak_b*IL6R*(IL6>0)*PP.JAK*(PP.JAK>0);
    # IL6_v14: IL6R_JAK => IL6R + JAK; k_il6r_jak_ub*IL10R_JAK*(IL10R_JAK>0);
    # IL6_v15: IL6R_JAK + IL6 => IL6_R_JAK; k_il6rjack_il6_b*IL6R_JAK*(IL6R_JAK>0)*IL6*(IL6>0) 
    # IL6_v16: IL6_R_JAK => IL6R_JAK + IL6; k_il6rjack_il6_ub*IL6_R_JAK*(IL6_R_JAK>0);
    # IL6_v17: IL6_R_JAK => pIL6_R_JAK; k_il6rjack_a*IL6_R_JAK*(IL6_R_JAK>0)
    # IL6_v18: pIL6_R_JAK => IL6_R_JAK; k_il6rjack_ua*pIL6_R_JAK*(pIL6_R_JAK>0);
    # IL6_v19: pIL6_R_JAK + PP.STAT3 => pIL6_R_JAK_STAT3; k_il6complex_b*pIL6_R_JAK*(pIL6_R_JAK>0)*PP.STAT3*(PP.STAT3>0);
    # IL6_v110: pIL6_R_JAK_STAT3 => PP.pSTAT3 + pIL6_R_JAK; k_stat3_a*pIL6_R_JAK_STAT3*(pIL6_R_JAK_STAT3>0);
    
    IL6_v12: IL6 => deg; k_il6_d*IL6;
    IL6_v13: $IL6R + IL6 => IL6_R; k_il6r_b*IL6*(IL6>0) 
    IL6_v14: IL6_R => $IL6R + IL6; k_il6r_ub*IL6_R*(IL6_R>0);
    IL6_v15: IL6_R => pIL6_R; k_il6r_a*IL6_R*(IL6_R>0)
    IL6_v16: pIL6_R => IL6_R; k_il6r_da*pIL6_R*(pIL6_R>0);
    IL6_v17: $STAT3_s + pIL6_R => PP.pSTAT3; (F_stat3_a-1)*(F_stat3_a>1)
    
    ## PI3K/Akt Pathway
    PP.PI3K => PP.pPI3K; (F_pi3k_a-1)*(F_pi3k_a>1);
    
    // IL6 transcriptional stimulation 
    IL6_v21: $IL6_prod + PP.NFKB_n => IL6_m + PP.NFKB_n; k_il6_p*F_nfkb_il6_p;
    IL6_v22: IL6_m => IL6; k_il6m_il6*IL6_m*(IL6_m>0);
    IL6_v23: IL6 => deg; k_il6_d*IL6;
    IL6_v24: IL6_m => deg; k_il6m_d*IL6_m;

    // Variables
    ## Known ##
    IL6 = 0; # no extracellular IL6 initially
    k_il6_d = 0.693/(1*60); # 1h half life
    k_il6m_d = 0.693/(1*60); # 1h half life
     
    STAT3_s = 1 #fixed ones
    IL6R = 1
    
    ## unknown ## 
    IL6_R = 100
    pIL6_R = 100
    k_il6r_b = 0.01
    k_il6r_ub = 0.01
    k_il6r_a = 0.01
    k_il6r_da = 0.01
    k_stat3_a = 0.01
    kd_stat3_a = 100
    o_stat3_a = 0
    k_pi3k_a = 0.01
    kd_pi3k_a = 100
    o_pi3k_a = 0
    
    IL6_m = 5
    k_il6m_il6 = .1
    k_il6_p = 5
    k_nfkb_il6_p = 1
    kd_nfkb_il6_p = 100

    IL6_R_0 = IL6_R;
    IL6_m_0 = IL6_m;
    NFKB_n_0 = NFKB_n;
    IL6_0 = IL6;
    pSTAT3_0 = PP.pSTAT3;
    pIL6_R_0 = pIL6_R
    at (time > 0): IL6_R_0 = IL6_R;
    at (time > 0): IL6_m_0 = IL6_m;
    at (time > 0): IL6_0 = IL6;
    at (time > 0): NFKB_n_0 = NFKB_n;
    at (time > 0): pSTAT3_0 = PP.pSTAT3;
    at (time > 0): pIL6_R_0 = pIL6_R
    
    // assignements
    ee = 0.0001
    nIL6_m := IL6_m/(IL6_m_0+ee);
    nIL6 := IL6/(IL6_0+ee);
    nNFKB_n := NFKB_n/(NFKB_n_0+ee);
    npSTAT3 := PP.pSTAT3/(pSTAT3_0+ee);
    npIL6_R := pIL6_R/(pIL6_R_0+ee);
    F_stat3_a := 1+k_stat3_a*(pIL6_R-pIL6_R_0)/(pIL6_R-pIL6_R_0 + kd_stat3_a) - o_stat3_a
    F_pi3k_a := 1+k_pi3k_a*(pIL6_R-pIL6_R_0)/(pIL6_R-pIL6_R_0 + kd_pi3k_a) - o_pi3k_a
    F_nfkb_il6_p := 1+k_nfkb_il6_p*PP.NFKB_n/(PP.NFKB_n+kd_nfkb_il6_p);
end

"""

IL6_model = te.loada(IL6_model_str)
# Zhao_model = te.loadSBMLModel(dirs.dir_Zhao_model)
# species_IDs = Zhao_model.getFloatingSpeciesIds()
LPS_model = te.loadSBMLModel(dirs.dir_LPS_model)
species_IDs = LPS_model.getFloatingSpeciesIds()
IL6_model_m = tools.assign_surrogate_names(IL6_model,species_IDs)
IL6_model_m.exportToSBML(dirs.dir_IL6_model)


