
import tellurium as te
import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
from tools import dirs, tools

IL8_model_str = """
import "LPS_sbml.xml";
model IL8_model()
    PP: LPS_model();
    compartment comp1;
    IFNGR is PP.IFNGR;
    IL4R is PP.IL4R;
    IRAK4 is PP.IRAK4;
    aTRAF6 is PP.aTRAF6;
    NFKB_n is PP.NFKB_n;
    // adjustments to zhao's model:
    ### IL8 upregulates IRAK recruitment/production ###
    PP.v246: $PP.irak4_prod + PP.m93 + PP.LPS + IL8 => PP.IRAK4 + PP.m93 + PP.LPS + IL8; PP.k246*PP.irak4_prod*(1 - PP.m93/(PP.m93 + PP.ka246))* PP.F_LPS_irak + F_il8_irak*(F_il8_irak>0);
    ### IL8 upregulates IFNGR production ###
    # PP.v80: $PP.IFNGR_prod + IL8_R => PP.IFNGR + IL8_R; PP.k80*PP.IFNGR_prod*F_il8_ifngr ;
    ### IL8 downregulates IL4R production ###
    #PP.v119: $PP.il4r_prod + PP.pSTAT3D_n + IL8_R => PP.IL4R + PP.pSTAT3D_n + IL8_R; PP.k119*PP.il4r_prod*(0.05 + PP.pSTAT3D_n/(PP.pSTAT3D_n + PP.ka119))*F_il8_il4p;    
    
    // new reactions
    ### IL8 receptor production and degradation ###
    IL8_v1: $IL8R_prod => IL8R; k_il8r_p;
    IL8_v2: IL8R => deg; k_il8r_d*IL8R;
    ### IL8/IL8R binding ###
    IL8_v3: IL8 + IL8R => IL8_R; k_il8_b*IL8*(IL8>0)*IL8R*(IL8R>0);
    IL8_v32: IL8_R => deg; k_il8r_d*IL8_R;
    IL8_v4: $IL8_prod + PP.NFKB_n => IL8_m + PP.NFKB_n; k_il8_p*F_nfkb_il8_p;
    IL8_v5: IL8_m => IL8; k_il8m_il8*IL8_m*(IL8_m>0);
    IL8_v6: IL8 => deg; k_il8_d*IL8;
    IL8_v7: IL8_m => deg; k_il8m_d*IL8_m;
    
    // Variables:
    IL8_prod = 1;

    IL8 = 0; # no extracellular IL8 initially
    IL8R = 100000;
    k_il8_d = 0.693/(4*60); # 4h half life
    k_il8m_d = 0.693/(10*60); # 10 half life


IL8_m = 1.00664143378666;
k_il8m_il8 = 0.7055534274159563;
k_il8_p = 7.7392148114967085;
kd_nfkb_il8_p = 9060.08392469455;
IL8_R = 9.170081838497477;
k_il8_b = 0.14605373485078582;
k_il8r_d = 0.07574267625162434;
k_il8r_p = 7.9394693179698415;
kd_il8_irak_p = 832.0720706453722;
k_il8_irak_p = 8931.088667122709
    
    kd_il8_ifngr_p = 1; 
    k_il8_ifngr_p = 1;
    kd_il8_il4r_p = 1;
    n_il8_il4r_p = 1;
   
    
    IL8R_0 = IL8R;
    IL8_R_0 = IL8_R;
    IL8_m_0 = IL8_m;
    IFNGR_0 = IFNGR;
    IL4R_0 = IL4R;
    IRAK4_0 = IRAK4;
    aTRAF6_0 = aTRAF6;
    NFKB_n_0 = NFKB_n;
    IL8_0 = IL8;
    at (time > 0): IFNGR_0 = IFNGR;
    at (time > 0): IL4R_0 = IL4R;
    at (time > 0): IL8R_0 = IL8R;
    at (time > 0): IL8_R_0 = IL8_R;
    at (time > 0): IL8_m_0 = IL8_m;
    at (time > 0): IL8_0 = IL8;
    at (time > 0): IRAK4_0 = IRAK4;
    at (time > 0): aTRAF6_0 = aTRAF6; 
    at (time > 0): NFKB_n_0 = NFKB_n;
    
    // assignements
    nIL8_m := IL8_m/IL8_m_0;
    nIL8 := IL8/IL8_0;
    nIL8R := IL8R/IL8R_0;
    nIFNGR := IFNGR/IFNGR_0;
    nIL4R := IL4R/IL4R_0;
    nIRAK4 := IRAK4/IRAK4_0;
    naTRAF6 := aTRAF6/aTRAF6_0;
    nNFKB_n := NFKB_n/NFKB_n_0;
    nIL8_R := IL8_R/IL8_R_0;
    F_il8_irak := 1+ k_il8_irak_p*(IL8_R-IL8_R_0)/((IL8_R-IL8_R_0)+kd_il8_irak_p);
    F_il8_ifngr := 1+ k_il8_ifngr_p*(IL8_R-IL8_R_0)/((IL8_R-IL8_R_0)+kd_il8_ifngr_p)
    F_il8_il4p := ((IL8_R_0+kd_il8_il4r_p)/(IL8_R+kd_il8_il4r_p))^n_il8_il4r_p
    F_nfkb_il8_p := PP.NFKB_n/(PP.NFKB_n+kd_nfkb_il8_p);
end

"""

IL8_model = te.loada(IL8_model_str)
# Zhao_model = te.loadSBMLModel(dirs.dir_Zhao_model)
# species_IDs = Zhao_model.getFloatingSpeciesIds()
LPS_model = te.loadSBMLModel(dirs.dir_LPS_model)
species_IDs = LPS_model.getFloatingSpeciesIds()
IL8_model_m = tools.assign_surrogate_names(IL8_model,species_IDs)
IL8_model_m.exportToSBML(dirs.dir_IL8_model)


