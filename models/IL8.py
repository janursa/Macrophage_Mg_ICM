
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
    IL1b is PP.IL1b;
    IL10 is PP.IL10;
    TNFa is PP.TNFa;
    IRAK4 is PP.IRAK4;
    aTRAF6 is PP.aTRAF6;
    NFKB_n is PP.NFKB_n;
    v226 is PP.v226;
    v221 is PP.v221;
    k226 is PP.k226;
    IKB is PP.IKB;
    Ikb_prod is PP.Ikb_prod;
    k221 is PP.k221;
    // adjustments to zhao's model:
    ### IL8 upregulates IRAK recruitment/production ###
    PP.v246: $PP.irak4_prod + PP.m93 + PP.LPS + IL8 => PP.IRAK4 + PP.m93 + PP.LPS + IL8; PP.k246*PP.irak4_prod*(1 - PP.m93/(PP.m93 + PP.ka246))* PP.F_LPS_irak + F_il8_irak;
    # PP.v246: $PP.irak4_prod + PP.m93 => PP.IRAK4 + PP.m93; PP.k246*PP.irak4_prod*(1 - PP.m93/(PP.m93 + PP.ka246))* F_il8_irak;
    ### IL8 upregulates IFNGR production ###
    #PP.v80: $PP.IFNGR_prod + IL8_R => PP.IFNGR + IL8_R; PP.k80*PP.IFNGR_prod*F_il8_ifngr ;
    ### IL8 downregulates IL4R production ###
    #PP.v119: $PP.il4r_prod + PP.pSTAT3D_n + IL8_R => PP.IL4R + PP.pSTAT3D_n + IL8_R; PP.k119*PP.il4r_prod*(0.05 + PP.pSTAT3D_n/(PP.pSTAT3D_n + PP.ka119))*F_il8_il4p;    
    
    // new reactions
    ### IL8 receptor production and degradation ###
    v1: $IL8R_prod -> IL8R; k_il8r_p;
    v2: IL8R -> deg; k_il8r_d*IL8R;
    ### IL8/IL8R binding ###
    v3: IL8 + IL8R -> IL8_R; k_il8_b*IL8*IL8R - k_il8_ub*IL8_R;
    v4: $IL8_prod + PP.NFKB_n -> IL8 + PP.NFKB_n; k_il8_p*F_nfkb_il8_p;
    v5: IL8 -> deg; k_il8_d*IL8;
    
    // Variables:
    IL8_prod = 1;


IL8 = 7.408568221151385; # pg/ml
k_il8_p = 5.970462870815936;
kd_nfkb_il8_p = 2085.514347894692;
k_il8_d = 0.37595294602525076;
IL8R = 4800.422519863175;
IL8_R = 7.689186472021902;
k_il8_b = 0.0011844646079401944;
k_il8_ub = 0.7266871178144181;
k_il8r_p = 917.4566246712063;
k_il8r_d = 0.19039877492259133;
kd_il8_irak_p = 33.35345100816778 
    
    kd_il8_ifngr_p = 1; 
    n_il8_ifngr_p = 1;
    kd_il8_il4r_p = 1;
    n_il8_il4r_p = 1;
   
    
    IL8R_0 = IL8R;
    IL8_R_0 = IL8_R;
    IL8_0 = IL8;
    IFNGR_0 = IFNGR;
    IL4R_0 = IL4R;
    IL1b_0 = IL1b;
    IL10_0 = IL10;
    TNFa_0 = TNFa
    IRAK4_0 = IRAK4;
    aTRAF6_0 = aTRAF6;
    NFKB_n_0 = NFKB_n;
    IKB_0 = IKB;
    at (time > 0): IFNGR_0 = IFNGR;
    at (time > 0): IL4R_0 = IL4R;
    at (time > 0): IL8R_0 = IL8R;
    at (time > 0): IL8_R_0 = IL8_R;
    at (time > 0): IL8_0 = IL8;
    at (time > 0): IL1b_0 = IL1b;
    at (time > 0): IL10_0 = IL10;
    at (time > 0): TNFa_0 = TNFa;
    at (time > 0): IRAK4_0 = IRAK4;
    at (time > 0): aTRAF6_0 = aTRAF6; 
    at (time > 0): NFKB_n_0 = NFKB_n;
    at (time > 0): IKB_0 = IKB;
    
    // assignements
    nIL8 := IL8/IL8_0;
    nIL8R := IL8R/IL8R_0;
    nIFNGR := IFNGR/IFNGR_0;
    nIL4R := IL4R/IL4R_0;
    nIL1b := IL1b/IL1b_0;
    nIL10 := IL10/IL10_0;
    nTNFa := TNFa/TNFa_0;
    nIRAK4 := IRAK4/IRAK4_0;
    naTRAF6 := aTRAF6/aTRAF6_0;
    nNFKB_n := NFKB_n/NFKB_n_0;
    nIKB := IKB/IKB_0;
    F_il8_irak := ((IL8_R+kd_il8_irak_p)/(IL8_R_0+kd_il8_irak_p));
    F_il8_ifngr := ((IL8_R+kd_il8_ifngr_p)/(IL8_R_0+kd_il8_ifngr_p))^n_il8_ifngr_p
    F_il8_il4p := ((IL8_R_0+kd_il8_il4r_p)/(IL8_R+kd_il8_il4r_p))^n_il8_il4r_p
    F_nfkb_il8_p := PP.NFKB_n/(PP.NFKB_n+kd_nfkb_il8_p);
end

"""

IL8_model = te.loada(IL8_model_str)
Zhao_model = te.loadSBMLModel(dirs.dir_Zhao_model)
species_IDs = Zhao_model.getFloatingSpeciesIds()
IL8_model_m = tools.assign_surrogate_names(IL8_model,species_IDs)
IL8_model_m.exportToSBML(dirs.dir_IL8_model)


