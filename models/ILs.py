
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
model ILs_model()
    PP: LPS_model();
    compartment comp1;
    IFNGR is PP.IFNGR;
    IL4R is PP.IL4R;
    IRAK4 is PP.IRAK4;
    aTRAF6 is PP.aTRAF6;
    NFKB_n is PP.NFKB_n;
    // IL6 signaling pathway 
    ## IL6/R/JAK/STAT3
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

    // IL6 Variables
    ## Known ##
    IL6 = 0; # no extracellular IL6 initially
    k_il6_d = 0.693/(1*60); # 1h half life
    k_il6m_d = 0.693/(1*60); # 1h half life
     
    STAT3_s = 1 #fixed ones
    IL6R = 1
    
    ## unknown ## 
     IL6_R = 7.748742282127936;
     pIL6_R = 1.2311883376917345;
     k_il6r_b = 0.5904512813489826;
     k_il6r_ub = 0.43759281734498223;
     k_il6r_a = 0.5796630495201404;
     k_il6r_da = 0.2265100029645159;
     k_stat3_a = 628.1246136436464;
     kd_stat3_a = 227738.9303508934;
     o_stat3_a = 0.020407706587385555;
     k_pi3k_a = 438.02503363015575;
     kd_pi3k_a = 8083.748439722462;
     o_pi3k_a = 0.617622547347842;
     IL6_m = 2.5612346988978407;
     k_il6m_il6 = 0.725319111642058;
     k_il6_p = 0.025787551075836745;
     k_nfkb_il6_p = 155.00165490551376;
     kd_nfkb_il6_p = 50055.041687676334

    
    // IL8 signaling pathway 
    IL8_v11: $IL8R + IL8 => IL8_R; k_il8r_b*IL8*(IL8>0);
    IL8_v12: IL8_R => $IL8R + IL8; k_il8r_ub*IL8_R*(IL8_R>0);
    IL8_v13: IL8_R => pIL8_R; k_il8r_a*IL8_R*(IL8_R>0);
    IL8_v14: pIL8_R => IL8_R; k_il8r_da*pIL8_R*(pIL8_R>0);
    ### IL8 upregulates IRAK recruitment/production ###
    IL8_v15: $PP.irak4_prod + IL8 => IRAK4 + IL8;  (F_il8_irak-1)*(F_il8_irak>1);

    // IL8 transcriptional stimulation 
    IL8_v21: $IL8_prod + PP.NFKB_n => IL8_m + PP.NFKB_n; k_il8_p*F_nfkb_il8_p*(F_nfkb_il8_p>0);
    IL8_v22: IL8_m => IL8; k_il8m_il8*IL8_m*(IL8_m>0);
    IL8_v23: IL8 => deg; k_il8_d*IL8;
    IL8_v24: IL8_m => deg; k_il8m_d*IL8_m;
    
    // IL8 Variables:
    IL8_prod = 1;
    IL8 = 0; # no extracellular IL8 initially
    k_il8_d = 0.693/(4*60); # 4h half life
    k_il8m_d = 0.693/(10*60); # 10 half life

    
IL8_R = 6.195009596795672;
pIL8_R = 1.7473583035488112;
k_il8r_b = 0.516826345402611;
k_il8r_ub = 0.9852906413788729;
k_il8r_a = 0.9302971850170958;
k_il8r_da = 0.948106328130202;
kd_il8_irak_p = 52589.193719531635;
k_il8_irak_p = 16416.675889591734;
o_il8_irak_p = 0.09546118942554127;
IL8_m = 0.10265149156962394;
k_il8m_il8 = 0.11412357371927684;
k_il8_p = 0.12098375907390846;
kd_nfkb_il8_p = 8274.262288279933

    // IL8 assignements
    nIL8_m := IL8_m/IL8_m_0;
    nIL8 := IL8/IL8_0;
    nIFNGR := IFNGR/IFNGR_0;
    nIL4R := IL4R/IL4R_0;
    nIRAK4 := IRAK4/IRAK4_0;
    naTRAF6 := aTRAF6/aTRAF6_0;
    nNFKB_n := NFKB_n/NFKB_n_0;
    npIL8_R := pIL8_R/pIL8_R_0;
    nIL6 := IL6/IL6_0;
    F_il8_irak := 1+ k_il8_irak_p*(pIL8_R-pIL8_R_0)/((pIL8_R-pIL8_R_0)+kd_il8_irak_p) - o_il8_irak_p;
    F_nfkb_il8_p := PP.NFKB_n/(PP.NFKB_n+kd_nfkb_il8_p);
    
    // IL6 assignements
    ee = 0.0001
    nIL6_m := IL6_m/(IL6_m_0+ee);
    nIL6 := IL6/(IL6_0+ee);
    npSTAT3 := PP.pSTAT3/(pSTAT3_0+ee);
    npIL6_R := pIL6_R/(pIL6_R_0+ee);
    F_stat3_a := 1+k_stat3_a*(pIL6_R-pIL6_R_0)/(pIL6_R-pIL6_R_0 + kd_stat3_a) - o_stat3_a
    F_pi3k_a := 1+k_pi3k_a*(pIL6_R-pIL6_R_0)/(pIL6_R-pIL6_R_0 + kd_pi3k_a) - o_pi3k_a
    F_nfkb_il6_p := 1+k_nfkb_il6_p*PP.NFKB_n/(PP.NFKB_n+kd_nfkb_il6_p);
    
    // initial conditions, both IL6 and IL8
    IL8_m_0 = IL8_m;
    IFNGR_0 = IFNGR;
    IL4R_0 = IL4R;
    IRAK4_0 = IRAK4;
    aTRAF6_0 = aTRAF6;
    NFKB_n_0 = NFKB_n;
    IL8_0 = IL8;
    pIL8_R_0 = pIL8_R;
    IL6_R_0 = IL6_R;
    IL6_m_0 = IL6_m;
    IL6_0 = IL6;
    pSTAT3_0 = PP.pSTAT3;
    pIL6_R_0 = pIL6_R
    at (time > 0): IFNGR_0 = IFNGR;
    at (time > 0): IL4R_0 = IL4R;
    at (time > 0): IL8_m_0 = IL8_m;
    at (time > 0): IL8_0 = IL8;
    at (time > 0): IRAK4_0 = IRAK4;
    at (time > 0): aTRAF6_0 = aTRAF6; 
    at (time > 0): NFKB_n_0 = NFKB_n;
    at (time > 0): pIL8_R_0 = pIL8_R;
    at (time > 0): IL6_R_0 = IL6_R;
    at (time > 0): IL6_m_0 = IL6_m;
    at (time > 0): IL6_0 = IL6;
    at (time > 0): pSTAT3_0 = PP.pSTAT3;
    at (time > 0): pIL6_R_0 = pIL6_R
end

"""

ILs_model = te.loada(ILs_model_str)
# Zhao_model = te.loadSBMLModel(dirs.dir_Zhao_model)
# species_IDs = Zhao_model.getFloatingSpeciesIds()
LPS_model = te.loadSBMLModel(dirs.dir_LPS_model)
species_IDs = LPS_model.getFloatingSpeciesIds()
ILs_model_m = tools.assign_surrogate_names(ILs_model,species_IDs)
ILs_model_m.exportToSBML(dirs.dir_ILs_model)


