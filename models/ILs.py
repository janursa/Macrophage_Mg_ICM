
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
    AP1_n is PP.AP1_n;
    // IL6 signaling pathway 
    ## IL6/R/JAK/STAT3
    IL6_v12: IL6 => deg; k_il6_d*IL6;
    IL6_v13: $IL6R_JACK + IL6 => IL6_R_JACK; k_il6r_b*IL6*(IL6>0) 
    IL6_v14: IL6_R_JACK => $IL6R_JACK + IL6; k_il6r_ub*IL6_R_JACK*(IL6_R_JACK>0);
    IL6_v15: IL6_R_JACK => pIL6_R_JACK; k_il6r_a*IL6_R_JACK*(IL6_R_JACK>0)
    IL6_v16: pIL6_R_JACK => IL6_R_JACK; k_il6r_da*pIL6_R_JACK*(pIL6_R_JACK>0);
    IL6_v17: $STAT3_s + pIL6_R_JACK => PP.pSTAT3 + IL6_R_JACK; (F_stat3_a-1)*(F_stat3_a>1)
    ## PI3K/Akt Pathway
    IL6_v18: PP.PI3K +IL6_R_JACK => PP.pPI3K + IL6_R_JACK; (F_pi3k_a-1)*(F_pi3k_a>1);
    
    // IL6 transcriptional stimulation 
    IL6_v21: $IL6_prod + PP.NFKB_n => IL6_m + PP.NFKB_n; k_il6_p+(F_nfkb_il6_p-1)*(F_nfkb_il6_p>1);
    IL6_v22: IL6_m => IL6; k_il6m_il6*IL6_m*(IL6_m>0);
    IL6_v23: IL6 => deg; k_il6_d*IL6;
    IL6_v24: IL6_m => deg; k_il6m_d*IL6_m;

    // IL6 Variables
    ## Known ##
    k_il6_d = 0.693/(1*60); # 1h half life
    k_il6m_d = 0.693/(1*60); # 1h half life
     
    STAT3_s = 1 #fixed ones
    IL6R_JACK = 1
    
    ## unknown ## 
     IL6 = 3.896510900860572;
     IL6_R_JACK = 1.90211958021639;
     pIL6_R_JACK = 7.151487103680696;
     k_il6r_b = 0.22041975560303922;
     k_il6r_ub = 0.28784102911843257;
     k_il6r_a = 0.8991179649849886;
     k_il6r_da = 0.36458082030419275;
     k_stat3_a = 19.09548272953805;
     kd_stat3_a = 16.108585912326816;
     o_stat3_a = 0.04027301267027811;
     k_pi3k_a = 285.0102124852435;
     kd_pi3k_a = 183.71671477158088;
     o_pi3k_a = 0.3393357854058401;
     IL6_m = 3.079490302967387;
     k_il6m_il6 = 0.10687115145921056;
     k_il6_p = 0.11548377600305315
     
     k_nfkb_il6_p = 155.00165490551376;
     kd_nfkb_il6_p = 50055.041687676334
     
    // IL6 assignements
    ee = 0.0001
    nIL6_m := IL6_m/(IL6_m_0+ee);
    nIL6 := IL6/(IL6_0+ee);
    npSTAT3 := PP.pSTAT3/(pSTAT3_0+ee);
    F_stat3_a := 1+k_stat3_a*(pIL6_R_JACK-pIL6_R_JACK_0)/(pIL6_R_JACK-pIL6_R_JACK_0 + kd_stat3_a) - o_stat3_a
    F_pi3k_a := 1+k_pi3k_a*(pIL6_R_JACK-pIL6_R_JACK_0)/(pIL6_R_JACK-pIL6_R_JACK_0 + kd_pi3k_a) - o_pi3k_a
    F_nfkb_il6_p := 1+k_nfkb_il6_p*(PP.NFKB_n-NFKB_n_0)/(PP.NFKB_n-NFKB_n_0+kd_nfkb_il6_p);
    

    // IL8 transcriptional stimulation 
    IL8_v11: $IL8_prod + PP.NFKB_n + PP.AP1_n => IL8_m + PP.NFKB_n + PP.AP1_n; k_il8_p+(F_nfkb_il8_p-1)*(F_nfkb_il8_p>1)+(F_ap1_il8_p-1)*(F_ap1_il8_p>1);
    IL8_v12: IL8_m => IL8; k_il8m_il8*IL8_m*(IL8_m>0);
    IL8_v13: IL8 => deg; k_il8_d*IL8;
    IL8_v14: IL8_m => deg; k_il8m_d*IL8_m;
    
    // IL8 signaling pathway 
    IL8_v21: $IL8R + IL8 => IL8_R; k_il8r_b*IL8*(IL8>0);
    IL8_v22: IL8_R => $IL8R + IL8; k_il8r_ub*IL8_R*(IL8_R>0);
    IL8_v23: IL8_R => pIL8_R; k_il8r_a*IL8_R*(IL8_R>0);
    IL8_v24: pIL8_R => IL8_R; k_il8r_da*pIL8_R*(pIL8_R>0);
    ### IL8 upregulates IRAK recruitment/production ###
    IL8_v25: $PP.irak4_prod + IL8 => IRAK4 + IL8;  (F_il8_irak-1)*(F_il8_irak>1);
    ### IL8 activates Rho GTPase
    IL8_26: $RHO + pIL8_R => aRHO + pIL8_R; (F_rho_a-1)*(F_rho_a>1);
    IL8_27: PP.PI3K + aRHO => PP.pPI3K + aRHO; (F_rho_pi3k_a-1)*(F_rho_pi3k_a>1);
    IL8_28: $STAT3_s + aRHO => PP.pSTAT3 + aRHO; (F_rho_stat3_a-1)*(F_rho_stat3_a>1)
    IL8_28: PP.JNK + aRHO => PP.pJNK + aRHO; (F_rho_jnk_a-1)*(F_rho_jnk_a>1)
    
    // IL8 Variables:
    IL8_prod = 1;
    k_il8_d = 0.693/(4*60); # 4h half life
    k_il8m_d = 0.693/(10*60); # 10 half life
    IL8R = 1
    RHO = 1
    aRHO = 0
    
 IL8 = 0; 
 IL8_m = 0;
 k_il8_p = 1
 k_il8m_il8 = 0.11412357371927684;
 IL8_R = 7.955153855553573;
 pIL8_R = 0.17994138285450756;
 k_il8r_b = 0.05553178676819742;
 k_il8r_ub = 0.8587736241604975;
 k_il8r_a = 0.6195023644260487;
 k_il8r_da = 0.9204202223990089;
 kd_il8_irak_p = 26044.044416341196;
 k_il8_irak_p = 42848.4639373502;
 o_il8_irak_p = 0.10546420546500945;
 k_rho_a = 39.496623324783286;
 kd_rho_a = 75501.38120137055;
 o_rho_a = 0.16221064736700586;
 k_rho_pi3k_a = 95741.02954720988;
 kd_rho_pi3k_a = 1557.6419724575026;
 o_rho_pi3k_a = 0.03610475397019314;
 k_rho_stat3_a = 21946.417237773316;
 kd_rho_stat3_a = 37308.68663289963;
 o_rho_stat3_a = 0.24260054422200417;
 k_rho_jnk_a = 83917.9235530874;
 kd_rho_jnk_a = 22770.275011467245;
 o_rho_jnk_a = 0.140118655310842
    
    
    k_nfkb_il8_p = 100
    kd_nfkb_il8_p = 100
    o_nfkb_il8_p = 0
    k_ap1_il8_p = 100
    kd_ap1_il8_p = 100
    o_ap1_il8_p = 0
    
    

    // IL8 assignements
    nIL8_m := (IL8_m+ee)/(IL8_m_0+ee);
    nIL8 := (IL8+ee)/(IL8_0+ee);
    nIFNGR := (IFNGR+ee)/(IFNGR_0+ee);
    nIL4R := (IL4R+ee)/(IL4R_0+ee);
    nIRAK4 := (IRAK4+ee)/(IRAK4_0+ee);
    naTRAF6 := (aTRAF6+ee)/(aTRAF6_0+ee);
    nNFKB_n := (NFKB_n+ee)/(NFKB_n_0+ee);
    npIL8_R := (pIL8_R+ee)/(pIL8_R_0+ee);
    nIL6 := (IL6+ee)/(IL6_0+ee);
    F_il8_irak := 1+ k_il8_irak_p*(pIL8_R-pIL8_R_0)/((pIL8_R-pIL8_R_0)+kd_il8_irak_p) - o_il8_irak_p;
    F_nfkb_il8_p := 1+ k_nfkb_il8_p*(NFKB_n-NFKB_n_0)/(NFKB_n-NFKB_n_0+kd_nfkb_il8_p) - o_nfkb_il8_p;
    F_ap1_il8_p := 1+k_ap1_il8_p*(AP1_n-AP1_n_0)/(AP1_n-AP1_n_0+kd_ap1_il8_p) - o_ap1_il8_p
    
    F_rho_a = 1
    F_rho_pi3k_a = 1
    F_rho_stat3_a = 1
    F_rho_jnk_a = 1
    
    # F_rho_a := 1+ k_rho_a*(pIL8_R-pIL8_R_0)/((pIL8_R-pIL8_R_0)+kd_rho_a) - o_rho_a
    # F_rho_pi3k_a := 1+ k_rho_pi3k_a*(aRHO)/(aRHO+kd_rho_pi3k_a) - o_rho_pi3k_a
    # F_rho_stat3_a := 1 + k_rho_stat3_a*aRHO/(aRHO+kd_rho_stat3_a) - o_rho_stat3_a
    # F_rho_jnk_a := 1 + k_rho_jnk_a*aRHO/(aRHO+kd_rho_jnk_a) - o_rho_jnk_a
    // initial conditions, both IL6 and IL8
    IL8_m_0 = IL8_m;
    IFNGR_0 = IFNGR;
    IL4R_0 = IL4R;
    IRAK4_0 = IRAK4;
    aTRAF6_0 = aTRAF6;
    NFKB_n_0 = NFKB_n;
    IL8_0 = IL8;
    pIL8_R_0 = pIL8_R;
    aRHO_0 = aRHO
    AP1_n_0 = AP1_n
    IL6_m_0 = IL6_m;
    IL6_0 = IL6;
    pSTAT3_0 = PP.pSTAT3;
    pIL6_R_JACK_0 = pIL6_R_JACK
    at (time > 0): IFNGR_0 = IFNGR;
    at (time > 0): IL4R_0 = IL4R;
    at (time > 0): IL8_m_0 = IL8_m;
    at (time > 0): IL8_0 = IL8;
    at (time > 0): IRAK4_0 = IRAK4;
    at (time > 0): aTRAF6_0 = aTRAF6; 
    at (time > 0): NFKB_n_0 = NFKB_n;
    at (time > 0): pIL8_R_0 = pIL8_R;
    at (time > 0): IL6_m_0 = IL6_m;
    at (time > 0): IL6_0 = IL6;
    at (time > 0): pSTAT3_0 = PP.pSTAT3;
    at (time > 0): pIL6_R_JACK_0 = pIL6_R_JACK
end

"""

ILs_model = te.loada(ILs_model_str)
# Zhao_model = te.loadSBMLModel(dirs.dir_Zhao_model)
# species_IDs = Zhao_model.getFloatingSpeciesIds()
LPS_model = te.loadSBMLModel(dirs.dir_LPS_model)
species_IDs = LPS_model.getFloatingSpeciesIds()
ILs_model_m = tools.assign_surrogate_names(ILs_model,species_IDs)
ILs_model_m.exportToSBML(dirs.dir_ILs_model)


