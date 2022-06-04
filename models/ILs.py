
import tellurium as te
import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
from tools import dirs, common

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
    pSTAT3 is PP.pSTAT3;
    // IL6 signaling pathway 
    ## IL6/R/JAK/STAT3
    IL6_v12: IL6 => deg; k_il6_d*IL6;
    IL6_v13: $IL6R_JACK + IL6 => IL6_R_JACK; k_il6r_b*IL6*(IL6>0) 
    IL6_v14: IL6_R_JACK => $IL6R_JACK + IL6; k_il6r_ub*IL6_R_JACK*(IL6_R_JACK>0);
    IL6_v15: IL6_R_JACK => pIL6_R_JACK; k_il6r_a*IL6_R_JACK*(IL6_R_JACK>0)
    IL6_v16: pIL6_R_JACK => IL6_R_JACK; k_il6r_da*pIL6_R_JACK*(pIL6_R_JACK>0);
    IL6_v17: $STAT3_s + pIL6_R_JACK => PP.pSTAT3 + IL6_R_JACK; (F_il6_stat3_a-1)*(F_il6_stat3_a>1)
    ## PI3K/Akt Pathway
    IL6_v18: PP.PI3K +IL6_R_JACK => PP.pPI3K + IL6_R_JACK; (F_il6_pi3k_a-1)*(F_il6_pi3k_a>1);
    
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
     IL6 = 7.734409288209127;
     IL6_R_JACK = 3.6276941315489175;
     pIL6_R_JACK = 9.229260900735166;
     k_il6r_b = 0.7128317112588906;
     k_il6r_ub = 0.4651946733609629;
     k_il6r_a = 0.6118072709651382;
     k_il6r_da = 0.523238216165728;
     IL6_m = 1.8347868135057839;
     k_il6m_il6 = 0.21085307745239146;
     k_il6_p = 0.15588824779440813;
     k_il6_stat3_a = 20.36297930994367;
     kd_il6_stat3_a = 243.32105159256753;
     o_il6_stat3_a = 0.09585643406063904;
     k_il6_pi3k_a = 298.2198274470974;
     kd_il6_pi3k_a = 474.9074578169857;
     o_il6_pi3k_a = 0.6477359906173604;
     k_nfkb_il6_p = 64.48561570146519;
     kd_nfkb_il6_p = 54997.2624657116;
     o_nfkb_il6_p = 0.4201709508109893
     
    // IL6 assignements
    ee = 0.0001
    nIL6_m := IL6_m/(IL6_m_0+ee);
    nIL6 := IL6/(IL6_0+ee);
    npSTAT3 := pSTAT3/(pSTAT3_0+ee);
    F_il6_stat3_a := 1+k_il6_stat3_a*(pIL6_R_JACK-pIL6_R_JACK_0)/(pIL6_R_JACK-pIL6_R_JACK_0 + kd_il6_stat3_a) - o_il6_stat3_a
    F_il6_pi3k_a := 1+k_il6_pi3k_a*(pIL6_R_JACK-pIL6_R_JACK_0)/(pIL6_R_JACK-pIL6_R_JACK_0 + kd_il6_pi3k_a) - o_il6_pi3k_a
    F_nfkb_il6_p := 1+k_nfkb_il6_p*(PP.NFKB_n-NFKB_n_0)/(PP.NFKB_n-NFKB_n_0+kd_nfkb_il6_p) - o_nfkb_il6_p;
    
#####------------------ IL8 ------------------------------------## 
    // IL8 transcriptional stimulation 
    IL8_v11: $IL8_prod + PP.NFKB_n + PP.AP1_n => IL8_m + PP.NFKB_n + PP.AP1_n; k_il8_p+(F_nfkb_il8_p-1)*(F_nfkb_il8_p>1)+(F_ap1_il8_p-1)*(F_ap1_il8_p>1);
    # IL8_v11: $IL8_prod + PP.NFKB_n + PP.AP1_n => IL8_m + PP.NFKB_n + PP.AP1_n; k_il8_p;
    IL8_v12: IL8_m => IL8; k_il8m_il8*IL8_m*(IL8_m>0);
    IL8_v13: IL8 => ; k_il8_d*IL8;
    IL8_v14: IL8_m => ; k_il8m_d*IL8_m;
    
    # // IL8 signaling pathway 
    IL8_v21: $IL8R + IL8 => IL8_R; k_il8r_b*IL8*(IL8>0);
    IL8_v22: IL8_R => $IL8R + IL8; k_il8r_ub*IL8_R*(IL8_R>0);
    IL8_v23: IL8_R => pIL8_R; k_il8r_a*IL8_R*(IL8_R>0);
    IL8_v24: pIL8_R => IL8_R; k_il8r_da*pIL8_R*(pIL8_R>0);
    ### IL8 upregulates IRAK recruitment/production ###
    IL8_v25: $PP.irak4_prod + pIL8_R => IRAK4 + pIL8_R;  (F_il8_irak-1)*(F_il8_irak>1);
    ### IL8 activates Rho GTPase
    IL8_26: $RHO + pIL8_R => aRHO + pIL8_R; (F_rho_a-1)*(F_rho_a>1);
    PP.IKB_NFKB + aRHO -> NFKB + IKB + aRHO; (F_rho_nfkb_a-1)*(F_rho_nfkb_a>1)
    IL8_27: PP.PI3K + aRHO => PP.pPI3K + aRHO; (F_rho_pi3k_a-1)*(F_rho_pi3k_a>1);
    IL8_28: $STAT3_s + aRHO => pSTAT3 + aRHO; (F_rho_stat3_a-1)*(F_rho_stat3_a>1)
       
    // IL8 Variables:
    IL8_prod = 1;
    k_il8_d = 0.693/(4*60); # 4h half life
    k_il8m_d = 0.693/(10*60); # 10 half life
    IL8R = 1
    RHO = 1
    aRHO = 0
    
IL8 = 8.812673844671163;
IL8_m = 0.021007864853922698;
k_il8_p = 0.013385842315493335;
k_il8m_il8 = 0.6198133323503698;
IL8_R = 4.82893203252714;
pIL8_R = 2.894467456365114;
k_il8r_b = 0.9662534764186338;
k_il8r_ub = 0.6196709975901828;
k_il8r_a = 0.6851200290039197;
k_il8r_da = 0.8813486044040026;
kd_il8_irak_p = 75192.46845579546;
k_il8_irak_p = 12538.604430541294;
o_il8_irak_p = 0.40763272189839656

    k_nfkb_il8_p = 125.12279150140773;
    kd_nfkb_il8_p = 75848.72994092123;
    o_nfkb_il8_p = 0.1649642604709619;
    k_ap1_il8_p = 7.6818306966347905;
    kd_ap1_il8_p = 60042.60509605664;
    o_ap1_il8_p = 0.24217754216712167

  


k_rho_a = 29.052917002477976;
kd_rho_a = 8129.6230662774815;
o_rho_a = 0.9276414555739282;
     k_rho_nfkb_a = 77972.8218479906;
     kd_rho_nfkb_a = 10575.98479567533;
     o_rho_nfkb_a = 0.0758851472315869;
 

 k_rho_pi3k_a = 95741.02954720988;
 kd_rho_pi3k_a = 1557.6419724575026;
 o_rho_pi3k_a = 0.03610475397019314;
 k_rho_stat3_a = 21946.417237773316;
 kd_rho_stat3_a = 37308.68663289963;
 o_rho_stat3_a = 0.24260054422200417;
 k_rho_jnk_a = 83917.9235530874;
 kd_rho_jnk_a = 22770.275011467245;
 o_rho_jnk_a = 0.140118655310842
    
    
    
    
    

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


    
    F_rho_a := 1+ k_rho_a*(pIL8_R-pIL8_R_0)/((pIL8_R-pIL8_R_0)+kd_rho_a) - o_rho_a
    F_rho_pi3k_a := 1+ k_rho_pi3k_a*(aRHO)/(aRHO+kd_rho_pi3k_a) - o_rho_pi3k_a
    F_rho_stat3_a := 1 + k_rho_stat3_a*aRHO/(aRHO+kd_rho_stat3_a) - o_rho_stat3_a
    F_rho_nfkb_a := 1 + k_rho_nfkb_a*aRHO/(aRHO+kd_rho_nfkb_a) - o_rho_nfkb_a
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
    pSTAT3_0 = pSTAT3;
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
    at (time > 0): pSTAT3_0 = pSTAT3;
    at (time > 0): pIL6_R_JACK_0 = pIL6_R_JACK
end

"""

ILs_model = te.loada(ILs_model_str)
# Zhao_model = te.loadSBMLModel(dirs.dir_Zhao_model)
# species_IDs = Zhao_model.getFloatingSpeciesIds()
LPS_model = te.loadSBMLModel(dirs.dir_LPS_model)
species_IDs = LPS_model.getFloatingSpeciesIds()
ILs_model_m = common.assign_surrogate_names(ILs_model,species_IDs)
ILs_model_m.exportToSBML(dirs.dir_ILs_model)


