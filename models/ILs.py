
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
    IFNGR is PP.IFNGR;
    IL4R is PP.IL4R;
    IRAK4 is PP.IRAK4;
    aTRAF6 is PP.aTRAF6;
    NFKB_n is PP.NFKB_n;
    AP1_n is PP.AP1_n;
    NFKB_n_0 is PP.NFKB_n_0;
    pCREB is PP.pCREB
    // Global equations
    _pSTAT3 => ; k_pstat3_d*_pSTAT3; 
    2 _pSTAT3 -> PP.pSTAT3D; k_pstat3_b*_pSTAT3^2 - k_pstat3_ub*PP.pSTAT3D;
    _pPI3K => ; k_ppi3k_d*_pPI3K
    PP.AKT1 + _pPI3K -> PP.pAKT1 + _pPI3K; K_ppi3k_a*PP.AKT1*_pPI3K - K_ppi3k_da*PP.pAKT1
   
    // Global variables 
    _pSTAT3 = 10 #local definition of _pSTAT3
    k_pstat3_d = 0.1
    k_pstat3_b = 0.1
    k_pstat3_ub = 0.1 
    _pPI3K = 10
    k_ppi3k_d = 0.1
    K_ppi3k_a = 0.1
    K_ppi3k_da = 0.1
    
     // IL6 signaling pathway 
    # IL6/R/JAK/STAT3
    IL6_v11: $IL6R_JACK + IL6 -> IL6_R_JACK; k_il6r_b*IL6*(IL6>0) - k_il6r_ub*IL6_R_JACK*(IL6_R_JACK>0)
    IL6_v12: IL6_R_JACK -> pIL6_R_JACK; k_il6r_a*IL6_R_JACK*(IL6_R_JACK>0) - k_il6r_da*pIL6_R_JACK*(pIL6_R_JACK>0)
    IL6_v13: $PP.STAT3 + pIL6_R_JACK => _pSTAT3 + pIL6_R_JACK; (F_il6_stat3_a-1)*(F_il6_stat3_a>1)*$PP.STAT3
    ## PI3K/Akt Pathway
    IL6_v14: $PP.PI3K +pIL6_R_JACK => _pPI3K + pIL6_R_JACK; (F_il6_pi3k_a-1)*(F_il6_pi3k_a>1)*PP.PI3K;
    
    // IL6 transcriptional stimulation 
    IL6_v21: $IL6_prod + PP.NFKB_n + PP.AP1_n + PP.pCREB => IL6_m + PP.NFKB_n + PP.AP1_n + PP.pCREB; k_il6_p+(F_nfkb_il6_p-1)*(F_nfkb_il6_p > 1)+(F_ap1_il6_p-1)*(F_ap1_il6_p>1)+(F_creb_il6_p-1)*(F_creb_il6_p>1);
    IL6_v22: IL6_m => IL6; k_il6m_il6*IL6_m*(IL6_m>0);
    IL6_v23: IL6 => ; k_il6_d*IL6;
    IL6_v24: IL6_m => ; k_il6m_d*IL6_m;

    // IL6 Variables
    ## Known ##
    k_il6_d = 0.693/(1*60); # 1h half life
    k_il6m_d = 0.693/(1*60); # 1h half life
    IL6R_JACK = 1
    ## unknown ## 
IL6 = 5.718738013069269;
IL6_R_JACK = 3.846369232849832;
pIL6_R_JACK = 5.287250355164202;
k_il6r_b = 0.0018293173288517206;
k_il6r_ub = 0.5354092623374117;
k_il6r_a = 0.003700948281175842;
k_il6r_da = 0.12753921902121917;
IL6_m = 6.283443153954376;
k_il6m_il6 = 0.16257285413784056;
k_il6_p = 0.13728047784941566;
     
     k_il6_stat3_a = 33.37388538700593;
     kd_il6_stat3_a = 6751.838577010346;
     o_il6_stat3_a = 0.10996297477107037;
     k_il6_pi3k_a = 814.8304465923425;
     kd_il6_pi3k_a = 8185.1658880053255;
     o_il6_pi3k_a = 0.3479259508892852;
     k_nfkb_il6_p = 56.92130543958485;
     kd_nfkb_il6_p = 56676.22713103319;
     o_nfkb_il6_p = 0.5580940352873225;
     
     k_ap1_il6_p = 11.398946355325336;
     kd_ap1_il6_p = 84278.15598566113;
     o_ap1_il6_p = 0.5956424887786982;
     k_creb_il6_p = 27.56655805407555;
     kd_creb_il6_p = 74668.6241501626;
     o_creb_il6_p = 0.09579108504083
    
    // IL6 assignements
    ee = 0.00001
    nIL6_m := IL6_m/(IL6_m_0+ee);
    nIL6 := IL6/(IL6_0+ee);
    npSTAT3 := _pSTAT3/(_pSTAT3_0+ee);
    npPI3K := _pPI3K/(_pPI3K_0+ee);
    diff_il6 := (pIL6_R_JACK - pIL6_R_JACK_0)*(pIL6_R_JACK-pIL6_R_JACK_0)
    F_il6_stat3_a := 1+k_il6_stat3_a * diff_il6/(diff_il6 + kd_il6_stat3_a) - o_il6_stat3_a
    F_il6_pi3k_a := 1+k_il6_pi3k_a * diff_il6/(diff_il6 + kd_il6_pi3k_a) - o_il6_pi3k_a
    diff_nfkb_s := (NFKB_n - NFKB_n_0)*(NFKB_n > NFKB_n_0)
    F_nfkb_il6_p := 1+k_nfkb_il6_p*diff_nfkb_s/(diff_nfkb_s + kd_nfkb_il6_p) - o_nfkb_il6_p;
    diff_ap1_s := (AP1_n - AP1_n_0)*(AP1_n>AP1_n_0)
    F_ap1_il6_p := 1+k_ap1_il6_p*diff_ap1_s/(diff_ap1_s + kd_ap1_il6_p) - o_ap1_il6_p;
    F_creb_il6_p := 1+k_creb_il6_p*PP.pCREB/(PP.pCREB + kd_creb_il6_p) - o_creb_il6_p;
    
#####------------------ IL8 ------------------------------------## 
    // IL8 transcriptional stimulation 
    IL8_v11: $IL8_prod + PP.NFKB_n + PP.AP1_n => IL8_m + PP.NFKB_n + PP.AP1_n; k_il8_p+(F_nfkb_il8_p-1)*(F_nfkb_il8_p>1)+(F_ap1_il8_p-1)*(F_ap1_il8_p>1);
    IL8_v12: IL8_m => IL8; k_il8m_il8*IL8_m*(IL8_m>0);
    IL8_v13: IL8 => ; k_il8_d*IL8;
    IL8_v14: IL8_m => ; k_il8m_d*IL8_m;
    
    // IL8 signaling pathway 
    IL8_v21: $IL8R + IL8 => IL8_R; k_il8r_b*IL8*(IL8>0) - k_il8r_ub*IL8_R*(IL8_R>0);
    IL8_v22: IL8_R => pIL8_R; k_il8r_a*IL8_R*(IL8_R>0) - k_il8r_da*pIL8_R*(pIL8_R>0);
    ### IL8 upregulates IRAK recruitment/production ###
    IL8_v23: $PP.irak4_prod + pIL8_R => PP.IRAK4 + pIL8_R;  (F_il8_irak_p-1)*(F_il8_irak_p>1);
    ### IL8 activates Rho GTPase
    IL8_24: $RHO + pIL8_R => aRHO + pIL8_R; k_rho_a*pIL8_R/(pIL8_R+kd_rho_a) ;
    IL8_25: aRHO => $RHO; k_rho_da*aRHO;
    IL8_26: PP.IKB_NFKB + aRHO => PP.NFKB + PP.IKB + aRHO; (F_rho_nfkb_a-1)*(F_rho_nfkb_a>1)*PP.IKB_NFKB
    IL8_27: $PP.PI3K + aRHO => _pPI3K + aRHO; (F_rho_pi3k_a-1)*(F_rho_pi3k_a>1)*PP.PI3K;
    IL8_28: $PP.STAT3 + aRHO => _pSTAT3 + aRHO; (F_rho_stat3_a-1)*(F_rho_stat3_a>1)*PP.STAT3
       
    // IL8 Variables:
    IL8_prod = 1;
    k_il8_d = 0.693/(4*60); # 4h half life
    k_il8m_d = 0.693/(10*60); # 10 half life
    IL8R = 1
    RHO = 1
    aRHO = 0
    
IL8 = 4.724288120656182;
IL8_m = 1.2138313631672166;
k_il8_p = 0.024397471013379857;
k_il8m_il8 = 0.14596703636802016;
IL8_R = 2.310405727945194;
pIL8_R = 6.666584229891921;
k_il8r_b = 0.6167265055192227;
k_il8r_ub = 0.19347421128048325;
k_il8r_a = 0.014547178907135416;
k_il8r_da = 0.46775013717966707;
  

    
k_nfkb_il8_p = 1.2380664811373663;
kd_nfkb_il8_p = 768.6185511351218;
o_nfkb_il8_p = 0.8408046136068019;
k_ap1_il8_p = 1.2774381921529425;
kd_ap1_il8_p = 786.3889374437415;
o_ap1_il8_p = 0.9967236511541081;
  
k_rho_a = 1.1167951937654834;
kd_rho_a = 2456.098246911869;
k_rho_da = 0.6171084280662366;

    k_il8_irak_p = 97560.46809730714;
    kd_il8_irak_p = 1003.5070641167404;
    o_il8_irak_p = 0.7770775119227259;
    k_rho_nfkb_a = 17585.928627588208;
    kd_rho_nfkb_a = 57935.08687938105;
    o_rho_nfkb_a = 6.648843175971475e-06
    
    

  


 

k_rho_pi3k_a = 1022.0687993213766;
kd_rho_pi3k_a = 98995.62566590015;
o_rho_pi3k_a = 0.004795232546251698;
k_rho_stat3_a = 6393.972581541672;
kd_rho_stat3_a = 83267.05020810368;
o_rho_stat3_a = 0.07719829702965487;
k_nfkb_il6_p = 3.0720042242873364;
kd_nfkb_il6_p = 66321.6855500037;
o_nfkb_il6_p = 0.09076515545524551
    
    
    
    
    

    // IL8 assignements
    nIL8_m := (IL8_m+ee)/(IL8_m_0+ee);
    nIL8 := (IL8+ee)/(IL8_0+ee);
    nIFNGR := (IFNGR+ee)/(IFNGR_0);
    nIL4R := (IL4R+ee)/(IL4R_0);
    nIRAK4 := (PP.IRAK4)/(IRAK4_0);
    naTRAF6 := (aTRAF6+ee)/(aTRAF6_0);
    npIL8_R := (pIL8_R+ee)/(pIL8_R_0+ee);
    nIL6 := (IL6+ee)/(IL6_0+ee);
    diff_il8_s := (pIL8_R - pIL8_R_0)*(pIL8_R>pIL8_R_0)
    F_il8_irak_p := 1+ k_il8_irak_p * diff_il8_s/(diff_il8_s + kd_il8_irak_p) - o_il8_irak_p;
    F_nfkb_il8_p := 1+ k_nfkb_il8_p * NFKB_n/(NFKB_n + kd_nfkb_il8_p) - o_nfkb_il8_p;
    F_ap1_il8_p := 1+k_ap1_il8_p*AP1_n/(AP1_n + kd_ap1_il8_p) - o_ap1_il8_p


    F_rho_pi3k_a := 1+ k_rho_pi3k_a*(aRHO)/(aRHO+kd_rho_pi3k_a) - o_rho_pi3k_a
    F_rho_stat3_a := 1 + k_rho_stat3_a*aRHO/(aRHO+kd_rho_stat3_a) - o_rho_stat3_a
    F_rho_nfkb_a := 1 + k_rho_nfkb_a*aRHO/(aRHO+kd_rho_nfkb_a) - o_rho_nfkb_a
    // initial conditions, both IL6 and IL8
    IL8_m_0 = IL8_m;
    IFNGR_0 = IFNGR;
    IL4R_0 = IL4R;
    IRAK4_0 = PP.IRAK4;
    aTRAF6_0 = aTRAF6;
    IL8_0 = IL8;
    pIL8_R_0 = pIL8_R;
    aRHO_0 = aRHO
    AP1_n_0 = AP1_n
    IL6_m_0 = IL6_m;
    IL6_0 = IL6;
    _pSTAT3_0 = _pSTAT3;
    _pPI3K_0 = _pPI3K;
    pIL6_R_JACK_0 = pIL6_R_JACK
    AP1_n_0 = PP.AP1_n
    pCREB_0 = PP.pCREB
    at (time > 0): IFNGR_0 = IFNGR;
    at (time > 0): IL4R_0 = IL4R;
    at (time > 0): IL8_m_0 = IL8_m;
    at (time > 0): IL8_0 = IL8;
    at (time > 0): aTRAF6_0 = aTRAF6; 
    at (time > 0): pIL8_R_0 = pIL8_R;
    at (time > 0): IL6_m_0 = IL6_m;
    at (time > 0): IL6_0 = IL6;
    at (time > 0): _pSTAT3_0 = _pSTAT3;
    at (time > 0): pIL6_R_JACK_0 = pIL6_R_JACK
    at (time > 0): _pPI3K_0 = _pPI3K 
end

"""

ILs_model = te.loada(ILs_model_str)
# Zhao_model = te.loadSBMLModel(dirs.dir_Zhao_model)
# species_IDs = Zhao_model.getFloatingSpeciesIds()
LPS_model = te.loadSBMLModel(dirs.dir_LPS_model)
species_IDs = LPS_model.getFloatingSpeciesIds()
ILs_model_m = common.assign_surrogate_names(ILs_model,species_IDs)
ILs_model_m.exportToSBML(dirs.dir_ILs_model)


