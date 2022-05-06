
import tellurium as te
import sys
sys.path.insert(0,'/Users/matin/Downloads/testProjs/intracellular_M')
from tools import dirs

M1_model_str = """
model Mg_model()
  compartment comp1;
  species $Mg_e in comp1, Mg_f in comp1, TRPM in comp1, ATP in comp1, Mg_ATP in comp1;
  species $Mg in comp1, Mg_IM in comp1, TRPM_n in comp1, $cons in comp1, $ATP_prod in comp1;
  species $TRPM_prod in comp1, $nMg in comp1, $nMg_f in comp1;
  species $IM in comp1, M7CK in comp1, M7CK_n in comp1, $deg in comp1;
  species $H3S10 in comp1, pH3S10 in comp1, $nMg_IM in comp1, $nMg_ATP in comp1;
  species  h3s10_prod in comp1;


  // surragate names

  // Assignments: normalization
    Mg := Mg_f + Mg_ATP + Mg_IM;
    nMg := Mg/Mg_0;
    nMg_f := Mg_f/Mg_f_0;
    nMg_IM := Mg_IM/Mg_IM_0;
    nMg_ATP := Mg_ATP/Mg_ATP_0;
    npH3S10 := pH3S10/pH3S10_0;
    nTRPM := TRPM/TRPM_0;
    nTRPM_n := TRPM_n/TRPM_n_0;
    nM7CK_n := M7CK_n/M7CK_n_0;
   // Assignment Rules: functions
   F_trpm_i := (TRPM)/(TRPM_0);
    F_trpm_i := 1;
    F_mg_i := ((Mg_0+kd_mg_i)/(Mg+kd_mg_i))^n_mg_i;

  // Reactions:
   v1: $Mg_e + TRPM => Mg_f + TRPM; k_m_i*Mg_e*(Mg_e > 0)*F_mg_i*F_trpm_i - k_m_e*Mg_f*(Mg_f > 0);
   v2: $ATP_prod + Mg_f => ATP + Mg_f; k_atp_p0 + k_atp_pm*Mg_f;
   v3: ATP => $cons; k_atp_c*ATP;

   v4: Mg_f +  ATP -> Mg_ATP; k_matp_b*Mg_f^n_matp_b*(Mg_f > 0)*ATP*(ATP > 0) - k_matp_ub*Mg_ATP*(Mg_ATP > 0);
   v5: Mg_f +  $IM -> Mg_IM; k_mim_b*Mg_f^n_mim_b*(Mg_f > 0)*IM*(IM > 0) - k_mim_ub*Mg_IM*(Mg_IM > 0);

   v6: $TRPM_prod + Mg_ATP => TRPM_n + Mg_ATP; k_ntrpm_p*TRPM_prod*(Mg_ATP)^n_ntrpm_p;
    v7: TRPM_n => TRPM; k_ntrpm_t*TRPM_n^n_ntrpm_t;
    v8: TRPM => M7CK; k_trpm_s*TRPM;
    v9: M7CK -> M7CK_n; k_m7ck_t1*M7CK^n_m7ck_t*(M7CK > 0) - k_m7ck_t2*M7CK_n*(M7CK_n > 0);
    v10: M7CK => $deg; k_m7ck_deg*M7CK*(M7CK>0);
    v11: $h3s10_prod -> H3S10; k_h3s10_p;
    v12: H3S10 -> deg; k_h3s10_d*H3S10;
    v13: H3S10 + M7CK_n -> pH3S10 + M7CK_n; k_h3s10_a*H3S10*M7CK_n^n_h3s10_a*(H3S10 > 0) - k_h3s10_da*pH3S10*(pH3S10 > 0);

  // Species initializations:
    Mg_e = 0.8;
    Mg_f = 0.8;
    TRPM = 1;
    ATP = 24.348;
    Mg_ATP  =  9.013479998793457;
     Mg_IM  =  2.810066255278443;
    TRPM_n = 1;
    cons = 0;
    ATP_prod = 1;
    TRPM_prod = 1;
    h3s10_prod = 1;
    IM = 88.6098019188817;
    IM_prod = 1;
    M7CK = 1;
    M7CK_n = 1;
    deg = 0;
    H3S10 = 1;
    pH3S10 = 1;
 
   # assignments initialization
   pH3S10_0 = pH3S10;
   Mg_0 = Mg_f + Mg_ATP + Mg_IM;
   Mg_f_0 = Mg_f;
   Mg_IM_0 = Mg_IM;
   TRPM_0 = TRPM; 
   Mg_ATP_0 = Mg_ATP;
   TRPM_n_0 = TRPM_n;
   M7CK_n_0 = M7CK_n;
   at (time > 0): Mg_f_0 = Mg_f;
   at (time > 0): pH3S10_0 = pH3S10;
   at (time > 0): Mg_0 = Mg_f + Mg_ATP + Mg_IM;
   at (time > 0): Mg_f_0 = Mg_f;
   at (time > 0): Mg_IM_0 = Mg_IM;
   at (time > 0): TRPM_0 = TRPM;
   at (time > 0): Mg_ATP_0 = Mg_ATP;
   at (time > 0): TRPM_n_0 = TRPM_n;
   at (time > 0):  M7CK_n_0 = M7CK_n;
  
  
  // Variable initializations:
      
      
     
      
     k_m_i  =  0.9771043536710443;
     k_m_e  =  0.6730171850345992;
     k_matp_b  =  0.5055075330283233;
     k_matp_ub  =  0.9464132182191956;
     k_mim_b  =  0.45017904044421053;
     k_mim_ub  =  0.479558321478773;
     n_matp_b  =  0.10607032768807301;
     n_mim_b  =  9.970539499072103;
     kd_mg_i  =  0.2408360732050383
     

  
  k_atp_p0 = 9.8145812614764;
  k_atp_c = 0.62216362042029;
  k_atp_pm = 1.13;
  n_mg_i = 1;


k_ntrpm_p  = 0.6664977332147098
n_ntrpm_p  = 1
k_ntrpm_t  = 0.7599026682313925
n_ntrpm_t  = 1
k_trpm_s  = 0.7578004751224269
k_m7ck_t1  = 0.5638119158774044
k_m7ck_t2  = 0.4542293364622845
n_m7ck_t  = 1
k_m7ck_deg  = 0.8500625313026188
  
  k_h3s10_p  = 7.035444348793543 ;
k_h3s10_a  = 0.9589420876901643 ;
k_h3s10_d  = 0.7221697589064646 ;
k_h3s10_da  = 0.807474039427079 ;
n_h3s10_a  = 1
  

end
"""
model_sbml = te.loada(M1_model_str)
model_sbml.exportToSBML(dirs.dir_M1_model)
