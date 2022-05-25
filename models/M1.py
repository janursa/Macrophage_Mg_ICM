
import tellurium as te
import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
from tools import dirs, tools

M1_model_str = """
model Mg_model()
  compartment comp1;
  species Mg_e in comp1, Mg_f in comp1, TRPM in comp1, ATP in comp1, Mg_ATP in comp1;
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
    nATP := ATP/ATP_0;
   // Assignment Rules: functions
    F_trpm_i := (TRPM)/(TRPM_0);
    F_trpm_i := 1;

  // Reactions:
    # v1: Mg_e + TRPM => Mg_f + TRPM; k_m_i*Mg_e*(Mg_e > 0)*F_trpm_i - k_m_e*Mg_f*(Mg_f > 0);
    M1_v1: Mg_e => Mg_f ; k_m_i*(Mg_e-Mg_f)*(Mg_e > Mg_f);
    M1_v2: $ATP_prod + Mg_f => ATP + Mg_f; k_atp_p0 + k_atp_pm*Mg_f;
    M1_v3: ATP => $cons; k_atp_c*ATP;
    
    M1_v4: Mg_f +  ATP -> Mg_ATP; k_matp_b*Mg_f*(Mg_f > 0)*ATP*(ATP > 0) - k_matp_ub*Mg_ATP*(Mg_ATP > 0);
    M1_v5: Mg_f +  $IM -> Mg_IM; k_mim_b*Mg_f*(Mg_f > 0)*IM*(IM > 0) - k_mim_ub*Mg_IM*(Mg_IM > 0);

    M1_v6: $TRPM_prod + Mg_ATP => TRPM_n + Mg_ATP; k_ntrpm_p*TRPM_prod*(Mg_ATP)^n_ntrpm_p;
    M1_v7: TRPM_n => TRPM; k_ntrpm_t*TRPM_n;
    M1_v8: TRPM => M7CK; k_trpm_s*TRPM;
    M1_v9: M7CK -> M7CK_n; k_m7ck_t1*M7CK*(M7CK > 0) - k_m7ck_t2*M7CK_n*(M7CK_n > 0);
    M1_v10: M7CK => $deg; k_m7ck_deg*M7CK*(M7CK>0);
    M1_v11: $h3s10_prod -> H3S10; k_h3s10_p;
    M1_v12: H3S10 -> deg; k_h3s10_d*H3S10;
    M1_v13: H3S10 + M7CK_n -> pH3S10 + M7CK_n; k_h3s10_a*H3S10*M7CK_n*(H3S10 > 0) - k_h3s10_da*pH3S10*(pH3S10 > 0);

  // Species initializations:
    Mg_e = 0.8;
    Mg_f = 0.8;
    TRPM = 1;
    IM = 1.0199739743020637;
    Mg_ATP = 17.57597150605035;
    Mg_IM = 0.12931885964186912;
    ATP = 4.4200492287922275;
    TRPM_n = 1;
    cons = 0;
    ATP_prod = 1;
    TRPM_prod = 1;
    h3s10_prod = 1;
    IM_prod = 1;
    M7CK = 1;
    M7CK_n = 1;
    deg = 0;
    H3S10 = 1;
    pH3S10 = 1;
 
   # assignments initialization
    
   Mg_e_0 = Mg_e;
   Mg_0 = Mg_f + Mg_ATP + Mg_IM;
   Mg_f_0 = Mg_f;
   Mg_IM_0 = Mg_IM;
   TRPM_0 = TRPM; 
   Mg_ATP_0 = Mg_ATP;
   TRPM_n_0 = TRPM_n;
   M7CK_n_0 = M7CK_n;
   pH3S10_0 = pH3S10;
   ATP_0 = ATP;
   at (time > 0): Mg_f_0 = Mg_f;
   at (time > 0): pH3S10_0 = pH3S10;
   at (time > 0): Mg_0 = Mg_f + Mg_ATP + Mg_IM;
   at (time > 0): Mg_f_0 = Mg_f;
   at (time > 0): Mg_IM_0 = Mg_IM;
   at (time > 0): TRPM_0 = TRPM;
   at (time > 0): Mg_ATP_0 = Mg_ATP;
   at (time > 0): TRPM_n_0 = TRPM_n;
   at (time > 0):  M7CK_n_0 = M7CK_n;
   at (time > 0):  Mg_e_0 = Mg_e;
   at (time > 0):  ATP_0 = ATP;
  
  
  // Variable initializations:
      
      
        k_atp_p0 = 0.002255124807943189;
        k_atp_pm = 8.491754273290441;
        k_atp_c = 0.590124477903334;
        k_m_i = 0.35411879710555166;
        k_matp_b = 0.8632245228883954;
        k_matp_ub = 0.45252482532151195;
        k_mim_b = 0.15981115222635706;
        k_mim_ub = 0.9897381814819478
     

    
k_ntrpm_p = 0.0010111311947780088;
n_ntrpm_p = 2.303172314228545;
k_ntrpm_t = 0.7410061040831892;
k_trpm_s = 0.9968437547905072;
k_m7ck_t1 = 0.9219270599665684;
k_m7ck_t2 = 0.7127451109467265;
k_m7ck_deg = 0.9623775964715906
  
k_h3s10_p  = 7.035444348793543;
k_h3s10_a  = 0.9589420876901643;
k_h3s10_d  = 0.7221697589064646;
k_h3s10_da  = 0.807474039427079;
n_h3s10_a  = 0.7377166985974215
  

end
"""
model_sbml = te.loada(M1_model_str)
model_sbml.exportToSBML(dirs.dir_M1_model)
