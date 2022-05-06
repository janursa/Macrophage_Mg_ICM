import tellurium as te
import sys
sys.path.insert(0,'/Users/matin/Downloads/testProjs/intracellular_M')
from tools import dirs, tools
combined_model_str = """
import "M1_sbml.xml";
import "IL8_sbml.xml";
model combined()
    M1: Mg_model();
    PP: IL8_model();
    compartment comp1;
    
    // surragate variables
    Mg_f_0 is M1.Mg_f_0;
    Mg_f is M1.Mg_f;
    pH3S10 is M1.pH3S10;
    pH3S10_0 is M1.pH3S10_0;
    Mg_e is M1.Mg_e;
    Mg is M1.Mg;
    nMg_ATP is M1.nMg_ATP;
    nTRPM is M1.nTRPM;
    nTRPM_n is M1.nTRPM_n;
    nM7CK_n is M1.nM7CK_n;
    nMg_f is M1.nMg_f;
    nMg is M1.nMg;
    npH3S10 is M1.npH3S10;
    ATP is M1.ATP;
    Mg_ATP is M1.Mg_ATP;
    k_atp_p0 is M1.k_atp_p0;
    k_atp_pm is M1.k_atp_pm;
    k_atp_c is M1.k_atp_c;
    IM is M1.IM;
    Mg_IM is M1.Mg_IM;
    
    // adjustments to previous model:
    ### Mg downregulate IKB degradation ###
    PP.v226: PP.IKB + M1.Mg_f => deg + M1.Mg_f; PP.k226*PP.IKB* F_mg_ikb_d;
    ### Mg pH3S10 upregulates IKB production ###
    PP.v221: $PP.Ikb_prod + PP.NFKB_n + M1.pH3S10 => PP.IKB + PP.NFKB_n + M1.pH3S10; PP.k221*PP.Ikb_prod*PP.NFKB_n*F_h3s10_ikb;
    ### IL8 production is regulated by pH3S10
    PP.v4: $PP.IL8_prod + PP.nNFKB_n + M1.pH3S10 -> PP.IL8 + PP.nNFKB_n + M1.pH3S10; PP.k_il8_p*PP.F_nfkb_il8_p* F_p3s10_il8_p;

    // new reactions
    
    
    // Variables:
    kd_mg_ikb_d = 88.948;
    n_mg_ikb_d = 23.43;
    
    kd_h3s10_ikb_p = 1.207;
    n_h3s10_ikb_p = 29.52;

    kd_h3s10_il8_p = 1;
    n_h3s10_il8_p = 1;

    
    
    // assignements
    F_mg_ikb_d := ((Mg_f_0+kd_mg_ikb_d)/(Mg_f+kd_mg_ikb_d))^n_mg_ikb_d;
    F_h3s10_ikb := ((pH3S10+kd_h3s10_ikb_p)/(pH3S10_0+kd_h3s10_ikb_p))^n_h3s10_ikb_p
    # F_p3s10_il8_p := ((M1.pH3S10+kd_h3s10_il8_p)/(M1.pH3S10_0+kd_h3s10_il8_p))^n_h3s10_il8_p;
    F_p3s10_il8_p := M1.pH3S10^n_h3s10_il8_p/(M1.pH3S10_0^n_h3s10_il8_p+kd_h3s10_il8_p);
end
"""

combined = te.loada(combined_model_str)


if True: # only labels of IL8 are replaced
    IL8_model = te.loadSBMLModel(dirs.dir_IL8_model)
    species_IDs = IL8_model.getFloatingSpeciesIds()
    combined_m = tools.assign_surrogate_names(combined,species_IDs)
    combined_m.exportToSBML(dirs.dir_model)
if False: # both labels of IL8 and M1 are replaced
    IL8_model = te.loadSBMLModel(dirs.dir_IL8_model)
    Mg_model = te.loadSBMLModel(dirs.dir_M1_model)
    species_IDs_1 = IL8_model.getFloatingSpeciesIds()
    #species_IDs_2 = Mg_model.getFloatingSpeciesIds()
    combined_m_1 = tools.assign_surrogate_names(combined,species_IDs_1)
    #combined_m_2 = tools.assign_surrogate_names(combined_m_1,species_IDs_2,prefix='M1_')
    #combined_m_2.exportToSBML(dirs.dir_model)
    combined_m_1.exportToSBML(dirs.dir_model)
    