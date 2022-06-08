import tellurium as te
import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
from tools import dirs, common
combined_model_str = """
import "M1_sbml.xml";
import "ILs_sbml.xml";
# import "LPS_sbml.xml";
model combined()
    M1: Mg_model();
    PP: ILs_model();
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
    nATP is M1.nATP;
    Mg_ATP is M1.Mg_ATP;
    k_atp_p0 is M1.k_atp_p0;
    k_atp_pm is M1.k_atp_pm;
    k_atp_c is M1.k_atp_c;
    IM is M1.IM;
    Mg_IM is M1.Mg_IM;
    Mg_ATP_0 is M1.Mg_ATP_0;
    
    
    // adjustments to previous model:
    ### Mg downregulate IKB degradation ###
    PP.v226: PP.IKB + M1.Mg_f => deg + M1.Mg_f; PP.k226 * PP.IKB * F_mg_ikb_d;
    ### Mg pH3S10 upregulates IKB production ###
    $PP.Ikb_prod + M1.pH3S10 => PP.IKB + M1.pH3S10; (F_h3s10_ikb-1)*(F_h3s10_ikb>1);
    ### IL8 production is regulated by pH3S10
    $PP.IL8_prod + M1.pH3S10 -> PP.IL8_m + M1.pH3S10; (F_p3s10_il8_p-1)*(F_p3s10_il8_p>1);
    // new reactions
    
    
    // Variables:
    kd_mg_ikb_d = 88.948;
    
    kd_h3s10_ikb_p = 1.207;
    k_h3s10_ikb_p = 29.52;

    kd_h3s10_il8_p = 1;
    k_h3s10_il8_p = 1;
    
    o_h3s10_ikb = 1
    o_h3s10_il8 = 0

    
    
    // assignements
    F_mg_ikb_d := kd_mg_ikb_d/(Mg_f-Mg_f_0+kd_mg_ikb_d);
    F_h3s10_ikb := 1+k_h3s10_ikb_p*(pH3S10-pH3S10_0)/((pH3S10-pH3S10_0)+kd_h3s10_ikb_p) - o_h3s10_ikb
    F_p3s10_il8_p := 1+k_h3s10_il8_p*(pH3S10-pH3S10_0)/((pH3S10-pH3S10_0)+kd_h3s10_il8_p) - o_h3s10_il8;
end
"""

combined = te.loada(combined_model_str)


if False: # only labels of IL8 are replaced
    ILs_model = te.loadSBMLModel(dirs.dir_ILs_model)
    species_IDs = ILs_model.getFloatingSpeciesIds()
    combined_m = common.assign_surrogate_names(combined,species_IDs)
    combined_m.exportToSBML(dirs.dir_model)
if True: # both labels of IL8 and M1 are replaced
    ILs_model = te.loadSBMLModel(dirs.dir_ILs_model)
    Mg_model = te.loadSBMLModel(dirs.dir_M1_model)
    species_IDs_1 = ILs_model.getFloatingSpeciesIds()
    # LPS_model = te.loadSBMLModel(dirs.dir_LPS_model)
    # Mg_model = te.loadSBMLModel(dirs.dir_M1_model)
    # species_IDs_1 = LPS_model.getFloatingSpeciesIds()
    #species_IDs_2 = Mg_model.getFloatingSpeciesIds()
    combined_m_1 = common.assign_surrogate_names(combined,species_IDs_1)
    #combined_m_2 = tools.assign_surrogate_names(combined_m_1,species_IDs_2,prefix='M1_')
    #combined_m_2.exportToSBML(dirs.dir_model)
    combined_m_1.exportToSBML(dirs.dir_model)
    
