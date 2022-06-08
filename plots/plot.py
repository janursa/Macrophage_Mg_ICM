import tellurium as te
import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
import sys
import numpy as np
import json
import matplotlib.pyplot as plt
import tellurium as te
plt.rcParams["font.family"] = "serif"
plt.style.use('seaborn-deep')
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
from data.observations import observations,t2m,select_obs
from models.params import fixed_params
from tools import dirs
from models.models import Macrophage
from plots import funcs 

target_package = 'ILs_2'
params = {**fixed_params}
def reload_params(params): # apply inferred params
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
    return params
# params = reload_params(params)
# params ={
#     "IL6": 8.441726867774754,
#     "IL6_R_JACK": 2.7002660692419136,
#     "pIL6_R_JACK": 3.5276321274421223,
#     "k_il6r_b": 0.0013278248494205114,
#     "k_il6r_ub": 0.6568592484751244,
#     "k_il6r_a": 0.018089666126617165,
#     "k_il6r_da": 0.7730647818263322,
#     "IL6_m": 3.0698231828018265,
#     "k_il6m_il6": 0.031416216121147256,
#     "k_il6_p": 0.1294112407915442,
#     "k_il6_stat3_a": 100000,
#     "kd_il6_stat3_a": 6486.436550851764,
#     "o_il6_stat3_a": 0.08432958682693681,
#     "k_il6_pi3k_a": 4331.081514377869,
#     "kd_il6_pi3k_a": 3392.142897763703,
#     "o_il6_pi3k_a": 0.09850423997068924,
#     "k_nfkb_il6_p": 912.3406319572266,
#     "kd_nfkb_il6_p": 12135.648768289873,
#     "o_nfkb_il6_p": 0.18340659356693279,
#     "k_ap1_il6_p": 473.75366659242803,
#     "kd_ap1_il6_p": 35471.14238828033,
#     "o_ap1_il6_p": 0.0631708369263767,
#     "IL8": 8.279515486677937,
#     "IL8_m": 9.438754980852947,
#     "k_il8_p": 0.5606364069221925,
#     "k_il8m_il8": 0.4910973278526436,
#     "IL8_R": 2.955364943276446,
#     "pIL8_R": 9.54426163654274,
#     "k_il8r_b": 0.6352444750409134,
#     "k_il8r_ub": 0.004243082189243741,
#     "k_il8r_a": 0.000542673717717268,
#     "k_il8r_da": 0.3461950078853718,
#     "kd_il8_irak_p": 19930.818813762686,
#     "k_il8_irak_p": 681694.8071137026,
#     "o_il8_irak_p": 0.5562148609288003,
#     "k_nfkb_il8_p": 1.076335514576897,
#     "kd_nfkb_il8_p": 793.6925052258482,
#     "o_nfkb_il8_p": 0.8367114238943978,
#     "k_ap1_il8_p": 2.342419953801823,
#     "kd_ap1_il8_p": 443.46183878108263,
#     "o_ap1_il8_p": 0.7780485495255755,
#     "k_rho_a": 9.947025761519058,
#     "kd_rho_a": 7659.4526077809915,
#     "k_rho_da": 0.5509902079943436,
#     "k_rho_nfkb_a": 32137.478732923297,
#     "kd_rho_nfkb_a": 31394.118065700502,
#     "o_rho_nfkb_a": 0.017073113000682538,
#     "k_rho_pi3k_a": 11677.916209752191,
#     "kd_rho_pi3k_a": 9817.371483185714,
#     "o_rho_pi3k_a": 0.022718007150504982,
#     "k_rho_stat3_a": 1778.6403685658934,
#     "kd_rho_stat3_a": 15410.829235952398,
#     "o_rho_stat3_a": 0.22120538822280533
# }

# params = {
#     "IL6": 8.207368038808369,
#     "IL6_R_JACK": 7.382108591939595,
#     "pIL6_R_JACK": 5.560866674872395,
#     "k_il6r_b": 0.002298941353970352,
#     "k_il6r_ub": 0.9882144525182053,
#     "k_il6r_a": 0.03277794458358041,
#     "k_il6r_da": 0.8723730926293621,
#     "IL6_m": 6.381986443077651,
#     "k_il6m_il6": 0.013277007612750302,
#     "k_il6_p": 0.15742142651960833,
#     "k_il6_stat3_a": 1963.1660609851292,
#     "kd_il6_stat3_a": 4911.542847755854,
#     "o_il6_stat3_a": 0.016574160820316652,
#     "k_il6_pi3k_a": 2521.645887466939,
#     "kd_il6_pi3k_a": 8107.674211697587,
#     "o_il6_pi3k_a": 0.04816035173026567,
#     "k_nfkb_il6_p": 313.3297337148234,
#     "kd_nfkb_il6_p": 77816.03818838934,
#     "o_nfkb_il6_p": 0.03538467594482242,
#     "k_ap1_il6_p": 907.0639817497131,
#     "kd_ap1_il6_p": 36805.76943191598,
#     "o_ap1_il6_p": 0.03546404423264837,
#     "IL8": 6.7579420504668555,
#     "IL8_m": 1.4033001183761207,
#     "k_il8_p": 0.11165470364522889,
#     "k_il8m_il8": 0.683398153712759,
#     "IL8_R": 3.7676648232224226,
#     "pIL8_R": 6.248804605926299,
#     "k_il8r_b": 0.14318054909835387,
#     "k_il8r_ub": 0.0023977389843861663,
#     "k_il8r_a": 0.000584870294840889,
#     "k_il8r_da": 0.40736133584619816,
#     "kd_il8_irak_p": 13552.136489492143,
#     "k_il8_irak_p": 549753.5262566409,
#     "o_il8_irak_p": 0.11195424837794171,
#     "k_nfkb_il8_p": 1.1946604913257488,
#     "kd_nfkb_il8_p": 379.39674667540424,
#     "o_nfkb_il8_p": 0.534849701679407,
#     "k_ap1_il8_p": 1.5367389280679618,
#     "kd_ap1_il8_p": 587.8842153676262,
#     "o_ap1_il8_p": 0.8372440286445306,
#     "k_rho_a": 19.35251181267097,
#     "kd_rho_a": 8750.580123329362,
#     "k_rho_da": 0.9345189057052146,
#     "k_rho_nfkb_a": 24924.32199104117,
#     "kd_rho_nfkb_a": 41345.398346111295,
#     "o_rho_nfkb_a": 0.0004653202895435471,
#     "k_rho_pi3k_a": 82032.6351804462,
#     "kd_rho_pi3k_a": 34020.02849984018,
#     "o_rho_pi3k_a": 0.0875868826740035,
#     "k_rho_stat3_a": 3815.8623298855327,
#     "kd_rho_stat3_a": 87534.62610511703,
#     "o_rho_stat3_a": 0.12116072144134876
# }

# params = {
#     "IL6": 5.4370102257715365,
#     "IL6_R_JACK": 2.463356323589118,
#     "pIL6_R_JACK": 1.849107258284274,
#     "k_il6r_b": 0.44075140849316885,
#     "k_il6r_ub": 0.6493205249598956,
#     "k_il6r_a": 0.4270556206063457,
#     "k_il6r_da": 0.22582593386762706,
#     "IL6_m": 4.576711588346838,
#     "k_il6m_il6": 0.3932979489856194,
#     "k_il6_p": 3.4710282869815017,
#     "k_il6_stat3_a": 1951.121971969158,
#     "kd_il6_stat3_a": 7318.074225489829,
#     "o_il6_stat3_a": 0.42793367800462123,
#     "k_il6_pi3k_a": 1266.540409509294,
#     "kd_il6_pi3k_a": 4598.969155571763,
#     "o_il6_pi3k_a": 0.4849854774932578,
#     "k_nfkb_il6_p": 33.93520028798429,
#     "kd_nfkb_il6_p": 72442.04587322613,
#     "o_nfkb_il6_p": 0.07627002255439408,
#     "k_ap1_il6_p": 240.69769424843963,
#     "kd_ap1_il6_p": 84333.84255485586,
#     "o_ap1_il6_p": 0.008153212009223143,
#     "k_creb_il6_p": 11.408246267173212,
#     "kd_creb_il6_p": 72487.83964114796,
#     "o_creb_il6_p": 0.5165947002014888,
#     "IL8": 9.007692694221326,
#     "IL8_m": 0.9209962714924638,
#     "k_il8_p": 0.01726296193977106,
#     "k_il8m_il8": 0.6214970639297708,
#     "IL8_R": 4.19077277865208,
#     "pIL8_R": 3.8213439139543124,
#     "k_il8r_b": 0.02161783977154752,
#     "k_il8r_ub": 0.496702850242793,
#     "k_il8r_a": 0.35558316753820685,
#     "k_il8r_da": 0.768923681728841,
#     "kd_il8_irak_p": 282522.7279195321,
#     "k_il8_irak_p": 911956.2708446574,
#     "o_il8_irak_p": 0.31236589163039247,
#     "k_nfkb_il8_p": 1.0381760913336393,
#     "kd_nfkb_il8_p": 979.5911484020489,
#     "o_nfkb_il8_p": 0.9159888394230121,
#     "k_ap1_il8_p": 1.2691172150609304,
#     "kd_ap1_il8_p": 456.9686461242685,
#     "o_ap1_il8_p": 0.9182426473795741,
#     "k_rho_a": 3.411447676886212,
#     "kd_rho_a": 6871.233986011801,
#     "k_rho_da": 0.4935009817946129,
#     "k_rho_nfkb_a": 18512.259719458943,
#     "kd_rho_nfkb_a": 71035.36196641353,
#     "o_rho_nfkb_a": 0.00044230546029566664,
#     "k_rho_pi3k_a": 8818.378265969557,
#     "kd_rho_pi3k_a": 99056.02852261194,
#     "o_rho_pi3k_a": 0.22113616205429887,
#     "k_rho_stat3_a": 2279.361710806342,
#     "kd_rho_stat3_a": 47182.69054735674,
#     "o_rho_stat3_a": 0.06661670957536459,
#     "_pSTAT3": 1.648340019534821,
#     "k_pstat3_d": 0.5780014457545202,
#     "k_pstat3_b": 0.19856751836254527,
#     "k_pstat3_ub": 0.7290493085817544,
#     "_pPI3K": 44.322216620555416,
#     "k_ppi3k_d": 0.6530846648225805,
#     "K_ppi3k_a": 0.7698986492411686,
#     "K_ppi3k_da": 0.1663185949643658
# }
print(params)
# params['k_il6_stat3_a'] = 10
# params['k_il8_irak_p'] = 10

print('t2m: {} '.format(t2m))
flags = [
    'M1',
    # 'LPS',
    # 'IL6',
    # 'IL8',
    # 'combined'
]

if 'M1' in flags : 
    model_t = 'combined'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('M1 is plotting')
    funcs.P1_eq(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    funcs.P11_plot(model_sbml=model_sbml,params=params,observations=observations)
    funcs.P12_plot (model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
if 'LPS' in flags : 
    model_t = 'LPS'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('LPS is plotting')
    funcs.LPS_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations) 
if 'IL6' in flags: 
    model_t = 'ILs'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('IL6 is plotting')
    fig = funcs.P2_IL6_IC_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        
    fig = funcs.P2_IL6_CYs_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        

if 'IL8' in flags: 
    model_t = 'ILs'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('IL8 is plotting')
    fig2 = funcs.P2_eq(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        
    fig2 = funcs.P2_ICs_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        
    # fig2 = funcs.P2_receptors_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    # fig2 = funcs.P23_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)

if 'combined' in flags : 
    model_t = 'combined'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('Combined is plotting')
    fig = funcs.P3_eq_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P3_IKB_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P3_NFKB_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P3_cytokines1_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P3_cytokines2_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P3_cytokines3_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P3_NFKB_14d_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)

    # _dir = os.path.join(dirs.dir_outputs,'plots','P3.png')