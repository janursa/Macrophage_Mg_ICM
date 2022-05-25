import tellurium as te
import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
from tools import dirs, tools
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
from tools import dirs, tools
from models.models import Macrophage
from plots import funcs 

params = {**fixed_params}
def reload_params(params): # apply inferred params
    target_package = 'IL8'
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
    return params
params = reload_params(params)
# params = {
#     "IL8_m": 6.571315008130573,
#     "k_il8_p": 6.837693370999659,
#     "kd_nfkb_il8_p": 6869.739497597992,
#     "IL8R": 84497.71197384532,
#     "IL8_R": 79822.75472764584,
#     "k_il8_b": 0.44794292479208003,
#     # "k_il8_ub": 0.5949379311724777,
#     # "k_il8r_p": 280.4238364787151,
#     # "k_il8r_d": 0.7158984163360824,
#     # "kd_il8_irak_p": 48486.502074111755,
#     # "k_il8_irak_p": 46550.06431383127
#     }
print(params)
# params['n_h3s10_il8_p'] = 4
print('t2m: {} '.format(t2m))
#params['IL8'] = 1
#params['IL8R'] = 100
flags = [
    # 'M1',
    # 'LPS',
    'IL8',
    # 'combined'
]

# params['kd_lps_irak_p'] = 50
# params['n_lps_irak_p'] = 1
# params['k_lps_d'] = .1
# params['k_lps_irak_p'] = 100000


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

    
if 'IL8' in flags: 
    model_t = 'IL8'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('IL8 is plotting')
    fig2 = funcs.P2_eq(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        
    fig2 = funcs.P21_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        
    # fig2 = funcs.P22_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig2 = funcs.P23_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
#     _dir = os.path.join(dirs.dir_outputs,'plots','P2_2.png')
    # plt.savefig(_dir)

if 'combined' in flags : 
    model_t = 'combined'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('Combined is plotting')
    fig = funcs.P31_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P32_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P33_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P34_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P35_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)

    # _dir = os.path.join(dirs.dir_outputs,'plots','P3.png')