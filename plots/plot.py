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
    target_package = 'LPS'
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
    return params
params = reload_params(params)
print(params)
# params['n_h3s10_il8_p'] = 4
print('t2m: {} '.format(t2m))
#params['IL8'] = 1
#params['IL8R'] = 100
flags = [
    # 'P1',
    # 'P2',
    'LPS',
    # 'P3'
]

# params['kd_lps_irak_p'] = 1
# params['n_lps_irak_p'] = 5

if 'P1' in flags : 
    model_t = 'combined'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('P1_3 is plotting')
    funcs.P1_eq(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    funcs.P11_plot(model_sbml=model_sbml,params=params,observations=observations)
    funcs.P12_plot (model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
if 'LPS' in flags : 
    model_t = 'LPS'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('LPS is plotting')
    funcs.LPS_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)

    
if 'P2' in flags: 
    model_t = 'combined'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('P2 is plotting')
    fig2 = funcs.P2_eq(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        
    fig2 = funcs.P21_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        
    fig2 = funcs.P22_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig2 = funcs.P23_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
#     _dir = os.path.join(dirs.dir_outputs,'plots','P2_2.png')
    # plt.savefig(_dir)

if 'P3' in flags : 
    model_t = 'combined'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('P3 is plotting')
    fig = funcs.P31_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P32_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P33_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P34_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig = funcs.P35_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)

    # _dir = os.path.join(dirs.dir_outputs,'plots','P3.png')