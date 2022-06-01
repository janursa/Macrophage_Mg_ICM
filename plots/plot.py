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
    target_package = 'ILs'
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
    return params
params = reload_params(params)
print(params)
# params['k_il8_p'] = 0.1
# params['o_il8_irak_p'] = 0

print('t2m: {} '.format(t2m))
flags = [
    # 'M1',
    # 'LPS',
    'IL6',
    'IL8',
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
if 'IL8' in flags: 
    model_t = 'ILs'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('IL8 is plotting')
    fig2 = funcs.P2_eq(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        
    fig2 = funcs.P2_ICs_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        
    # fig2 = funcs.P2_receptors_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    fig2 = funcs.P23_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
if 'IL6' in flags: 
    model_t = 'ILs'
    model_sbml = Macrophage.create_sbml_model(model_t)
    macrophage_obj = Macrophage(model_t = model_t)
    print('IL6 is plotting')
    fig = funcs.P2_IL6_IC_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        
    fig = funcs.P2_IL6_CYs_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)                        

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
    fig = funcs.P3_cytokines4_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)

    # _dir = os.path.join(dirs.dir_outputs,'plots','P3.png')