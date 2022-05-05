main_dir = '/Users/matin/Downloads/testProjs/intracellular_M'
import sys
sys.path.insert(0,main_dir)
import numpy as np
import json
import matplotlib.pyplot as plt
import tellurium as te
import os
plt.rcParams["font.family"] = "serif"
plt.style.use('seaborn-deep')
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
from data.observations import observations,t2m,select_obs
from models.params import fixed_params
from tools import dirs, tools
from models.models import Macrophage
from plots import funcs 
params = {**fixed_params}
if True: # apply inferred params
    target_package = 'P2'
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
print(params)
print('t2m: {} '.format(t2m))
flags = {
    'P1': False,
    'P4__': False,
    'P2': True
}
for key,value in flags.items():
    if key == 'P1' and value : 
        model_t = 'M1'
        model_sbml = Macrophage.models[model_t]
        macrophage_obj = Macrophage(model_t = model_t)
        params = {}
        print('P1_3 is plotting')
        funcs.P1_eq_plot(model_sbml=model_sbml,params={},observations=observations)
        plt.savefig(os.path.join(dirs.dir_outputs,'P11.png'))
        funcs.P1_qualitative_plot (model_sbml=model_sbml,params=params,observations=observations)
        plt.savefig(os.path.join(dirs.dir_outputs,'P12.png'))
        fig = funcs.P1_plot (model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
        plt.savefig(os.path.join(dirs.dir_outputs,'P13.png'))
    elif key == 'P4' and value : 
        dir_model = dirs.dir_model
        model_sbml = te.loadSBMLModel(dir_model)
        macrophage_obj = Macrophage(dir_model = dir_model)
        print('P4 is plotting')
        funcs.P4_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
    elif key == 'P2' and value : 
        model_t = 'IL8'
        model_sbml = Macrophage.models[model_t]
        macrophage_obj = Macrophage(model_t = model_t)
        print('P2 is plotting')
        funcs.P2_plot(model_sbml=model_sbml,model_macrophage=macrophage_obj,params=params,observations=observations)
        plt.savefig(os.path.join(dirs.dir_outputs,'P2.png'))
