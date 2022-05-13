#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 17:48:16 2022

@author: matin
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:42:32 2022

@author: matin
"""


main_dir = '/Users/matin/Downloads/testProjs/intracellular_M'
import sys
import os
import json
sys.path.insert(0,main_dir)
from data.observations import observations,t2m,select_obs,packages
from models.params import fixed_params
from tools import dirs, tools
from models.models import Macrophage
from plots import funcs 
import matplotlib
import matplotlib.pyplot as plt
import tellurium as te

params = {**fixed_params}
def reload_params(params): # apply inferred params
    target_package = 'P3'
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
    return params
params = reload_params(params)
model_t = 'combined'

model_sbml = Macrophage.create_sbml_model(model_t)

tags = ['IL8','F_p3s10_il8_p']
duration = 24*60
inputs = {}
ctr = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags)

fig = plt.figure(figsize=(len(tags)*2,2))
jj = 1
for tag in tags:
    ax = fig.add_subplot(1,len(tags),jj)
    tt = 0
    ax.plot(ctr['time'][tt:],ctr[tag][tt:],label = 'ctr')
    ax.set_title(tag)
    ax.legend()
    jj+=1
fig.tight_layout()