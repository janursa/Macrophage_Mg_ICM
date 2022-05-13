#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 15:09:44 2022

@author: matin
"""


import tellurium as te
import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
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
    target_package = 'LPS'
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
    return params
# params = reload_params(params)
model_t = 'LPS'

model_sbml = Macrophage.create_sbml_model(model_t)
print(model_sbml['TNFa'])
tags = ['IRAK4','TRAF6','aTRAF6','IKB','TNFa']
duration = 8*60
inputs = {}
ctr = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags)
inputs = {'LPS':50}
rr_8 = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags)

fig = plt.figure(figsize=(len(tags)*2,2))
jj = 1
for tag in tags:
    ax = fig.add_subplot(1,len(tags),jj)
    tt = 0
    ax.plot(ctr['time'][tt:],ctr[tag][tt:],label = 'ctr')
    ax.plot(rr_8['time'][tt:],rr_8[tag][tt:],label = 'LPS 50')
    # ax.plot(rr_dpt8['time'][tt:],rr_dpt8[tag][tt:],label = 'Mg 0.08')
    ax.set_title(tag)
    ax.legend()
    jj+=1
fig.tight_layout()