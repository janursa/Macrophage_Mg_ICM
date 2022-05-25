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
    target_package = 'IL8'
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
    return params
params = reload_params(params)
params['k_il8_ifngr_p'] = 10000
# params['kd_il8_ifngr_p'] = params['kd_il8_irak_p']
# params['kd_il8_irak_p'] = 10000
# params =  {'IL8_m': 7.034761864791667, 'k_il8_p': 9.465279786024439, 'kd_nfkb_il8_p': 5309.896779512799, 'IL8R': 100, 'IL8_R': 100, 'k_il8_b': 0.022726688427806785}

model_t = 'IL8'

# obj = Macrophage(model_t = 'IL8')
# obj.run(params=params,studies=select_obs(packages['IL8']))

model_sbml = Macrophage.create_sbml_model(model_t)


# tags = ['F_il8_irak','IRAK4','NFKB_n']
tags = ['pTAK1','IL1b','IL10','TNFa','IL4']
activation = False
duration = 72*60
inputs = {}
ctr = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=activation)

inputs = {'IL8':10000}
stim = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=activation)

fig = plt.figure(figsize=(len(tags)*2,2))
jj = 1
for tag in tags:
    ax = fig.add_subplot(1,len(tags),jj)
    tt = 0
    ax.plot(ctr['time'][tt:],ctr[tag][tt:],label = 'ctr')
    ax.plot(stim['time'][tt:],stim[tag][tt:],linestyle='--',label = 'stim')
    ax.set_title(tag)
    ax.legend()
    jj+=1
fig.tight_layout()