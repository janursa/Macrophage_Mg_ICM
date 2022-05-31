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
    target_package = 'IL6'
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
    return params
params = reload_params(params)

model_t = 'IL6'

# obj = Macrophage(model_t = model_t)
# print(obj.run(params=params,studies=select_obs(packages[model_t])))

model_sbml = Macrophage.create_sbml_model(model_t)
# params['o_stat3_a'] = 0
# params['k_stat3_a'] = 100

# tags = ['F_il8_irak','IRAK4','NFKB_n']
tags = ['IL1b','pPI3K','pAKT1','ARG1','TNFa','IL10']
# tags = ['F_h3s10_ikb','F_p3s10_il8_p']
activation = False
duration = 24*60
inputs = {}
ctr = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=activation)
# print(ctr['F_stat3_a'])
inputs = {'IL6':50*1000}
stim = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=activation)
# print(ctr['IL6_R_0'])
fig = plt.figure(figsize=(len(tags)*2,2))
jj = 1
for tag in tags:
    ax = fig.add_subplot(1,len(tags),jj)
    tt = 0
    ax.plot(ctr['time'][tt:],ctr[tag][tt:],label = 'ctr')
    ax.plot(stim['time'][tt:],stim[tag][tt:],linestyle='--',label = 'stim')
    ax.set_title(tag)
    ax.legend()
    # if tag == 'F_stat3_a':
    #     ax.set_ylim([0,2])
    jj+=1
fig.tight_layout()