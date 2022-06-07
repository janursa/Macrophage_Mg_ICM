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
from tools import dirs, common
from models.models import Macrophage
from plots import funcs 
import matplotlib
import matplotlib.pyplot as plt
import tellurium as te

model_t = 'ILs'
target_package = 'ILs_p3'
params = {**fixed_params}
def reload_params(params): # apply inferred params
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
    return params
params = reload_params(params)
# params['k_rho_stat3_a'] = 500
obj = Macrophage(model_t = model_t)
error_mean, error = obj.run(params=params,studies=select_obs(packages[target_package]))
# print(error)

model_sbml = Macrophage.create_sbml_model(model_t)


# tags = ['IL8','IL8_m','IL8_m_0','F_ap1_il8_p']
tags = ['NFKB_n','pSTAT3D_n','AP1_n','pAKT_t','IL1b','IL10','TNFa','IL6']
# tags = ['IL8','IL8_m','IL8_R','pIL8_R','F_il8_irak','IRAK4','pIL8_R_0']
# tags = ['F_h3s10_ikb','F_p3s10_il8_p']
activation = False
duration = 24*60
inputs = {}
ctr = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=activation)
# print(ctr['F_il8_irak'])
inputs = {'IL8':10*1000}
stim = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=activation)
# print(ctr['IL6_R_0'])
row = int(len(tags)/2)
col = 3
fig = plt.figure(figsize=(row*3,3*col))
jj = 1
for tag in tags:
    ax = fig.add_subplot(row,col,jj)
    tt = 0
    ax.plot(ctr['time'][tt:],ctr[tag][tt:],linewidth=3,label = 'ctr')
    ax.plot(stim['time'][tt:],stim[tag][tt:],linewidth=3,linestyle='--',label = 'stim')
    ax.set_title(tag,fontdict ={'size':20})
    ax.legend()
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(20)

    jj+=1
fig.tight_layout()