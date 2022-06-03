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
    target_package = 'ILs'
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
    return params
params = reload_params(params)

params ={'IL8': 0.05417800783813398, 'IL8_m': 9.010354898782708, 'k_il8_p': 43.80976911279544, 'k_il8m_il8': 0.7349143557052686, 'IL8_R': 3.9819106417736605, 'pIL8_R': 3.120912652647376, 'k_il8r_b': 0.4533389467185141, 'k_il8r_ub': 0.4610177027654007, 'k_il8r_a': 0.15800112940461108, 'k_il8r_da': 0.8935678273623348, 'kd_il8_irak_p': 55506.50900111589, 'k_il8_irak_p': 44497.45504584355, 'o_il8_irak_p': 0.7431874694708296, 'k_rho_a': 445.5189104749942, 'kd_rho_a': 18335.221110170503, 'o_rho_a': 0.6875494540296361, 'k_rho_pi3k_a': 62947.927598727205, 'kd_rho_pi3k_a': 7402.372139830644, 'o_rho_pi3k_a': 0.9613267006001454, 'k_rho_stat3_a': 59068.43972027341, 'kd_rho_stat3_a': 39817.44239652883, 'o_rho_stat3_a': 0.08343997652724, 'k_rho_jnk_a': 25958.808506109814, 'kd_rho_jnk_a': 8383.891699878732, 'o_rho_jnk_a': 0.501202636096848, 'k_nfkb_il6_p': 756.3263286777271, 'kd_nfkb_il6_p': 570.8704540351173}
model_t = 'ILs'

# obj = Macrophage(model_t = model_t)
# print(obj.run(params=params,studies=select_obs(packages[model_t])))

model_sbml = Macrophage.create_sbml_model(model_t)


# tags = ['IL8','IL8_m','IL8_m_0','F_ap1_il8_p']
tags = ['IL6','IL6_m','IL6_m_0','IL6_R_JACK','F_pi3k_a']
# tags = ['IL8','IL8_m','IL8_R','pIL8_R','F_il8_irak','IRAK4','pIL8_R_0']
# tags = ['F_h3s10_ikb','F_p3s10_il8_p']
activation = False
duration = 24*60
inputs = {}
ctr = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=activation)
# print(ctr['F_il8_irak'])
# inputs = {'IL8':100*1000}
# stim = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=activation)
# print(ctr['IL6_R_0'])
fig = plt.figure(figsize=(len(tags)*2,2))
jj = 1
for tag in tags:
    ax = fig.add_subplot(1,len(tags),jj)
    tt = 0
    ax.plot(ctr['time'][tt:],ctr[tag][tt:],label = 'ctr')
    # ax.plot(stim['time'][tt:],stim[tag][tt:],linestyle='--',label = 'stim')
    ax.set_title(tag)
    ax.legend()
    # if tag == 'F_stat3_a':
    #     ax.set_ylim([0,2])
    jj+=1
fig.tight_layout()