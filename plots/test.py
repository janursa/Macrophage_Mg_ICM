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
target_package = 'ILs_1'
params = {**fixed_params}
def reload_params(params): # apply inferred params
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
    return params
# params = reload_params(params)
params = {
    "IL6": 5.4370102257715365,
    "IL6_R_JACK": 2.463356323589118,
    "pIL6_R_JACK": 1.849107258284274,
    "k_il6r_b": 0.44075140849316885,
    "k_il6r_ub": 0.6493205249598956,
    "k_il6r_a": 0.4270556206063457,
    "k_il6r_da": 0.22582593386762706,
    "IL6_m": 4.576711588346838,
    "k_il6m_il6": 0.3932979489856194,
    "k_il6_p": 3.4710282869815017,
    "k_il6_stat3_a": 1951.121971969158,
    "kd_il6_stat3_a": 7318.074225489829,
    "o_il6_stat3_a": 0.42793367800462123,
    "k_il6_pi3k_a": 1266.540409509294,
    "kd_il6_pi3k_a": 4598.969155571763,
    "o_il6_pi3k_a": 0.4849854774932578,
    "k_nfkb_il6_p": 33.93520028798429,
    "kd_nfkb_il6_p": 72442.04587322613,
    "o_nfkb_il6_p": 0.07627002255439408,
    "k_ap1_il6_p": 240.69769424843963,
    "kd_ap1_il6_p": 84333.84255485586,
    "o_ap1_il6_p": 0.008153212009223143,
    "k_creb_il6_p": 11.408246267173212,
    "kd_creb_il6_p": 72487.83964114796,
    "o_creb_il6_p": 0.5165947002014888,
    "IL8": 9.007692694221326,
    "IL8_m": 0.9209962714924638,
    "k_il8_p": 0.01726296193977106,
    "k_il8m_il8": 0.6214970639297708,
    "IL8_R": 4.19077277865208,
    "pIL8_R": 3.8213439139543124,
    "k_il8r_b": 0.02161783977154752,
    "k_il8r_ub": 0.496702850242793,
    "k_il8r_a": 0.35558316753820685,
    "k_il8r_da": 0.768923681728841,
    "kd_il8_irak_p": 282522.7279195321,
    "k_il8_irak_p": 911956.2708446574,
    "o_il8_irak_p": 0.31236589163039247,
    "k_nfkb_il8_p": 1.0381760913336393,
    "kd_nfkb_il8_p": 979.5911484020489,
    "o_nfkb_il8_p": 0.9159888394230121,
    "k_ap1_il8_p": 1.2691172150609304,
    "kd_ap1_il8_p": 456.9686461242685,
    "o_ap1_il8_p": 0.9182426473795741,
    "k_rho_a": 3.411447676886212,
    "kd_rho_a": 6871.233986011801,
    "k_rho_da": 0.4935009817946129,
    "k_rho_nfkb_a": 18512.259719458943,
    "kd_rho_nfkb_a": 71035.36196641353,
    "o_rho_nfkb_a": 0.00044230546029566664,
    "k_rho_pi3k_a": 8818.378265969557,
    "kd_rho_pi3k_a": 99056.02852261194,
    "o_rho_pi3k_a": 0.22113616205429887,
    "k_rho_stat3_a": 2279.361710806342,
    "kd_rho_stat3_a": 47182.69054735674,
    "o_rho_stat3_a": 0.06661670957536459,
    "_pSTAT3": 1.648340019534821,
    "k_pstat3_d": 0.5780014457545202,
    "k_pstat3_b": 0.19856751836254527,
    "k_pstat3_ub": 0.7290493085817544,
    "_pPI3K": 44.322216620555416,
    "k_ppi3k_d": 0.6530846648225805,
    "K_ppi3k_a": 0.7698986492411686,
    "K_ppi3k_da": 0.1663185949643658
}

# obj = Macrophage(model_t = model_t)
# error_mean, error = obj.run(params=params,studies=select_obs(packages[target_package]))
# print(error_mean)
# F_nfkb_il6_p := 1+32.61753482998196*(PP.NFKB_n-NFKB_n_0)/(PP.NFKB_n-NFKB_n_0+209.0757011467067) - 0.494621418549415;

model_sbml = Macrophage.create_sbml_model(model_t)
# tags = ['pPI3K','IRAK4','pIRAK4','m146b','SOCS1','A20','IL1b_R','aTRAF6']
# tags = ['F_mg_ikb_d','IKB','IKB_0','NFKB_n']
# tags = ['IKB_0','IKB','nNFKB_n','pSTAT3D_n','AP1_n','pAKT_t','IL1b','IL10','TNFa','pJNK']
# tags = ['NFKB_n','TNFa','IL6','IL8','IL1b','IL10']
# tags = ['NFKB_n','TNFa','IL1b','IL10']
# tags = ['F_nfkb_il6_p','F_ap1_il6_p','F_creb_il6_p']
tags= ['pSTAT3','pPI3K','AKT1','pAKT1']
activation = False
duration = 24*60
inputs = {}
ctr = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=activation)
# print(ctr['NFKB_n'][-1])
inputs = {'IL8':100000}
stim = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=activation)
# print(stim['NFKB_n'][-1])
row = int(len(tags))
col = min([len(tags),3])
fig = plt.figure(figsize=(col*4,row*3))
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