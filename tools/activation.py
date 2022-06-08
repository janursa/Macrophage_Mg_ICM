import tellurium as te
import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import numpy as np
import tellurium as te
import json
import copy
from tools import dirs,common
from models.models import Macrophage

def activation_zhao():
    ## applies TNFa to activate Zhao's model
    zz_model = te.loadSBMLModel(dirs.dir_Zhao_model)
    # apply 10ng/ml TNFa for 12h
    inputs = {'TNFa':common.c_2_ac['TNFa']*10}
    species_IDs = zz_model.getFloatingSpeciesIds()
    Macrophage.run_sbml_model(model_sbml = zz_model,params = {**inputs},duration = 12*60, selections = ['time']+species_IDs)

    # store the stimulated values of the species
    activation_stimuli = {}
    for ID in species_IDs:
        activation_stimuli[ID] = zz_model[ID]
    with open(dirs.dir_activation_stimuli,'w') as ff:
        ff.write(json.dumps(activation_stimuli,indent=4))
        
    selection = 'NFKB_n'
    def scenario1():
        results = Macrophage.run_sbml_model(model_sbml=zz_model,params={**inputs},duration=36*60,selections = ['time',selection])
        return results
    def scenario2():
        rr1 = Macrophage.run_sbml_model(model_sbml=zz_model,params={**inputs},duration=12*60,selections = ['time',selection])
        rr2 = Macrophage.run_sbml_model(model_sbml=zz_model,params={**inputs},duration=24*60,selections = ['time',selection],activation=True)
        return rr1,rr2
    rr1 = scenario1()
    rr21,rr22 = scenario2()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    ax.plot(rr1['time'],rr1[selection],label='rr1')
    ax.plot(rr21['time'],rr21[selection],label='rr21')
    ax.plot(rr22['time']+max(rr21['time']),rr22[selection],label='rr22')
    ax.legend()
def activation_LPS():
    ## applies LPS to activate LPS model
    model = te.loadSBMLModel(dirs.dir_LPS_model)
    # apply LPS (1.0 μg/ml) for 24h
    inputs = {'LPS': 10*1000,
             'IFNG':50*common.c_2_ac['IFNG']}
    duration = 24*60
    species_IDs = model.getFloatingSpeciesIds()
    Macrophage.run_sbml_model(model_sbml = model,params = {**inputs},duration = duration, selections = ['time']+species_IDs)

    # store the stimulated values of the species
    activation_stimuli = {}
    for ID in species_IDs:
        activation_stimuli[ID] = model[ID]
    with open(dirs.dir_activation_stimuli,'w') as ff:
        ff.write(json.dumps(activation_stimuli,indent=4))
    selections = ['NFKB_n','TNFa','IL10']
    def scenario1():
        results = Macrophage.run_sbml_model(model_sbml=model,params={**inputs},duration=2*duration,selections = ['time']+selections)
        return results
    def scenario2():
        rr1 = Macrophage.run_sbml_model(model_sbml=model,params={**inputs},duration=duration,selections = ['time']+selections)
        rr2 = Macrophage.run_sbml_model(model_sbml=model,params={**inputs},duration=duration,selections = ['time']+selections,activation=True)
        return rr1,rr2
    rr1 = scenario1()
    rr21,rr22 = scenario2()
    fig = plt.figure()
    nn = len(selections)
    jj=1
    for selection in selections:
        ax = fig.add_subplot(1,nn,jj)
        ax.plot(rr1['time'],rr1[selection],label='rr1')
        ax.plot(rr21['time'],rr21[selection],label='rr21')
        ax.plot(rr22['time']+max(rr21['time']),rr22[selection],label='rr22')
        ax.legend()
        ax.set_title(selection)
        jj+=1
    fig.tight_layout()
