## functions
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
from tools import dirs
from models.models import Macrophage



def normalize(study_tag,target,sims,exps):
    ff = lambda ctr,vector: [i/ctr for i in vector]
    n_sims = ff(sims[0],sims)
    return n_sims,exps

def assign_surrogate_names(model,selections,prefix = 'PP_'):
    model_script  = model.getAntimony()
    index = model_script.find('end')
    model_script = model_script.replace(prefix,'')
    return te.loada(model_script)



class InvalidParams(Exception):
    def __init__(self,message = ''):
        pass

def calibrate(**args):
    cost_function = args['cost_function']
    free_params = args['free_params']
    bounds = list(free_params.values())
    maxiter = args['maxiter']
    workers = args['workers']
    callback = args['callback']

    results = differential_evolution(func = cost_function,bounds=bounds,maxiter=maxiter,workers=workers,
        disp=True,callback=callback)

    inferred_params = {}
    for key,value in zip(free_params.keys(),results.x):
        inferred_params[key] = value
    return inferred_params,results.fun



###----convert from concentration to copy number and vice versa ----#######
mws = { #molecular weights/ Da
    'IL1b' : 31000,
    'IFNG': 16879,
    'TNFa': 26000,
    'IL4': 20000,
    'V165a': 45000,
    'V165b': 45000,
    'IL10': 18000,
    'IL8': 8800 
} 
c_2_ac = {} # concentration to absolute copy
for key,value in mws.items():
    c_2_ac[key] = (6.022*10**23)/(value*10**9)/10**6
ac_2_c = {} # absolute copy number to concentration
for key,value in c_2_ac.items():
    ac_2_c[key] = 10**6*(value*10**9)/(6.022*10**23) #TODO
c_2_ac['O2'] = 120400000/21 # % to absolute copy for oxygen
ac_2_c['02'] = 21/120400000 # absolute copy to %
                            # 
                            

def indexing(real_time,time_vector):
    index = next(x[0] for x in enumerate(time_vector) if x[1] >= real_time)
    return index
def edit_matlab_model(input_file,output_file): # edit matlab sbml
    with open(input_file,'r') as file:
        content = file.read()
    lines = content.splitlines()
    id_name_map = {}
    for line in lines:
        id_i = line.find(' id')
        if id_i!=-1:
            name_i = line.find(' name')
            if name_i == -1:
                print(line) # this is an error
            def find_name(line):
                name_i = line.find('name')
                start_i = line[name_i:].find("\"") #find the first "
                a_start_line = line[name_i+start_i:]
                end_i = line[name_i+start_i+1:].find("\"")
                return line[name_i+start_i:name_i+start_i+end_i+2]
            name = find_name(line)
            def find_ID(line):
                start_i = line[id_i:].find("\"")
                end_i = line[start_i+id_i+1:].find("\"")
                return line[id_i+start_i:id_i+start_i+end_i+2]
            ID = find_ID(line)
            id_name_map.update({ID[1:-1]:name[1:-1]})
    for ID,name in id_name_map.items():
        new_name = name.replace("/","_")
        new_name = new_name.replace(" ","_")
        new_name = new_name.replace("-","_")
        content = content.replace(ID,new_name)
        content = content.replace(ID[1:-1],new_name[1:-1])

    with open(output_file,"w") as file:
        file.write(content)
def edit_M1_model():
    ##--- import sbml model and edits it to remove invalid ids----#####
    edit_matlab_model(input_file=dirs.dir_M1_matlab_model,output_file=dirs.dir_M1_model)
    model = te.loadSBMLModel(dirs.dir_M1_model)
    model_str = model.getAntimony()
    with open('models/M1_sbml_str.txt','w') as file:
        file.write(model_str)
def edit_zhao_model():
    import sys
    import tellurium as te
    from tools import dirs
    from tools.tools import edit_matlab_model
    edit_matlab_model(input_file=dirs.dir_Zhao_model_original,output_file=dirs.dir_Zhao_model)
    model = te.loadSBMLModel(dirs.dir_Zhao_model)
    model_str = model.getAntimony()
    with open('models/%s_str.txt'%'Zhao_sbml','w') as file:
        file.write(model_str)

def activation_zhao():
    ## applies TNFa to activate Zhao's model
    zz_model = te.loadSBMLModel(dirs.dir_Zhao_model)
    # apply 10ng/ml TNFa for 12h
    inputs = {'TNFa':tools.c_2_ac['TNFa']*10}
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
    inputs = {'LPS':1000}
    species_IDs = model.getFloatingSpeciesIds()
    Macrophage.run_sbml_model(model_sbml = model,params = {**inputs},duration = 24*60, selections = ['time']+species_IDs)

    # store the stimulated values of the species
    activation_stimuli = {}
    for ID in species_IDs:
        activation_stimuli[ID] = model[ID]
    with open(dirs.dir_activation_stimuli,'w') as ff:
        ff.write(json.dumps(activation_stimuli,indent=4))
    selections = ['NFKB_n','TNFa']
    def scenario1():
        results = Macrophage.run_sbml_model(model_sbml=model,params={**inputs},duration=48*60,selections = ['time']+selections)
        return results
    def scenario2():
        rr1 = Macrophage.run_sbml_model(model_sbml=model,params={**inputs},duration=24*60,selections = ['time']+selections)
        rr2 = Macrophage.run_sbml_model(model_sbml=model,params={**inputs},duration=24*60,selections = ['time']+selections,activation=True)
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

    
