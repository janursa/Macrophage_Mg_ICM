## functions
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import numpy as np
import tellurium as te
import json
import copy


def assign_surrogate_names(model,selections,prefix = 'PP_'):
    model_script  = model.getAntimony()
    index = model_script.find('end')
    model_script = model_script.replace(prefix,'')
    return te.loada(model_script)

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


