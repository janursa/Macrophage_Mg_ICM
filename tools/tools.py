## functions
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import numpy as np
import tellurium as te
import json

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

class Calibration_class:
    def __init__(self,model,fixed_params,free_params,studies):
        self.free_params = free_params
        self.fixed_params = fixed_params
        self.model = model
        self.studies = studies
    def cost_function(self,calib_params_values):
        calib_params = {}
        for key,value in zip(self.free_params.keys(),calib_params_values):
            calib_params[key] = value
        params = {**calib_params,**self.fixed_params}
        error = self.model.run(params=params,studies=self.studies)

        return error

    def optimize(self,n_proc,max_iters,disp,strategy,tol,callback):
        results = differential_evolution(func = self.cost_function,bounds=list(self.free_params.values()),maxiter=max_iters,workers=n_proc,strategy=strategy,
            disp=disp,tol=tol,callback=callback)

        inferred_params = {}
        for key,value in zip(self.free_params.keys(),results.x):
            inferred_params[key] = value
        return inferred_params,results.fun

def calibrate(model,fixed_params,free_params,studies,max_iters = 20, n_proc=1,disp=True,strategy='best1bin',tol=1e-4,callback=None):
    calib_obj = Calibration_class(model,fixed_params,free_params,studies)

    inferred_params,error = calib_obj.optimize(n_proc=n_proc,max_iters=max_iters,disp=disp,strategy=strategy,tol=tol,callback=callback)
    return inferred_params

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
                            
def run_model(model,params,target_keys,duration,study=''):
    model.reset()
#     model = te.loadSBMLModel(dir_model)
    for key,value in params.items():
        model[key] = value
    if study == 'S12_IkBa_mg':
        model.integrator.variable_step_size = True
        model.integrator.absolute_tolerance = 1e-3 
        model.integrator.relative_tolerance = 1e-3 
        try:
            results = model.simulate(0,duration,selections = ['TIME']+target_keys)
        except RuntimeError:
            raise RuntimeError('study : {} didnt converge'.formaat(study))
    elif study == 'Q21_IkBa':
        jj = duration
        while True:
            steps = jj
            try:
                results = model.simulate(start = 0,end = duration,steps =steps ,selections = ['TIME']+target_keys)
                break
            except RuntimeError:
                jj-=1
            if jj < 20:
                print('Invalid parameter set')
                raise InvalidParams('run model didnt converge')
    else:
        results = model.simulate(start = 0,end = duration,steps = int(duration/2) ,selections = ['TIME']+target_keys)

    return results
def indexing(real_time,time_vector):
    index = next(x[0] for x in enumerate(time_vector) if x[1] >= real_time)
    return index


