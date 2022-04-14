## functions
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import numpy as np
import tellurium as te
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
###------some plots--------#####
def plot(y,x=None,target=''):
    fig = plt.figure()
    fig.canvas.draw()
    ax = fig.add_subplot(1, 1, 1)
    if x == None:
        ax.plot(y)
    else:
        ax.plot(x,y)
    ax.set_title(target)
    plt.show()
def run_model(model,params,target_keys,duration=500):
    model.reset()
#     model = te.loadSBMLModel(dir_model)
    for key,value in params.items():
        model[key] = value

    results = model.simulate(0,duration,duration,selections = ['TIME']+target_keys)
    return results
def indexing(real_time,time_vector):
    index = next(x[0] for x in enumerate(time_vector) if x[1] >= real_time)
    return index
###------calibration algorithm---------###

class Calibrate_2:
    """
    Target should be given in this format:
    target = [{'target_key':{'time':[t0,t1,t2,...],
                 'values':[v1,v2,v3,...]}},..]
    """
    def __init__(self,model,target,free_params,max_iteration):
        self.model = model
        self.target = target
        self.free_params = free_params
        self.max_iteration = max_iteration
    def cost_function(self,params_values):
        ## create the parameter set for the given parameter values
        params = {}
        for key,value in zip(self.free_params.keys(),params_values):
            params[key] = value
        errors = []
        for target_item in self.target:
            results_item = target_item['results']
            inputs_item = target_item['inputs']
            for key,value in inputs_item.items():
                params[key] = value
            results = run_model(model=self.model, params = params, target_keys= list(results_item.keys()))
            tag_errors=[]
            for tag in list(results_item.keys()):
                n_abs_values= []
                for i in range(len(results_item[tag]['time'])):
                    index = indexing(results_item[tag]['time'][i],results['time'])
                    result_i = results[tag][index] 
                    target_i = results_item[tag]['values'][i]
                    abs_value = abs(result_i - target_i)
        #             mean_value = np.mean([results[tag][i], samples[tag][i]])
                    mean_value = target_i
                    n_abs_values.append(abs_value/mean_value)
                tag_error = np.mean(n_abs_values)
                tag_errors.append(tag_error)
            errors.append(np.mean(tag_errors))
        return np.mean(errors)
    def optimize(self):
        calib_obj = Calibration(free_params=self.free_params,cost_function=self.cost_function,max_iters=self.max_iteration)
        params = calib_obj.optimize()
        return params
