import tellurium as te
import numpy as np

## run function
def run(model,targets, duration):
    model.integrator.absolute_tolerance = 1e-9
    model.integrator.relatice_tolerance = 1e-9
    results = model.simulate(start = 0, end = duration,steps = duration,
                             selections = targets)
    return results


def average(stack_results):
    """
    Averages the stack results of all iterations
    """
    mean_results = {} # mean results for each target
    stack = [] 
    for tag in PARAMS.targets:
        for i in range(PARAMS.replica_n):
            stack.append(stack_results[i][tag])

        stack = np.array(stack)
        mean_results.update({tag:list(np.mean(stack,axis=0))})
        stack=[] 
    return mean_results
def initial_conditions(model,calib_params):
    keys = list(PARAMS.free_params.keys())
    for i in range(len(keys)):
        free_param_name = keys[i]
        model[free_param_name] = calib_params[i]


class PARAMS:
    targets = ['NFKB', 'pIKK', 'TAK1']
    duration = 1500
    free_params = {'Mg_e':[0,1000], 'NEMO_IKK':[0,1000],'k301':[0,1000],'k302':[0,1000],
                'k303':[0,1000],'k304':[0,1000],'k308':[0,1000]}

    replica_n = 1

Zhao_2021 = te.loadSBMLModel("Zhao_2021.xml") 
Mg_model = te.loadSBMLModel("Mg_M.xml") 




