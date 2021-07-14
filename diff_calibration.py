from scipy.optimize import differential_evolution

import numpy as np
import json
# import matplotlib.pyplot as plt
from models import *

# load the samples from the original model
with open('samples.json') as json_file:
	samples = json.load(json_file)
# Set-up hyperparameters
class SETTINGS:
	params = list(PARAMS.free_params.values())
	print(params)
	max_iters = 50

def run_replicas(model,calib_params):
    """
    run the model for a given number of replicas
    
    """
    stack_results = [] # store for each iteration
    for i in range(PARAMS.replica_n):
        model.reset()
        initial_conditions(model,calib_params)
        results_dict=run(model,targets=PARAMS.targets,duration=PARAMS.duration)
        stack_results.append(results_dict)
    return stack_results


def cost_function(calib_params):
	# print('\n',calib_params)
		
	# run the model
	stack_results = run_replicas(model = Mg_model,
								 calib_params=calib_params)
	mean_results = average(stack_results)

	# calculate the error for each tag by comparing the results to the original model
	tag_errors = []
	for tag in PARAMS.targets:
		abs_diff =np.mean(abs(np.array(mean_results[tag])-np.array(samples[tag])))
		means = [np.mean(samples[tag]),np.mean(mean_results[tag])]
		mean = np.mean(means)
		tag_error = abs_diff/mean
		# tag_error = np.mean(abs_diff)
		# print('tag {} abs_diff {} mean {} tag_error {}'.format(tag,abs_diff,mean,tag_error))
		tag_errors.append(tag_error)
	error = np.mean(tag_errors)
	# print('error',error)
	return error



def optimize():
	# Call instance of PSO
	results = differential_evolution(cost_function,bounds=SETTINGS.params,disp=True,maxiter=SETTINGS.max_iters)
	return results
results = optimize()
print(results)
np.savetxt('inferred_values.csv',results.x,delimiter=",")