import pyswarms as ps
from pyswarms.utils.functions import single_obj as fx
import numpy as np
import json
# import matplotlib.pyplot as plt
from models import *

# load the samples from the original model
with open('samples.json') as json_file:
    samples = json.load(json_file)

# Set-up hyperparameters
class SETTINGS:
	options = {'c1': 0.5, 'c2': 2, 'w':2}
	dimensions = len(PARAMS.free_params)
	bounds = ([0 for i in range(dimensions)],[1000 for i in range(dimensions)])
	n_particles = 100
	iters = 100

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
	costs = [] # cost values for each particle
	for j in range(SETTINGS.n_particles):
		# run the model
		stack_results = run_replicas(model = Mg_model,
								 calib_params=calib_params[j])
		mean_results = average(stack_results)
		# calculate the error for each tag by comparing the results to the original model
		tag_errors = []
		for tag in PARAMS.targets:
			abs_diff =abs(np.array(mean_results[tag])-np.array(samples[tag]))
			means = [np.mean(samples[tag]),np.mean(mean_results[tag])]
			mean = np.mean(means)
			tag_error = np.mean(abs_diff/mean)
			tag_errors.append(tag_error)
		costs.append(np.mean(tag_errors))

	return costs

def optimize():
	# Call instance of PSO
	optimizer = ps.single.GlobalBestPSO(n_particles=SETTINGS.n_particles, 
		dimensions=SETTINGS.dimensions, options=SETTINGS.options, bounds=SETTINGS.bounds)
	# optimizer = ps.single.GlobalBestPSO(n_particles=SETTINGS.n_particles, 
	# 	dimensions=SETTINGS.dimensions, options=SETTINGS.options)

	# Perform optimization
	cost, pos = optimizer.optimize(cost_function,iters=SETTINGS.iters)

	return cost,pos
cost,pos = optimize()
