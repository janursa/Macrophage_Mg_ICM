import pyswarms as ps
from pyswarms.utils.functions import single_obj as fx
import numpy as np
import json
# import matplotlib.pyplot as plt
from models import *
from observations import observations
from calibration import cost_function_study

# load the samples from the original model
with open('samples.json') as json_file:
    samples = json.load(json_file)

# Set-up hyperparameters
class SETTINGS:
	options = {'c1': 0.5, 'c2': 2, 'w':2}
	dimensions = len(PARAMS.free_params.keys())
	bounds = ([item[0] for item in PARAMS.free_params.values()],[item[1] for item in PARAMS.free_params.values()])
	n_particles = 100
	iters = 100


def categorize(calib_params): # divides the calib values between the model parameters and those for studies
	free_params_model_list = calib_params[0:PARAMS.free_params_model_n]
	inferred_params_model = {}
	for key,value in zip(PARAMS.free_params_model.keys(),free_params_model_list):
	    inferred_params_model[key] = value
	hyperparams_list = calib_params[PARAMS.free_params_model_n:PARAMS.free_params_model_n+PARAMS.hyperparams_n]
	hyperparams = {}
	for key,value in zip(PARAMS.hyperparams.keys(),hyperparams_list):
		hyperparams[key] = value
	free_params = {**inferred_params_model,**hyperparams}
	return free_params,inferred_params_model,hyperparams


def cost_function(calib_params_partcles):
	costs = [] # cost values for each particle
	for j in range(SETTINGS.n_particles):
		calib_params = calib_params_partcles[j]
		free_params,free_params_model,hyperparams = categorize(calib_params)

		def vs_Zhao(): # to compare the model's results with the Zhao
			reset(Mg_M,params=free_params_model)
			selections = ['TIME']+PARAMS.targets
			sim_results = Mg_M.simulate(0,PARAMS.duration,PARAMS.duration,selections=selections)
			# calculate the error for each tag by comparing the results to the original model
			tag_errors = []
			for tag in PARAMS.targets:
				abs_diff =np.mean([abs(i-j) for i,j in zip(sim_results[tag],samples[tag])])
				means = [np.mean(samples[tag]),np.mean(sim_results[tag])]
				mean = np.mean(means)
				tag_error = abs_diff/mean
				# tag_error = np.mean(abs_diff)
				# print('tag {} abs_diff {} mean {} tag_error {}'.format(tag,abs_diff,mean,tag_error))
			tag_errors.append(tag_error)
			return np.mean(tag_errors)
		def vs_observations(): # to compare the model's predictions with the observations
			studies_errors = []
			for study in observations['studies']:
				study_error = cost_function_study(study=study,free_params_model=free_params_model,hyperparams=hyperparams)
				studies_errors.append(study_error)

			return np.mean(studies_errors)

		# error1 = vs_Zhao()
		error2 = vs_observations()

		costs.append(error2)

	return costs

def optimize():
	# Call instance of PSO
	optimizer = ps.single.GlobalBestPSO(n_particles=SETTINGS.n_particles, 
		dimensions=SETTINGS.dimensions, options=SETTINGS.options, bounds=SETTINGS.bounds)
	# Perform optimization
	cost, pos = optimizer.optimize(cost_function,iters=SETTINGS.iters)

	return cost,pos
cost,pos = optimize()
inferred_params = {}
for key,value in zip(PARAMS.free_params.keys(),pos):
	inferred_params[key] = value
with open('inferred_params.json','w') as file:
	file.write(json.dumps(inferred_params, indent = 4))
