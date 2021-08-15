from scipy.optimize import differential_evolution

import numpy as np
import json
# import matplotlib.pyplot as plt
from models import *

#// load the samples from the original model
with open('samples.json') as json_file:
	samples = json.load(json_file)

#// load the observation
with open('observations.json') as file:
	observations = json.load(file)
#// Set-up hyperparameters
class SETTINGS:
	params = list(PARAMS.free_params.values())
	print(PARAMS.free_params)
	max_iters = 10

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

def categorize(calib_params): # divides the calib values between the model parameters and those for studies
	free_params_model_list = calib_params[0:PARAMS.free_params_model_n]
	inferred_params_model = {}
	for key,value in zip(PARAMS.free_params_model.keys(),free_params_model_list):
	    inferred_params_model[key] = value
	free_params_observations_list = calib_params[PARAMS.free_params_model_n:PARAMS.free_params_model_n+PARAMS.free_params_observations_n]
	free_params_observations = {}
	for key,value in zip(PARAMS.free_params_observations.keys(),free_params_observations_list):
		free_params_observations[key] = value
	free_params = {**inferred_params_model,**free_params_observations}
	return free_params,inferred_params_model,free_params_observations

def cost_function(calib_params):

	free_params,free_params_model,free_params_observations = categorize(calib_params)

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
			study_observations = observations[study]
			measurement_scheme = study_observations['measurement_scheme']
			simulation_period = study_observations['experiment_period']*60 # changing from hour to minute
			selections = ['TIME']+ list(measurement_scheme.keys())
			IDs_errors = []
			for ID in study_observations['IDs']:
				reset(Mg_M,free_params_model) #reset the model
				ID_inputs = study_observations[ID]['inputs']
				ID_observations = study_observations[ID]['expectations']
				## apply the boundary condition
				for key,value in ID_inputs.items():
					Mg_M[key] = value
				ID_results = Mg_M.simulate(0,simulation_period,simulation_period,selections=selections)
				tag_errors = []
				for key in ID_observations.keys():
					exp = ID_observations[key]['mean'] # the whole array
					sim = ID_results[key]
					
					sim_selected = [sim[i*60-1] for i in measurement_scheme[key]]

					sim_selected = [i*free_params_observations[study] for i in sim_selected]

					abs_diff = [np.abs(exp_item - sim_item) for exp_item,sim_item in zip(exp,sim_selected)]
					norm_abs_diff = np.mean(abs_diff)/((np.mean(exp)+np.mean(sim_selected))/2)

					tag_errors.append(norm_abs_diff)
				IDs_errors.append(np.mean(tag_errors))
			studies_errors.append(np.mean(IDs_errors))

		return np.mean(studies_errors)

	# error1 = vs_Zhao()
	error2 = vs_observations()

	# return error1
	return error2

def optimize():
	# Call instance of PSO
	results = differential_evolution(cost_function,bounds=SETTINGS.params,disp=True,maxiter=SETTINGS.max_iters)
	return results
results = optimize()
inferred_params = {}
for key,value in zip(PARAMS.free_params.keys(),results.x):
	inferred_params[key] = value
with open('inferred_params.json','w') as file:
	file.write(json.dumps(inferred_params,indent=4))
