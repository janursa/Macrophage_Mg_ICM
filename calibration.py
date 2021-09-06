from scipy.optimize import differential_evolution

import numpy as np
import json
# import matplotlib.pyplot as plt
from models import *

#// load the samples from the original model
with open('samples.json') as json_file:
	samples = json.load(json_file)

#// load the observation
from observations import observations
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
	hyperparams_list = calib_params[PARAMS.free_params_model_n:PARAMS.free_params_model_n+PARAMS.hyperparams_n]
	hyperparams = {}
	for key,value in zip(PARAMS.hyperparams.keys(),hyperparams_list):
		hyperparams[key] = value
	free_params = {**inferred_params_model,**hyperparams}
	return free_params,inferred_params_model,hyperparams

def cost_function_study(study,free_params_model,hyperparams,select_sim = True):
	study_observations = observations[study]
	measurement_scheme = study_observations['measurement_scheme']
	simulation_period = study_observations['experiment_period'] # changing from hour to minute
	selections = ['TIME']+ list(measurement_scheme.keys())
	IDs_results = {}
	for ID in study_observations['IDs']:
		reset(Mg_M,free_params_model) #reset the model
		ID_inputs = study_observations[ID]['inputs']
		ID_observations = study_observations[ID]['expectations']
		## apply the boundary condition
		for key,value in ID_inputs.items():
			if key == 'Mg_e':
				Mg_M[key] = value*free_params_model['Mg_copy']
			else:
				Mg_M[key] = value
		ID_results = Mg_M.simulate(0,simulation_period,simulation_period,selections=selections)
		
		ID_results_selected = {}
		for key in measurement_scheme.keys():
			sim = ID_results[key]
			if select_sim == True:
				sim_selected = [sim[i-1] for i in measurement_scheme[key]]
			else:
				sim_selected = sim
			# sim_selected = [i*hyperparams[study] for i in sim_selected] #TODO check this out later
			ID_results_selected[key] = sim_selected
		IDs_results[ID] = ID_results_selected


	if study == 'Quao_2021_TRPM':
		#// normalize the outputs
		base_ID = 'Mg_.8'
		tag = 'TRPM'
		# print(IDs_results)
		IDs_results_norm = {}
		for ID,data in IDs_results.items():
			IDs_results_norm.update({ID:{}})
			n_values=[]
			for i in range(len(data[tag])):
				if IDs_results[base_ID][tag][i] != 0:
					n_value = data[tag][i]/IDs_results[base_ID][tag][i]
				else:
					raise ValueError('Denominator is zero')
					n_value = None #TODO: this is not correct
				n_values.append(n_value)
			# print(n_values)

			IDs_results_norm[ID].update({tag:n_values})
		IDs_results = IDs_results_norm

	#// now calculate the error
	IDs_errors = {}
	for ID in study_observations['IDs']:
		tag_errors = []
		ID_results = IDs_results[ID] 
		for key in measurement_scheme.keys():
			exp = ID_observations[key]['mean'] # the whole array
			sim = ID_results[key]
			# print('ID {} sim {} exp {}'.format(ID,sim,exp))
			# if key == 'Mg':
			# 	sim_selected = [i*hyperparams['Mg_normalization_f'] for i in sim_selected]
			abs_diff = [np.abs(exp_item - sim_item) for exp_item,sim_item in zip(exp,sim)]
			norm_abs_diff = np.mean(abs_diff)/((np.mean(exp)+np.mean(sim))/2)

			tag_errors.append(norm_abs_diff)
		IDs_errors[ID] = np.mean(tag_errors)
	# print(IDs_errors)
	# exit(2)
	return np.mean(list(IDs_errors.values())),IDs_results
def cost_function(calib_params):

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
			try:
				study_error,_ = cost_function_study(study=study,free_params_model=free_params_model,hyperparams=hyperparams)
			except:
				study_error = 1
			studies_errors.append(study_error)

		return np.mean(studies_errors)

	# error1 = vs_Zhao()
	error2 = vs_observations()

	# return error1
	return error2

def optimize():
	# Call instance of PSO
	results = differential_evolution(cost_function,bounds=SETTINGS.params,disp=True,maxiter=SETTINGS.max_iters)
	return results
if __name__ == '__main__':
	results = optimize()
	inferred_params = {}
	for key,value in zip(PARAMS.free_params.keys(),results.x):
		inferred_params[key] = value
	with open('inferred_params.json','w') as file:
		file.write(json.dumps(inferred_params,indent=4))
