from scipy.optimize import differential_evolution
import numpy as np
import json
# from models import *
# from observations import observations


#// load the samples from the original model
with open('samples.json') as json_file:
	samples = json.load(json_file)

class CALIBRATION:
	max_iters = 20
	free_params_model = { 
						  'k301_1':[0,1], # extra to intra transportation of Mg
						  'k301_2':[0,1], # saturation coeff of Mg intracellular
						  'k302':[0,1],# degradation rate of Mg

						  'k303':[0,1000], # production rate of IKK
						  'ikk_0':[0,1],#default production coeff

						  'k311_1':[0,100], # Mg-based production coeff of TRPM
						  'k311_2':[0,100], # fixed production coeff of TRPM
						  'n300':[0,10], # n as the power of Mg
						  'k312':[0,.5], # saturation coeff of TRPM
						  # 'k313':[0,1], # degradation of TRPM
						  'TRPM':[0,10000], #initial concentration of TRPM
						  'TRPM_0':[0,1], #default production coeff

						  # 'k314':[0,1], # production coeff M7CKs
						  # 'k315':[0,1000], # saturation coeff of M7CKs
						  # 'M7CKs':[0,10000],

						  # 'k316':[0,1], #  coeff of M7CKs nuclear translocation
						  # 'k317':[0,1], # saturation coeff of M7CKs_n
						  # 'k318':[0,1], # coeff of M7CKs cytoplasm translocation
						  # 'k319':[0,100], # saturation coeff of M7CKs
						  # 'M7CKs_n':[0,10000], 

						  # 'k320':[0,100], # coeff of H3S10 activation
						  # 'k321':[0,100], # saturation coeff of H3S10
						  # 'k322':[0,1], # degradation rate of H3S10
						  # 'k323':[0,1], # phosphorylatio rate of H3S10
						  # 'k324':[0,100], # saturation rate of pH3S10
						  # 'H3S10':[0,10000],
						  # 'pH3S10': [0,10000],

						  # 'k325':[0,100], # production rate of IL8
						  # 'k326':[0,100], # saturation rate of IL8
						  # 'k327':[0,1], # degradation of IL8
						  # 'k328': [0,1], # intra to extra transportation of IL8
						  # 'k329': [0,1], # saturation in extracellular IL8
						  # 'k330': [0,1], # extra to intra transportation of IL8
						  # 'k331': [0,1], # saturation in intracellular IL8
						  # 'IL8':[0,10000]
						  }

	free_params_model_n = len(free_params_model.keys())

	hyperparams = {
	# 'Quao_2021_TRPM':[0,1000],
					# 'Mg_normalization_f':[0,1000], # to scale the models' outputs to the standrdized format of the observations
					# 'Quao_2021_Mg':[0,1]
					}
	hyperparams_n = len(hyperparams.keys())

	free_params = {**free_params_model, **hyperparams}

	def __init__(self,model):
		self.model = model
		
	def categorize(self,calib_params): # divides the calib values between the model parameters and those for studies
		free_params_model_list = calib_params[0:self.free_params_model_n]
		inferred_params_model = {}
		for key,value in zip(self.free_params_model.keys(),free_params_model_list):
			inferred_params_model[key] = value
		hyperparams_list = calib_params[self.free_params_model_n:self.free_params_model_n+self.hyperparams_n]
		hyperparams = {}
		for key,value in zip(self.hyperparams.keys(),hyperparams_list):
			hyperparams[key] = value
		free_params = {**inferred_params_model,**hyperparams}
		return free_params,inferred_params_model,hyperparams
	@staticmethod
	def indexing(real_time,time_vector):
		index = next(x[0] for x in enumerate(time_vector) if x[1] >= real_time)
		return index
	def cost_function_study(self,study,free_params_model,hyperparams=None,select_sim = True):
		study_observations = observations[study]
		measurement_scheme = study_observations['measurement_scheme']
		simulation_period = study_observations['experiment_period'] # changing from hour to minute
		selections = ['TIME']+ list(measurement_scheme.keys())
		IDs_results = {}
		for ID in study_observations['IDs']:
			self.model.reset(free_params_model) #reset the model
			ID_inputs = study_observations[ID]['inputs']
			ID_observations = study_observations[ID]['expectations']
			## apply the boundary condition
			for key,value in ID_inputs.items():
				self.model.set(key,value)

			ID_results = self.model.simulate(0,simulation_period,selections=selections)
			time_vector = ID_results['time']
			
			ID_results_selected = {'time':time_vector}
			for key in measurement_scheme.keys():
				sim = ID_results[key]
				if select_sim == True:
					selected_indices = [self.indexing(i,time_vector) for i in measurement_scheme[key]]
					sim_selected = [sim[i] for i in selected_indices]
				else:
					sim_selected = sim
				# sim_selected = [i*hyperparams[study] for i in sim_selected] #TODO check this out later
				ID_results_selected[key] = sim_selected
			# print('originl ',ID_results)
			# print('selected ',ID_results_selected)
			IDs_results[ID] = ID_results_selected


		if study == 'Quao_2021_TRPM':
			#// normalize the outputs
			base_ID = 'Mg_.8'
			tag = 'TRPM'
			IDs_results_norm = {}
			for ID,data in IDs_results.items():
				IDs_results_norm.update({ID:{'time':data['time']}})
				n_values=[]
				for i in range(len(data[tag])):
					if IDs_results[base_ID][tag][i] != 0:
						n_value = data[tag][i]/IDs_results[base_ID][tag][i]
					else:
						raise ValueError('Denominator is zero')
						n_value = None #TODO: this is not correct
					n_values.append(n_value)

				IDs_results_norm[ID].update({tag:n_values})
			IDs_results = IDs_results_norm
		if study == 'Quao_2021_Mg':
			#// normalize the outputs
			base_index = 0
			tag = 'Mg'
			IDs_results_norm = {}
			for ID,data in IDs_results.items():
				IDs_results_norm.update({ID:{'time':data['time']}})
				n_values=[]
				for i in range(len(data[tag])):
					if data[tag][base_index] != 0:
						n_value = (data[tag][i]/data[tag][base_index])*100
					else:
						raise ValueError('Denominator is zero')
						n_value = None #TODO: this is not correct
					n_values.append(n_value)

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
				# if key == 'Mg':
				# 	sim_selected = [i*hyperparams['Mg_normalization_f'] for i in sim_selected]
				abs_diff = [np.abs(exp_item - sim_item) for exp_item,sim_item in zip(exp,sim)]
				norm_abs_diff = np.mean(abs_diff)/((np.mean(exp)+np.mean(sim))/2)
				tag_errors.append(norm_abs_diff)
			IDs_errors[ID] = np.mean(tag_errors)
		return np.mean(list(IDs_errors.values())),IDs_results
	def cost_function(self,calib_params):

		free_params,free_params_model,hyperparams = self.categorize(calib_params)

		def vs_Zhao(): # to compare the model's results with the Zhao
			
			from sample import targets, duration
			self.model.reset(params=free_params_model)
			selections = ['TIME']+targets
			sim_results = self.model.simulate(0,duration,selections=selections)
			exp_time_vector = samples['time']
			sim_time_vector = sim_results['time']
			# calculate the error for each tag by comparing the results to the original model
			check_points = [i for i in range(duration)]
			exp_indices = [self.indexing(i,exp_time_vector) for i in check_points]
			sim_indices = [self.indexing(i,sim_time_vector) for i in check_points]

			tag_errors = []
			for tag in targets:
				n_abs_values= []

				for exp_index,sim_index in zip(exp_indices,sim_indices):
					abs_value = abs(sim_results[tag][sim_index] - samples[tag][exp_index])
					mean_value = np.mean([sim_results[tag][sim_index], samples[tag][exp_index]])
					n_abs_values.append(abs_value/mean_value)
					
				tag_error = np.mean(n_abs_values)
			tag_errors.append(tag_error)
			return np.mean(tag_errors)
		def vs_observations(): # to compare the model's predictions with the observations
			studies_errors = []
			for study in observations['studies']:
				try:
					study_error,_ = self.cost_function_study(study=study,free_params_model=free_params_model,hyperparams=hyperparams)
				except:
					study_error = 1
				studies_errors.append(study_error)

			return np.mean(studies_errors)

		# error1 = vs_Zhao()
		error2 = vs_observations()
		# error = np.mean([error1,error2])
		# return error1
		# print(calib_params)
		return error2
		# return error

	def optimize(self):
		# Call instance of PSO
		results = differential_evolution(self.cost_function,
			bounds=list(self.free_params.values()),
			disp=True,
			maxiter=self.max_iters)
		inferred_params_model = {}
		for key,value in zip(self.free_params_model.keys(),results.x):
			inferred_params_model[key] = value
		with open(os.path.join(output_dir,'inferred_params_model.json'),'w') as file:
			file.write(json.dumps(inferred_params_model,indent=4))

		inferred_hyperparameters = {}
		for key,value in zip(self.hyperparams.keys(),results.x):
			inferred_hyperparameters[key] = value
		with open(os.path.join(output_dir,'inferred_hyperparameters.json'),'w') as file:
			file.write(json.dumps(inferred_hyperparameters,indent=4))

if __name__ == '__main__':
	model = MG_MODEL()
	obj = CALIBRATION(model)
	obj.optimize()
