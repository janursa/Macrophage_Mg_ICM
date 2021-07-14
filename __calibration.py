from PSO import PSO
import numpy as np
import json
# import matplotlib.pyplot as plt
from models import Mg_model, PARAMS,run

# load the samples from the original model
with open('samples.json') as json_file:
    samples = json.load(json_file)

# Set-up hyperparameters
class SETTINGS:
	options = {'c1': 2, 'c2': 0, 'w':.1}
	dimensions = len(PARAMS.free_params)
	n_particles = 20
	iters = 100


def cost_function(calib_params):
	costs = [] # cost values for each particle
	for j in range(SETTINGS.n_particles):
		# set the parameters of the model
		for i in range(len(PARAMS.free_params)):
			free_param_name = PARAMS.free_params[i]
			Mg_model[free_param_name] = calib_params[j,i]
		# run the model
		results = run(Mg_model,selections = PARAMS.targets, duration = PARAMS.duration)
		# calculate the error for each tag by comparing the results to the original model
		tag_errors = []
		for tag in PARAMS.targets:
			print('\n results',results[tag])
			print('samples',samples[tag])
			tag_error = np.mean(abs(results[tag]-samples[tag]))
			print('error',tag_error)
			tag_errors.append(tag_error)
		# final error is the mean of all tags
		costs.append(np.mean(tag_errors))

	return costs

def optimize():

	particles = np.random.uniform(-5, 5, (SETTINGS.n_particles, SETTINGS.dimensions))
	velocities = (np.random.random((SETTINGS.n_particles, SETTINGS.dimensions)) - 0.5) / 10

	pso_1 = PSO(particles.copy(), velocities.copy(), cost_function, 
		w=SETTINGS.options['w'], c_1=SETTINGS.options['c1'], c_2=SETTINGS.options['c2'], 
		max_iter=100,auto_coef=False)
	while pso_1.next():
		print(pso_1.iter)

optimize()