
from pathlib import Path
import os
import sys
import copy
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
dir_output = dir_file
import matplotlib.pyplot as plt
from data.observations import observations
from tools.tools import run_model
from plots.plots_classes import Plot_bar

tick_font_size = 20
def adjust_font(ax):
	for label in (ax.get_xticklabels() + ax.get_yticklabels()):
		label.set_fontname('Times New Roman')
		label.set_fontsize(tick_font_size)

def plot_eq_mg(model,params):
	study = 'eq_mg'
	obs = observations[study]
	target_tags = list(obs['measurement_scheme'].keys())
	extra_targets = []
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']

	figsize = (12,6)
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model,params = params_copy,target_keys=target_tags+extra_targets,duration=duration)
	ii = 1

	fig = plt.figure(figsize=figsize)
	colors = ['r','g','b','purple']
	ii = 0
	for target_tag in target_tags: # one subgraph for each target_tag
		ax = fig.add_subplot(1,4,ii+1)
		obs_xx = obs['measurement_scheme'][target_tag]
		obs_yy = obs[ID]['expectations'][target_tag]['mean']
		x = results['time']
		ax.plot(x,results[target_tag],color=colors[ii],label = 'S:%s'%target_tag)
		ax.plot(obs_xx,obs_yy,linestyle='--',color=colors[ii], label = 'E:%s'%target_tag)
		ii+=1
		ax.legend()
		if target_tag == 'Mg':
			ax.set_ylim(0,20)
		elif target_tag == 'Mg_f':
			ax.set_ylim(0,2)
		elif target_tag == 'Mg_ATP_n' or target_tag == 'Mg_IM_n':
			ax.set_ylim(0,2)
	ax.set_title(study)

	fig = plt.figure(figsize=(4,4))

	for target_tag in extra_targets: # one subgraph for each target_tag
		x = results['time']
		ax.plot(x,results[target_tag],label = 'S:%s'%target_tag)
		ax.plot(obs_xx,obs_yy,linestyle='--',color=colors[ii], label = 'E:%s'%target_tag)
		ii+=1
	ax.legend()
	ax.set_title(study)
	return fig
def plot_Q21_eq_trpm(model_sbml,params):
	study = 'Q21_eq_trpm'
	obs = observations[study]
	target_tags = list(obs['measurement_scheme'].keys())
	extra_targets = []
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']

	figsize = (12,4)
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model_sbml,params = params_copy,target_keys=target_tags+extra_targets,duration=duration)
	ii = 1

	fig = plt.figure(figsize=figsize)
	colors = ['r','g','b','purple']
	ii = 0
	for target_tag in target_tags: # one subgraph for each target_tag
		ax = fig.add_subplot(1,3,ii+1)
		obs_xx = obs['measurement_scheme'][target_tag]
		obs_yy = obs[ID]['expectations'][target_tag]['mean']
		x = results['time']
		ax.plot(x,results[target_tag],color=colors[ii],label = 'S:%s'%target_tag)
		ax.plot(obs_xx,obs_yy,linestyle='--',color=colors[ii], label = 'E:%s'%target_tag)
		ii+=1
		ax.legend()
	ax.set_title(study)

	return fig
def plot_Q21_eq_h3s10(model_sbml,params):
	study = 'Q21_eq_h3s10'
	obs = observations[study]
	target_tags = list(obs['measurement_scheme'].keys())
	extra_targets = []
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']

	figsize = (12,4)
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model_sbml,params = params_copy,target_keys=target_tags+extra_targets,duration=duration)
	ii = 1

	fig = plt.figure(figsize=figsize)
	colors = ['r','g','b','purple']
	ii = 0
	for target_tag in target_tags: # one subgraph for each target_tag
		ax = fig.add_subplot(1,3,ii+1)
		obs_xx = obs['measurement_scheme'][target_tag]
		obs_yy = obs[ID]['expectations'][target_tag]['mean']
		x = results['time']
		ax.plot(x,results[target_tag],color=colors[ii],label = 'S:%s'%target_tag)
		ax.plot(obs_xx,obs_yy,linestyle='--',color=colors[ii], label = 'E:%s'%target_tag)
		ii+=1
		ax.legend()
	ax.set_ylim(0,2)
	ax.set_title(study)

	return fig
def plot_R05_mg_f_n(model,params):
	study = 'R05_mg_f_n'
	obs = observations[study]
	target_tags = list(obs['measurement_scheme'].keys())
	extra_plot_tags = ['Mg','Mg_ATP']
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']

	figsize = (8,4)
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model,params = params_copy,target_keys=target_tags+extra_plot_tags,duration=duration)

	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(1,2,1)
	target_tag = target_tags[0]
		
	x = results['time']
	ax.plot(x,results[target_tag],label = 'S: %s'%target_tag)

	obs_xx = obs['measurement_scheme'][target_tag]
	obs_yy = obs[ID]['expectations'][target_tag]['mean']
	ax.plot(obs_xx,obs_yy,linestyle='--',color='r', label = 'E: %s'%target_tag)
	ax.legend()
	adjust_font(ax)
	##-----extra plot----#
	ax = fig.add_subplot(1,2,2)
	for tag in extra_plot_tags:
		ax.plot(x,results[tag],label = 'S: %s'%tag)
	ax.legend()
	
	ax.set_title(study)
	return fig
def plot_Q21_nTRPM(model_sbml,macrophage_model,params):
	study = 'Q21_nTRPM'
	obs = observations[study]

	##-- figure 1
	results = macrophage_model.simulate_study(params= params,study=study)
	plot_obj = Plot_bar(study=study,observations=observations,errors={},destination = dir_output)
	fig = plot_obj.plot(simulation_results=results)
	##-- figure 2
	if False:
		figsize = (8,4)
		fig2 = plt.figure(figsize=figsize)
		ax = fig2.add_subplot(1,2,1)
		IDs = obs['IDs']
		target_tags = ['Mg_ATP']
		duration = obs['experiment_period']
		for ID in IDs:
			inputs = obs[ID]['inputs']
			params_copy = copy.deepcopy(params)
			for key,value in inputs.items():
				params_copy[key] = value
			results = run_model(model_sbml,params = params_copy,target_keys=target_tags,duration=duration)
			for tag in target_tags:
				ax.plot(results['time'],results[tag], label = '%s: %s'%(ID,tag))
		ax.legend()
		ax = fig2.add_subplot(1,2,2)
		IDs = obs['IDs']
		target_tags = ['nTRPM']
		duration = obs['experiment_period']
		for ID in IDs:
			inputs = obs[ID]['inputs']
			params_copy = copy.deepcopy(params)
			for key,value in inputs.items():
				params_copy[key] = value
			results = run_model(model_sbml,params = params_copy,target_keys=target_tags,duration=duration)
			for tag in target_tags:
				ax.plot(results['time'],results[tag], label = '%s: %s'%(ID,tag))
		ax.legend()
		
	# return fig,fig2
def plot_Q21_TRPM(model_sbml,macrophage_model,params):
	study = 'Q21_TRPM'
	obs = observations[study]

	# figure 1
	results = macrophage_model.simulate_study(params= params,study=study)
	plot_obj = Plot_bar(study=study,observations=observations,errors={},destination = dir_output)
	fig = plot_obj.plot(simulation_results=results)
	if False:
		# figure 2
		figsize = (4,4)
		fig2 = plt.figure(figsize=figsize)
		ax = fig2.add_subplot(1,1,1)
		IDs = obs['IDs']
		target_tags = ['TRPM']
		duration = obs['experiment_period']
		for ID in IDs:
			inputs = obs[ID]['inputs']
			params_copy = copy.deepcopy(params)
			for key,value in inputs.items():
				params_copy[key] = value
			results = run_model(model_sbml,params = params_copy,target_keys=target_tags,duration=duration)
			for tag in target_tags:
				ax.plot(results['time'],results[tag], label = '%s: %s'%(ID,tag))
		ax.legend()
		
	# return fig,fig2
def plot_Q21_Mg(model,params):
	study = 'Q21_Mg'
	obs = observations[study]
	target_tags = list(obs['measurement_scheme'].keys())
	extra_plot_tags = ['Mg','Mg_ATP']
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']

	figsize = (8,4)
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model,params = params_copy,target_keys=target_tags+extra_plot_tags,duration=duration)

	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(1,2,1)
	target_tag = target_tags[0]
		
	x = results['time']
	ax.plot(x,results[target_tag],label = 'S: %s'%target_tag)

	obs_xx = obs['measurement_scheme'][target_tag]
	obs_yy = obs[ID]['expectations'][target_tag]['mean']
	ax.plot(obs_xx,obs_yy,linestyle='--',color='r', label = 'E: %s'%target_tag)
	ax.legend()
	adjust_font(ax)
	##-----extra plot----#
	ax = fig.add_subplot(1,2,2)
	for tag in extra_plot_tags:
		ax.plot(x,results[tag],label = 'S: %s'%tag)
	ax.legend()
	
	ax.set_title(study)
	return fig
def plot_Q21_nM7CK(model_sbml,macrophage_model,params):
	study = 'Q21_nM7CK'
	obs = observations[study]

	# figure 1
	results = macrophage_model.simulate_study(params= params,study=study)
	plot_obj = Plot_bar(study=study,observations=observations,errors={},destination = dir_output)
	fig = plot_obj.plot(simulation_results=results)
	if False:
		# figure 2
		figsize = (6,6)
		fig2 = plt.figure(figsize=figsize)
		ax = fig2.add_subplot(1,1,1)
		IDs = obs['IDs']
		target_tags = ['nM7CK']
		duration = obs['experiment_period']
		for ID in IDs:
			inputs = obs[ID]['inputs']
			params_copy = copy.deepcopy(params)
			for key,value in inputs.items():
				params_copy[key] = value
			results = run_model(model_sbml,params = params_copy,target_keys=target_tags,duration=duration)
			for tag in target_tags:
				ax.plot(results['time'],results[tag], label = '%s: %s'%(ID,tag))
		ax.legend()
		
	# return fig,fig2
def plot_Q21_H3S10(model_sbml,macrophage_model,params):
	study = 'Q21_H3S10'
	obs = observations[study]

	# figure 1
	results = macrophage_model.simulate_study(params= params,study=study)
	plot_obj = Plot_bar(study=study,observations=observations,errors={},destination = dir_output)
	fig = plot_obj.plot(simulation_results=results)
	# figure 2
	figsize = (6,6)
	fig2 = plt.figure(figsize=figsize)
	ax = fig2.add_subplot(1,1,1)
	IDs = obs['IDs']
	target_tags = ['pH3S10']
	duration = obs['experiment_period']
	for ID in IDs:
		inputs = obs[ID]['inputs']
		params_copy = copy.deepcopy(params)
		for key,value in inputs.items():
			params_copy[key] = value
		results = run_model(model_sbml,params = params_copy,target_keys=target_tags,duration=duration)
		for tag in target_tags:
			ax.plot(results['time'],results[tag], label = '%s: %s'%(ID,tag))
	ax.legend()
		
	return fig,fig2
def plot_R05_mg_n(model,params):
	study = 'R05_mg_n'
	obs = observations[study]
	target_tags = list(obs['measurement_scheme'].keys())
	extra_plot_tags = ['Mg','Mg_ATP']
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']

	figsize = (8,4)
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model,params = params_copy,target_keys=target_tags+extra_plot_tags,duration=duration)
	ii = 1

	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(1,2,1)
	for target_tag in target_tags: # one subgraph for each target_tag
		
		x = results['time']
		ax.plot(x,results[target_tag],label = 'S: %s'%target_tag)

		obs_xx = obs['measurement_scheme'][target_tag]
		obs_yy = obs[ID]['expectations'][target_tag]['mean']
		ax.plot(obs_xx,obs_yy,linestyle='--',color='r', label = 'E: %s'%target_tag)
		ax.legend()
		
		ii+=1
	##-----extra plot----#
	ax = fig.add_subplot(1,2,2)
	for tag in extra_plot_tags:
		ax.plot(x,results[tag],label = 'S: %s'%tag)
	ax.legend()
	ax.set_title(study)
	return fig
def plot_S12_mg(Mg,obs):
	figsize = (4,4)
	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(1,1,1)
	obs_xx = obs['measurement_scheme']['Mg_n']
	obs_yy = obs['Mg_2a5']['expectations']['Mg_n']['mean']
	Mg_0 = Mg[0]
	obs_yy_scaled = [i*Mg_0 for i in obs_yy] 
	x = range(0,len(Mg))
	ax.plot(x,Mg,label = 'S')
	ax.plot(obs_xx,obs_yy_scaled,linestyle='--',color='r', label = 'E')
	ax.legend()
	ax.set_title('Norm Mg2+')
	return fig