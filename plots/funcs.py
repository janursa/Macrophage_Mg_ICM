from pathlib import Path
import os
import sys
import copy
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
dir_output = dir_file
import matplotlib.pyplot as plt
import numpy as np
from tools.tools import run_model
from plots.graphs import graph_bar,graph_line,PSL

labels = {
	'Mg_f':'free Mg$^{2+}$',
	'Mg_f_n':'free Mg$^{2+}$',
	'Mg':'total Mg$^{2+}$',
	'Mg_n':'total Mg$^{2+}$',

	'TRPM': 'TRPM',
	'nTRPM': 'Nuclear TRPM',
	'nM7CK': 'Nuclear M7CK',
	'pH3S10': 'Phos H3S10',
}

def P4_plot(model_sbml,model_macrophage,params,observations): # plotting sim vs exp measurements for P4
	figsize = (4,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	nn = 2
	jj = 1

	ax = fig.add_subplot(1,nn,jj)

	plot_S12_IkBa_mg(ax=ax,model_sbml=model_sbml,model_macrophage=model_macrophage,params=params,observations=observations)
	jj+=1

	ax = fig.add_subplot(1,nn,jj)

	plot_Q21_IkBa(ax=ax,model_sbml=model_sbml,model_macrophage=model_macrophage,params=params,observations=observations)
	jj+=1
	
	fig.tight_layout()

def plot_Q21_IkBa(ax,model_sbml,model_macrophage,params,observations):
	study_tag = 'Q21_IkBa'
	obs = observations[study_tag]
	target = list(observations[study_tag]['measurement_scheme'].keys())[0]
	
	sim_results = model_macrophage.simulate_study(study_tag = study_tag,params= params,study=observations[study_tag])

	plot_obj = graph_bar(study=study_tag,observations=observations,errors={},destination = dir_output)
	plot_obj.plot(ax = ax,simulation_results=sim_results)

	### figure 2###
	tag = 'IKB'
	input_tag = 'Mg_e'
	input_values = [0.08,0.5,8] 
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	for value in input_values:
		params_copy = {**params,input_tag:value}
		results = run_model(model=model_sbml,params =params_copy,duration = 180,target_keys=[tag])
		xx = results['time']
		ax.plot(xx, results[tag], label = '{}'.format(value))
	ax.legend()

def plot_S12_IkBa_mg(ax,model_sbml,model_macrophage,params,observations):
	study_tag = 'S12_IkBa_mg'
	obs = observations[study_tag]
	target = list(observations[study_tag]['measurement_scheme'].keys())[0]
	
	sim_results = model_macrophage.simulate_study(study_tag=study_tag,params= params,study=observations[study_tag])

	plot_obj = graph_bar(study=study_tag,observations=observations,errors={},destination = dir_output)
	plot_obj.plot(ax = ax,simulation_results=sim_results)

	### figure 2###
	tag = 'IKB'
	input_tag = 'Mg_e'
	input_values = [0.5,2.5] 
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	for value in input_values:
		params_copy = {**params,input_tag:value}
		results = run_model(model=model_sbml,params =params_copy,duration = 180,target_keys=[tag])
		xx = results['time']
		ax.plot(xx, results[tag], label = '{}'.format(value))
	ax.legend()

def P1_3_eq_plot(model_sbml,params,observations): # equalibrium for P1 to P3
	figsize = (12,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	nn = 3
	jj = 1

	ax1 = fig.add_subplot(1,nn,jj)
	ax2 = fig.add_subplot(1,nn,jj+1)
	plot_eq_mg(ax1=ax1,ax2=ax2,model_sbml=model_sbml,params=params,observations=observations)
	jj+=2
	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_eq(ax=ax,model_sbml = model_sbml ,params = params,observations	=observations)
	jj+=1

	fig.tight_layout()

def P1_3_qualitative_plot(model_sbml,params,observations):
	figsize = (12,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	nn = 3
	jj = 1

	ax = fig.add_subplot(1,nn,jj)
	plot_R05_mg_f_n(ax=ax,model_sbml=model_sbml,params=params,observations=observations)
	
	jj+=1
	ax = fig.add_subplot(1,nn,jj)
	plot_R05_mg_n(ax=ax,model_sbml=model_sbml,params=params,observations=observations)
	jj+=1

	fig.tight_layout()
def P1_3_plot(model_sbml,model_macrophage,params,observations): # plotting sim vs exp measurements for P1 to p3
	figsize = (20,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	nn = 5
	jj = 1

	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_Mg(ax=ax,model_sbml=model_sbml,params=params,observations=observations)
	jj+=1
	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_nTRPM(ax=ax,model_sbml = model_sbml,model_macrophage = model_macrophage ,params = params,observations=observations)
	jj+=1
	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_TRPM(ax=ax,model_sbml = model_sbml,model_macrophage = model_macrophage ,params = params,observations=observations)
	jj+=1
	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_nM7CK(ax=ax,model_sbml = model_sbml,model_macrophage = model_macrophage ,params = params,observations=observations)
	jj+=1
	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_H3S10(ax=ax,model_sbml = model_sbml,model_macrophage = model_macrophage ,params = params,observations=observations)
	jj+=1

	fig.tight_layout()



def plot_eq_mg(ax1,ax2,model_sbml,params,observations):
	axs = [ax1,ax2]
	study = 'eq_mg'
	obs = observations[study]
	target_tags = list(obs['measurement_scheme'].keys())
	extra_targets = []
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model_sbml,params = params_copy,target_keys=target_tags+extra_targets,duration=duration)
	ii = 0
	for target_tag in target_tags: # one subgraph for each target_tag
		if target_tag == 'Mg_ATP_n':
			continue
		obs_xx = obs['measurement_scheme'][target_tag]
		obs_yy = obs[ID]['expectations'][target_tag]['mean']
		x = results['time']
		sims = results[target_tag]
		exps = obs_yy
		axs[ii].plot(obs_xx,exps,linestyle='--',color='r', label = 'Expectation')
		axs[ii].plot(x,sims,color='black',linewidth = PSL.line_width,label = '%s'%labels[target_tag])
		axs[ii].legend()
		ii+=1

	graph_line.postprocess(ax=axs[0],study=study+'_f')
	graph_line.postprocess(ax=axs[1],study=study)

def plot_Q21_eq(ax,model_sbml,params,observations):
	study1,study2 = 'Q21_eq_trpm','Q21_eq_h3s10'
	target_tags1 = list(observations[study1]['measurement_scheme'].keys())

	target_tags2 = list(observations[study2]['measurement_scheme'].keys())

	target_tags = target_tags1 + target_tags2

	duration = observations[study1]['experiment_period']
	IDs = observations[study1]['IDs']
	ID = IDs[0] #only 1 ID
	inputs = observations[study1][ID]['inputs']
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model_sbml,params = params_copy,target_keys=target_tags,duration=duration)

	ii = 0

	ax.plot(results['time'],[1 for i in range(len(results['time']))],linestyle='--',color='r', label = 'Expectation')
	for target_tag in target_tags: # one subgraph for each target_tag
		x = results['time']
		ax.plot(x,results[target_tag],linewidth=PSL.line_width, color=graph_line.colors[ii],label = '%s'%labels[target_tag])
		ii+=1
	
	graph_line.postprocess(ax=ax,study='Q21_eq')

def plot_Q21_eq_h3s10(ax,model_sbml,params,observations):
	study = 'Q21_eq_h3s10'
	obs = observations[study]
	target_tags = list(obs['measurement_scheme'].keys())
	extra_targets = []
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model_sbml,params = params_copy,target_keys=target_tags+extra_targets,duration=duration)

	ii = 0
	for target_tag in target_tags: # one subgraph for each target_tag
		obs_xx = obs['measurement_scheme'][target_tag]
		obs_yy = obs[ID]['expectations'][target_tag]['mean']
		x = results['time']
		ax.plot(x,results[target_tag],color=graph_line.colors[ii],linewidth=PSL.line_width, label = 'S: %s'%labels[target_tag])
		ax.plot(obs_xx,obs_yy,linestyle='--',color=graph_line.colors[ii], label = 'E: %s'%labels[target_tag])
		ii+=1
		ax.legend()
	graph_line.postprocess(ax=ax,study=study)


def plot_R05_mg_f_n(ax,model_sbml,params,observations):
	study = 'R05_mg_f_n'
	obs = observations[study]
	target_tags = list(obs['measurement_scheme'].keys())
	extra_plot_tags = ['Mg','Mg_ATP']
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']

	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model_sbml,params = params_copy,target_keys=target_tags+extra_plot_tags,duration=duration)

	target_tag = target_tags[0]
		
	x = results['time']
	ax.plot(x,results[target_tag],color = 'black',linewidth=PSL.line_width,label = 'S')

	obs_xx = obs['measurement_scheme'][target_tag]
	obs_yy = obs[ID]['expectations'][target_tag]['mean']
	ax.plot(obs_xx,obs_yy,linestyle='--',color = 'r', label = 'E')
	ax.legend()
	
	graph_line.postprocess(ax=ax,study=study)
	##-----extra plot----#

def plot_Q21_nTRPM(ax,model_sbml,model_macrophage,params,observations):
	study_tag = 'Q21_nTRPM'
	
	target = list(observations[study_tag]['measurement_scheme'].keys())[0]
	##-- figure 1
	sim_results = model_macrophage.simulate_study(params= params,study=observations[study_tag])

	plot_obj = graph_bar(study=study_tag,observations=observations,errors={},destination = dir_output)
	plot_obj.plot(ax = ax,simulation_results=sim_results)
	# exp_target_results,sim_target_results = sort(study=study,sim_results=sim_results)
	# # _,processed_detailed_errors_sorted = sort(study=study,sim_results= processed_detailed_errors)
	# x_exps,x_sims = bar_positions(study=study,sim_results= sim_results)

	# sims = [item[0] for item in sim_target_results[target]]
	# ax.bar(x=x_sims,height=sims,width = PSB.bar_width, label = "S", 
	# 		facecolor = graph.colors[0],
	# 		 edgecolor="black", yerr =  0,
	# 		 error_kw = dict(capsize= PSB.error_bar_width))
	# exp_values = [exp_target_results[target][i]['mean'] for i in range(len(exp_target_results[target]))]
	# exp_std = [exp_target_results[target][i]['std'] for i in range(len(exp_target_results[target]))]
	# exp_values = [item[0] for item in exp_values]
	# exp_std = [item[0] for item in exp_std]
	# ax.bar(x=x_exps,height=exp_values,width = PSB.bar_width, label = 'E', 
	# 			facecolor = graph.colors[1],hatch=r'\\\\',
	# 			 edgecolor="black", yerr =  exp_std,
	# 			 error_kw = dict(capsize= PSB.error_bar_width))
	# x_ticks = [(i+j)/2 for i,j in zip(x_sims,x_exps)]
	# graph.postprocess(ax=ax,study=study,graph_t='bar')
	# processed_detailed_errors_sorted_values = [item[0] for item in processed_detailed_errors_sorted[target]]
	# ax.bar(x=x_sim,height=sim_values,width = self.bar_width, label = "S", 
	# 		facecolor = self.colors[0],
	# 		 edgecolor="black", yerr =  processed_detailed_errors_sorted_values,
	# 		 error_kw = dict(capsize= self.error_bar_width))
	
	# plot_obj = graph_bar(study=study,observations=observations,errors={},destination = dir_output)
	# fig = plot_obj.plot(simulation_results=results)
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
def plot_Q21_TRPM(ax,model_sbml,model_macrophage,params,observations):
	study_tag = 'Q21_TRPM'
	# figure 1
	results = model_macrophage.simulate_study(params= params,study=observations[ study_tag])
	plot_obj = graph_bar(study=study_tag,observations=observations,errors={},destination = dir_output)
	fig = plot_obj.plot(ax=ax,simulation_results=results)
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
def plot_Q21_Mg(ax,model_sbml,params,observations):
	study = 'Q21_Mg'
	obs = observations[study]
	target_tags = list(obs['measurement_scheme'].keys())
	extra_plot_tags = ['Mg','Mg_ATP']
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']

	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model_sbml,params = params_copy,target_keys=target_tags+extra_plot_tags,duration=duration)

	target_tag = target_tags[0]
		
	x = results['time']
	ax.plot(x,results[target_tag],color='black',linewidth=PSL.line_width, label = 'S')

	obs_xx = obs['measurement_scheme'][target_tag]
	obs_yy = obs[ID]['expectations'][target_tag]['mean']
	ax.plot(obs_xx,obs_yy,linestyle='--',color='r', linewidth=PSL.line_width, label = 'E')
	ax.legend()
	graph_line.postprocess(ax,study)
	
def plot_Q21_nM7CK(ax,model_sbml,model_macrophage,params,observations):
	study_tag = 'Q21_nM7CK'
	# figure 1
	results = model_macrophage.simulate_study(params= params,study=observations[study_tag])
	plot_obj = graph_bar(study=study_tag,observations=observations,errors={},destination = dir_output)
	fig = plot_obj.plot(ax=ax,simulation_results=results)
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
def plot_Q21_H3S10(ax,model_sbml,model_macrophage,params,observations):
	study_tag = 'Q21_H3S10'
	# figure 1
	results = model_macrophage.simulate_study(params= params,study=observations[study_tag])
	plot_obj = graph_bar(study=study_tag,observations=observations,errors={},destination = dir_output)
	fig = plot_obj.plot(ax=ax,simulation_results=results)
	# figure 2
	if False:
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
		
def plot_R05_mg_n(ax,model_sbml,params,observations):
	study = 'R05_mg_n'
	obs = observations[study]
	target_tags = list(obs['measurement_scheme'].keys())
	extra_plot_tags = ['Mg','Mg_ATP']
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model_sbml,params = params_copy,target_keys=target_tags+extra_plot_tags,duration=duration)
	ii = 0
	target_tag = target_tags[0]
		
	x = results['time']
	ax.plot(x,results[target_tag],color='black',linewidth=PSL.line_width,label = 'S')

	obs_xx = obs['measurement_scheme'][target_tag]
	obs_yy = obs[ID]['expectations'][target_tag]['mean']
	ax.plot(obs_xx,obs_yy,linestyle='--',color='r', label = 'E')
	ax.legend()
	
	graph_line.postprocess(ax=ax,study=study)
	##-----extra plot----#
	# ax = fig.add_subplot(1,2,2)
	# for tag in extra_plot_tags:
	# 	ax.plot(x,results[tag],label = 'S: %s'%tag)
	# ax.legend()
	# return fig
def plot_S12_mg(Mg,observations):
	figsize = (4,4)
	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(1,1,1)
	obs_xx = observations['measurement_scheme']['Mg_n']
	obs_yy = observations['Mg_2a5']['expectations']['Mg_n']['mean']
	Mg_0 = Mg[0]
	obs_yy_scaled = [i*Mg_0 for i in obs_yy] 
	x = range(0,len(Mg))
	ax.plot(x,Mg,label = 'S')
	ax.plot(obs_xx,obs_yy_scaled,linestyle='--',color='r', label = 'E')
	ax.legend()
	ax.set_title('Norm Mg2+')
	return fig