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
from models.models import Macrophage
from plots.plotTools import plotTools, Specs,labels
from data.observations import t2m

def P21_plot(model_sbml,model_macrophage,params,observations):
	figsize = (20,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	obs = observations

	width_ratios= [1,1,1,1.5]
	jj = 1
	fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)

	g_size = 4
	study_tag,targets = 'eq_IL8',['nIL8','nIL8R']
	ID = obs[study_tag]['IDs'][0]
	plotTools.run_plot_line_multi_target(ax=axes[0],study_tag=study_tag,targets=targets,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])

	g_size = 4
	study_tag,target = 'M05_IT','nIRAK4'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=axes[1],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)

	g_size = 4
	study_tag,target = 'M05_IT','naTRAF6'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=axes[2],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)

	g_size = 6
	study_tag,target = 'M05_NFKB','nNFKB_n'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=axes[3],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
		
	fig.tight_layout()
def P22_plot(model_sbml,model_macrophage,params,observations):
	figsize = (21,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	obs = observations
	nn = 5
	jj = 1
	fig, axes = plt.subplots(1, nn, gridspec_kw={'width_ratios': [1.5,1.5,1.5,1.5,1]},figsize=figsize)

	g_size = 6
	study_tag,target = 'M18','nIFNGR'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=axes[0],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)

	g_size = 6
	study_tag,target = 'M18','nIL4R'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=axes[1],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)

	g_size = 6
	study_tag,target = 'M18','nIL1b'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=axes[2],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)

	g_size = 6
	study_tag,target = 'M18','nIL10'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=axes[3],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)

	# g_size = 5
	# study_tag,target = 'Q21_Mg_IL8','IL8'
	# IDs = obs[study_tag]['IDs']
	# plotTools.run_plot_bar(ax=ax0,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar2',IDs=IDs)

	# g_size = 4
	# study_tag,targets = 'eq_IL8',['nIL8','nIL8R']
	# ID = obs[study_tag]['IDs'][0]
	# plotTools.run_plot_line_multi_target(ax=axes[4],study_tag=study_tag,targets=targets,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])

	
	fig.tight_layout()


def plot_Q21_Mg_IL8(ax,model_sbml,model_macrophage,params,observations):
	study_tag = 'Q21_Mg_IL8'
	target = 'IL8'
	obs = observations[study_tag]
	
	sim_results = model_macrophage.simulate_study(study_tag = study_tag,params= params,study=observations[study_tag])

	IDs = observations[study_tag]['IDs']
	IDs = ['0.08 mM','8 mM']

	# plotTools.run_plot_bar(ax=ax,model = model_macrophage,study_tag=study_tag,target=target,sims=sim_results,study=observations[study_tag],plot_t= 'bar2', IDs=IDs)


	### figure 2###
	target_keys = ['IL8']
	input_tag = 'Mg_e'
	input_values = [0.08,0.5,8] 
	fig = plt.figure()
	jj = 1
	for tag in target_keys:
		ax = fig.add_subplot(1,2,jj)
		for value in input_values:
			params_copy = {**params,input_tag:value}
			results = run_model(model=model_sbml,params =params_copy,duration = 180,target_keys=[tag])
			xx = results['time']
			ax.plot(xx, results[tag], label = '{}: {}'.format(tag,value))
		ax.legend()
		jj+=1
	fig.tight_layout()

def P4_plot(model_sbml,model_macrophage,params,observations): # plotting sim vs exp measurements for P4
	figsize = (6,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	obs = observations
	nn = 2
	jj = 1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'S12_IkBa_mg','IKB'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
	jj+=1


	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'Q21_IkBa','IKB'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
	jj+=1
	
	fig.tight_layout()

def P1_eq_plot(model_sbml,params,observations): # equalibrium for P1 to P3
	figsize = (12,4)
	fig = plt.figure(figsize=figsize)
	obs = observations
	fig.canvas.draw()
	nn = 3
	jj = 1

	ax1 = fig.add_subplot(1,nn,jj)
	ax2 = fig.add_subplot(1,nn,jj+1)
	# print(model_sbml['Mg_e'])
	plot_eq_mg(ax1=ax1,ax2=ax2,model_sbml=model_sbml,params=params,observations=observations)
	
	
	jj+=2
	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_eq(ax=ax,model_sbml = model_sbml ,params = params,observations	=observations)
	jj+=1

	fig.tight_layout()

def P1_qualitative_plot(model_sbml,params,observations):
	figsize = (12,4)
	obs = observations
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	nn = 3
	jj = 1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'R05_mg_f_n','nMg_f'
	ID = obs[study_tag]['IDs'][0]
	plotTools.run_plot_line(ax=ax,study_tag=study_tag,target=target,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
	jj+=1
	
	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'R05_mg_n','nMg'
	ID = obs[study_tag]['IDs'][0]
	plotTools.run_plot_line(ax=ax,study_tag=study_tag,target=target,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
	
	jj+=1

	fig.tight_layout()
def P1_plot(model_sbml,model_macrophage,params,observations): # plotting sim vs exp measurements for P1 to p3
	fig_tag = 'P1_3'
	figsize = (20,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	obs = observations
	nn = 5
	jj = 1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target,ID = 'Q21_Mg','nMg_f','Mg_8'
	plotTools.run_plot_line(ax=ax,study_tag=study_tag,target=target,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
	jj+=1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'Q21_nTRPM','nTRPM'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
	jj+=1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'Q21_TRPM','TRPM'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
	jj+=1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'Q21_nM7CK','nM7CK'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
	jj+=1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'Q21_H3S10','pH3S10'
	IDs = obs[study_tag]['IDs']
	plotTools.run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
	jj+=1

	fig.tight_layout()
	plt.savefig(os.path.join(dir_output, fig_tag+'.svg'),bbox_inches='tight')






def plot_eq_mg(ax1,ax2,model_sbml,params,observations):
	axs = [ax1,ax2]
	study = 'eq_mg'
	obs = observations[study]
	targets = list(obs['measurement_scheme'].keys())
	extra_targets = []
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']

	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	model_sbml.reset()
	# print(model_sbml['Mg_f'])
	results =Macrophage.run_sbml_model(model_sbml=model_sbml,params = params_copy,selections=targets+extra_targets,duration=duration)
	
	# print(model_sbml['Mg_e'])

	ii = 0
	study1,study2 = study+'_f',study
	specs1,specs2 = Specs(study1),Specs(study2)
	for target in targets: # one subgraph for each target
		if target == 'nMg_ATP':
			continue
		obs_xx = obs['measurement_scheme'][target]
		obs_yy = obs[ID]['expectations'][target]['mean']
		x = results['time']
		sims = results[target]
		exps = obs_yy
		axs[ii].plot(obs_xx,exps,linestyle='--',color='r', label = 'Expectation')
		axs[ii].plot(x,sims,color='black',linewidth = specs1.line_width,label = '%s'%labels[target])
		axs[ii].legend()
		ii+=1
	plotTools.ax_postprocess(ax=axs[0],study_tag=study1,specs=specs1,sims=results)
	plotTools.ax_postprocess(ax=axs[1],study_tag=study2,specs=specs2,sims=results)

def plot_Q21_eq(ax,model_sbml,params,observations):
	study_tag = 'Q21_eq'
	study1,study2 = 'Q21_eq_trpm','Q21_eq_h3s10'
	targets1 = list(observations[study1]['measurement_scheme'].keys())

	targets2 = list(observations[study2]['measurement_scheme'].keys())

	targets = targets1 + targets2

	duration = observations[study1]['experiment_period']
	IDs = observations[study1]['IDs']
	ID = IDs[0] #only 1 ID
	inputs = observations[study1][ID]['inputs']
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results =Macrophage.run_sbml_model(model_sbml,params = params_copy,selections=targets,duration=duration)

	ii = 0

	specs = Specs(study_tag)

	ax.plot(results['time'],[1 for i in range(len(results['time']))],linestyle='--',color='r', label = 'Expectation')
	for target in targets: # one subgraph for each target
		x = results['time']
		ax.plot(x,results[target],linewidth=specs.line_width, color=specs.colors[ii],label = '%s'%labels[target])
		ii+=1
	
	plotTools.ax_postprocess(ax=ax,study_tag=study_tag,specs=specs,sims=results)

def plot_Q21_eq_h3s10(ax,model_sbml,params,observations):
	study = 'Q21_eq_h3s10'
	obs = observations[study]
	targets = list(obs['measurement_scheme'].keys())
	extra_targets = []
	duration = obs['experiment_period']
	IDs = obs['IDs']
	ID = IDs[0] #only 1 ID
	inputs = obs[ID]['inputs']
	
	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	results = run_model(model_sbml,params = params_copy,target_keys=targets+extra_targets,duration=duration)

	ii = 0
	for target in targets: # one subgraph for each target
		obs_xx = obs['measurement_scheme'][target]
		obs_yy = obs[ID]['expectations'][target]['mean']
		x = results['time']
		ax.plot(x,results[target],color=Graph_line.colors[ii],linewidth=PSL.line_width, label = 'S: %s'%labels[target])
		ax.plot(obs_xx,obs_yy,linestyle='--',color=Graph_line.colors[ii], label = 'E: %s'%labels[target])
		ii+=1
		ax.legend()
	Graph_line.postprocess(ax=ax,study=study)


def plot_S12_mg(Mg,observations):
	figsize = (4,4)
	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(1,1,1)
	obs_xx = observations['measurement_scheme']['nMg']
	obs_yy = observations['Mg_2a5']['expectations']['nMg']['mean']
	Mg_0 = Mg[0]
	obs_yy_scaled = [i*Mg_0 for i in obs_yy] 
	x = range(0,len(Mg))
	ax.plot(x,Mg,label = 'S')
	ax.plot(obs_xx,obs_yy_scaled,linestyle='--',color='r', label = 'E')
	ax.legend()
	ax.set_title('Norm Mg2+')
	return fig