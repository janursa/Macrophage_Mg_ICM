### from plots classes ####
import matplotlib.pyplot as plt
from pathlib import Path
import os
import sys
import copy
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
from data.observations import observations,t2m
plt.rcParams["font.family"] = "serif"
plt.style.use('seaborn-deep')
# plt.style.use('Agg')

plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

font_type = 'Arial'

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


def sort(study, sim_results): #sorts the results for bar plot
	measurement_scheme = observations[study]['measurement_scheme']
	exp_target_results = {}
	for target in measurement_scheme:
		exp_target_results[target] = []
	for target in measurement_scheme:
		for ID in observations[study]['IDs']:
			ID_observations = observations[study][ID]['expectations']
			exp_target_results[target].append(ID_observations[target])
	sim_target_results = {}
	for target in measurement_scheme:
		sim_target_results[target] = []
	for target in measurement_scheme:
		for ID,ID_result in sim_results.items():
			sim_target_results[target].append(ID_result[target])
	return exp_target_results,sim_target_results

def bar_positions(study, sim_results): # determines bar positions for bar plot
	measurement_scheme = observations[study]['measurement_scheme']
	delta = .12
	for i in range(len(measurement_scheme)):
		x_exp =[float(j) + delta for j in range(len(sim_results.keys()))]
		x_sim =[float(j) - delta for j in range(len(sim_results.keys()))]

	return x_exp,x_sim



class graph_line:
	colors = ['peru','violet','lime' ,  'yellowgreen',  'skyblue']
	def determine_xlabel(study):
		label = 'Time (h)'
		return label

	def determine_ylabel(study):
		label = 'Relative Concentration \n (To initial condition)'
		if study == 'eq_mg' or study == 'eq_mg_f':
			label = 'Concentration (ng/ml)'
		return label

	def determine_xlim(study):

		if study == '':
			pass
		else:
			raise ValueError('Define')
		return lim
	def determine_ylim(study):

		if study == 'eq_mg':
			lim = [18,19]
		elif study == 'eq_mg_f':
			lim = [0,1]
		else:
			raise ValueError('Define')
		return lim
	def xticks(ax,study):
		x_ticks = ax.get_xticks()
		if study == 'Q21_Mg':
			x_ticks_ad = np.array([0,1,2,3])*60/t2m
		else:
			x_ticks_ad = x_ticks[1:-1]
			
		x_ticks_labels = [int(i*t2m/60) for i in x_ticks_ad]
		
		ax.set_xticks(ticks = x_ticks_ad)
		ax.set_xticklabels(x_ticks_labels)
	def yticks(ax,study):
		
		if study == 'eq_mg_f':
			y_ticks_ad = [0,0.5,1]
			y_ticks_label = y_ticks_ad
		elif study == 'eq_mg':
			y_ticks_ad = [18+0.5*i for i in range(0,3)]
			y_ticks_label = y_ticks_ad
		else:
			raise ValueError('define')
			
		ax.set_yticks(ticks = y_ticks_ad)
		ax.set_yticklabels(y_ticks_label)

	def axes(ax,study):
		tick_font_size = PSL.tick_font_size
		for label in (ax.get_xticklabels() + ax.get_yticklabels()):
			label.set_fontname(font_type)
			label.set_fontsize(tick_font_size)
		ax.set_ylabel(graph_line.determine_ylabel(study=study),fontdict ={'family':font_type,'size':PSL.title_font_size})
		ax.set_xlabel(graph_line.determine_xlabel(study=study),fontdict ={'family':font_type,'size':PSL.title_font_size})
		
		graph_line.xticks(ax,study)
		try:
			graph_line.yticks(ax,study)
		except:
			pass

		try:
			ax.set_xlim(graph_line.determine_xlim(study))
		except:
			pass

		try:
			ax.set_ylim(graph_line.determine_ylim(study))
		except:
			pass
	def legend(ax,study):
		legend_location = (1,1.1)
		if study == 'R05_mg_f_n':
			legend_location = (1.1,.7)
		elif study == 'eq_mg' or study == 'eq_mg_f':
			legend_location = (1,1.1)
		elif study == 'R05_mg_n':
			legend_location = (1.1,1.1)
		elif study == 'Q21_Mg':
			legend_location = (1.1,.8)
		elif study == 'Q21_eq':
			legend_location = (1.14,.75)

		elif study == 'Q21_nTRPM' or study == 'Q21_TRPM' or study == 'Q21_nM7CK' or study == 'Q21_H3S10':
			legend_location = (0.7,1.1)
		ncol = 1
		ax.legend(bbox_to_anchor=legend_location,loc = 'upper right', borderaxespad=2,prop={ 'family':font_type,'size':PSL.legend_font_size},ncol=ncol)

	def title(ax,study):
		title = ''
		if study == 'eq_mg':
			# title = 'Intracellular Mg$^{2+}$'
			title = 'Equilibrium'
			pass
		elif study == 'eq_mg_f':
			# title = 'free Mg$^{2+}$'
			pass
		elif study == 'Q21_eq_trpm':
			# title = 'Equilibrium'
			pass
		elif study == 'R05_mg_f_n':
			title = 'Free Mg$^{2+}$ ions'
		elif study == 'R05_mg_n':
			title = 'Total Mg$^{2+}$ ions'
		elif study == 'Q21_Mg':
			title = 'Free Mg$^{2+}$ ions'

		ax.set_title(title,fontdict ={'family':font_type,'size':PSL.title_font_size},pad = 10)

	def postprocess(ax,study):
		graph_line.axes(ax,study)
		graph_line.legend(ax,study)
		graph_line.title(ax,study)
	
class PSL: #plot specification line
	legend_font_size = 15
	tick_font_size = 21
	title_font_size = 21 
	R2_font_size  = 18
	line_width = 3



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
from data.observations import observations,t2m
from tools.tools import run_model
from plots.plots_classes import graph_bar


def P1_3_eq_plot(model_sbml,params): # equalibrium for P1 to P3
	figsize = (12,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	nn = 3
	jj = 1

	ax1 = fig.add_subplot(1,nn,jj)
	ax2 = fig.add_subplot(1,nn,jj+1)
	plot_eq_mg(ax1=ax1,ax2=ax2,model_sbml=model_sbml,params=params)
	jj+=2
	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_eq(ax=ax,model_sbml = model_sbml ,params = params)
	jj+=1

	fig.tight_layout()

def P1_3_qualitative_plot(model_sbml,params):
	figsize = (12,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	nn = 3
	jj = 1

	ax = fig.add_subplot(1,nn,jj)
	plot_R05_mg_f_n(ax=ax,model_sbml=model_sbml,params=params)
	
	jj+=1
	ax = fig.add_subplot(1,nn,jj)
	plot_R05_mg_n(ax=ax,model_sbml=model_sbml,params=params)
	jj+=1

	fig.tight_layout()
def P1_3_plot(model_sbml,model_macrophage,params):
	figsize = (20,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	nn = 5
	jj = 1

	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_Mg(ax=ax,model_sbml=model_sbml,params=params)
	jj+=1
	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_nTRPM(ax=ax,model_sbml = model_sbml,macrophage_model = model_macrophage ,params = params)
	jj+=1
	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_TRPM(ax=ax,model_sbml = model_sbml,macrophage_model = model_macrophage ,params = params)
	jj+=1
	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_nM7CK(ax=ax,model_sbml = model_sbml,macrophage_model = model_macrophage ,params = params)
	jj+=1
	ax = fig.add_subplot(1,nn,jj)
	plot_Q21_H3S10(ax=ax,model_sbml = model_sbml,macrophage_model = model_macrophage ,params = params)
	jj+=1

	fig.tight_layout()



def plot_eq_mg(ax1,ax2,model_sbml,params):
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

def plot_Q21_eq(ax,model_sbml,params):
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

def plot_Q21_eq_h3s10(ax,model_sbml,params):
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


def plot_R05_mg_f_n(ax,model_sbml,params):
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

def plot_Q21_nTRPM(ax,model_sbml,macrophage_model,params):
	study = 'Q21_nTRPM'
	obs = observations[study]
	target = list(observations[study]['measurement_scheme'].keys())[0]
	##-- figure 1
	sim_results = macrophage_model.simulate_study(params= params,study=study)

	plot_obj = graph_bar(study=study,observations=observations,errors={},destination = dir_output)
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
def plot_Q21_TRPM(ax,model_sbml,macrophage_model,params):
	study = 'Q21_TRPM'
	obs = observations[study]

	# figure 1
	results = macrophage_model.simulate_study(params= params,study=study)
	plot_obj = graph_bar(study=study,observations=observations,errors={},destination = dir_output)
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
def plot_Q21_Mg(ax,model_sbml,params):
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
	
def plot_Q21_nM7CK(ax,model_sbml,macrophage_model,params):
	study = 'Q21_nM7CK'
	obs = observations[study]

	# figure 1
	results = macrophage_model.simulate_study(params= params,study=study)
	plot_obj = graph_bar(study=study,observations=observations,errors={},destination = dir_output)
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
def plot_Q21_H3S10(ax,model_sbml,macrophage_model,params):
	study = 'Q21_H3S10'
	obs = observations[study]

	# figure 1
	results = macrophage_model.simulate_study(params= params,study=study)
	plot_obj = graph_bar(study=study,observations=observations,errors={},destination = dir_output)
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
		
def plot_R05_mg_n(ax,model_sbml,params):
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