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
from plots import plots
from data.observations import t2m


labels = {
	'Mg_f':'free Mg$^{2+}$',
	'Mg_f_n':'free Mg$^{2+}$',
	'Mg':'total Mg$^{2+}$',
	'Mg_n':'total Mg$^{2+}$',

	'TRPM': 'TRPM',
	'nTRPM': 'Nuclear TRPM',
	'nM7CK': 'Nuclear M7CK',
	'pH3S10': 'Phos H3S10',

	'IL8_n': 'norm IL8',
	'IL8R_n': 'norm IL8R'
}
class process_data:
	def sort(sims,exps,target,plot_t):
		sims_sorted = []
		exps_sorted = []
		for ID in exps['IDs']:
			ID_exps = exps[ID]['expectations'][target]
			ID_sims = sims[ID][target]
			exps_sorted.append(ID_exps['mean'])
			sims_sorted.append(ID_sims)
		if plot_t == 'bar1':
			sims_sorted = [item[0] for item in sims_sorted]
			exps_sorted = [item[0] for item in exps_sorted]
		
		return sims_sorted,exps_sorted

	def normalize(study_tag,sims,exps):
		sims_n = []
		exps_n = []
		flag_sim = False
		flag_exp = False
		if study_tag=='Q21_nTRPM' or study_tag=='Q21_TRPM' or study_tag=='Q21_nM7CK' or study_tag=='Q21_H3S10':
			ctr_sim = sims[1]
			flag_sim = True
		elif study_tag == 'S12_IkBa_mg':
			ctr_sim = sims[0]
			flag_sim = True
		elif study_tag == 'Q21_IkBa':
			ctr_sim = sims[0]
			ctr_exp = exps[0]
			flag_sim = True
			flag_exp = True
		

		if flag_sim:
			for item in sims:
				sims_n.append(item/ctr_sim)
		else:
			sims_n = sims

		if flag_exp:
			for item in exps:
				exps_n.append(item/ctr_exp)
		else:
			exps_n = exps

		return sims_n,exps_n


class Specs:
	def __init__(self,study_tag):
		self.line_width = 3
		self.bar_width = .2
		self.error_bar_width = 5
		self.colors = ['lime' , 'violet', 'yellowgreen', 'peru', 'skyblue']
		self.legend_font_size = 15
		self.tick_font_size = 21
		self.title_font_size = 21 
		self.delta = .12
		self.R2_font_size  = 18
		self.D = 1.1 # the length in which all Mg dosages are plotted in a certain time point
		if study_tag == 'Q21_Mg_IL8':
			self.bar_width = .18
			self.delta = .099
			self.D = 1.3
		
		
	@staticmethod
	def determine_title(study_tag):
		label = study_tag
		if study_tag == 'Q21_nTRPM':
			label = 'Nuclear TRPM'
		elif study_tag == 'Q21_TRPM':
			label = 'TRPM'
		elif study_tag == 'Q21_nM7CK':
			label = 'Nuclear M7CK'
		elif study_tag == 'Q21_H3S10':
			label = 'Phosphorylated H3S10'
		elif study_tag == 'Q21_Mg_IL8':
			label = 'IL8'
		elif study_tag == 'Q21_eq_trpm':
			pass
		elif study_tag == 'R05_mg_f_n':
			label = 'Free Mg$^{2+}$ ions'
		elif study_tag == 'R05_mg_n':
			label = 'Total Mg$^{2+}$ ions'
		elif study_tag == 'Q21_Mg':
			label = 'Free Mg$^{2+}$ ions'
		elif study_tag == 'eq_IL8':
			label = 'IL8 equalibrium'
		elif study_tag == 'eq_mg_f' or study_tag == 'Q21_eq' or study_tag == 'eq_mg':
			label = ''
		return label
	@staticmethod
	def determine_xlabel(study_tag):
		label = 'Time (h)'
		if study_tag == 'Q21_nTRPM' or study_tag == 'Q21_TRPM' or study_tag == 'Q21_nM7CK' or study_tag == 'Q21_H3S10':
			label = 'Mg conc. (ng/ml)'
		elif study_tag == 'Q21_Mg_IL8':
			label = ''
		elif study_tag == 'M18':
			label = 'IL8 (pg/ml)'
		return label
		
	@staticmethod
	def bar_positions(sims,delta,D,plot_t):
		if plot_t == 'bar1':
			x_exp =[float(j) + delta for j in range(len(sims))]
			x_sim =[float(j) - delta for j in range(len(sims))]
			xs = [[x_exp,x_sim]]
		elif plot_t == 'bar2':
			IDs_n = len(sims)
			d = D/(IDs_n+1) # the length allocation to a pair of exp-sim
			xs = [] # the location of bars sorted by Mg count
			for i in range(IDs_n):
				x_exp =[(float(j)-D/2) + d*(i+1) - delta for j in range(len(sims[0]))]
				x_sim =[(float(j)-D/2) + d*(i+1) + delta for j in range(len(sims[0]))]
				xs.append([x_exp,x_sim])
		return xs
	@staticmethod
	def determine_graph_size(study_tag):
		graph_size = [3,3]
		if study_tag == 'Q21_nTRPM' or study_tag == 'Q21_TRPM' or study_tag == 'Q21_nM7CK' or study_tag == 'Q21_H3S10':
			graph_size = [3,3]
		elif study_tag == 'Q21_Mg_IL8':
			graph_size = [4,3]
		return graph_size
		
	@staticmethod
	def determine_ylabel(study_tag):
		label = 'Relative Concentration \n (To initial quantity)'
		if study_tag == 'eq_mg' or study_tag == 'eq_mg_f':
			label = 'Concentration (mM)'
		elif study_tag == 'Q21_nTRPM' or study_tag == 'Q21_TRPM' or study_tag == 'Q21_nM7CK' or study_tag == 'Q21_H3S10':
			label = 'Relative intensity \n (To control)'
		elif study_tag == 'Q21_Mg_IL8':
			label = 'Concentration (pg/ml)'
		return label
	@staticmethod
	def determine_xlim(study_tag):
		if study_tag == '':
			pass
		else:
			raise ValueError('Define')
		return lim
	@staticmethod
	def determine_ylim(study_tag):
		if study_tag == 'eq_mg':
			lim = [18,19]
		elif study_tag == 'eq_mg_f':
			lim = [0,1]
		elif study_tag == 'eq_IL8':
			lim = [0.5,1.2]
		else:
			raise ValueError('Define')
		return lim
	@staticmethod
	def determine_xticks(ax,study_tag):
		ticks = ax.get_xticks()
		if study_tag == 'Q21_nTRPM' or study_tag == 'Q21_TRPM' or study_tag == 'Q21_nM7CK' or study_tag == 'Q21_H3S10':
			adj_ticks = [0,1,2]
			adj_labels = ['0.08','0.8','8']
		elif study_tag == 'Q21_Mg_IL8':
			adj_ticks = [0,1]
			adj_labels = ['6h','72h']
		elif study_tag == 'Q21_Mg':
			adj_ticks = np.array([0,1,2,3])*60/t2m	
			adj_labels = [int(i*t2m/60) for i in adj_ticks]
		elif study_tag == 'M18':
			adj_ticks = [0,1]
			adj_labels = ['Ctr','0.01']
		else:
			adj_ticks = ticks[1:-1]
			adj_labels = [int(i*t2m/60) for i in adj_ticks]
			
		return adj_ticks,adj_labels

	@staticmethod
	def determine_yticks(study_tag):
		if study_tag == 'eq_mg_f':
			y_ticks_ad = [0,0.5,1]
			y_ticks_label = y_ticks_ad
		elif study_tag == 'eq_mg':
			y_ticks_ad = [18+0.5*i for i in range(0,3)]
			y_ticks_label = y_ticks_ad
		else:
			raise ValueError('define')
			
		return y_ticks_ad,y_ticks_label
	@staticmethod
	def determine_legend(ax,study_tag):
		position = (1,1.1)
		if study_tag == 'R05_mg_f_n':
			position = (1.1,.7)
		elif study_tag == 'eq_mg' or study_tag == 'eq_mg_f':
			position = (1,1.1)
		elif study_tag == 'R05_mg_n':
			position = (1.1,1.1)
		elif study_tag == 'Q21_Mg':
			position = (1.1,.8)
		elif study_tag == 'Q21_eq':
			position = (1.14,.75)
		elif study_tag == 'Q21_nTRPM' or study_tag == 'Q21_TRPM' or study_tag == 'Q21_nM7CK' or study_tag == 'Q21_H3S10':
			position = (0.7,1.1)
		elif study_tag == 'Q21_Mg_IL8':
			position = (1.14,.75)
		elif study_tag == 'eq_IL8':
			position = (1.1,.7)
		ncol = 1
		return position,ncol

	
	

def run_plot_bar(ax,model,params,study_tag,target,study,plot_t,IDs=[]): 
	sims = model.simulate_study(study_tag=study_tag,params= params,study=study)

	sims,exps = process_data.sort(sims=sims,exps=study,target = target,plot_t=plot_t)
	sims,exps = process_data.normalize(study_tag = study_tag, sims= sims, exps=exps)
	specs = Specs(study_tag)
	xs = specs.bar_positions(sims,specs.delta,specs.D,plot_t=plot_t)
	for i in range(len(xs)):
		if plot_t == 'bar1':
			labels = ['S','E']
			sims_ID = sims
			exps_ID = exps
		elif plot_t == 'bar2':
			labels = ['S-%s'%IDs[i],'E-%s'%IDs[i]]
			sims_ID = sims[i]
			exps_ID = exps[i]
		plots.tools.plot_bar(ax=ax,specs=specs,x_exp=xs[i][0],x_sim=xs[i][1],sims=sims_ID,exps=exps_ID,labels=labels, plot_i = i)
	
	plots.tools.ax_postprocess(ax=ax,study_tag=study_tag,specs=specs,sims=sims)

def run_plot_line(ax,study_tag,model_sbml,params,target,ID,study):
	duration = study['experiment_period']
	inputs = study[ID]['inputs']
	specs = Specs(study_tag)

	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	sims = run_model(model_sbml,params = params_copy,target_keys=['TIME',target],duration=duration)
	
	x = sims['time']
	ax.plot(x,sims[target],color='black',linewidth=specs.line_width, label = 'S')
	obs_xx = study['measurement_scheme'][target]
	obs_yy = study[ID]['expectations'][target]['mean']
	ax.plot(obs_xx,obs_yy,linestyle='--',color='r', linewidth=specs.line_width, label = 'E')
	ax.legend()
	plots.tools.ax_postprocess(ax,study_tag=study_tag,sims=sims,specs=specs)

def run_plot_line_multi_target(ax,study_tag,model_sbml,params,targets,ID,study):
	duration = study['experiment_period']
	inputs = study[ID]['inputs']
	specs = Specs(study_tag)

	params_copy = copy.deepcopy(params)
	for key,value in inputs.items():
		params_copy[key] = value
	sims = run_model(model_sbml,params = params_copy,target_keys=['TIME']+targets,duration=duration)

	x = sims['time']

	obs_yy = [1 for i in range(len(x))]
	ax.plot(x,obs_yy,linestyle='--',color='r', linewidth=specs.line_width, label = 'Expectation')

	jj=0
	for target in targets:
		ax.plot(x,sims[target],color=specs.colors[jj],linewidth=specs.line_width, label = 'S: {}'.format(labels[target]))
		jj+=1
	ax.legend()
	plots.tools.ax_postprocess(ax,study_tag=study_tag,sims=sims,specs=specs)


	
def P5_plot(model_sbml,model_macrophage,params,observations): # plotting sim vs exp measurements for P4
	figsize = (17,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	obs = observations
	nn = 4
	jj = 1
	fig, (ax0, ax1, ax2,ax3) = plt.subplots(1, nn, gridspec_kw={'width_ratios': [1.25, 1,1,1]},figsize=figsize)


	study_tag,target = 'Q21_Mg_IL8','IL8'
	IDs = obs[study_tag]['IDs']
	run_plot_bar(ax=ax0,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar2',IDs=IDs)



	study_tag,targets = 'eq_IL8',['IL8_n','IL8R_n']
	ID = obs[study_tag]['IDs'][0]
	run_plot_line_multi_target(ax=ax1,study_tag=study_tag,targets=targets,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])

	study_tag,target = 'M18','nIFNGR'
	IDs = obs[study_tag]['IDs']
	run_plot_bar(ax=ax2,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)

	study_tag,target = 'M18','nIL4R'
	IDs = obs[study_tag]['IDs']
	run_plot_bar(ax=ax3,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)

	
	fig.tight_layout()

def plot_Q21_Mg_IL8(ax,model_sbml,model_macrophage,params,observations):
	study_tag = 'Q21_Mg_IL8'
	target = 'IL8'
	obs = observations[study_tag]
	
	sim_results = model_macrophage.simulate_study(study_tag = study_tag,params= params,study=observations[study_tag])

	IDs = observations[study_tag]['IDs']
	IDs = ['0.08 mM','8 mM']

	# run_plot_bar(ax=ax,model = model_macrophage,study_tag=study_tag,target=target,sims=sim_results,study=observations[study_tag],plot_t= 'bar2', IDs=IDs)


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
	run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
	jj+=1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'Q21_IkBa','IKB'
	IDs = obs[study_tag]['IDs']
	run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
	jj+=1
	
	fig.tight_layout()

def P1_3_eq_plot(model_sbml,params,observations): # equalibrium for P1 to P3
	figsize = (12,4)
	fig = plt.figure(figsize=figsize)
	obs = observations
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
	obs = observations
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	nn = 3
	jj = 1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'R05_mg_f_n','Mg_f_n'
	ID = obs[study_tag]['IDs'][0]
	run_plot_line(ax=ax,study_tag=study_tag,target=target,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
	jj+=1
	
	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'R05_mg_n','Mg_n'
	ID = obs[study_tag]['IDs'][0]
	run_plot_line(ax=ax,study_tag=study_tag,target=target,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
	
	jj+=1

	fig.tight_layout()
def P1_3_plot(model_sbml,model_macrophage,params,observations): # plotting sim vs exp measurements for P1 to p3
	fig_tag = 'P1_3'
	figsize = (20,4)
	fig = plt.figure(figsize=figsize)
	fig.canvas.draw()
	obs = observations
	nn = 5
	jj = 1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target,ID = 'Q21_Mg','Mg_f_n','Mg_8'
	run_plot_line(ax=ax,study_tag=study_tag,target=target,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
	jj+=1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'Q21_nTRPM','nTRPM'
	IDs = obs[study_tag]['IDs']
	run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
	jj+=1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'Q21_TRPM','TRPM'
	IDs = obs[study_tag]['IDs']
	run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
	jj+=1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'Q21_nM7CK','nM7CK'
	IDs = obs[study_tag]['IDs']
	run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
	jj+=1

	ax = fig.add_subplot(1,nn,jj)
	study_tag,target = 'Q21_H3S10','pH3S10'
	IDs = obs[study_tag]['IDs']
	run_plot_bar(ax=ax,model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
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
	results = run_model(model_sbml,params = params_copy,target_keys=targets+extra_targets,duration=duration)
	ii = 0
	study1,study2 = study+'_f',study
	specs1,specs2 = Specs(study1),Specs(study2)
	for target in targets: # one subgraph for each target
		if target == 'Mg_ATP_n':
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
	plots.tools.ax_postprocess(ax=axs[0],study_tag=study1,specs=specs1,sims=results)
	plots.tools.ax_postprocess(ax=axs[1],study_tag=study2,specs=specs2,sims=results)

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
	results = run_model(model_sbml,params = params_copy,target_keys=targets,duration=duration)

	ii = 0

	specs = Specs(study_tag)

	ax.plot(results['time'],[1 for i in range(len(results['time']))],linestyle='--',color='r', label = 'Expectation')
	for target in targets: # one subgraph for each target
		x = results['time']
		ax.plot(x,results[target],linewidth=specs.line_width, color=specs.colors[ii],label = '%s'%labels[target])
		ii+=1
	
	plots.tools.ax_postprocess(ax=ax,study_tag=study_tag,specs=specs,sims=results)

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