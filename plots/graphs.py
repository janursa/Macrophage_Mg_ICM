import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
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




class graph_bar:
	"""
	Plots the results of a study by allocting a figure for each target and a bar for each ID
	"""
	def __init__(self,study,observations,errors,destination=''):
		self.measurement_scheme = observations[study]['measurement_scheme']
		self.study = study
		self.observations = observations
		self.colors = ['lime' , 'violet', 'yellowgreen', 'peru', 'skyblue']
		self.errors = errors
		self.destination = destination

		if study == 'Q21_nTRPM' or study == 'Q21_TRPM' or study == 'Q21_nM7CK' or study == 'Q21_H3S10':
			self.graph_size = [3,3]
		self.bar_width = .2
		self.error_bar_width = 5
		self.colors = ['lime' , 'violet', 'yellowgreen', 'peru', 'skyblue']
		self.legend_font_size = 18
		self.tick_font_size = 21
		self.title_font_size = 21 
		self.delta = .12
		self.legend_location=(1.5,1.2)
		self.R2_font_size  = 18

		# else:
		# 	raise ValueError('input not defined')
	@staticmethod
	def determine_yrange(study,target):
		if study == 'Qiao_2021_IL8_IL1b':
			yrange_value = [0,60]
		return yrange_value
	@staticmethod
	def determine_xrange(study,target):
		
		raise ValueError('define')

	@staticmethod
	def p_values_positions(study,target,xx):
		if study == 'Qiao_2021_IL8_IL1b':
			# print(study,xx)
			Xs = xx[1]+.08
			Ys = 50
		return Xs,Ys

	def draw_R2(self,ax,study,target):
		value = round(1-self.errors[target],2)
		if study == 'Qiao_2021_Mg':
			xy = (55,22)
		
		else:
			raise ValueError('draw_R2 is not defined for {}'.format(study))

		ax.annotate('$\\mathdefault{\\overline{R^2} = %.2f } $'%value, xy=xy,fontsize = self.R2_font_size)

	@staticmethod
	def draw_p_values(ax,study,target,Xs,Ys):
		if study == 'Qiao_2021_IL8_IL1b':
			barplot_annotate_stars(significance=0.001,
					center=Xs,height=Ys)
		
	@staticmethod
	def determine_ylabel(study,target):
		label = ''
		if study == 'Q21_nTRPM' or study == 'Q21_TRPM' or study == 'Q21_nM7CK' or study == 'Q21_H3S10':

			label = 'Relative intensity \n (To control)'
		
		return label
	@staticmethod
	def determine_title(study,target):
		label = study
		if study == 'Q21_nTRPM':
			label = 'Nuclear TRPM'
		elif study == 'Q21_TRPM':
			label = 'TRPM'
		elif study == 'Q21_nM7CK':
			label = 'Nuclear M7CK'
		elif study == 'Q21_H3S10':
			label = 'Phosphorylated H3S10'

		return label
	@staticmethod
	def determine_xlabel(study,target):
		label = ''
		if study == 'Q21_nTRPM' or study == 'Q21_TRPM' or study == 'Q21_nM7CK' or study == 'Q21_H3S10':
			label = 'Mg conc. (ng/ml)'
		return label
	@staticmethod
	def normalize_exp_values(study,exps):
		exps_n = []
		if study=='Q21_IkBa':
			ctr = exps[0]
		else:
			return exps
		for item in exps:
			exps_n.append(item/ctr)
		return exps_n
	@staticmethod
	def normalize_sim_values(study,sims):
		sims_n = []
		if study=='Q21_nTRPM' or study=='Q21_TRPM' or study=='Q21_nM7CK' or study=='Q21_H3S10':
			ctr = sims[1]
		elif study == 'S12_IkBa_mg':
			ctr = sims[0]
		elif study == 'Q21_IkBa':
			ctr = sims[0]
		else:
			raise ValueError('define')
		for item in sims:
			sims_n.append(item/ctr)
		return sims_n
	
	def plot(self,ax,simulation_results,processed_detailed_errors={}):
		##/ sort out based on the measurement_scheme
		IDs = list(simulation_results.keys())
		x_labels = self.adjust_x_label(IDs)
		exp_target_results,sim_target_results = self.sort(simulation_results)
		_,processed_detailed_errors_sorted = self.sort(processed_detailed_errors)
		x_exp,x_sim,base_x = self.bar_positions(simulation_results)

		##/ plot for each target
		target_n = len(self.measurement_scheme)

		target = list(self.measurement_scheme.keys())[0]

		# for target,ii in zip(self.measurement_scheme.keys(),range(target_n)):
		sim_values = [item[0] for item in sim_target_results[target]]
		sim_values = graph_bar.normalize_sim_values(study = self.study, sims= sim_values)
		processed_detailed_errors_sorted_values = [item[0] for item in processed_detailed_errors_sorted[target]]
		# ax.bar(x=x_sim,height=sim_values,width = self.bar_width, label = "S", 
		# 		facecolor = self.colors[0],
		# 		 edgecolor="black", yerr =  processed_detailed_errors_sorted_values,
		# 		 error_kw = dict(capsize= self.error_bar_width))
		ax.bar(x=x_sim,height=sim_values,width = self.bar_width, label = "S", 
				facecolor = self.colors[0],
				 edgecolor="black", yerr =  0,
				 error_kw = dict(capsize= self.error_bar_width))
		exp_values = [exp_target_results[target][i]['mean'] for i in range(len(exp_target_results[target]))]
		exp_std = [exp_target_results[target][i]['std'] for i in range(len(exp_target_results[target]))]
		exp_values = [item[0] for item in exp_values]
		exp_values = graph_bar.normalize_exp_values(study = self.study, exps= exp_values)
		exp_std = [item[0] for item in exp_std]

		ax.bar(x=x_exp,height=exp_values,width = self.bar_width, label = 'E', 
				facecolor = self.colors[1],hatch=r'\\\\',
				 edgecolor="black", yerr =  exp_std,
				 error_kw = dict(capsize= self.error_bar_width))
		# if ii == 0: # legend only for the first target
		# 	if self.study == 'Qiao_2021_IL8_IL1b' or self.study == 'Qiao_2021_IL1b':
		# 		pass
		# 	else:
		# 		if self.study == 'Chen_2018' or self.study == 'Valles_2020_TNFa' or self.study == 'Valles_2020_IL10':
		# 			ncol = 1
		# 		else:
		# 			ncol = 2
		graph_line.legend(ax=ax,study=self.study)
		x_ticks = [(i+j)/2 for i,j in zip(x_sim,x_exp)]
		x_ticks_adj,x_labels_adj = self.add_adjustements(ax=ax,study=self.study,target=target,base_x=base_x,x_ticks=x_ticks,x_labels=x_labels)
		self.finalize_and_save(ax=ax,target=target,exp_xx=x_ticks_adj,x_ticks=x_ticks_adj,x_labels=x_labels_adj)
		# return fig
	def add_adjustements(self,ax,study,target,base_x,x_ticks,x_labels):
		x_ticks_adj = copy.deepcopy(x_ticks)
		x_labels_adj = copy.deepcopy(x_labels)
		if study == 'Chen_2018':
			x_ticks_adj.insert(0,base_x)
			x_labels_adj.insert(0,'ctr')
			if target == 'ALP':
				ax.bar(x=base_x,height=1,width = self.bar_width, label = 'Exp', 
					facecolor = self.colors[1],hatch=r'\\\\',
					 edgecolor="black", yerr =  0,
					 error_kw = dict(capsize= self.error_bar_width))
			elif target == 'ARS':
				ax.bar(x=base_x,height=0.1,width = self.bar_width, label = 'Exp', 
					facecolor = self.colors[1],hatch=r'\\\\',
					 edgecolor="black", yerr =  0,
					 error_kw = dict(capsize= self.error_bar_width))
		return x_ticks_adj,x_labels_adj
	def finalize_and_save(self,ax,target,x_ticks,x_labels,exp_xx):
		try:
			ax.set_ylim(self.determine_yrange(target = target, study = self.study))
		except:
			pass
		try:
			ax.set_xlim(self.determine_xrange(target = target, study = self.study))
		except:
			pass
		
		if self.study == 'Ber_2016':
			ax.set_xticks(ticks = [])
			ax.set_xticklabels([])
		else:
			ax.set_xticks(ticks = x_ticks)
			ax.set_xticklabels(x_labels)

		for label in (ax.get_xticklabels() + ax.get_yticklabels()):
			label.set_fontname(font_type)
			label.set_fontsize(self.tick_font_size)
		# unit_value = graph_bar.unit(study = self.study,target=target)
		ax.set_ylabel(graph_bar.determine_ylabel(study=self.study,target = target),fontdict ={'family':font_type,'size':self.title_font_size})
		ax.set_xlabel(graph_bar.determine_xlabel(study=self.study,target = target),fontdict ={'family':font_type,'size':self.title_font_size})

		ax.set_title(graph_bar.determine_title(study=self.study,target = target),fontdict ={'family':font_type,'size':self.title_font_size})
		
		# Xs,Ys = graph_bar.p_values_positions(study=self.study,target=target,xx = exp_xx)	
		# graph_bar.draw_p_values(ax,study=self.study,target=target,Xs=Xs,Ys=Ys)


		# self.draw_R2(ax,study=self.study,target=target)
		
		if self.study == 'Qiao_2021_Mg':
			ax.get_xaxis().set_major_formatter(
				matplotlib.ticker.FuncFormatter(lambda x, p: format(int(int(x)/24), ',')))
		plt.savefig(os.path.join(self.destination, self.study+'.svg'),bbox_inches='tight')

	def sort(self,sim_results):
		exp_target_results = {}
		for target in self.measurement_scheme:
			exp_target_results[target] = []
		for target in self.measurement_scheme:
			for ID in self.observations[self.study]['IDs']:
				ID_observations = self.observations[self.study][ID]['expectations']
				exp_target_results[target].append(ID_observations[target])
		sim_target_results = {}
		for target in self.measurement_scheme:
			sim_target_results[target] = []
		for target in self.measurement_scheme:
			for ID,ID_result in sim_results.items():
				sim_target_results[target].append(ID_result[target])
		return exp_target_results,sim_target_results
	def bar_positions(self,sim_results):
		for i in range(len(self.measurement_scheme)):
			x_exp =[float(j) + self.delta for j in range(len(sim_results.keys()))]
			x_sim =[float(j) - self.delta for j in range(len(sim_results.keys()))]
		base = None

		return x_exp,x_sim,base

	def adjust_x_label(self,labels):
		adj_labels = []
		for label in labels:
			if label == 'Mg_.08':
				adj_labels.append('0.08')
			elif label == 'Mg_.8':
				adj_labels.append('0.8')
			elif label == 'Mg_8':
				adj_labels.append('8')
			elif label == 'ctr':
				adj_labels.append('0.5')
			elif label == 'Mg_2dot5':
				adj_labels.append('2.5')
			else:
				raise ValueError('not defined')
		return adj_labels
class graph_bar_2(graph_bar):
	"""
	This one has more days (time points). Previous one had only the final time points
	"""
	def __init__(self,**args):
		super().__init__(**args)
	def bar_positions(self,study,IDs,checkpoints):
		IDs_n = len(IDs)
		d = self.D/(IDs_n+1) # the length allocation to a pair of exp-sim
		xs = [] # the location of bars sorted by Mg count
		for i in range(IDs_n):
			x_exp =[(float(j)-self.D/2) + d*(i+1) - self.delta for j in checkpoints]
			x_sim =[(float(j)-self.D/2) + d*(i+1) + self.delta for j in checkpoints]
			xs.append([x_exp,x_sim])
		return xs
