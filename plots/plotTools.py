import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import os
import sys
import copy
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
from data.observations import observations,t2m
from models.models import Macrophage
plt.rcParams["font.family"] = "serif"
plt.style.use('seaborn-deep')
# plt.style.use('Agg')

plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]


font_type = 'Arial'

labels = {
		'Mg_f':'free Mg$^{2+}$',
		'nMg_f':'free Mg$^{2+}$',
		'Mg':'total Mg$^{2+}$',
		'nMg':'total Mg$^{2+}$',

		'TRPM': 'TRPM',
		'nTRPM': 'Nuclear TRPM',
		'nM7CK': 'Nuclear M7CK',
		'pH3S10': 'Phos H3S10',

		'nIL8': 'norm IL8',
		'nIL8R': 'norm IL8R'
	}

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
	def determine_title(study_tag,target=''):
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
		elif study_tag == 'M18':
			label = target
		elif study_tag == 'M05_IT':
			label = target
		elif study_tag == 'M05_NFKB':
			label = target
		return label
	@staticmethod
	def determine_xlabel(study_tag):
		label = 'Time (h)'
		if study_tag == 'Q21_nTRPM' or study_tag == 'Q21_TRPM' or study_tag == 'Q21_nM7CK' or study_tag == 'Q21_H3S10':
			label = 'Mg conc. (ng/ml)'
		elif study_tag == 'Q21_Mg_IL8':
			label = ''
		elif study_tag == 'M18' or study_tag == 'M05_IT' or study_tag == 'M05_NFKB':
			label = 'IL8 (ng/ml)'
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
		raise ValueError('this should be depricated')
		return graph_size
		
	@staticmethod
	def determine_ylabel(study_tag):
		label = 'Relative quantity \n (To control)'
		if study_tag == 'eq_mg' or study_tag == 'eq_mg_f':
			label = 'Concentration (mM)'
		elif study_tag == 'Q21_Mg_IL8':
			label = 'Concentration (pg/ml)'
		elif study_tag =='eq_IL8':
			label = 'Relative quantity \n (To initial value)'

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
			adj_ticks = [0,1,2,3,4]
			adj_labels = ['Ctr','0.01','0.1','1','10']
		elif study_tag == 'M05_IT':
			adj_ticks = [0,1]
			adj_labels = ['Ctr','100']
		elif study_tag == 'M05_NFKB':
			adj_ticks = [0,1,2,3,4,5]
			adj_labels = ['Ctr','1','10','10','100','1000']
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
		elif study_tag == 'M18':
			position = (.7,1)
		elif study_tag == 'M05_IT':
			position = (.5,1)
		elif study_tag == 'M05_NFKB':
			position = (.4,1)
		ncol = 1
		return position,ncol



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



class plotTools:
	
	@staticmethod
	def plot_bar(ax,specs,x_exp,x_sim,sims,exps,labels,plot_i):
		ax.bar(x=x_sim,height=sims,width = specs.bar_width, label = labels[0], 
				facecolor = specs.colors[plot_i*2],
				 edgecolor="black", yerr =  0,
				 error_kw = dict(capsize= specs.error_bar_width))

		ax.bar(x=x_exp,height=exps,width = specs.bar_width, label = labels[1], 
				facecolor = specs.colors[plot_i*2+1],hatch=r'\\\\',
				 edgecolor="black", yerr =  0,
				 error_kw = dict(capsize= specs.error_bar_width))
	# def plot_1(self,ax,x_exp,x_sim,sims,exps):
	# 	ax.bar(x=x_sim,height=sims,width = self.specs.bar_width, label = "S", 
	# 			facecolor = self.specs.colors[0],
	# 			 edgecolor="black", yerr =  0,
	# 			 error_kw = dict(capsize= self.specs.error_bar_width))

	# 	ax.bar(x=x_exp,height=exps,width = self.specs.bar_width, label = 'E', 
	# 			facecolor = self.specs.colors[1],hatch=r'\\\\',
	# 			 edgecolor="black", yerr =  0,
	# 			 error_kw = dict(capsize= self.specs.error_bar_width))
	
	@staticmethod
	def ax_postprocess(ax,study_tag,sims,specs,target=''):
		x_ticks = [i for i in range(len(sims))]
		try:
			y_ticks_ad,y_ticks_label = specs.determine_yticks(study_tag=study_tag)
			ax.set_yticks(ticks = y_ticks_ad)
			ax.set_yticklabels(y_ticks_label)
		except ValueError:
			pass

		
		x_ticks_ad,x_ticks_labels_adj = specs.determine_xticks(ax=ax,study_tag=study_tag)
		ax.set_xticks(ticks = x_ticks_ad)
		ax.set_xticklabels(x_ticks_labels_adj)
		
		

		try:
			ax.set_xlim(specs.determine_xlim(study_tag))
		except ValueError:
			pass

		try:
			ax.set_ylim(specs.determine_ylim(study_tag))
		except ValueError:
			pass

		position,ncol =specs.determine_legend(ax=ax,study_tag=study_tag)
		ax.legend(bbox_to_anchor=position,loc = 'upper right', borderaxespad=2,prop={'family':font_type,'size':specs.legend_font_size},ncol=ncol)

		for label in (ax.get_xticklabels() + ax.get_yticklabels()):
			label.set_fontname(font_type)
			label.set_fontsize(specs.tick_font_size)
		ax.set_ylabel(specs.determine_ylabel(study_tag=study_tag),fontdict ={'family':font_type,'size':specs.title_font_size})
		ax.set_xlabel(specs.determine_xlabel(study_tag=study_tag),fontdict ={'family':font_type,'size':specs.title_font_size})
		ax.set_title(specs.determine_title(study_tag=study_tag,target=target),fontdict ={'family':font_type,'size':specs.title_font_size})
		

	@staticmethod
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
			plotTools.plot_bar(ax=ax,specs=specs,x_exp=xs[i][0],x_sim=xs[i][1],sims=sims_ID,exps=exps_ID,labels=labels, plot_i = i)
		
		plotTools.ax_postprocess(ax=ax,study_tag=study_tag,specs=specs,sims=sims,target=target)

	@staticmethod
	def run_plot_line(ax,study_tag,model_sbml,params,target,ID,study):
		duration = study['experiment_period']
		inputs = study[ID]['inputs']
		specs = Specs(study_tag)

		params_copy = copy.deepcopy(params)
		for key,value in inputs.items():
			params_copy[key] = value
		sims = Macrophage.run_sbml_model(model_sbml=model_sbml,params = params_copy,selections=['TIME',target],duration=duration)
		
		x = sims['time']
		ax.plot(x,sims[target],color='black',linewidth=specs.line_width, label = 'S')
		obs_xx = study['measurement_scheme'][target]
		obs_yy = study[ID]['expectations'][target]['mean']
		ax.plot(obs_xx,obs_yy,linestyle='--',color='r', linewidth=specs.line_width, label = 'E')
		ax.legend()
		plotTools.ax_postprocess(ax,study_tag=study_tag,sims=sims,specs=specs)

	@staticmethod
	def run_plot_line_multi_target(ax,study_tag,model_sbml,params,targets,ID,study):
		duration = study['experiment_period']
		inputs = study[ID]['inputs']
		specs = Specs(study_tag)

		params_copy = copy.deepcopy(params)
		for key,value in inputs.items():
			params_copy[key] = value
		
		sims = Macrophage.run_sbml_model(model_sbml = model_sbml,params = params_copy,selections=['TIME']+targets,duration=duration)

		x = sims['time']

		obs_yy = [1 for i in range(len(x))]
		ax.plot(x,obs_yy,linestyle='--',color='r', linewidth=specs.line_width, label = 'Expectation')

		jj=0
		for target in targets:
			ax.plot(x,sims[target],color=specs.colors[jj],linewidth=specs.line_width, label = 'S: {}'.format(labels[target]))
			jj+=1
		ax.legend()
		plotTools.ax_postprocess(ax,study_tag=study_tag,sims=sims,specs=specs)


