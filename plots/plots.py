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




class tools:
	
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
	def ax_postprocess(ax,study_tag,sims,specs):
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
		ax.set_title(specs.determine_title(study_tag=study_tag),fontdict ={'family':font_type,'size':specs.title_font_size})
		
		
