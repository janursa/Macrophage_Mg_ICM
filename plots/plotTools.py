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
from tools import tools

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Arial"] + plt.rcParams["font.serif"]

font_type = 'Arial'

labels = {
        'Mg_f':r'free $\mathdefault{Mg^{2+}}$',
        'nMg_f':r'free $\mathdefault{Mg^{2+}}$',
        'Mg':r'total $\mathdefault{Mg^{2+}}$',
        'nMg':r'total $\mathdefault{Mg^{2+}}$',
        'nMg_ATP': r'$\mathdefault{Mg.ATP}$',
        'nATP': r'$\mathdefault{ATP}$',
        
        'nIKB': 'Cytosolic IkBa',
        'nNFKB_n': 'Nuclear NFkB',
        'nTRPM': 'Cytosolic TRPM',
        'nTRPM_n': 'Nuclear TRPM',
        'nM7CK_n': 'Nuclear M7CK',
        'npH3S10': 'Activated H3S10',

        'IL8': 'IL8',
        'nIL8_m': 'IL8 mRNA',
        'IL8_R': 'IL8/IL8R complex',
        'IL8R': 'IL8R',
        'nIL8': 'IL8',
        'nIL10': 'IL10',
        'nTNFa': 'TNFa',
        'nIL1b': 'IL1b',
        'nIFNGR': r'IFN-Y receptor',
        'IFNGR':r'IFN-Y receptor',
        'nIL4R': r'IL4 receptor',
        'nIL6':'IL6',
        'F_il8_irak':'F_il8_irak',

        'nIRAK4':'Cytosolic IRAK',
        'naTRAF6':'Activated cytosolic TRAF6',
        'npSTAT3': 'Phos STAT3',
        'npIL6_R': 'Phos IL6/R',
        'F_stat3_a': 'F_stat3_a',
        'F_pi3k_a':'F_pi3k_a',

        'F_h3s10_ikb': 'H3S10 effect on IkBa',
        'F_p3s10_il8_p': 'H3S10 effect on IL8'
    }

class Specs:
    def __init__(self,study_tag):
        self.line_width = 3
        self.bar_width = .2
        self.error_bar_width = 5
        self.colors = ['lime' , 'violet', 'yellowgreen', 'peru', 'skyblue','b']
        self.legend_font_size = 15
        self.tick_font_size = 21
        self.title_font_size = 21
        self.delta = .12
        self.R2_font_size  = 18
        self.D = 1.1 # the length in which all Mg dosages are plotted in a certain time point
        self.pvalue_size = 21
        if  study_tag == 'M05_NFKBn' or study_tag == 'M18':
            self.bar_width = .3
            self.delta = .17
            self.D = 1.3
        
        
    @staticmethod
    def determine_title(study_tag,target='',duration=None):
        prefix = study_tag.split('_')[0]
        if study_tag == 'eq_mg' or study_tag == 'Q21_eq' or study_tag == 'eq_IL8' or study_tag == 'eq_combined' or study_tag == 'eq_IL6':
            # label = '(A) Equalibrium \n (no stimulus)'
            label = ''
        elif study_tag == 'R05_nMg_f' or study_tag == 'Q21_Mg':
            label = prefix+': '+labels[target]
        
        elif study_tag == 'Q21_14d':
            label = prefix+': '+labels[target]+' (%s'%(int(duration/60/24))+'d)'
        else:
            label = study_tag+'\n'+prefix+': '+labels[target]+' (%s'%(int(duration/60))+'h)'
        return label
    @staticmethod
    def determine_xlabel(study_tag):
        label = 'Time (h)'
        if study_tag == 'Q21_M1'  or study_tag == 'S12_IKBa_mg' or study_tag == 'S12_NFKBn_mg' or study_tag == 'Q21_IkBa_6h' or study_tag == 'Q21_IkBa_72h' or study_tag == 'Q21_NFKBn_72h' or study_tag == 'Q21_Mg_IL8' or study_tag == 'Q21_Mg_IL8' or study_tag == 'Q21_14d' or study_tag == 'Z19_IKB_NFKB' or study_tag == 'Z19_IL10' or study_tag == 'Q21_cytokines_72h' or study_tag == 'F18_cytokines' or study_tag == 'B20_NFKBn' or study_tag == 'B20_TNFa' :
            label = 'Mg2+ ions (mM)'
        elif  study_tag == 'Q21_IkBa':
            label = ''
        elif study_tag == 'M18' or study_tag == 'M05_IT' or study_tag == 'M05_NFKBn':
            label = 'IL8 (ng/ml)'
        elif study_tag == 'B17' or study_tag == 'F17' or study_tag == 'N03' or study_tag == 'F14':
            label = 'IL6 (ng/ml)'
        
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
                x_exp =[(float(j)-D/2) + d*(i+1) + delta for j in range(len(sims[0]))]
                x_sim =[(float(j)-D/2) + d*(i+1) - delta for j in range(len(sims[0]))]
                xs.append([x_exp,x_sim])
        return xs
        
    @staticmethod
    def determine_ylabel(study_tag):
        label = 'Relative quantity \n (To control)'
        if study_tag =='eq_IL8' or study_tag =='Q21_Mg' or study_tag == 'eq_mg' or study_tag == 'Q21_eq' or study_tag =='R05_nMg_f':
            label = 'Relative quantity \n (To initial value)'
        elif study_tag =='Q21_IkBa' or study_tag =='Q21_IkBa_6h' or study_tag =='Q21_IkBa_72h' :
            label = "% Input"

        return label
    @staticmethod
    def determine_xlim(study_tag,target=''):
        if study_tag == '':
            pass
        else:
            raise ValueError('Define')
        return lim
    @staticmethod
    def determine_ylim(study_tag,target,max_y):
        offset = .5
        if study_tag == 'eq_mg':
            lim = [0,1.35]
        elif study_tag == 'Q21_eq':
            lim = [-.1,1.5]
        elif study_tag == 'eq_IL6':
            if target == 'F_stat3_a':
                lim = [0,2]
            else:
                raise ValueError('Define')
        else:
            if max_y == None:
                raise ValueError('Define')
            lim = [0,max_y+offset]
        
        
        return lim
    @staticmethod
    def determine_xticks(ax,study_tag):
        ticks = ax.get_xticks()
        if study_tag == 'Q21_M1' or study_tag == 'Q21_14d' or study_tag == 'Q21_cytokines_72h' or study_tag == 'Q21_Mg_IL8':
            adj_ticks = [0,1]
            adj_labels = ['ctr','8']
        elif study_tag == 'S12_IKBa_mg' or  study_tag == 'S12_NFKBn_mg':
            adj_ticks = [0,1]
            adj_labels = ['ctr','2.5']
        elif study_tag == 'Q21_IkBa':
            adj_ticks = [0,1]
            adj_labels = ['6h','72h']
        elif study_tag == 'Q21_Mg':
            adj_ticks = np.array([0,1,2,3])*60/t2m
            adj_labels = [int(i*t2m/60) for i in adj_ticks]
        elif study_tag == 'M18':
            adj_ticks = [0,1,2,3,4]
            adj_labels = ['ctr','0.01','0.1','1','10']
        elif study_tag == 'M05_IT':
            adj_ticks = [0,1]
            adj_labels = ['ctr','100']
        elif study_tag == 'M05_NFKBn':
            adj_ticks = [0,1,2,3,4]
            adj_labels = ['ctr','1','10','10','100']
        elif study_tag == 'Q21_IkBa_6h' or study_tag == 'Q21_IkBa_72h':
            adj_ticks = [0,1]
            adj_labels = ['0.08','8']
        elif study_tag == 'Z19_IKB_NFKB' or study_tag == 'Z19_IL10':
            adj_ticks = [0,1,2]
            adj_labels = ['ctr','4','20']
        elif study_tag == 'F18_cytokines' or study_tag == 'B20_NFKBn' or study_tag == 'B20_TNFa':
            adj_ticks = [0,1]
            adj_labels = ['ctr','5']
        elif study_tag == 'S12_LPS':
            adj_ticks = [int(ii*60) for ii in [0,1,2,8]]
            adj_labels = [0,1,2,8]
        elif study_tag == 'B17':
            adj_ticks = [0,1]
            adj_labels = ['ctr','60']
        elif study_tag == 'F17':
            adj_ticks = [0,1,2,3]
            adj_labels = ['ctr','50','100','200']
        elif study_tag == 'N03':
            adj_ticks = [0,1,2,3]
            adj_labels = ['ctr','0.1','10','100']
        elif study_tag == 'F14':
            adj_ticks = [0,1]
            adj_labels = ['ctr','10']

        else:
            # pass
            adj_ticks = ticks[1:-1]
            adj_labels = [int(i*t2m/60) for i in adj_ticks]
            
        return adj_ticks,adj_labels

    @staticmethod
    def determine_yticks(study_tag):
#        if study_tag == '':
##            y_ticks_ad = [18+0.5*i for i in range(0,3)]
##            y_ticks_label = y_ticks_ad
#
#        else:
        raise ValueError('define')
            
        return y_ticks_ad,y_ticks_label
    @staticmethod
    def determine_legend(ax,study_tag):
        position = (.6,1.1)
        if study_tag == 'R05_nMg_f':
            position = (1.1,1)
        elif study_tag == 'R05_mg_n':
            position = (1.1,1.1)
        elif study_tag == 'Q21_Mg':
            position = (1.1,1)
        elif study_tag == 'Q21_eq':
            position = (1.1,.65)
        elif study_tag == 'eq_mg':
            position = (1.05,.7)
        elif study_tag == 'Q21_M1' or study_tag == 'Q21_cytokines_72h':
            position = (0.7,1.1)
        elif study_tag == 'eq_IL8':
            position = (1.1,.7)
        elif study_tag == 'M18':
            position = (.3,1.1)
        elif study_tag == 'M05_IT':
            position = (.5,1)
        elif study_tag == 'M05_NFKBn':
            position = (.4,1)
        ncol = 1
        return position,ncol



class process_data:
    def sort(sims,exps,target,plot_t):
        sims_sorted = []
        exps_sorted = []
        exps_std_sorted = []
        pvalues = []
        for ID in exps['IDs']:
            ID_exps = exps[ID]['expectations'][target]
            ID_sims = sims[ID][target]
            exps_sorted.append(ID_exps['mean'])
            exps_std_sorted.append(ID_exps['std'])
            pvalues.append(ID_exps['pvalue'])
            sims_sorted.append(ID_sims)
        if plot_t == 'bar1':
            sims_sorted = [item[0] for item in sims_sorted]
            exps_sorted = [item[0] for item in exps_sorted]
            exps_std_sorted = [item[0] for item in exps_std_sorted]
            pvalues = [item[0] for item in pvalues]
        
        return sims_sorted,exps_sorted,exps_std_sorted,pvalues





class plotTools:
    @staticmethod
    def barplot_annotate_stars(ax,significance, center, height, specs):

        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        if significance == None:
            return
        if type(significance) == str:
            text = significance
        else:
            text = ''
            p = .05
    
            while significance < p:
                text += '*'
                p /= 10.

        kwargs = dict(ha='center', va='center')
        kwargs['fontsize'] = specs.pvalue_size
        mid = (center, height)

        ax.text(*mid, text, **kwargs)
        
        
    @staticmethod
    def draw_p_values(ax,pvalues,x_exp,exps,exps_std,specs):
        offset = .1
        for i in range(len(pvalues)):
            x = x_exp[i]
            y = exps[i]+exps_std[i]+offset
            pvalue = pvalues[i]
            plotTools.barplot_annotate_stars(ax=ax,significance = pvalue,
                        center = x, height = y, specs = specs)
    @staticmethod
    def plot_bar(ax,specs,x_exp,x_sim,sims,exps,exps_std,labels,plot_i):
        ax.bar(x=x_sim,height=sims,width = specs.bar_width, label = labels[0],
                facecolor = specs.colors[plot_i*2],
                 edgecolor="black", yerr =  0,
                 error_kw = dict(capsize= specs.error_bar_width))

        ax.bar(x=x_exp,height=exps,width = specs.bar_width, label = labels[1],
                facecolor = specs.colors[plot_i*2+1],hatch=r'\\\\',
                 edgecolor="black", yerr =  exps_std,
                 error_kw = dict(capsize= specs.error_bar_width))
        
    # def plot_1(self,ax,x_exp,x_sim,sims,exps):
    #   ax.bar(x=x_sim,height=sims,width = self.specs.bar_width, label = "S",
    #           facecolor = self.specs.colors[0],
    #            edgecolor="black", yerr =  0,
    #            error_kw = dict(capsize= self.specs.error_bar_width))

    #   ax.bar(x=x_exp,height=exps,width = self.specs.bar_width, label = 'E',
    #           facecolor = self.specs.colors[1],hatch=r'\\\\',
    #            edgecolor="black", yerr =  0,
    #            error_kw = dict(capsize= self.specs.error_bar_width))
    
    @staticmethod
    def ax_postprocess(ax,study_tag,sims,specs,duration,target='',max_y=None):
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
            ax.set_xlim(specs.determine_xlim(study_tag,target))
        except ValueError:
            pass

        try:
            ax.set_ylim(specs.determine_ylim(study_tag,target,max_y))
        except ValueError:
            pass

        position,ncol =specs.determine_legend(ax=ax,study_tag=study_tag)
        ax.legend(bbox_to_anchor=position,loc = 'upper right', borderaxespad=2,prop={'family':font_type,'size':specs.legend_font_size},ncol=ncol)

        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontname(font_type)
            label.set_fontsize(specs.tick_font_size)
        ax.set_ylabel(specs.determine_ylabel(study_tag=study_tag),fontdict ={'family':font_type,'size':specs.title_font_size})
        ax.set_xlabel(specs.determine_xlabel(study_tag=study_tag),fontdict ={'family':font_type,'size':specs.title_font_size})
        ax.set_title(specs.determine_title(study_tag=study_tag,duration=duration,target=target),fontdict ={'family':font_type,'size':specs.title_font_size})
        

    @staticmethod
    def run_plot_bar(ax,model,params,study_tag,target,study,plot_t,IDs=[]):
        """ run and plot bar
        """
        sims = model.simulate_study(study_tag=study_tag,params= params,study=study)
#        print('--',study_tag,target,sims)
        sims,exps,exps_std,pvalues = process_data.sort(sims=sims,exps=study,target = target,plot_t=plot_t)
        sims,exps = tools.normalize(study_tag = study_tag,target=target, sims= sims, exps=exps)
        specs = Specs(study_tag)
        xs = specs.bar_positions(sims,specs.delta,specs.D,plot_t=plot_t)
        max_y = 0 # to determine ylim
        for i in range(len(xs)):
            if plot_t == 'bar1':
                labels = ['S','E']
                sims_ID = sims
                exps_ID = exps
                exps_std_ID = exps_std
            elif plot_t == 'bar2':
                labels = ['S-%s'%IDs[i],'E-%s'%IDs[i]]
                sims_ID = sims[i]
                exps_ID = exps[i]
                exps_std_ID = exps_std[i]
            plotTools.plot_bar(ax=ax,specs=specs,x_exp=xs[i][0],x_sim=xs[i][1],sims=sims_ID,exps=exps_ID,exps_std=exps_std_ID,labels=labels, plot_i = i)
            plotTools.draw_p_values(ax=ax,specs=specs,x_exp=xs[i][0], exps=exps_ID,exps_std=exps_std_ID,pvalues=pvalues)
            
            exp_y = max(np.array(exps_ID) + np.array(exps_std_ID)) # to find out the maximum space required for the y range
            sim_y = max(sims_ID)
            
            if exp_y>max_y:
                max_y = exp_y
            if sim_y>max_y:
                max_y = sim_y
        duration = study['duration']
        plotTools.ax_postprocess(ax=ax,study_tag=study_tag,specs=specs,sims=sims,duration=duration,target=target,max_y=max_y)

    @staticmethod
    def run_plot_line(ax,study_tag,model_sbml,params,target,ID,study):
        duration = study['duration']
        inputs = study[ID]['inputs']
        specs = Specs(study_tag)
        params_copy = copy.deepcopy(params)
        for key,value in inputs.items():
            params_copy[key] = value
        sims = Macrophage.run_sbml_model(model_sbml=model_sbml,params = params_copy,selections=['TIME',target],duration=duration)
#        sims = Macrophage.run_sbml_model_recursive(model_sbml=model_sbml,params = params_copy,selections=['TIME',target],duration=duration)
        ii_0 = 0 # start ploting from this index onward
        x = sims['time'][ii_0:]
        ax.plot(x,sims[target][ii_0:],color='black',linewidth=specs.line_width, label = 'S')
        obs_xx = study['selections'][target]
        obs_yy = study[ID]['expectations'][target]['mean']
        ax.scatter(obs_xx,obs_yy,color='r', linewidth=specs.line_width, label = 'E')
        ax.legend()
        plotTools.ax_postprocess(ax,study_tag=study_tag,sims=sims,duration=duration,specs=specs,target=target)

    @staticmethod
    def run_plot_line_multi_target(ax,study_tag,model_sbml,params,targets,ID,study):
        duration = study['duration']
        inputs = study[ID]['inputs']
        specs = Specs(study_tag)

        params_copy = copy.deepcopy(params)
        for key,value in inputs.items():
            params_copy[key] = value
        sims_raw = Macrophage.run_sbml_model(model_sbml = model_sbml,params = params_copy,selections=['TIME']+targets,duration=duration)
        sims = {}
        for target in targets:
            sims[target] = sims_raw[target]
        # if study_tag == 'eq_IL8':
        #     # sims['IL8_R'] = np.array(sims_raw['IL8_R'])/model_sbml['IL8_R_0']
        #     sims['IL8_m'] = np.array(sims_raw['IL8_m'])/model_sbml['IL8_m_0']
        #     sims['IFNGR'] = np.array(sims_raw['IFNGR'])/model_sbml['IFNGR_0']
        i_0 = 0
        x = sims_raw['time'][i_0:]

#        obs_yy = [1 for i in range(len(x))]
#        ax.plot(x,obs_yy,linestyle='--',color='r', linewidth=specs.line_width, label = 'Expectation')
        jj=0
        for target in targets:
            label = target
            if target in list(labels.keys()):
                label = labels[target]
            ax.plot(x,sims[target][i_0:],color=specs.colors[jj],linewidth=specs.line_width, label = 'S: {}'.format(label))
            jj+=1
        ax.legend()
        plotTools.ax_postprocess(ax,study_tag=study_tag,sims=sims,duration=duration,specs=specs)


