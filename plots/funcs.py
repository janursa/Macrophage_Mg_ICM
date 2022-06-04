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

def LPS_plot(model_sbml,model_macrophage,params,observations):
    figsize = (12,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)
    
    g_size = 4
    study_tag,target = 'S12_LPS','nIKB'
    ID = obs[study_tag]['IDs'][0]
    plotTools.run_plot_line(ax=axes[jj],study_tag=study_tag,target=target,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
    jj+=1

    g_size = 4
    study_tag,target = 'S12_LPS','nTNFa'
    ID = obs[study_tag]['IDs'][0]
    plotTools.run_plot_line(ax=axes[jj],study_tag=study_tag,target=target,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
    jj+=1

    g_size = 4
    study_tag,target = 'B20_LPS','nTNFa'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    fig.tight_layout()

def P2_eq(model_sbml,model_macrophage,params,observations):
    figsize = (16,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1,1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)

    g_size = 4
    study_tag,targets = 'eq_IL8', ['nIL8_m','F_il8_irak']
    ID = obs[study_tag]['IDs'][0]
    plotTools.run_plot_line_multi_target(ax=axes[jj],study_tag=study_tag,targets=targets,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
    jj+=1

    g_size = 4
    study_tag,targets = 'eq_IL8', ['F_nfkb_il8_p','F_ap1_il8_p']
    ID = obs[study_tag]['IDs'][0]
    plotTools.run_plot_line_multi_target(ax=axes[jj],study_tag=study_tag,targets=targets,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
    jj+=1

    g_size = 4
    study_tag,targets = 'eq_IL8', ['F_rho_nfkb_a','F_rho_stat3_a','F_rho_pi3k_a']
    ID = obs[study_tag]['IDs'][0]
    plotTools.run_plot_line_multi_target(ax=axes[jj],study_tag=study_tag,targets=targets,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
    jj+=1

    g_size = 4
    study_tag,targets = 'eq_IL8',list(observations['eq_IL8']['selections'].keys())
    ID = obs[study_tag]['IDs'][0]
    plotTools.run_plot_line_multi_target(ax=axes[jj],study_tag=study_tag,targets=targets,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
    jj+=1
    

    
    fig.tight_layout()

def P2_ICs_plot(model_sbml,model_macrophage,params,observations):
    figsize = (18,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1,2]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)

    
    g_size = 4
    study_tag,target = 'M05_IT','nIRAK4'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    g_size = 4
    study_tag,target = 'M05_IT','naTRAF6'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    g_size = 8
    study_tag,target = 'M05_NFKBn','nNFKB_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    fig.tight_layout()

def P2_IL6_IC_plot(model_sbml,model_macrophage,params,observations):
    # figsize = (28,4)
    figsize = (10,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1.5]
    # width_ratios= [1,1,1,1.25]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)

    g_size = 4
    study_tag,targets = 'eq_IL6',list(observations['eq_IL6']['selections'].keys())
    ID = obs[study_tag]['IDs'][0]
    plotTools.run_plot_line_multi_target(ax=axes[jj],study_tag=study_tag,targets=targets,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
    jj+=1

    # g_size = 4
    # study_tag,target = 'B17','npSTAT3'
    # IDs = obs[study_tag]['IDs']
    # plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    # jj+=1

    # axes[jj].set_title('B17: nNFKB_n')
    # g_size = 4
    # study_tag,target = 'B17','nNFKB_n'
    # IDs = obs[study_tag]['IDs']
    # plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    # jj+=1

    g_size = 6
    study_tag,target = 'F17','npSTAT3'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    

    fig.tight_layout()

def P2_IL6_CYs_plot(model_sbml,model_macrophage,params,observations):
    figsize = (14,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1.5,1,1.5]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)

    g_size = 6
    study_tag,target = 'F17','nIL10'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    g_size = 4
    study_tag,target = 'F14','nIL10'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1


    g_size = 6
    study_tag,target = 'N03','nTNFa'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    fig.tight_layout()

def P2_receptors_plot(model_sbml,model_macrophage,params,observations):
    figsize = (12,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations
    
    
    width_ratios = [1.5,1.5]
    nn = len(width_ratios)
    fig, axes = plt.subplots(1, nn, gridspec_kw={'width_ratios': width_ratios},figsize=figsize)
    jj = 0
    
    g_size = 6
    study_tag,target = 'M18','nIFNGR'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    g_size = 6
    study_tag,target = 'M18','nIL4R'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    fig.tight_layout()
    
def P23_plot(model_sbml,model_macrophage,params,observations):
    figsize = (24,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations
    
    width_ratios = [1.5,1.5,1.5,1.5]
    nn = len(width_ratios)
    fig, axes = plt.subplots(1, nn, gridspec_kw={'width_ratios': width_ratios},figsize=figsize)
    jj = 0
    g_size = 6
    study_tag,targets = 'M18',list(observations['M18']['selections'].keys())
    IDs = obs[study_tag]['IDs']
    for target in targets:
        plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
        jj+=1
#
#     g_size = 6
#     study_tag,target = 'M18','nIL10'
#     IDs = obs[study_tag]['IDs']
#     plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
#     jj+=1
# #
#     g_size = 6
#     study_tag,target = 'M18','nTNFa'
#     IDs = obs[study_tag]['IDs']
#     plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
#     jj+=1

#     g_size = 6
#     study_tag,target = 'M18','nIL6'
#     IDs = obs[study_tag]['IDs']
#     plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
#     jj+=1
    
    fig.tight_layout()


def P3_eq_plot(model_sbml,model_macrophage,params,observations): # equalibrium for P1 to P3
    figsize = (8,4)
    fig = plt.figure(figsize=figsize)
    obs = observations
    fig.canvas.draw()
    width_ratios= [1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)


    study_tag,targets = 'eq_combined',['F_h3s10_ikb','F_p3s10_il8_p']
    ID = obs[study_tag]['IDs'][0]
    plotTools.run_plot_line_multi_target(ax=axes[jj],study_tag=study_tag,targets=targets,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
    jj+=1


    fig.tight_layout()

def P3_IKB_plot(model_sbml,model_macrophage,params,observations): # Mg regulates IKB/NFKB
    figsize = (9,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1.25]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)


    size = 4
    study_tag,target = 'S12_IKBa_mg','nIKB'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    size = 5
    study_tag,target = 'Z19_IKB_NFKB','nIKB'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
   

    
    fig.tight_layout()
    
def P3_NFKB_plot(model_sbml,model_macrophage,params,observations): # Mg regulates NFKB
    figsize = (13,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1.25,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)
    
    size = 4
    study_tag,target = 'S12_NFKBn_mg','nNFKB_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    size = 5
    study_tag,target = 'Z19_IKB_NFKB','nNFKB_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    size = 4
    study_tag,target = 'B20_NFKBn','nNFKB_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    fig.tight_layout()
       

def P3_cytokines1_plot(model_sbml,model_macrophage,params,observations): # Mg regulates IL8,IL10
    figsize = (17,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1.25,1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)


    size = 4
    study_tag,target = 'Q21_Mg_IL8','nIL8'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    size = 5
    study_tag,target = 'Z19_IL10','nIL10'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    size = 4
    study_tag,target = 'Q21_cytokines_72h','nIL10'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
#
    size = 4
    study_tag,target = 'F18_cytokines','nIL10'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    
    
#
    
    fig.tight_layout()

def P3_cytokines2_plot(model_sbml,model_macrophage,params,observations): # Mg regulates TNF
    figsize = (12,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)
    
    
    size = 4
    study_tag,target = 'Q21_cytokines_72h','nTNFa'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
        
    size = 4
    study_tag,target = 'F18_cytokines','nTNFa'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    size = 4
    study_tag,target = 'B20_TNFa','nTNFa'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    fig.tight_layout()
    
def P3_cytokines3_plot(model_sbml,model_macrophage,params,observations): # Mg regulates nIL1b
    figsize = (8,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)

    size = 4
    study_tag,target = 'Q21_cytokines_72h','nIL1b'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    
    size = 4
    study_tag,target = 'F18_cytokines','nIL1b'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    fig.tight_layout()

def P3_cytokines4_plot(model_sbml,model_macrophage,params,observations): # Mg regulates IL8
    figsize = (8,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)

    size = 4
    study_tag,target = 'Q21_14d','nIKB'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    size = 4
    study_tag,target = 'Q21_14d','nNFKB_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    fig.tight_layout()


def P1_eq(model_sbml,model_macrophage,params,observations): # equalibrium for P1 to P3
    figsize = (8,4)
    fig = plt.figure(figsize=figsize)
    obs = observations
    fig.canvas.draw()
    width_ratios= [1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)
    
    plot_eq_mg(ax=axes[jj],model_sbml=model_sbml,params=params,observations=observations)
    jj+=1
    
    plot_Q21_eq(ax=axes[jj],model_sbml = model_sbml ,params = params,observations =observations)
    jj+=1


    fig.tight_layout()

def P11_plot(model_sbml,params,observations): # equalibrium for P1 to P3
    figsize = (8,4)
    fig = plt.figure(figsize=figsize)
    obs = observations
    fig.canvas.draw()
    width_ratios= [1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)
    
    
    study_tag,target = 'R05_nMg_f','nMg_f'
    ID = obs[study_tag]['IDs'][0]
    plotTools.run_plot_line(ax=axes[jj],study_tag=study_tag,target=target,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
    jj+=1
    
    study_tag,target,ID = 'Q21_Mg','nMg_f','Mg_8'
    plotTools.run_plot_line(ax=axes[jj],study_tag=study_tag,target=target,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
    jj+=1

    fig.tight_layout()

def P12_plot(model_sbml,model_macrophage,params,observations): # plotting sim vs exp measurements for P1 to p3
    fig_tag = 'P12'
    figsize = (20,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1,1,1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)

    study_tag,target = 'Q21_M1','nATP'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    study_tag,target = 'Q21_M1','nTRPM_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    study_tag,target = 'Q21_M1','nTRPM'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    study_tag,target = 'Q21_M1','nM7CK_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    study_tag,target = 'Q21_M1','npH3S10'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    fig.tight_layout()
    # plt.savefig(os.path.join(dir_output, fig_tag+'.svg'),bbox_inches='tight')



def plot_eq_mg(ax,model_sbml,params,observations):

    study_tag = 'eq_mg'
    obs = observations[study_tag]
    targets = ['nMg_f', 'nMg', 'nMg_ATP']
    duration = obs['duration']
    IDs = obs['IDs']
    ID = IDs[0] #only 1 ID
    inputs = obs[ID]['inputs']
    specs = Specs(study_tag)

    
    params_copy = copy.deepcopy(params)
    for key,value in inputs.items():
        params_copy[key] = value

    selections = ['nMg_f', 'Mg','Mg_ATP' ]
    sims_raw = Macrophage.run_sbml_model(model_sbml=model_sbml,params = params_copy,selections=selections,duration=duration)
    sims = {}
    for selection in selections:
        sims[selection] = sims_raw[selection]
    sims['nMg'] = np.array(sims['Mg'])/18.5;
    sims['nMg_ATP'] = np.array(sims['Mg_ATP'])/model_sbml['Mg_ATP_0'];
    
    i_0 = 1
    x = sims_raw['time'][i_0:]

#    obs_yy = [1 for i in range(len(x))]
#    ax.plot(x,obs_yy,linestyle='--',color='r', linewidth=specs.line_width, label = 'Expectation')
    
    jj=0
    for target in targets:
        ax.plot(x,sims[target][i_0:],color=specs.colors[jj],linewidth=specs.line_width, label = r'S: {}'.format(labels[target]))
        
        jj+=1
    ax.legend()
    plotTools.ax_postprocess(ax,study_tag=study_tag,sims=sims,duration=duration,specs=specs)

    
def plot_Q21_eq(ax,model_sbml,params,observations):
    study_tag = 'Q21_eq'
    study1,study2 = 'Q21_eq_trpm','Q21_eq_h3s10'
    targets1 = list(observations[study1]['selections'].keys())

    targets2 = list(observations[study2]['selections'].keys())

    targets = targets1 + targets2

    duration = observations[study1]['duration']
    IDs = observations[study1]['IDs']
    ID = IDs[0] #only 1 ID
    inputs = observations[study1][ID]['inputs']
    
    params_copy = copy.deepcopy(params)
    for key,value in inputs.items():
        params_copy[key] = value
    results =Macrophage.run_sbml_model(model_sbml,params = params_copy,selections=targets,duration=duration)
#    results =Macrophage.run_sbml_model_recursive(model_sbml,params = params_copy,selections=targets,duration=duration)

    ii = 0

    specs = Specs(study_tag)

#    ax.plot(results['time'],[1 for i in range(len(results['time']))],linestyle='--',color='r', label = 'Expectation')
    for target in targets: # one subgraph for each target
        x = results['time']
        ax.plot(x,results[target],linewidth=specs.line_width, color=specs.colors[ii],label = '%s'%labels[target])
        ii+=1
    
    plotTools.ax_postprocess(ax=ax,study_tag=study_tag,specs=specs,duration=duration,sims=results)

