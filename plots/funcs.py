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
    figsize = (22,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1,1,2]
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

    g_size = 8
    study_tag,target = 'M05_NFKBn','nNFKB_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[3],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
        
    fig.tight_layout()
def P22_plot(model_sbml,model_macrophage,params,observations):
    figsize = (12,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations
    
    
    width_ratios = [1.5,1.5]
    nn = len(width_ratios)
    fig, axes = plt.subplots(1, nn, gridspec_kw={'width_ratios': width_ratios},figsize=figsize)
    jj = 0
#    g_size = 6
#    study_tag,target = 'M18','nIFNGR'
#    IDs = obs[study_tag]['IDs']
#    plotTools.run_plot_bar(ax=axes[0],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
#
#    g_size = 6
#    study_tag,target = 'M18','nIL4R'
#    IDs = obs[study_tag]['IDs']
#    plotTools.run_plot_bar(ax=axes[1],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)

    g_size = 6
    study_tag,target = 'M18','nIL1b'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    g_size = 6
    study_tag,target = 'M18','nIL10'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    # g_size = 4
    # study_tag,targets = 'eq_IL8',['nIL8','nIL8R']
    # ID = obs[study_tag]['IDs'][0]
    # plotTools.run_plot_line_multi_target(ax=axes[4],study_tag=study_tag,targets=targets,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])

    
    fig.tight_layout()



def P31_plot(model_sbml,model_macrophage,params,observations): # Mg regulates IKB/NFKB
    figsize = (12,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)


    size = 4
    study_tag,target = 'S12_IKBa_mg','nIKB'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
#    size = 4
#    study_tag,target = 'Q21_IkBa_6h','IKB'
#    IDs = obs[study_tag]['IDs']
#    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
#    jj+=1
    
#    size = 4
#    study_tag,target = 'Q21_IkBa_72h','IKB'
#    IDs = obs[study_tag]['IDs']
#    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
#    jj+=1

#    size = 5
#    study_tag,target = 'Q21_IkBa','IKB'
#    IDs = obs[study_tag]['IDs']
#    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar2',IDs=IDs)
#    jj+=1

    size = 4
    study_tag,target = 'S12_NFKBn_mg','nNFKB_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
    size = 4
    study_tag,target = 'Q21_NFKBn_72h','nNFKB_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    fig.tight_layout()
def P32_plot(model_sbml,model_macrophage,params,observations): # Mg regulates IL8
    figsize = (8,4)
    fig = plt.figure(figsize=figsize)
    fig.canvas.draw()
    obs = observations

    width_ratios= [1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)


    size = 4
    study_tag,target = 'Q21_Mg_IL8','nIL8'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1
    
#    
    
    fig.tight_layout()

def P1_eq_plot(model_sbml,params,observations): # equalibrium for P1 to P3
    figsize = (9,4)
    fig = plt.figure(figsize=figsize)
    obs = observations
    fig.canvas.draw()
    nn = 2
    jj = 1

    ax = fig.add_subplot(1,nn,jj)
    jj+=1
    plot_eq_mg(ax=ax,model_sbml=model_sbml,params=params,observations=observations)
    

    ax = fig.add_subplot(1,nn,jj)
    plot_Q21_eq(ax=ax,model_sbml = model_sbml ,params = params,observations =observations)
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
    study_tag,target = 'R05_nMg_f','nMg_f'
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

    width_ratios= [1,1,1,1,1]
    jj = 0
    fig, axes = plt.subplots(1, len(width_ratios), gridspec_kw={'width_ratios': width_ratios},figsize=figsize)


    study_tag,target,ID = 'Q21_Mg','nMg_f','Mg_8'
    plotTools.run_plot_line(ax=axes[jj],study_tag=study_tag,target=target,ID=ID,model_sbml=model_sbml,params=params,study=observations[study_tag])
    jj+=1

    study_tag,target = 'Q21_TRPMn','nTRPM_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    study_tag,target = 'Q21_TRPM','nTRPM'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    study_tag,target = 'Q21_M7CKn','nM7CK_n'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    study_tag,target = 'Q21_H3S10','npH3S10'
    IDs = obs[study_tag]['IDs']
    plotTools.run_plot_bar(ax=axes[jj],model=model_macrophage,params=params,study_tag=study_tag,target=target,study=obs[study_tag],plot_t='bar1',IDs=IDs)
    jj+=1

    fig.tight_layout()
    plt.savefig(os.path.join(dir_output, fig_tag+'.svg'),bbox_inches='tight')






def plot_eq_mg(ax,model_sbml,params,observations):

    study_tag = 'eq_mg'
    obs = observations[study_tag]
    targets = ['nMg_f', 'nMg', 'nMg_ATP']
    duration = obs['experiment_period']
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

    obs_yy = [1 for i in range(len(x))]
    ax.plot(x,obs_yy,linestyle='--',color='r', linewidth=specs.line_width, label = 'Expectation')
    
    jj=0
    for target in targets:
        ax.plot(x,sims[target][i_0:],color=specs.colors[jj],linewidth=specs.line_width, label = r'S: {}'.format(labels[target]))
        
        jj+=1
    ax.legend()
    plotTools.ax_postprocess(ax,study_tag=study_tag,sims=sims,specs=specs)

    
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
#    results =Macrophage.run_sbml_model_recursive(model_sbml,params = params_copy,selections=targets,duration=duration)

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
