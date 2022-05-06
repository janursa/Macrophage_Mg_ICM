main_dir = '/Users/matin/Downloads/testProjs/intracellular_M'
import sys
import os
import json
sys.path.insert(0,main_dir)
from data.observations import observations,t2m,select_obs
from models.params import fixed_params
from tools import dirs, tools
from models.models import Macrophage
from plots import funcs 
import matplotlib
import matplotlib.pyplot as plt
import tellurium as te

params = {**fixed_params}
if True: # apply inferred params
    target_package = 'P12'
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
# params['n_h3s10_il8_p'] =5
# inferred_params = {
#     "kd_mg_ikb_d": 0.6372241221592411,
#     "n_mg_ikb_d": 8.644881763454306,
#     "kd_h3s10_il8_p": 0.02857981446166838,
#     "n_h3s10_il8_p": 5,
#     "kd_h3s10_ikb_p": 4.782521786248381,
#     "n_h3s10_ikb_p": 24.29685438635782
# }
# params = {**params,**inferred_params}
# params['n_h3s10_il8_p'] = 4
model_t = 'M1'

model_sbml = Macrophage.create_sbml_model(model_t)

tags = ['IL8','NFKB_n','pH3S10']
duration = 200*60
inputs = {}
ctr = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=True)
inputs = {'Mg_e':8}
rr_8 = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=True)
inputs = {'Mg_e':0.08}
rr_08 = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=duration,params={**params,**inputs},selections=['time']+tags,activation=True)


fig = plt.figure()
ax = fig.add_subplot(1,3,1)
tt = 0
ax.plot(ctr['time'][tt:],ctr['IL8'][tt:],label = 'ctr')
ax.plot(rr_8['time'][tt:],rr_8['IL8'][tt:],label = 'Mg 8')
ax.plot(rr_08['time'][tt:],rr_08['IL8'][tt:],label = 'Mg 0.08')
ax.set_title('IL8')
ax.legend()
ax = fig.add_subplot(1,3,2)
tt = 0
ax.plot(ctr['time'][tt:],ctr['NFKB_n'][tt:],label = 'ctr')
ax.plot(rr_8['time'][tt:],rr_8['NFKB_n'][tt:],label = 'Mg 8')
ax.plot(rr_08['time'][tt:],rr_08['NFKB_n'][tt:],label = 'Mg 0.08')
ax.set_title('NFKB_n')
ax.legend()
ax = fig.add_subplot(1,3,3)
tt = 0
ax.plot(ctr['time'][tt:],ctr['pH3S10'][tt:],label = 'ctr')
ax.plot(rr_8['time'][tt:],rr_8['pH3S10'][tt:],label = 'Mg 8')
ax.plot(rr_08['time'][tt:],rr_08['pH3S10'][tt:],label = 'Mg 0.08')
ax.set_title('pH3S10')
ax.legend()

fig.tight_layout()
# _dir = os.path.join(dirs.dir_outputs,'plots','test.png')
# fig.show()
#plt.savefig(_dir)

