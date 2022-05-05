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
params = {**fixed_params}
if True: # apply inferred params
    target_package = 'P21'
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(target_package)),'r') as file:
        inferred_params = json.load(file)
    params = {**params,**inferred_params}
model_sbml = Macrophage.models['IL8']
# tags = ['IL4R_0','IL4R','nIL4R','IL8','F_il8_il4p','F_il8_ifngr','IL8_R','IFNGR','nIL1b','TNFa']
# tags = ['IFNGR_0','IFNGR','nIFNGR']
tags = ['IRAK4','aTRAF6','NFKB_n']
inputs = {}
rr1 = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=1500,params={**params,**inputs},selections=['time']+tags,activation=True)
inputs = {'IL8':100*1000}
rr2 = Macrophage.run_sbml_model(model_sbml=model_sbml,duration=1500,params={**params,**inputs},selections=['time']+tags,activation=True)

# 
# print(tags[0],rr2[tags[0]][0:2])
# print(tags[1],rr2[tags[1]][0:2])
# print(tags[2],rr2[tags[2]][0:2])
# 
# print(rr2['F_il8_il4p'][0:-3])

fig = plt.figure()
ax = fig.add_subplot(1,2,1)
ax.plot(rr1['time'],rr1['IRAK4'],label = 'ctr')
ax.plot(rr2['time'],rr2['IRAK4'],label = '100')
ax.set_title('IRAK4')
ax.legend()
ax = fig.add_subplot(1,2,2)
ax.plot(rr1['time'],rr1['aTRAF6'],label = 'ctr')
ax.plot(rr2['time'],rr2['aTRAF6'],label = '100')
ax.set_title('aTRAF6')
ax.legend()
fig.tight_layout()
plt.savefig('test.png')

