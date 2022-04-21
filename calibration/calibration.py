import os
import sys
main_dir = '/Users/matin/Downloads/testProjs/intracellular_M'
sys.path.insert(0,main_dir)
from tools import dirs
from tools.tools import calibrate
import json
from models.params import free_params_p, fixed_params
from data.observations import observations,packages
from models.models import Macrophage

target_package = 'P1'
free_params = free_params_p[target_package]
studies = {}
for study in packages[target_package]:
    studies[study] = observations[study]
print(free_params)
print(studies.keys())

# model_obj = Macrophage(dir_model=dirs.dir_model)
model_obj = Macrophage(dir_model=dirs.dir_M1_model)

def callback(xk, convergence):
    params = {**fixed_params}
    keys = list(free_params.keys())
    for ii in range(len(keys)):
        params[keys[ii]] = xk[ii]
    error = model_obj.run(params = params,studies=studies)
    if  error < 0.01:
        return True
    
class Strategies:
    best1bin = 'best1bin'
    rand1exp = 'rand1exp'
inferred_params = calibrate(model = model_obj,fixed_params = fixed_params, free_params=free_params, studies = studies, n_proc=1,disp=True,max_iters=200,strategy=Strategies.best1bin,callback=callback)
with open(dirs.dir_calib_output,'w') as f:
    f.write(json.dumps(inferred_params))
os.system('say "Hey Matin, calibration is done, come back"')