import os
import sys
from pathlib import Path

dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
from tools import dirs
from tools.tools import calibrate
import json
from models.params import free_params_p, fixed_params
from data.observations import observations,packages,select_obs
from models.models import Macrophage
# %load_ext autoreload
# %autoreload

class settings:
    model_t = 'IL8'
    target_package = 'P21'
    free_params = free_params_p[target_package]
    studies = select_obs(packages[target_package])
print(settings.free_params)
print(settings.studies.keys())

model = Macrophage(settings.model_t)

def output(inferred_params):
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(settings.target_package)),'w') as f:
        f.write(json.dumps(inferred_params,indent=4))
def callback(xk, convergence):
    params = {}
    keys = list(settings.free_params.keys())
    for ii in range(len(keys)):
        params[keys[ii]] = xk[ii]
    output(params)
    params = {**params,**fixed_params}
    error = model.run(params = params,studies=settings.studies)
    if  error < 0.01:
        return True

def cost_function(calib_params_values):
    calib_params = {}
    for key,value in zip(settings.free_params.keys(),calib_params_values):
        calib_params[key] = value
    params = {**calib_params,**fixed_params}
    error = model.run(params=params,studies=settings.studies)

    return error
    
class Strategies:
    best1bin = 'best1bin'
    rand1exp = 'rand1exp'
if __name__ == '__main__':
    workers = sys.argv[1]
    inferred_params,_ = calib_obj = calibrate(cost_function=cost_function, workers=workers,maxiter=100,callback=callback,free_params=settings.free_params)

    output(inferred_params)
os.system('say "Hey Matin, calibration is done, come back"')