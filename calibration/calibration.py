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
    model_t = 'M1'
    target_package = 'P12'
    free_params = free_params_p[target_package]
    studies = select_obs(packages[target_package])
print(settings.free_params)
print(settings.studies.keys())

model = Macrophage(settings.model_t)

def output(inferred_params,error):
    with open(os.path.join(dirs.dir_outputs,'inferred_params_{}.json'.format(settings.target_package)),'w') as f:
        f.write(json.dumps(inferred_params,indent=4))
    with open(os.path.join(dirs.dir_outputs,'error_{}.json'.format(settings.target_package)),'w') as f:
        f.write(str(error))
def callback(xk, convergence):
    params = {}
    keys = list(settings.free_params.keys())
    for ii in range(len(keys)):
        params[keys[ii]] = xk[ii]
    _params = params
    params = {**params,**fixed_params}
    error = model.run(params = params,studies=settings.studies)
    output(_params,error)
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
    inferred_params,error = calib_obj = calibrate(cost_function=cost_function, workers=workers,maxiter=500,callback=callback,free_params=settings.free_params)

    output(inferred_params,error)
os.system('say "Hey Matin, calibration is done, come back"')
