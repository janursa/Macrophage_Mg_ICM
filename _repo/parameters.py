import sys
import pathlib
import os
current_file = pathlib.Path(__file__).parent.absolute()
sys.path.insert(0,current_file)

from observations import observations
fixed_params = {
    
}
free_params_all = {
    
}

def specifications(study):
    
    if study == 'Quao_2021_Mg':
        studies = ['Quao_2021_Mg']
        candidate = {}
    else:
        raise ValueError('specifications of {} not defined.'.format(study))
    free_params  = {}
    for key in candidate:
        free_params[key] = free_params_all[key]
    obs = {'studies':studies}
    for study in studies:
        obs[study] = observations[study]

    return obs,free_params





