import numpy as np
import matplotlib.pyplot as plt
import json
from tools import dirs
import tellurium as te
from tools import dirs
from data.observations import observations
from tools.dirs import dir_model
from models.models import Macrophage
from models.params import fixed_params
from plots.plots import P1_3_eq_plot, P1_3_qualitative_plot, P1_3_plot

model_sbml = te.loadSBMLModel(dir_model)

with open(dirs.dir_calib_output,'r') as file:
    inferred_params = json.load(file)
# inferred_params['k_ntrpm_p0'] = .1
# inferred_params['k_ntrpm_pm'] = 4

params = {**inferred_params,**fixed_params}
# print(params)
# macrophage_obj = Macrophage(observations=observations)
# macrophage_obj.run(params)
flags = {
    'P1_3': True,
}
for key,value in flags.items():
    if key == 'P1_3' and value : 
        print('P1_3 is plotting')
        P1_3_eq_plot(model_sbml=model_sbml,params=params)
#         P1_3_qualitative_plot (model_sbml=model,params=params)
#         P1_3_plot (model_sbml=model,model_macrophage=macrophage_obj,params=params)
