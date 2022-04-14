import tellurium as te
import numpy as np
import json
import pathlib as pl
import os
file_dir = pl.Path(__file__).parent.absolute()
output_dir = os.path.join(file_dir,"results")
original_model_dir = os.path.join(file_dir,"Zhao_2021.xml")

Zhao_2021 = te.loadSBMLModel(original_model_dir)
targets = ['NFKB', 'pIKK', 'TAK1']
# targets = ['pIKK']
duration = 1500
stack_results = [] # store for each iteration


results=Zhao_2021.simulate(0,duration,selections=['TIME']+targets)
# mean_results_no_tag = np.array([mean_results[tag] for tag in mean_results.keys()])
results_dict = {'time':list(results['time'])} 
for key in targets:
    results_dict[key] = list(results[key])

with open('samples.json', 'w', encoding='utf-8') as f:
    json.dump(results_dict, f, ensure_ascii=False, indent=4)