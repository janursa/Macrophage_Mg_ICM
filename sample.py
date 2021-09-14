import tellurium as te
import numpy as np
import json

Zhao_2021 = te.loadSBMLModel("Zhao_2021.xml")
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