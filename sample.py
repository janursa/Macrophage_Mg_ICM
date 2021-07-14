from models import Zhao_2021, PARAMS,average,run
import numpy as np
import json

stack_results = [] # store for each iteration
for i in range(PARAMS.replica_n):
    Zhao_2021.reset()
    results_dict=run(Zhao_2021,targets=PARAMS.targets,duration=PARAMS.duration)
    stack_results.append(results_dict)
mean_results = average(stack_results)

# mean_results_no_tag = np.array([mean_results[tag] for tag in mean_results.keys()])

with open('samples.json', 'w', encoding='utf-8') as f:
    json.dump(mean_results, f, ensure_ascii=False, indent=4)