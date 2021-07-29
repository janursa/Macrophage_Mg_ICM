from models import Zhao_model, PARAMS,batch_run
import numpy as np
import json
import matplotlib.pyplot as plt
import tellurium as te
import time

means = []
begin= time.time()
for i in range(100):
    mean_results = batch_run(model=Zhao_model,targets=PARAMS.targets, duration=PARAMS.duration,replica_n=PARAMS.replica_n)
    means.append(mean_results['IKK'][-1])
# print(means)
end=time.time()
plt.plot(means)
plt.savefig('stoch.svg')
print('took ',end-begin)
# means = np.array(means)
# stds = []
# for i in range(PARAMS.duration):
#     stds.append(np.std(means[:,i]))
# print(stds)
# print('mean std',np.mean(stds))
# mean_results_no_tag = np.array([mean_results[tag] for tag in mean_results.keys()])

# with open('samples.json', 'w', encoding='utf-8') as f:
#     json.dump(mean_results, f, ensure_ascii=False, indent=4)