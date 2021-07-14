import pathlib
current_file = pathlib.Path(__file__).parent.absolute()
import os
import sys
ABC_path = os.path.join(current_file,'..','ABayesianC')
sys.path.insert(1,ABC_path)
from ABayesianC import tools
import json

from models import *

# load the samples from the original model
with open('samples.json') as json_file:
	samples = json.load(json_file)
def initial_conditions(model,calib_params):
    for key,value in calib_params.items():
        model[key] = value

def run_replicas(model,calib_params):
    """
    run the model for a given number of replicas
    
    """
    stack_results = [] # store for each iteration
    for i in range(PARAMS.replica_n):
        model.reset()
        initial_conditions(model,calib_params)
        results_dict=run(model,targets=PARAMS.targets,duration=PARAMS.duration)
        stack_results.append(results_dict)
    return stack_results

class Optimize(object):
	"""docstring for ClassName"""
	def __init__(self, calib_params):
		self.calib_params = calib_params
	def run(self):
		stack_results = run_replicas(model = Mg_model,
									 calib_params=self.calib_params)
		mean_results = average(stack_results)

		# calculate the error for each tag by comparing the results to the original model
		tag_errors = []
		for tag in PARAMS.targets:
			abs_diff =abs(np.array(mean_results[tag])-np.array(samples[tag]))
			means = [np.mean(samples[tag]),np.mean(mean_results[tag])]
			mean = np.mean(means)
			tag_error = np.mean(abs_diff/mean)
			# print('tag {} abs_diff {} mean {} tag_error {}'.format(tag,abs_diff,mean,tag_error))
			tag_errors.append(tag_error)
		error = np.mean(tag_errors)
		# print('error',error)
		return error

settings = {
	"MPI_flag": True,
	"sample_n": 100000,
	"top_n": 50,
    "replica_n": 1,
	"output_path": "ABC",
	"test": False,
	"model":Optimize
}


working_dir = os.getcwd()
output_dir = os.path.join(working_dir,settings["output_path"])
try:
	os.makedirs(output_dir)
except:
	pass
sys.path.insert(1,output_dir)

free_params = {'NEMO_IKK':[0,1000],'k301':[0,1000],'k302':[0,1000],'k303':[0,1000],
'k304':[0,1000],'k305':[0,1000],'k306':[0,1000],'k307':[0,1000],'k308':[0,1000]}


if __name__ == "__main__":
	obj = tools.ABC(settings=settings,free_params=free_params)
	obj.sample()
	obj.run()
	obj.postprocessing()
	# obj.run_tests()