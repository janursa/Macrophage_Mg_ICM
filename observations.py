observations = {
	'studies':['Quao_2021'],
	'Quao_2021':{# the effect of different Mg ion concentrations on IL8 and IL1b
		'exposure_period':None, #when this is none, the cells are exposed to the stimuli the whole time
		'culture_volume':None, #ml
		'experiment_period':72, # hours
		'measurement_scheme':{
			# 'IL8': [6,12,24,48,72], #hours
			'IL1b':[6,12,24,48,72]
			# 'IL1b':[24,48,72]
		},
		# 'IDs': ['Mg_.08','Mg_8'],
		'IDs': ['Mg_8'],
		'Mg_.08':{
			'inputs':{
						"Mg_e_mM": 0.08 # mM
					},
			"expectations": {			
				# "IL8": {
				# 	'mean':[10000,19000,26000,29500,32000], #pg/ml
				# 	'std': [0,0,0,0,0]
				# },
				"IL1b": {
					# 'mean':[1700,2950,2600,3000,3400], #pg/ml
					# 'std': [200,200,450,150,300]
					'mean':[3000,3400] #pg/ml
				},
			}
		},
		'Mg_8':{
			'inputs':{
						"Mg_e_mM": 8 # mM
					},
			"expectations": {			
				# "IL8": {
				# 	'mean':[16000,21000,30000,38000,37000], #pg/ml
				# 	'std': [0,3000,0,1000,1000]
				# },
				"IL1b": {
					'mean':[1950,1200,1000,1700,1950], #pg/ml
					# 'std': [750,350,0,200,250]
					# 'mean':[1000,1700,1950] #pg/ml
				},
			}	
		}
	},
	"scale": 1,
}
import json
with open('observations.json','w') as file:
	file.write(json.dumps(observations,indent=4))
