observations = {
	'studies':['Quao_2021_Mg'],
	"Quao_2021_Mg": { 
        "exposure_period": None,
        "culture_volume": None,
        "experiment_period": 13, #minutes
        "measurement_scheme": {
            "Mg": [0,6,8,12]
        },
        "IDs": ["Mg_8"],
        "Mg_8": {
            "inputs": {
                "Mg_e": 8
            },
            "expectations": {
                "Mg": {
                    "mean": [100,145,125,125],
                    "std": [0]
                }
            }
        }
    },
	'Quao_2021_TRPM':{# the effect of different Mg ion concentrations on TRPM
		'exposure_period':None, #when this is none, the cells are exposed to the stimuli the whole time
		'culture_volume':None, #ml
		'experiment_period':72*60, # minutes
		'measurement_scheme':{
			'TRPM': [72*60]
		},
		'IDs': ['Mg_.08','Mg_.8','Mg_8'],
		# 'IDs': ['Mg_.08','Mg_.8'],
		'Mg_.08':{
			'inputs':{
						"Mg_e": 0.08 # mM
					},
			"expectations": {			
				"TRPM": {
					'mean':[1.1], #pg/ml
					'std': [.1]
				}
			}
		},
		'Mg_.8':{
			'inputs':{
						"Mg_e": 0.8 # mM
					},
			"expectations": {			
				"TRPM": {
					'mean':[1], #pg/ml
					'std': [.01]
				}
			}
		},
		'Mg_8':{
			'inputs':{
						"Mg_e": 8 # mM
					},
			"expectations": {			
				"TRPM": {
					'mean': [2.5], #pg/ml
					'std': [0.7]
				}
			}
		}

	},
	'Quao_2021_1':{# the effect of different Mg ion concentrations on IL8 and IL1b
		'exposure_period':None, #when this is none, the cells are exposed to the stimuli the whole time
		'culture_volume':None, #ml
		'experiment_period':72*60, # hours
		'measurement_scheme':{
			# 'IL8': [6*60,12*60,24*60,48*60,72*60] #hours
			# 'IL8': [72] #hours
			# 'IL1b':[6,12,24,48,72]
			# 'IL1b':[24,48,72]
		},
		'IDs': ['Mg_.08','Mg_8'],
		'Mg_.08':{
			'inputs':{
						"Mg_e": 0.08 # mM
					},
			"expectations": {			
				"IL8": {
					# 'mean':[10000,19000,26000,29500,32000], #pg/ml
					# 'std': [0,0,0,0,0]
					'mean':[32000], #pg/ml
				},
				"IL1b": {
					# 'mean':[1700,2950,2600,3000,3400], #pg/ml
					# 'std': [200,200,450,150,300]
					'mean':[3000,3400] #pg/ml
				},
			}
		},
		'Mg_8':{
			'inputs':{
						"Mg_e": 8 # mM
					},
			"expectations": {			
				"IL8": {
					# 'mean':[16000,21000,30000,38000,37000], #pg/ml
					# 'std': [0,3000,0,1000,1000]
					'mean':[37000], #pg/ml
				},
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
