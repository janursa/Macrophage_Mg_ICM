t2m = 1 # relationship between each simulation time and minute

range_3h_10mStep = list(range(int(10/t2m),3*int(60/t2m),int(10/t2m)))
range_24h_60mStep = list(range(0,24*int(60/t2m),int(60/t2m)))
range_48h_60mStep = list(range(12*int(60/t2m),48*int(60/t2m),int(60/t2m)))
range_12h_60mStep = list(range(3*int(60/t2m),12*int(60/t2m),int(60/t2m))) # span of 3 to 12 hs


packages = {
	# 'P1' : ['eq_mg','R05_mg_f_n','R05_mg_n','Q21_Mg'], # mg entry and equalibrium
	'P1' : ['eq_mg','R05_mg_f_n','R05_mg_n'], # mg entry and equalibrium
	'P2' : ['Q21_nTRPM','Q21_TRPM','Q21_nM7CK','Q21_eq_trpm'], #mg affects TRPM 
	'P3' : ['Q21_H3S10','Q21_eq_h3s10'], # m7ck regulates H3S10
	'P4' : ['S12_IkBa_mg','Q21_IkBa'], # Mg regulate IkBa 
	# 'P5' : ['Q21_Mg_IL8'], # IL8 
	'P5' : ['eq_IL8'], # IL8 
}


observations = {	
	'S12_IkBa_mg': { # mg influences IkBa level
        "IDs": ['ctr',"Mg_2dot5"],
        "experiment_period":int(3*60/t2m),
        "measurement_scheme": {
            "IKB": [int(3*60/t2m)]
        },
        "ctr": {
            "inputs": {
                "Mg_e": 0.5 ,
            },
            "expectations": {
                "IKB":{'mean':[1],
                		'std':[0]} #normalized format
            }
        },
        "Mg_2dot5": {
            "inputs": {
                "Mg_e": 2.5 ,
            },
            "expectations": {
                "IKB":{'mean':[1.32],
                		'std':[0]} #normalized format
            }
        }
    },
    'Q21_IkBa':{# the effect of different Mg ion concentrations on TRPM_n
		# 'experiment_period':int(72*60/t2m), # minutes
		'experiment_period':int(6*60/t2m), # minutes
		'measurement_scheme':{
			# 'ZZ_IKB': [int(6*60/t2m),int(72*60/t2m)]
			'IKB': [int(6*60/t2m)]
		},
		'IDs': ['Mg_.08','Mg_8'],
		'Mg_.08':{
			'inputs':{
						"Mg_e": 0.08 # mM
					},
			"expectations": {			
				"IKB": {
					# 'mean':[0.05,1.7], 
					'mean':[0.05], 
					# 'std': [0,.1]
					'std': [0]
				}
			}
		},
		'Mg_8':{
			'inputs':{
						"Mg_e": 8 # mM
					},
			"expectations": {			
				"IKB": {
					# 'mean': [0.7,1.67],
					'mean': [0.7], 
					# 'std': [.1,.25]
					'std': [.1]
				}
			}
		}
	},

	'R05_mg_n': { # mg extrusion
        "IDs": ["Mg_02"],
        "experiment_period":int(48*60/t2m),
        "measurement_scheme": {
            "Mg_n": [0]+range_48h_60mStep
        },
        "Mg_02": {
            "inputs": {
                "Mg_e": 0.02,
            },
            "expectations": {
                "Mg_n":{'mean':[1]+[.5 for i in range_48h_60mStep]} #normalized format
            }
        }
    },
    'S12_mg': { # total mg goes almost 4 times for 2.5 mM
        "IDs": ["Mg_2a5"],
        "experiment_period":1*60,
        "measurement_scheme": {
            "Mg_n": [0,1*60]
        },
        "Mg_2a5": {
            "inputs": {
            	'Mg_f' : 0.4,
                "Mg_e": 2.5,
            },
            "expectations": {
                "Mg_n":{'mean':[1,3.6]} #normalized format
            }
        }
    },
    'R05_mg_f_n': { # normalized free mg for an increase in extracellular mg
        "IDs": ["Mg_19"],
        "experiment_period":int(12*60/t2m),
        "measurement_scheme": {
            "Mg_f_n": [0]+range_12h_60mStep,
            # "Mg": [0]+range_12h_60mStep
        },
        "Mg_19": {
            "inputs": {
                "Mg_e": 19,
            },
            "expectations": {
                "Mg_f_n":{'mean':[1]+[1.4 for i in range_12h_60mStep]}, #normalized format
                # "Mg":{'mean':[0,100]} #normalized format
            }
        }
    },
	'eq_mg': { # mg should stay in its equalibrium when there is no stimuli
        "IDs": ["Mg_0dot5"],
        "experiment_period":int(24*60/t2m),
        "measurement_scheme": {
            "Mg_f": range_24h_60mStep,
            'Mg': range_24h_60mStep,
            'Mg_ATP_n': range_24h_60mStep,
            # 'Mg_IM_n': range_24h_60mStep
        },
        "Mg_0dot5": {
            "inputs": {
                "Mg_e": 0.5,
            },
            "expectations": {
                "Mg_f":{'mean':[0.5 for i in range_24h_60mStep]}, 
                'Mg': {'mean':[18.5 for i in range_24h_60mStep]},
                'Mg_ATP_n': {'mean':[1 for i in range_24h_60mStep]},
                # 'Mg_IM_n': {'mean':[1 for i in range_24h_60mStep]}
            }
        }
    },
	'Q21_Mg': { # measurements on the intracellular Mg concentration for a given external Mg
        "exposure_time": None,
        "culture_volume": None,
        'experiment_period':int(3*60/t2m),
        "measurement_scheme": {
            "Mg_f_n": [0]+range_3h_10mStep # keep the stable concentration until 3 hours
        },
        "IDs": ["Mg_8"],
        "Mg_8": {
            "inputs": {
                "Mg_e": 8
            },
            "expectations": {
                "Mg_f_n": {
                    "mean": [1]+[1.25 for i in range_3h_10mStep],
                }
            }
        }
    },
    'Q21_eq_trpm':{# equalibrium in trpm related block
		'experiment_period':int(24*60/t2m), # minutes
		'measurement_scheme':{
			'TRPM': range_24h_60mStep,
			'nTRPM': range_24h_60mStep,
			'nM7CK': range_24h_60mStep
		},
		'IDs': ['Mg_.5'],
		'Mg_.5':{
			'inputs':{
						"Mg_e": 0.5 # mM
					},
			"expectations": {			
				"TRPM": {'mean':[1 for i in range_24h_60mStep]},
				"nTRPM": {'mean':[1 for i in range_24h_60mStep]},
				"nM7CK": {'mean':[1 for i in range_24h_60mStep]}
			}
		}
	},
	'Q21_eq_h3s10':{# equalibrium in h3s10 in related block
		'experiment_period':int(24*60/t2m), # minutes
		'measurement_scheme':{
			'pH3S10': range_24h_60mStep,
		},
		'IDs': ['Mg_.5'],
		'Mg_.5':{
			'inputs':{
						"Mg_e": 0.5 # mM
					},
			"expectations": {			
				"pH3S10": {'mean':[1 for i in range_24h_60mStep]}
			}
		}
	},
	'Q21_nTRPM':{# the effect of different Mg ion concentrations on TRPM_n
		'experiment_period':int(72*60/t2m), # minutes
		'measurement_scheme':{
			'nTRPM': [int(72*60/t2m)]
		},
		'IDs': ['Mg_.08','Mg_.8','Mg_8'],
		'Mg_.08':{
			'inputs':{
						"Mg_e": 0.08 # mM
					},
			"expectations": {			
				"nTRPM": {
					'mean':[1.1], 
					'std': [.1]
				}
			}
		},
		'Mg_.8':{
			'inputs':{
						"Mg_e": 0.8 # mM
					},
			"expectations": {			
				"nTRPM": {
					'mean':[1], 
					'std': [.01]
				}
			}
		},
		'Mg_8':{
			'inputs':{
						"Mg_e": 8 # mM
					},
			"expectations": {			
				"nTRPM": {
					'mean': [2.5], 
					'std': [0.7]
				}
			}
		}
	},
	'Q21_TRPM':{# the effect of different Mg ion concentrations on TRPM
		'experiment_period':int(72*60/t2m), # minutes
		'measurement_scheme':{
			'TRPM': [int(72*60/t2m)]
		},
		'IDs': ['Mg_.08','Mg_.8','Mg_8'],
		'Mg_.08':{
			'inputs':{
						"Mg_e": 0.08 # mM
					},
			"expectations": {			
				"TRPM": {
					'mean':[0.85], 
					'std': [0.15]
				}
			}
		},
		'Mg_.8':{
			'inputs':{
						"Mg_e": 0.8 # mM
					},
			"expectations": {			
				"TRPM": {
					'mean':[1], 
					'std': [0]
				}
			}
		},
		'Mg_8':{
			'inputs':{
						"Mg_e": 8 # mM
					},
			"expectations": {			
				"TRPM": {
					'mean': [1.6], 
					'std': [0.3]
				}
			}
		}
	},
	'Q21_nM7CK':{# the effect of different Mg ion concentrations on M7CKs_n
		'experiment_period':int(72*60/t2m), # minutes
		'measurement_scheme':{
			'nM7CK': [int(72*60/t2m)]
		},
		'IDs': ['Mg_.08','Mg_.8','Mg_8'],
		'Mg_.08':{
			'inputs':{
						"Mg_e": 0.08 # mM
					},
			"expectations": {			
				"nM7CK": {
					'mean':[0.5], 
					'std': [0.1]
				}
			}
		},
		'Mg_.8':{
			'inputs':{
						"Mg_e": 0.8 # mM
					},
			"expectations": {			
				"nM7CK": {
					'mean':[1], 
					'std': [0]
				}
			}
		},
		'Mg_8':{
			'inputs':{
						"Mg_e": 8 # mM
					},
			"expectations": {			
				"nM7CK": {
					'mean': [2.45],
					'std': [0.6]
				}
			}
		}
	},
	'Q21_H3S10':{# the effect of different Mg ion concentrations on M7CKs_n
		'experiment_period':int(72*60/t2m), # minutes
		'measurement_scheme':{
			'pH3S10': [int(72*60/t2m)]
		},
		'IDs': ['Mg_.08','Mg_.8','Mg_8'],
		'Mg_.08':{
			'inputs':{
						"Mg_e": 0.08 # mM
					},
			"expectations": {			
				"pH3S10": {
					'mean':[0.35], 
					'std': [0.17]
				}
			}
		},
		'Mg_.8':{
			'inputs':{
						"Mg_e": 0.8 # mM
					},
			"expectations": {			
				"pH3S10": {
					'mean':[1], 
					'std': [0]
				}
			}
		},
		'Mg_8':{
			'inputs':{
						"Mg_e": 8 # mM
					},
			"expectations": {			
				"pH3S10": {
					'mean': [2.3], 
					'std': [0.7]
				}
			}
		}
	},
	'eq_IL8':{# equalibrium of IL8
		'experiment_period':24*int(60/t2m), # hours
		'measurement_scheme':{
			'IL8_n':  range_24h_60mStep
		},
		'IDs': ['Mg_0dot5'],
		'Mg_0dot5':{
			'inputs':{
						"Mg_e": 0.5 # mM
					},
			"expectations": {			
				"IL8_n": {
					'mean':[1 for i in range_24h_60mStep]
				}
			}
		}
	},
	'Q21_Mg_IL8':{# the effect of different Mg ion concentrations on IL8 and IL1b
		'experiment_period':72*int(60/t2m), # hours
		'measurement_scheme':{
			'IL8': [6*int(60/t2m),72*int(60/t2m)] #hours
		},
		'IDs': ['Mg_.08','Mg_8'],
		'Mg_.08':{
			'inputs':{
						"Mg_e": 0.08 # mM
					},
			"expectations": {			
				"IL8": {
					'mean':[10000,32000], 
					'std': [0,0]
				},
				# "IL1b": {
				# 	'mean':[1700,2950,2600,3000,3400], 
				# 	'std': [200,200,450,150,300]
				# 	# 'mean':[3000,3400] 
				# },
			}
		},
		'Mg_8':{
			'inputs':{
						"Mg_e": 8 # mM
					},
			"expectations": {			
				"IL8": {
					'mean':[16000,37000],
					'std': [0,0]
					# 'mean':[37000], 
				},
				# "IL1b": {
				# 	'mean':[1950,1200,1000,1700,1950], 
				# 	'std': [750,350,0,200,250]
				# 	# 'mean':[1000,1700,1950] 
				# },
			}	
		}
	}
}
def select_obs(studies):
	obs = {}
	for study in studies:
		obs[study] = observations[study]
	return obs
import json
with open('observations.json','w') as file:
	file.write(json.dumps(observations,indent=4))
