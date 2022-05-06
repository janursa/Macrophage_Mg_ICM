t2m = 1 # relationship between each simulation time and minute

range_3h_10mStep = list(range(int(10/t2m),3*int(60/t2m),int(10/t2m)))
range_24h_60mStep = list(range(60,24*int(60/t2m),int(60/t2m)))
range_48h_60mStep = list(range(12*int(60/t2m),48*int(60/t2m),int(60/t2m)))
range_12h_60mStep = list(range(3*int(60/t2m),12*int(60/t2m),int(60/t2m))) # span of 3 to 12 hs


packages = {
	'P11' : ['eq_mg','R05_nMg_f','R05_mg_n','Q21_Mg'], # mg entry and equalibrium
	'P12' : ['Q21_TRPMn','Q21_TRPM','Q21_M7CKn','Q21_eq_trpm'],

	'P13' : ['Q21_H3S10','Q21_eq_h3s10'], 
	'P1': ['eq_mg','R05_nMg_f','R05_mg_n']+['Q21_TRPMn','Q21_TRPM','Q21_M7CKn','Q21_eq_trpm']+['eq_mg','R05_nMg_f','R05_mg_n','Q21_Mg'], # M1 model
	
	'P21' : ['eq_IL8','M05_IT'],
	#'M18','M05_NFKBn'

	'P3' : ['S12_IKBa_mg','S12_NFKBn_mg','Q21_NFKBn_72h','Q21_Mg_IL8'], # Mg regulate
    #Q21_Mg_IL1b

}


observations = {	
	'M05_IT': { # IL8 regulate IRAK4 and TRAF6
        'IDs': ['ctr','100'],
        'activation': True,
        'experiment_period':int(2*60/t2m),
        'measurement_scheme': {
            'nIRAK4': [int(2*60/t2m)],
            'naTRAF6': [int(2*60/t2m)],
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nIRAK4':{'mean':[1],
                		'std':[0]
                		}, #normalized format
                'naTRAF6':{'mean':[1],
                		'std':[0]
                		}, #normalized format
            }
        },
       	'100': {
            'inputs': {
            	'IL8': 100*1000
            },
            'expectations': {
                'nIRAK4':{'mean':[3.895],
                		'std':[0]
                		}, #normalized format
                'naTRAF6':{'mean':[2.587],
                		'std':[0]
                		}, #normalized format
            }
        },
        
    },
    'M05_NFKBn': { # IL8 regulate NFKB
        'IDs': ['ctr','0dot1','1','10','100','1000'],
        'activation': True,
        'experiment_period':int(2*60/t2m),
        'measurement_scheme': {
            'nNFKB_n': [int(2*60/t2m)],
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nNFKB_n':{'mean':[1],
                		'std':[0]
                		}, #normalized format
            }
        },
        '0dot1': {
            'inputs': {
            	'IL8': 1*1000
            },
            'expectations': {
                'nNFKB_n':{'mean':[1.1],
                		'std':[0]
                		}, #normalized format
            }
        },
        '1': {
            'inputs': {
            	'IL8': 1*1000
            },
            'expectations': {
                'nNFKB_n':{'mean':[1.2],
                		'std':[0]
                		}, #normalized format
            }
        },
       	'10': {
            'inputs': {
            	'IL8': 10*1000
            },
            'expectations': {
                'nNFKB_n':{'mean':[2.4],
                		'std':[0]
                		}, #normalized format
            }
        },
        '100': {
            'inputs': {
            	'IL8': 100*1000
            },
            'expectations': {
                'nNFKB_n':{'mean':[4.2],
                		'std':[0]
                		}, #normalized format
            }
        },
       	'1000': {
            'inputs': {
            	'IL8': 1000*1000
            },
            'expectations': {
                'nNFKB_n':{'mean':[4.7],
                		'std':[0]
                		}, #normalized format
            }
        },
    },
	'M18': { # IL8 regulate IL4R, IFNGR, and IL-1b
        'IDs': ['ctr','0dot01','0dot1','1','10'],
        # 'IDs': ['ctr','10'],
        'activation': True,
        'experiment_period':int(24*60/t2m),
        'measurement_scheme': {
            'nIFNGR': [int(24*60/t2m)],
            'nIL4R': [int(24*60/t2m)],
            'nIL1b': [int(24*60/t2m)],
            'nIL10': [int(24*60/t2m)],
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nIFNGR':{'mean':[2212/2212],
                		'std':[0]}, #normalized format
                'nIL4R':{'mean':[6251/6251],
                		'std':[0]},
                'nIL1b':{'mean':[390.0/390.0],
                		'std':[0]},
                'nIL10':{'mean':[126.8/126.8],
                		'std':[0]},
            }
        },
       	'0dot01': {
            'inputs': {
                'IL8': 10 ,#pg/ml
            },
            'expectations': {
                'nIFNGR':{'mean':[4587/2212],
                		'std':[0]}, #normalized format
                'nIL4R':{'mean':[5462/6251],
                		'std':[0]},
                'nIL1b':{'mean':[434.4/390.0],
                		'std':[0]},
                'nIL10':{'mean':[137.9/126.8],
                		'std':[0]},
            }
        },
        '0dot1': {
            'inputs': {
                'IL8': 100 ,
            },
            'expectations': {
                'nIFNGR':{'mean':[4564/2212],
                		'std':[0]}, #normalized format
                'nIL4R':{'mean':[5359/6251],
                		'std':[0]},
                'nIL1b':{'mean':[466.2/390.0],
                		'std':[0]},
                'nIL10':{'mean':[161.7/126.8],
                		'std':[0]},
            }
        },
       	'1': {
            'inputs': {
                'IL8': 1000 ,
            },
            'expectations': {
                'nIFNGR':{'mean':[4438/2212],
                		'std':[0]}, #normalized format
                'nIL4R':{'mean':[5367/6251],
                		'std':[0]},
                'nIL1b':{'mean':[498.3/390.0],
                		'std':[0]},
                'nIL10':{'mean':[150/126.8],
                		'std':[0]},
            }
        },
       	'10': {
            'inputs': {
                'IL8': 10000 ,
            },
            'expectations': {
                'nIFNGR':{'mean':[4668/2212],
                		'std':[0]}, #normalized format
                'nIL4R':{'mean':[5323/6251],
                		'std':[0]},
                'nIL1b':{'mean':[516/390.0],
                		'std':[0]},
                'nIL10':{'mean':[127/126.8],
                		'std':[0]},
            }
        },
        
    },
    'S12_NFKBn_mg': { # mg influences NFkB level
        'IDs': ['ctr','Mg_2dot5'],
        'activation': True,
        'experiment_period':int(3*60/t2m),
        'measurement_scheme': {
            'nNFKB_n': [int(3*60/t2m)]
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nNFKB_n':{'mean':[1],
                		'std':[0]} #normalized format
            }
        },
        'Mg_2dot5': {
            'inputs': {
             'Mg_e': 2.5
            },
            'expectations': {
                'nNFKB_n':{'mean':[0.7],
                		'std':[0]} #normalized format
            }
        },
    },
	'S12_IKBa_mg': { # mg influences IkBa level
        'IDs': ['ctr','Mg_2dot5'],
        'activation': False,
        'experiment_period':int(3*60/t2m),
        'measurement_scheme': {
            'nIKB': [int(3*60/t2m)]
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nIKB':{'mean':[1],
                		'std':[0]} #normalized format
            }
        },
        'Mg_2dot5': {
            'inputs': {
                'Mg_e': 2.5 ,
            },
            'expectations': {
                'nIKB':{'mean':[1.32],
                		'std':[0]} #normalized format
            }
        }
    },
    'Q21_IkBa':{# the effect of different Mg ion concentrations on IKB
		'experiment_period':int(72*60/t2m), # minutes
		'activation':True,
		'measurement_scheme':{
			'IKB': [int(6*60/t2m),int(72*60/t2m)]
		},
		'IDs': ['Mg_.08','Mg_8'],
		'Mg_.08':{
			'inputs':{
						'Mg_e': 0.08 # mM
					},
			'expectations': {
				'IKB': {
					'mean':[0.05,1.7], 
					'std': [0,.1]
				}
			}
		},
		'Mg_8':{
			'inputs':{
						'Mg_e': 8 # mM
					},
			'expectations': {
				'IKB': {
					'mean': [0.7,1.67], 
					'std': [.1,.25]
				}
			}
		}
	},
    'Q21_IkBa_6h':{# the effect of different Mg ion concentrations on IKB
        'experiment_period':int(6*60/t2m), # minutes
        'activation':True,
        'measurement_scheme':{
            'IKB': [int(6*60/t2m)]
        },
        'IDs': ['Mg_.08','Mg_8'],
        'Mg_.08':{
            'inputs':{
                        'Mg_e': 0.08 # mM
                    },
            'expectations': {
                'IKB': {
                    'mean':[0.05],
                    'std': [0]
                }
            }
        },
        'Mg_8':{
            'inputs':{
                        'Mg_e': 8 # mM
                    },
            'expectations': {
                'IKB': {
                    'mean': [0.7],
                    'std': [.1]
                }
            }
        }
    },
    'Q21_IkBa_72h':{# the effect of different Mg ion concentrations on IKB
        'experiment_period':int(72*60/t2m), # minutes
        'activation':True,
        'measurement_scheme':{
            'IKB': [int(72*60/t2m)]
        },
        'IDs': ['Mg_.08','Mg_8'],
        'Mg_.08':{
            'inputs':{
                        'Mg_e': 0.08 # mM
                    },
            'expectations': {
                'IKB': {
                    'mean':[1.7],
                    'std': [0.1]
                }
            }
        },
        'Mg_8':{
            'inputs':{
                        'Mg_e': 8 # mM
                    },
            'expectations': {
                'IKB': {
                    'mean': [1.67],
                    'std': [.25]
                }
            }
        }
    },
    'Q21_NFKBn_72h':{# the effect of different Mg ion concentrations on NFKB after 72h
        'experiment_period':int(72*60/t2m), # minutes
        'activation':True,
        'measurement_scheme':{
            'nNFKB_n': [int(72*60/t2m)]
        },
        'IDs': ['ctr','Mg_8'],
        'ctr':{
            'inputs':{
                       
                    },
            'expectations': {
                'nNFKB_n': {
                    'mean':[1],
                    'std': [0]
                }
            }
        },
        'Mg_8':{
            'inputs':{
                        'Mg_e': 8 # mM
                    },
            'expectations': {
                'nNFKB_n': {
                    'mean': [1],
                    'std': [0]
                }
            }
        }
    },

	'R05_mg_n': { # mg extrusion
        'IDs': ['Mg_02'],
        'activation':False,
        'experiment_period':int(48*60/t2m),
        'measurement_scheme': {
            'nMg': [0]+range_48h_60mStep
        },
        'Mg_02': {
            'inputs': {
                'Mg_e': 0.02,
            },
            'expectations': {
                'nMg':{'mean':[1]+[.5 for i in range_48h_60mStep]} #normalized format
            }
        }
    },
    'S12_mg': { # total mg goes almost 4 times for 2.5 mM
        'IDs': ['Mg_2a5'],
        'experiment_period':1*60,
        'activation': False,
        'measurement_scheme': {
            'nMg': [0,1*60]
        },
        'Mg_2a5': {
            'inputs': {
            	'Mg_f' : 0.4,
                'Mg_e': 2.5,
            },
            'expectations': {
                'nMg':{'mean':[1,3.6]} #normalized format
            }
        }
    },
    'R05_nMg_f': { # normalized free mg for an increase in extracellular mg
        'IDs': ['Mg_19'],
        'activation':False,
        'experiment_period':int(12*60/t2m),
        'measurement_scheme': {
            'nMg_f': [0]+range_12h_60mStep,
            # 'Mg': [0]+range_12h_60mStep
        },
        'Mg_19': {
            'inputs': {
                'Mg_e': 19,
            },
            'expectations': {
                'nMg_f':{'mean':[1]+[1.4 for i in range_12h_60mStep]}, #normalized format
                # 'Mg':{'mean':[0,100]} #normalized format
            }
        }
    },
	'eq_mg': { # mg should stay in its equalibrium when there is no stimuli
        'IDs': ['ctr'],
        'activation': False,
        'experiment_period':int(24*60/t2m),
        'measurement_scheme': {
            'Mg_f': range_24h_60mStep,
            'Mg': range_24h_60mStep,
            'nMg_ATP': range_24h_60mStep,
            # 'Mg_IM_n': range_24h_60mStep
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'Mg_f':{'mean':[0.8 for i in range_24h_60mStep]},
                'Mg': {'mean':[18.5 for i in range_24h_60mStep]},
                'nMg_ATP': {'mean':[1 for i in range_24h_60mStep]},
                # 'Mg_IM_n': {'mean':[1 for i in range_24h_60mStep]}
            }
        }
    },
	'Q21_Mg': { # measurements on the intracellular Mg concentration for a given external Mg
        'experiment_period':int(3*60/t2m),
        'activation': False,
        'measurement_scheme': {
            'nMg_f': [0]+range_3h_10mStep # keep the stable concentration until 3 hours
        },
        'IDs': ['Mg_8'],
        'Mg_8': {
            'inputs': {
                'Mg_e': 8
            },
            'expectations': {
                'nMg_f': {
                    'mean': [1]+[1.25 for i in range_3h_10mStep],
                }
            }
        }
    },
    'Q21_eq_trpm':{# equalibrium in trpm related block
		'experiment_period':int(24*60/t2m), # minutes
		'activation':False,
		'measurement_scheme':{
			'nTRPM': range_24h_60mStep,
			'nTRPM_n': range_24h_60mStep,
			'nM7CK_n': range_24h_60mStep
		},
		'IDs': ['ctr'],
		'ctr':{
			'inputs':{
					},
			'expectations': {
				'nTRPM': {'mean':[1 for i in range_24h_60mStep]},
				'nTRPM_n': {'mean':[1 for i in range_24h_60mStep]},
				'nM7CK_n': {'mean':[1 for i in range_24h_60mStep]}
			}
		}
	},
	'Q21_eq_h3s10':{# equalibrium in h3s10 in related block
		'experiment_period':int(24*60/t2m), # minutes
		'activation':False,
		'measurement_scheme':{
			'npH3S10': range_24h_60mStep,
		},
		'IDs': ['ctr'],
		'ctr':{
			'inputs':{
					},
			'expectations': {
				'npH3S10': {'mean':[1 for i in range_24h_60mStep]}
			}
		}
	},
	'Q21_TRPMn':{# the effect of different Mg ion concentrations on TRPM_n
		'experiment_period':int(72*60/t2m), # minutes
		'activation':False,
		'measurement_scheme':{
			'nTRPM_n': [int(72*60/t2m)]
		},
        'IDs': ['Mg_.8','Mg_8'],
		'Mg_.8':{
			'inputs':{
						'Mg_e': 0.8 # mM
					},
			'expectations': {
				'nTRPM_n': {
					'mean':[1], 
					'std': [.01]
				}
			}
		},
		'Mg_8':{
			'inputs':{
						'Mg_e': 8 # mM
					},
			'expectations': {
				'nTRPM_n': {
					'mean': [2.5], 
					'std': [0.7]
				}
			}
		}
	},
	'Q21_TRPM':{# the effect of different Mg ion concentrations on TRPM
		'experiment_period':int(72*60/t2m), # minutes
		'activation': False,
		'measurement_scheme':{
			'nTRPM': [int(72*60/t2m)]
		},
        'IDs': ['Mg_.8','Mg_8'],
		'Mg_.8':{
			'inputs':{
						'Mg_e': 0.8 # mM
					},
			'expectations': {
				'nTRPM': {
					'mean':[1], 
					'std': [0]
				}
			}
		},
		'Mg_8':{
			'inputs':{
						'Mg_e': 8 # mM
					},
			'expectations': {
				'nTRPM': {
					'mean': [1.6], 
					'std': [0.3]
				}
			}
		}
	},
	'Q21_M7CKn':{# the effect of different Mg ion concentrations on M7CKs_n
		'experiment_period':int(72*60/t2m), # minutes
		'activation': False,
		'measurement_scheme':{
			'nM7CK_n': [int(72*60/t2m)]
		},
        'IDs': ['Mg_.8','Mg_8'],
		'Mg_.8':{
			'inputs':{
						'Mg_e': 0.8 # mM
					},
			'expectations': {
				'nM7CK_n': {
					'mean':[1], 
					'std': [0]
				}
			}
		},
		'Mg_8':{
			'inputs':{
						'Mg_e': 8 # mM
					},
			'expectations': {
				'nM7CK_n': {
					'mean': [2.45],
					'std': [0.6]
				}
			}
		}
	},
	'Q21_H3S10':{# the effect of different Mg ion concentrations on M7CKs_n
		'experiment_period':int(72*60/t2m), # minutes
		'activation': False,
		'measurement_scheme':{
			'npH3S10': [int(72*60/t2m)]
		},
#		'IDs': ['Mg_.08','Mg_.8','Mg_8'],
        'IDs': ['Mg_.8','Mg_8'],
		'Mg_.08':{
			'inputs':{
						'Mg_e': 0.08 # mM
					},
			'expectations': {
				'npH3S10': {
					'mean':[0.35], 
					'std': [0.17]
				}
			}
		},
		'Mg_.8':{
			'inputs':{
						'Mg_e': 0.8 # mM
					},
			'expectations': {
				'npH3S10': {
					'mean':[1], 
					'std': [0]
				}
			}
		},
		'Mg_8':{
			'inputs':{
						'Mg_e': 8 # mM
					},
			'expectations': {
				'npH3S10': {
					'mean': [2.3], 
					'std': [0.7]
				}
			}
		}
	},
	'eq_IL8':{# equalibrium of IL8
		'experiment_period':24*int(60/t2m), # hours
		'activation': False,
		'measurement_scheme':{
			'nIL8':  range_24h_60mStep,
			'nIL8R':  range_24h_60mStep
		},
		'IDs': ['ctr'],
		'ctr':{
			'inputs':{
						
					},
			'expectations': {
				'nIL8': {
					'mean':[1 for i in range_24h_60mStep]
				},
				'nIL8R': {
					'mean':[1 for i in range_24h_60mStep]
				}
			}
		}
	},
	
 'Q21_Mg_IL8':{# the effect of different Mg ion concentrations on IL8
        'experiment_period':72*int(60/t2m), # hours
        'activation': True,
        'cellType': 'THP1',
        'measurement_scheme':{
            'nIL8': [72*int(60/t2m)]
        },
        'IDs': ['ctr','Mg_8'],
        'ctr':{
            'inputs':{
                       
                    },
            'expectations': {
                'nIL8': {
                    'mean':[1],
                    'std': [.05]
                },
            }
        },
        'Mg_8':{
            'inputs':{
                        'Mg_e': 8 # mM
                    },
            'expectations': {
                'nIL8': {
                    'mean':[2.7],
                    'std': [0.7]
                },
            }
        }
    },
    'Q21_Mg_IL1b':{# the effect of different Mg ion concentrations on IL8
        'experiment_period':72*int(60/t2m), # hours
        'activation': True,
        'cellType': 'THP1',
        'measurement_scheme':{
            'nIL1b': [72*int(60/t2m)]
        },
        'IDs': ['ctr','Mg_8'],
        'ctr':{
            'inputs':{
                       
                    },
            'expectations': {
                'nIL1b': {
                    'mean':[1],
                    'std': [.02]
                },
            }
        },
        'Mg_8':{
            'inputs':{
                        'Mg_e': 8 # mM
                    },
            'expectations': {
                'nIL1b': {
                    'mean':[.5],
                    'std': [0.3]
                },
            }
        }
    },
    
 
    
}
def select_obs(studies_keys):
	obs = {}
	for study in studies_keys:
		obs[study] = observations[study]
	return obs
#import json
#with open('observations.json','w') as file:
#	file.write(json.dumps(observations,indent=4))
