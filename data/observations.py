import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
from tools import tools

t2m = 1 # relationship between each simulation time and minute

range_3h_10mStep = list(range(int(10/t2m),3*int(60/t2m),int(10/t2m)))
range_24h_60mStep = list(range(60,24*int(60/t2m),int(60/t2m)))
range_48h_60mStep = list(range(12*int(60/t2m),48*int(60/t2m),int(60/t2m)))
range_12h_60mStep = list(range(3*int(60/t2m),12*int(60/t2m),int(60/t2m))) # span of 3 to 12 hs


packages = {
	'M11' : ['eq_mg','R05_nMg_f','Q21_Mg'], # mg entry and equalibrium
	'M12' : ['Q21_M1','Q21_eq_trpm','Q21_eq_h3s10'],
	'M1': [],
	# R05_mg_n (mg extrusion)
	'IL8' : ['eq_IL8','M05_IT','M05_NFKBn'],
    # 'IL8' : ['eq_IL8','M05_IT'],
    # 'IL8' : ['M05_IT'],
	#'M18','M05_NFKBn'
    #Q21_Mg_IL1b
    #Q21_NFKBn_72h is not quantitative
    #B20_NFKBn
    ##F18_cytokines
    ##Z19_IL10
    #S12_LPS
    #Q21_14d
    'LPS':['S12_LPS','B20_LPS'],
    #B20_LPS
    # 'combined' : ['S12_IKBa_mg','Z19_IKB_NFKB','S12_NFKBn_mg','Q21_Mg_IL8','B20_NFKBn','eq_combined'], # Mg regulate
    'combined' : ['S12_IKBa_mg','S12_NFKBn_mg','B20_NFKBn','Q21_Mg_IL8','eq_combined'], # Mg regulate
    # 'IL6':['F17']
    'IL6':['eq_IL6','F17','N03','F14']
    #'F17'
    # ,'B17'
    # 'eq_IL6'

}


observations = {
    'F14': { # IL6 stimulates 
        'IDs': ['ctr','IL6_10'],
        'activation': False,
        'cellType':'HM',
        'duration':int(24*60/t2m),
        'selections': {
            'nIL10': [int(24*60/t2m)],
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nIL10':{'mean':[40/40],
                        'std':[15/40],
                        'pvalue':[None]},
            }
        },
        'IL6_10': {
            'inputs': {
                'IL6': 10*1000
            },
            'expectations': {
                'nIL10':{'mean':[75/40],
                        'std':[35/40],
                        'pvalue':[None]},
            }
        }
    },
    'N03': { # IL6 stimulates nTNFa
        'IDs': ['ctr','IL6_0dot1','IL6_10','IL6_100'],
        'activation': False,
        'cellType':'HM',
        'duration':int(24*60/t2m),
        'selections': {
            'nTNFa': [int(24*60/t2m)],
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nTNFa':{'mean':[100/100],
                        'std':[0],
                        'pvalue':[None]},
            }
        },
        'IL6_0dot1': {
            'inputs': {
                'IL6': 0.1*1000
            },
            'expectations': {
                'nTNFa':{'mean':[50/100],
                        'std':[0],
                        'pvalue':[None]},
            }
        },
        'IL6_10': {
            'inputs': {
                'IL6': 10*1000
            },
            'expectations': {
                'nTNFa':{'mean':[55/100],
                        'std':[0],
                        'pvalue':[None]},
            }
        },
        'IL6_100': {
            'inputs': {
                'IL6': 100*1000
            },
            'expectations': {
                'nTNFa':{'mean':[55/100],
                        'std':[0],
                        'pvalue':[None]},
            }
        },
    },
    'eq_IL6':{ # eq 
        'IDs': ['ctr'],
        'activation': False,
        'duration':int(24*60/t2m),
        'selections': {
            # 'npIL6_R': range_24h_60mStep,
            'F_stat3_a': range_24h_60mStep,
            'F_pi3k_a': range_24h_60mStep

        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'npIL6_R':{'mean':[1 for i in range_24h_60mStep]}, #normalized format
                'F_stat3_a': {'mean':[1 for i in range_24h_60mStep]},
                'F_pi3k_a': {'mean':[1 for i in range_24h_60mStep]},
            }
        },
    },
    'F17': { # IL6 stimulates STAT3 
        'IDs': ['ctr','IL6_50','IL6_100','IL6_200'],
        'activation': False,
        'cellType':'PBMCs',
        'duration':int(24*60/t2m),
        'selections': {
            'npSTAT3': [int(24*60/t2m)],
            'nIL10': [int(24*60/t2m)],
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'npSTAT3':{'mean':[0.25/0.25],
                        'std':[0.05/0.25],
                        'pvalue':[None]},
                'nIL10':{'mean':[0.5/0.5],
                        'std':[0/0.5],
                        'pvalue':[None]}
            }
        },
        'IL6_50': {
            'inputs': {
             'IL6': 50*1000
            },
            'expectations': {
                'npSTAT3':{'mean':[0.32/0.25],
                        'std':[0.05/0.25],
                        'pvalue':[0.01]},
                'nIL10':{'mean':[0.8/0.5],
                        'std':[0.1/0.5],
                        'pvalue':[0.01]}
            }
        },
        'IL6_100': {
            'inputs': {
             'IL6': 100*1000
            },
            'expectations': {
                'npSTAT3':{'mean':[0.38/0.25],
                        'std':[0.04/0.25],
                        'pvalue':[0.01]},
                'nIL10':{'mean':[1/0.5],
                        'std':[0.2/0.5],
                        'pvalue':[0.001]}
            }
        },
        'IL6_200': {
            'inputs': {
             'IL6': 200*1000
            },
            'expectations': {                
                'npSTAT3':{'mean':[0.4/0.25],
                        'std':[0.03/0.25],
                        'pvalue':[0.001]},
                'nIL10':{'mean':[1.7/0.5],
                        'std':[0.3/0.5],
                        'pvalue':[0.001]}
            }
        },
    },
    'B17': { # IL6 stimulates STAT3 and NFKB
        'IDs': ['ctr','stim'],
        'activation': False,
        'cellType':'THP1',
        'duration':int(2*60/t2m),
        'selections': {
            'npSTAT3': [int(2*60/t2m)],
            # 'nNFKB_n': [int(2*60/t2m)]
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'npSTAT3':{'mean':[1],
                        'std':[0],
                        'pvalue':[None]}, #normalized format
                'nNFKB_n':{'mean':[1],
                        'std':[0],
                        'pvalue':[None]}
            }
        },
        'stim': {
            'inputs': {
             'IL6': 60*1000
            },
            'expectations': {
                'npSTAT3':{'mean':[40],
                        'std':[2],
                        'pvalue':[0.001]}, #normalized format
                'nNFKB_n':{'mean':[4],
                        'std':[0.5],
                        'pvalue':[0.01]}
            }
        },
    },
    'eq_combined':{ # eq 
        'IDs': ['ctr'],
        'activation': False,
        'duration':int(24*60/t2m),
        'selections': {
            'F_h3s10_ikb': range_24h_60mStep,
            'F_p3s10_il8_p': range_24h_60mStep,

        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'F_h3s10_ikb':{'mean':[1 for i in range_24h_60mStep]}, #normalized format
                'F_p3s10_il8_p':{'mean':[1 for i in range_24h_60mStep]}
            }
        },
    },
    'B20_LPS': { # LPS stimulates cytokine production
        'IDs': ['ctr','stim'],
        'activation': False,
        'duration':int(13*24*60/t2m),
        'selections': {
            'nTNFa': [int(13*24*60/t2m)]
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nTNFa':{'mean':[1],
                        'std':[0],
                        'pvalue':[None]} #normalized format
            }
        },
        'stim': {
            'inputs': {
             'LPS': 10*1000,
             'IFNG':50*tools.c_2_ac['IFNG']
            },
            'expectations': {
                'nTNFa':{'mean':[200],
                        'std':[0],
                        'pvalue':[None]} #normalized format
            }
        },
    },
    'S12_LPS': { # mg influences NFkB level
        'IDs': ['LPS'],
        'activation': False,
        'duration':int(24*60/t2m),
        'selections': {
            'nIKB': [int(ii*60/t2m) for ii in [0,1,2,8]],
            'nTNFa': [int(ii*60/t2m) for ii in [0,1,2,4,8,24]]
        },
        'LPS': {
            'inputs': {
             'LPS': 50
            },
            'expectations': {
                'nIKB':{'mean':[1,15,8,4],
                        'std':[0,0,0,0],
                        'pvalue':['ns','ns','ns','ns']},
                'nTNFa':{'mean':[1,27,8,7,4,5],
                        'std':[0,0,0,0,0,0],
                        'pvalue':['ns','ns','ns','ns','ns','ns']}
                         #normalized format
            }
        },
    },
    'B20_TNFa': { # Mg regulate NFKB
        'IDs': ['ctr','Mg'],
        'activation': False,
        'duration':int(72*60/t2m),
        'selections': {
            'nTNFa': [int(72*60/t2m)],
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nTNFa':{'mean':[150/150],
                        'std':[abs(250 - 100)/1.35/150],
                        'pvalue':[None]
                        }, #normalized format
            }
        },
        'Mg': {
            'inputs': {
                'Mg_e': 5
            },
            'expectations': {
                'nTNFa':{'mean':[80/150],
                        'std':[abs(120 - 75)/1.35/150],
                        'pvalue':[0.01]
                        }, #normalized format
            }
        }
    },
    'B20_NFKBn': { # Mg regulate NFKB
        'IDs': ['ctr','Mg'],
        'activation': False,
        'duration':int(24*60/t2m),
        'selections': {
            'nNFKB_n': [int(24*60/t2m)],
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nNFKB_n':{'mean':[1],
                        'std':[0.8],
                        'pvalue':[None]
                        }, #normalized format
            }
        },
        'Mg': {
            'inputs': {
                'Mg_e': 5
            },
            'expectations': {
                'nNFKB_n':{'mean':[0.6],
                        'std':[0.4],
                        'pvalue':[0.01]
                        }, #normalized format
            }
        }
    },
    'F18_cytokines':{
        'duration':int(72*60/t2m), # minutes
            'activation':False,
            'selections':{
                'nIL10': [int(72*60/t2m)],
                'nTNFa': [int(72*60/t2m)],
                'nIL1b': [int(72*60/t2m)],
            },
            'IDs': ['ctr','Mg_5'],
            'ctr':{
                'inputs':{
                        },
                'expectations': {
                    'nIL10': {
                        'mean':[17/17],
                        'std': [5/17],
                        'pvalue':[None]
                    },
                    'nTNFa': {
                        'mean':[16000/16000],
                        'std': [1000/16000],
                        'pvalue':[None]
                    },
                    'nIL1b': {
                        'mean':[125/125],
                        'std': [10/125],
                        'pvalue':[None]
                    }
                }
            },
            'Mg_5':{
                'inputs':{
                            'Mg_e': 8 # mM
                        },
                'expectations': {
                    'nIL10': {
                        'mean':[25/17],
                        'std': [4/17],
                        'pvalue':[None]
                    },
                    'nTNFa': {
                        'mean':[10000/16000],
                        'std': [800/16000],
                        'pvalue':[None]
                    },
                    'nIL1b': {
                        'mean':[85/125],
                        'std': [5/125],
                        'pvalue':[0.01]
                    }
                }
            }
            
    },
    'Q21_cytokines_72h':{
        'duration':int(72*60/t2m), # minutes
            'activation':False,
            'cellType':'BMM',
            'selections':{
                'nIL10': [int(72*60/t2m)],
                'nTNFa': [int(72*60/t2m)],
                'nIL1b': [int(72*60/t2m)],
            },
            'IDs': ['ctr','Mg_8'],
            'ctr':{
                'inputs':{
                        },
                'expectations': {
                    'nIL10': {
                        'mean':[1],
                        'std': [0.05],
                        'pvalue':[None]
                    },
                    'nTNFa': {
                        'mean':[1],
                        'std': [0.1],
                        'pvalue':[None]
                    },
                    'nIL1b': {
                        'mean':[1],
                        'std': [0.05],
                        'pvalue':[None]
                    }
                }
            },
            'Mg_8':{
                'inputs':{
                            'Mg_e': 8 # mM
                        },
                'expectations': {
                    'nIL10': {
                        'mean':[2.5],
                        'std': [0.3],
                        'pvalue':[0.01]
                    },
                    'nTNFa': {
                        'mean':[0.35],
                        'std': [0.05],
                        'pvalue':[0.01]
                    },
                    'nIL1b': {
                        'mean':[.2],
                        'std': [0.04],
                        'pvalue':[0.01]
                    }
                }
            }
        
    },
    'Z19_IL10':{
    'duration':int(6*60/t2m), # minutes
        'activation':False,
        'selections':{
            'nIL10': [int(6*60/t2m)],
        },
        'IDs': ['ctr','Mg_4','Mg_20'],
        'ctr':{
            'inputs':{
                    },
            'expectations': {
                'nIL10': {
                    'mean':[2.8/2.8],
                    'std': [.2/2.5],
                    'pvalue':[None]
                }
            }
        },
        'Mg_4':{
            'inputs':{
                        'Mg_e': 4 # mM
                    },
            'expectations': {
                'nIL10': {
                    'mean': [4.9/2.8],
                    'std': [.3/2.8],
                    'pvalue':[0.01]
                }
            }
        },
        'Mg_20':{
            'inputs':{
                        'Mg_e': 20 # mM
                    },
            'expectations': {
                'nIL10': {
                    'mean': [3.7/2.8],
                    'std': [.2/2.8],
                    'pvalue':[0.01]
                }
            }
        }
    },
    'Z19_IKB_NFKB':{
        'duration':int(6*60/t2m), # minutes
            'activation':False,
            'selections':{
                'nIKB': [int(6*60/t2m)],
                'nNFKB_n': [int(6*60/t2m)]
            },
            'IDs': ['ctr','Mg_4','Mg_20'],
            'ctr':{
                'inputs':{
                        },
                'expectations': {
                    'nIKB': {
                        'mean':[0.4/0.4],
                        'std': [0.03],
                        'pvalue':[None]
                    },
                    'nNFKB_n': {
                        'mean':[0.2/0.2],
                        'std': [0.03],
                        'pvalue':[None]
                    }
                }
            },
            'Mg_4':{
                'inputs':{
                            'Mg_e': 4 # mM
                        },
                'expectations': {
                    'nIKB': {
                        'mean': [0.8/0.4],
                        'std': [0.05],
                        'pvalue':[0.01]
                    },
                    'nNFKB_n': {
                        'mean':[0.18/0.2],
                        'std': [0.05],
                        'pvalue':[None]
                    }
                }
            },
            'Mg_20':{
                'inputs':{
                            'Mg_e': 20 # mM
                        },
                'expectations': {
                    'nIKB': {
                        'mean': [0.55/0.4],
                        'std': [0.04],
                        'pvalue':[0.01]
                    },
                    'nNFKB_n': {
                        'mean':[0.19/0.2],
                        'std': [0.08],
                        'pvalue':[None]
                    }
                }
            }
        },
    
    'Q21_14d':{
        'duration':int(14*24*60/t2m), # minutes
        'activation':False,
        'selections':{
            'nIKB': [int(14*24*60/t2m)],
            'nNFKB_n': [int(14*24*60/t2m)]
        },
        'IDs': ['ctr','Mg_8'],
        'ctr':{
            'inputs':{
                    },
            'expectations': {
                'nIKB': {
                    'mean':[1],
                    'std': [0],
                    'pvalue':[None]
                },
                'nNFKB_n': {
                    'mean':[1],
                    'std': [0],
                    'pvalue':[None]
                }
            }
        },
        'Mg_8':{
            'inputs':{
                        'Mg_e': 8 # mM
                    },
            'expectations': {
                'nIKB': {
                    'mean': [1.25],
                    'std': [0.2],
                    'pvalue':[0.01]
                },
                'nNFKB_n': {
                    'mean':[1.45],
                    'std': [0.05],
                    'pvalue':[0.01]
                }
            }
        }
    },
	'M05_IT': { # IL8 regulate IRAK4 and TRAF6
        'IDs': ['ctr','100'],
        'activation': False,
        'duration':int(2*60/t2m),
        'selections': {
            'nIRAK4': [int(2*60/t2m)],
            'naTRAF6': [int(2*60/t2m)],
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nIRAK4':{'mean':[1],
                		'std':[0],
                        'pvalue':[None]
                		}, #normalized format
                'naTRAF6':{'mean':[1],
                		'std':[0],
                        'pvalue':[None]
                		}, #normalized format
            }
        },
       	'100': {
            'inputs': {
            	'IL8': 100*1000
            },
            'expectations': {
                'nIRAK4':{'mean':[3.895],
                		'std':[0],
                        'pvalue':['ns']
                		}, #normalized format
                'naTRAF6':{'mean':[2.587],
                		'std':[0],
                        'pvalue':['ns']
                		}, #normalized format
            }
        },
        
    },
    'M05_NFKBn': { # IL8 regulate NFKB
        'IDs': ['ctr','0dot1','1','10','100','1000'],
        'activation': False,
        'duration':int(2*60/t2m),
        'selections': {
            'nNFKB_n': [int(2*60/t2m)],
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nNFKB_n':{'mean':[1],
                		'std':[0],
                        'pvalue':[None]
                		}, #normalized format
            }
        },
        '0dot1': {
            'inputs': {
            	'IL8': 1*1000
            },
            'expectations': {
                'nNFKB_n':{'mean':[1.1],
                		'std':[0],
                        'pvalue':['ns']
                		}, #normalized format
            }
        },
        '1': {
            'inputs': {
            	'IL8': 1*1000
            },
            'expectations': {
                'nNFKB_n':{'mean':[1.2],
                		'std':[0],
                        'pvalue':['ns']
                		}, #normalized format
            }
        },
       	'10': {
            'inputs': {
            	'IL8': 10*1000
            },
            'expectations': {
                'nNFKB_n':{'mean':[2.4],
                		'std':[0],
                        'pvalue':['ns']
                		}, #normalized format
            }
        },
        '100': {
            'inputs': {
            	'IL8': 100*1000
            },
            'expectations': {
                'nNFKB_n':{'mean':[4.2],
                		'std':[0],
                        'pvalue':['ns']
                		}, #normalized format
            }
        },
       	'1000': {
            'inputs': {
            	'IL8': 1000*1000
            },
            'expectations': {
                'nNFKB_n':{'mean':[4.7],
                		'std':[0],
                        'pvalue':['ns']
                		}, #normalized format
            }
        },
    },
	'M18': { # IL8 regulate IL4R, IFNGR, and IL-1b
        'IDs': ['ctr','0dot01','0dot1','1','10'],
        # 'IDs': ['ctr','10'],
        'activation': False,
        'duration':int(24*60/t2m),
        'selections': {
            # 'nIFNGR': [int(24*60/t2m)],
            # # 'nIL4R': [int(24*60/t2m)],
            'nIL1b': [int(24*60/t2m)],
            'nIL10': [int(24*60/t2m)],
            'nTNFa': [int(24*60/t2m)],
            'nIL6': [int(24*60/t2m)],
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nIFNGR':{'mean':[2212/2212],
                		'std':[abs(1913 - 3145)/1.35/2212],
                        'pvalue':[None]}, #normalized format
                'nIL4R':{'mean':[6251/6251],
                		'std':[abs(4517-8834)/1.35/6251],
                        'pvalue':[0.01]},
                'nIL1b':{'mean':[390.0/390.0],
                		'std':[269/390.0],
                        'pvalue':[None]},
                'nIL10':{'mean':[126/126],
                		'std':[abs(95.9-152.9)/1.35/126],
                        'pvalue':[None]},
                'nTNFa':{'mean':[1114.3/1114.3],
                    'std':[abs(722.2-1586.4)/1.35/1114],
                    'pvalue':[None]},
                'nIL6':{'mean':[14300/14300],
                    'std':[abs(9382-16584)/1.35/14300],
                    'pvalue':[None]}
            }
        },
       	'0dot01': {
            'inputs': {
                'IL8': 10 ,#pg/ml
            },
            'expectations': {
                'nIFNGR':{'mean':[4587/2212],
                		'std':[abs(3452-4625)/1.35/2212],
                        'pvalue':[0.01]}, #normalized format
                'nIL4R':{'mean':[5462/6251],
                		'std':[abs(3181-7475)/1.35/6251],
                        'pvalue':[0.01]},
                'nIL1b':{'mean':[434.4/390.0],
                		'std':[354/390.0],
                        'pvalue':[None]},
                'nIL10':{'mean':[137.9/126.8],
                		'std':[abs(103.1-195.0)/1.35/126],
                        'pvalue':[None]},
                'nTNFa':{'mean':[1223.3/1114.3],
                        'std':[abs(555.0-1636.0)/1.35/1114],
                    'pvalue':[None]},
                'nIL6':{'mean':[17850/14300],
                    'std':[abs(11850-21251)/1.35/14300],
                    'pvalue':[0.01]}
            }
        },
        '0dot1': {
            'inputs': {
                'IL8': 100 ,
            },
            'expectations': {
                'nIFNGR':{'mean':[4564/2212],
                		'std':[abs(3785-4595)/1.35/2212],
                        'pvalue':[0.01]}, #normalized format
                'nIL4R':{'mean':[5359/6251],
                		'std':[abs(3206-7697)/1.35/6251],
                        'pvalue':[0.01]},
                'nIL1b':{'mean':[466.2/390.0],
                		'std':[321/390.0],
                        'pvalue':[None]},
                'nIL10':{'mean':[161.7/126.8],
                		'std':[abs(96.8-182.0)/1.35/126],
                        'pvalue':[None]},
                'nTNFa':{'mean':[1199.9/1114.3],
                    'std':[abs(806.6-1666.4)/1.35/1114],
                    'pvalue':[None]},
                'nIL6':{'mean':[16182/14300],
                    'std':[abs(10776-20545)/1.35/14300],
                    'pvalue':[None]}
                
            }
        },
       	'1': {
            'inputs': {
                'IL8': 1000 ,
            },
            'expectations': {
                'nIFNGR':{'mean':[4438/2212],
                		'std':[abs(3248-4578)/1.35/2212],
                        'pvalue':[0.01]}, #normalized format
                'nIL4R':{'mean':[5367/6251],
                		'std':[abs(3136-7625)/1.35/6251],
                        'pvalue':[0.01]},
                'nIL1b':{'mean':[498.3/390.0],
                		'std':[274/390.0],
                        'pvalue':[None]},
                'nIL10':{'mean':[150/126.8],
                		'std':[abs(91.1-166.4)/1.35/126],
                        'pvalue':[None]},
                'nTNFa':{'mean':[1138.9/1114.3],
                        'std':[abs(780.3-1723.7)/1.35/1114],
                        'pvalue':[None]},
                'nIL6':{'mean':[15120/14300],
                    'std':[abs(12115-19545)/1.35/14300],
                    'pvalue':[0.01]}
            }
        },
       	'10': {
            'inputs': {
                'IL8': 10000 ,
            },
            'expectations': {
                'nIFNGR':{'mean':[4668/2212],
                		'std':[abs(3589-4685)/1.35/2212],
                        'pvalue':[0.01]}, #normalized format
                'nIL4R':{'mean':[5323/6251],
                		'std':[abs(3155-7715)/1.35/6251],
                        'pvalue':[0.01]},
                'nIL1b':{'mean':[516/390.0],
                		'std':[258/390.0],
                        'pvalue':[0.01]},
                'nIL10':{'mean':[127/126.8],
                		'std':[abs(90.2-151.7)/1.35/126],
                        'pvalue':[None]},
                'nTNFa':{'mean':[1136.1/1114.3],
                        'std':[abs(665.5-1346.7)/1.35/1114],
                        'pvalue':[None]},
                'nIL6':{'mean':[12554/14300],
                    'std':[abs(10212-16777)/1.35/14300],
                    'pvalue':[None]}
            }
        },
        
    },
    'S12_NFKBn_mg': { # mg influences NFkB level
        'IDs': ['ctr','Mg_2dot5'],
        'activation': False,
        'duration':int(3*60/t2m),
        'selections': {
            'nNFKB_n': [int(3*60/t2m)]
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nNFKB_n':{'mean':[1],
                		'std':[0],
                        'pvalue':[None]} #normalized format
            }
        },
        'Mg_2dot5': {
            'inputs': {
             'Mg_e': 2.5
            },
            'expectations': {
                'nNFKB_n':{'mean':[0.7],
                		'std':[0],
                        'pvalue':['ns']} #normalized format
            }
        },
    },
	'S12_IKBa_mg': { # mg influences IkBa level
        'IDs': ['ctr','Mg_2dot5'],
        'activation': False,
        'duration':int(3*60/t2m),
        'selections': {
            'nIKB': [int(3*60/t2m)]
        },
        'ctr': {
            'inputs': {
            },
            'expectations': {
                'nIKB':{'mean':[1],
                		'std':[0],
                        'pvalue':[None]} #normalized format
            }
        },
        'Mg_2dot5': {
            'inputs': {
                'Mg_e': 2.5 ,
            },
            'expectations': {
                'nIKB':{'mean':[1.32],
                		'std':[0],
                        'pvalue':['ns']} #normalized format
            }
        }
    },
    'Q21_IkBa':{# the effect of different Mg ion concentrations on IKB
		'duration':int(72*60/t2m), # minutes
		'activation':False,
		'selections':{
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
        'duration':int(6*60/t2m), # minutes
        'activation':False,
        'selections':{
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
        'duration':int(72*60/t2m), # minutes
        'activation':False,
        'selections':{
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
        'duration':int(72*60/t2m), # minutes
        'activation':False,
        'selections':{
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
        'duration':int(48*60/t2m),
        'selections': {
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
        'duration':1*60,
        'activation': False,
        'selections': {
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
        'duration':int(12*60/t2m),
        'selections': {
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
        'duration':int(24*60/t2m),
        'selections': {
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
        'duration':int(3*60/t2m),
        'activation': False,
        'selections': {
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
		'duration':int(24*60/t2m), # minutes
		'activation':False,
		'selections':{
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
		'duration':int(24*60/t2m), # minutes
		'activation':False,
		'selections':{
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
	'Q21_M1':{# the effect of different Mg ion concentrations on TRPM_n and the subsequent impact on H3P10
		'duration':int(72*60/t2m), # minutes
		'activation':False,
		'selections':{
			'nTRPM_n': [int(72*60/t2m)],
            'nTRPM': [int(72*60/t2m)],
            'nM7CK_n': [int(72*60/t2m)],
            'npH3S10': [int(72*60/t2m)],
            'nATP':[int(72*60/t2m)]

		},
        'IDs': ['ctr','Mg_8'],
		'ctr':{
			'inputs':{
						
					},
			'expectations': {
				'nTRPM_n': {
					'mean':[1], 
					'std': [.01],
                    'pvalue':[None]
				},
                'nTRPM': {
                    'mean':[1], 
                    'std': [0],
                    'pvalue':[None]
                },
                'nM7CK_n': {
                    'mean':[1], 
                    'std': [0],
                    'pvalue':[None]
                },
                'npH3S10': {
                    'mean':[1], 
                    'std': [0],
                    'pvalue':[None]
                },
                'nATP': {
                    'mean':[600/600], 
                    'std': [20/600],
                    'pvalue':[None]
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
					'std': [0.7],
                    'pvalue':[0.01]
				},
                'nTRPM': {
                    'mean': [1.6], 
                    'std': [0.3],
                    'pvalue':[0.01]
                },
                'nM7CK_n': {
                    'mean': [2.45],
                    'std': [0.6],
                    'pvalue':[0.01]
                },
                'npH3S10': {
                    'mean': [2.3], 
                    'std': [0.7],
                    'pvalue':[0.01]
                },
                'nATP': {
                    'mean':[750/600], 
                    'std': [40/600],
                    'pvalue':[0.01]
                }
			}
		}
	},
	
	'eq_IL8':{# equalibrium of IL8
		'duration':24*int(60/t2m), # hours
		'activation': False,
		'selections':{
			'nIL8_m':  range_24h_60mStep,
			'F_il8_irak':  range_24h_60mStep,
            # 'nIFNGR':  range_24h_60mStep
		},
		'IDs': ['ctr'],
		'ctr':{
			'inputs':{
						
					},
			'expectations': {
				'nIL8_m': {
					'mean':[1 for i in range_24h_60mStep]
				},
				'F_il8_irak': {
					'mean':[1 for i in range_24h_60mStep]
				},
			}
		}
	},
	
 'Q21_Mg_IL8':{# the effect of different Mg ion concentrations on IL8
        'duration':72*int(60/t2m), # hours
        'activation': False,
        'cellType': 'THP1',
        'selections':{
            'nIL8': [72*int(60/t2m)]
        },
        'IDs': ['ctr','Mg_8'],
        'ctr':{
            'inputs':{
                       
                    },
            'expectations': {
                'nIL8': {
                    'mean':[1],
                    'std': [.05],
                    'pvalue':[None]
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
                    'std': [0.7],
                    'pvalue':[0.01]
                },
            }
        }
    },
    'Q21_Mg_IL1b':{# the effect of different Mg ion concentrations on IL8
        'duration':72*int(60/t2m), # hours
        'activation': False,
        'cellType': 'THP1',
        'selections':{
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
