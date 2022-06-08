fixed_params_p = {
#     'P1':{
#            ### P12
#            "ATP": 4.4200492287922275,
#            "k_atp_p0": 0.002255124807943189,
#            "k_atp_pm": 8.491754273290441,
#            "k_atp_c": 0.590124477903334,
#            "IM": 1.0199739743020637,
#            "Mg_ATP": 17.57597150605035,
#            "Mg_IM": 0.12931885964186912,
#            "k_m_i": 0.35411879710555166,
#            "k_matp_b": 0.8632245228883954,
#            "k_matp_ub": 0.45252482532151195,
#            "k_mim_b": 0.15981115222635706,
#            "k_mim_ub": 0.9897381814819478
            ### P12
#            "k_ntrpm_p": 0.0010111311947780088,
#            "n_ntrpm_p": 2.303172314228545,
#            "k_ntrpm_t": 0.7410061040831892,
#            "k_trpm_s": 0.9968437547905072,
#            "k_m7ck_t1": 0.9219270599665684,
#            "k_m7ck_t2": 0.7127451109467265,
#            "k_m7ck_deg": 0.9623775964715906
            ### P13
#            "k_h3s10_p": 7.035444348793543,
#            "k_h3s10_a": 0.9589420876901643,
#            "k_h3s10_d": 0.7221697589064646,
#            "k_h3s10_da": 0.807474039427079,
#            "n_h3s10_a": 0.7377166985974215
#     },
#     'P2':{
#             "IL8": 7.408568221151385,
#             "k_il8_p": 5.970462870815936,
#             "kd_nfkb_il8_p": 2085.514347894692,
#             "k_il8_d": 0.37595294602525076,
#             "IL8R": 4800.422519863175,
#             "IL8_R": 7.689186472021902,
#             "k_il8_b": 0.0011844646079401944,
#             "k_il8_ub": 0.7266871178144181,
#             "k_il8r_p": 917.4566246712063,
#             "k_il8r_d": 0.19039877492259133,
#             "kd_il8_irak_p": 33.35345100816778
#         }
    
    
}
free_params_p = {  
    'M1':{
        ### P11
        'ATP':[.1,100],
        'k_atp_p0':[0.001,10],
        'k_atp_pm':[0.001,10],
        'k_atp_c':[.001,1],
        'IM':[.1,100],
        
        'Mg_ATP':[.1,100],
        'Mg_IM':[.1,100],
        'k_m_i':[0.001,1],
        'k_matp_b':[0.001,1],
        'k_matp_ub':[0.001,1],
        'k_mim_b':[0.001,1],
        'k_mim_ub':[0.001,1],
        ###P12
        'k_ntrpm_p':[0.001,1],
        'n_ntrpm_p':[.01,10],
        'k_ntrpm_t':[0.001,1],
        'k_trpm_s':[0.001,1],

        'k_m7ck_t1':[0.001,1],
        'k_m7ck_t2':[0.001,1],
        'k_m7ck_deg':[0.001,1],
        
        'k_h3s10_p':[0.01,10],
        'k_h3s10_a':[0.01,1],
        'k_h3s10_d':[0.01,1],
        'k_h3s10_da':[0.01,1],
    
    },
    'LPS':{
        'k_lps_a':[0,1],
        'k_lps_d':[0,1],
        'k_lpsa_d':[0,1],
        'k_lps':[.1,1000000],
        'kd_lps':[0.001,100000],
    },
    'IL6':{
        ## IL6
        'IL6' : [1,10], 
        'IL6_R_JACK' : [1,10], 
        'pIL6_R_JACK': [1,10], 
        'k_il6r_b' :[0,1],
        'k_il6r_ub':[0,1],
        'k_il6r_a':[0,1],
        'k_il6r_da':[0,1],
        'IL6_m' : [1,10],
        'k_il6m_il6' : [0,1],
        'k_il6_p' : [0,100],
        'k_il6_stat3_a':[1,100],
        'kd_il6_stat3_a':[100,10000],
        'o_il6_stat3_a':[0,1],
        'k_il6_pi3k_a':[1,1000],
        'kd_il6_pi3k_a':[1,10000],
        'o_il6_pi3k_a':[0,1],
        
        ## regulated IL6 production
        'k_nfkb_il6_p' : [1,100],
        'kd_nfkb_il6_p' : [1,100000],
        'o_nfkb_il6_p' : [0,1],
        'k_ap1_il6_p' : [1,100],
        'kd_ap1_il6_p' : [1,100000],
        'o_ap1_il6_p' : [0,1],
        'k_creb_il6_p' : [1,100],
        'kd_creb_il6_p' : [1,100000],
        'o_creb_il6_p' : [0,1],
    },
    'IL8':{
        ## IL8
        'IL8':[0,10],
        'IL8_m':[0,10],
        'k_il8_p':[0,100],
        'k_il8m_il8':[0,1],
        'IL8_R': [0.1,10],
        'pIL8_R': [0.1,10],
        'k_il8r_b': [0,1],
        'k_il8r_ub': [0,1],
        'k_il8r_a': [0,1],
        'k_il8r_da': [0,1],
        'kd_il8_irak_p': [1000,100000],
        'k_il8_irak_p': [1000,100000],
        'o_il8_irak_p': [0,1],
        ## regulated production of IL8
        'k_nfkb_il8_p':[1,1000],
        'kd_nfkb_il8_p':[1,1000],
        'o_nfkb_il8_p': [0,1],
        'k_ap1_il8_p': [1,1000],
        'kd_ap1_il8_p':[1,1000],
        'o_ap1_il8_p':[0,1],
        
        # 'k_rho_a':[1,1000],
        # 'kd_rho_a':[100,10000],
        # 'k_rho_da':[0,1],
        # 'k_rho_nfkb_a':[100,100000],
        # 'kd_rho_nfkb_a':[100,100000],
        # 'o_rho_nfkb_a':[0,1],
        # 'k_rho_nfkb_a':[1,1000],
        # 'kd_rho_nfkb_a':[1,10000],
        # 'o_rho_nfkb_a':[0,1],
        
    },
    # 'ILs':{
        
    #     'k_rho_pi3k_a':[1000,100000],
    #     'kd_rho_pi3k_a':[1000,100000],
    #     'o_rho_pi3k_a':[0,1],
    #     'k_rho_stat3_a':[1,10000],
    #     'kd_rho_stat3_a':[1,100000],
    #     'o_rho_stat3_a':[0,1],

    #     ## regulated IL6 production
    #     'k_nfkb_il6_p' : [1,100],
    #     'kd_nfkb_il6_p' : [1,100000],
    #     'o_nfkb_il6_p' : [0,1],
    #     'k_ap1_il6_p' : [1,100],
    #     'kd_ap1_il6_p' : [1,100000],
    #     'o_ap1_il6_p' : [0,1],
    #     'k_creb_il6_p' : [1,100],
    #     'kd_creb_il6_p' : [1,100000],
    #     'o_creb_il6_p' : [0,1],
    # },
    'ILs_1':{
        # 'IL6' : [1,10], 
        # 'IL6_R_JACK' : [1,10], 
        # 'pIL6_R_JACK': [1,10], 
        # 'k_il6r_b' :[0,1],
        # 'k_il6r_ub':[0,1],
        # 'k_il6r_a':[0,1],
        # 'k_il6r_da':[0,1],
        # 'IL6_m' : [1,10],
        # 'k_il6m_il6' : [0,1],
        # 'k_il6_p' : [0,100],
        # ## for IRAK and NFKB
        'k_il6_stat3_a':[1,10000],
        'kd_il6_stat3_a':[1,10000],
        'o_il6_stat3_a':[0,1],
        'k_il6_pi3k_a':[1,10000],
        'kd_il6_pi3k_a':[1,10000],
        'o_il6_pi3k_a':[0,1],
        
        ##  IL6 transcriptional factors
        'k_nfkb_il6_p' : [1,1000],
        'kd_nfkb_il6_p' : [1,100000],
        'o_nfkb_il6_p' : [0,1],
        'k_ap1_il6_p' : [1,1000],
        'kd_ap1_il6_p' : [1,100000],
        'o_ap1_il6_p' : [0,1],
        'k_creb_il6_p' : [1,100],
        'kd_creb_il6_p' : [1,100000],
        'o_creb_il6_p' : [0,1],  
        
        ## global factors
        '_pSTAT3': [1,100], 
        'k_pstat3_d': [0,1],
        'k_pstat3_b': [0,1],
        'k_pstat3_ub': [0,1], 
        '_pPI3K': [1,100],
        'k_ppi3k_d':[0,1],
        'K_ppi3k_a':[0,1],
        'K_ppi3k_da':[0,1]
    },
    'ILs_2':{
        ## IL8
        'IL8':[0,10],
        'IL8_m':[0,10],
        'k_il8_p':[0,100],
        'k_il8m_il8':[0,1],
        'IL8_R': [0.1,10],
        'pIL8_R': [0.1,10],
        'k_il8r_b': [0,1],
        'k_il8r_ub': [0,1],
        'k_il8r_a': [0,1],
        'k_il8r_da': [0,1],
        'kd_il8_irak_p': [1,1000000],
        'k_il8_irak_p': [1,1000000],
        'o_il8_irak_p': [0,1],

        ## regulated production of IL8
        'k_nfkb_il8_p':[1,1000],
        'kd_nfkb_il8_p':[1,1000],
        'o_nfkb_il8_p': [0,1],
        'k_ap1_il8_p': [1,1000],
        'kd_ap1_il8_p':[1,1000],
        'o_ap1_il8_p':[0,1],
        
        'k_rho_a':[1,1000],
        'kd_rho_a':[100,10000],
        'k_rho_da':[0,1],
        'k_rho_nfkb_a':[100,100000],
        'kd_rho_nfkb_a':[100,100000],
        'o_rho_nfkb_a':[0,1],
        'k_rho_pi3k_a':[1,100000],
        'kd_rho_pi3k_a':[1,100000],
        'o_rho_pi3k_a':[0,1],
        'k_rho_stat3_a':[1,10000],
        'kd_rho_stat3_a':[1,100000],
        'o_rho_stat3_a':[0,1],    
        ## global factors
        '_pSTAT3': [1,100], 
        'k_pstat3_d': [0,1],
        'k_pstat3_b': [0,1],
        'k_pstat3_ub': [0,1], 
        '_pPI3K': [1,100],
        'k_ppi3k_d':[0,1],
        'K_ppi3k_a':[0,1],
        'K_ppi3k_da':[0,1]
    },
    'ILs' :{
        # 'IL6' : [1,10], 
        # 'IL6_R_JACK' : [1,10], 
        # 'pIL6_R_JACK': [1,10], 
        # 'k_il6r_b' :[0,1],
        # 'k_il6r_ub':[0,1],
        # 'k_il6r_a':[0,1],
        # 'k_il6r_da':[0,1],
        # 'IL6_m' : [1,10],
        # 'k_il6m_il6' : [0,1],
        # 'k_il6_p' : [0,10],
        # ## for IRAK and NFKB
        'k_il6_stat3_a':[100,10000],
        'kd_il6_stat3_a':[100,10000],
        'o_il6_stat3_a':[0,1],
        'k_il6_pi3k_a':[100,10000],
        'kd_il6_pi3k_a':[100,10000],
        'o_il6_pi3k_a':[0,1],
        
        ##  IL6 transcriptional factors
        'k_nfkb_il6_p' : [10,1000],
        'kd_nfkb_il6_p' : [10000,100000],
        'o_nfkb_il6_p' : [0,1],
        'k_ap1_il6_p' : [100,1000],
        'kd_ap1_il6_p' : [1000,100000],
        'o_ap1_il6_p' : [0,1],
        'k_creb_il6_p' : [1,100],
        'kd_creb_il6_p' : [1,100000],
        'o_creb_il6_p' : [0,1], 
        ## IL8
        # 'IL8':[0,10],
        # 'IL8_m':[0,10],
        # 'k_il8_p':[0,100],
        # 'k_il8m_il8':[0,1],
        # 'IL8_R': [0.1,10],
        # 'pIL8_R': [0.1,10],
        # 'k_il8r_b': [0,1],
        # 'k_il8r_ub': [0,1],
        # 'k_il8r_a': [0,1],
        # 'k_il8r_da': [0,1],
        
        'kd_il8_irak_p': [10000,1000000],
        'k_il8_irak_p': [10000,1000000],
        'o_il8_irak_p': [0,1],

        ## regulated production of IL8
        'k_nfkb_il8_p':[1,1000],
        'kd_nfkb_il8_p':[1,1000],
        'o_nfkb_il8_p': [0,1],
        'k_ap1_il8_p': [1,1000],
        'kd_ap1_il8_p':[1,1000],
        'o_ap1_il8_p':[0,1],
        
        # 'k_rho_a':[1,1000],
        # 'kd_rho_a':[1000,10000],
        # 'k_rho_da':[0,1],
        # 
        'k_rho_nfkb_a':[1000,100000],
        'kd_rho_nfkb_a':[1000,100000],
        'o_rho_nfkb_a':[0,1],
        'k_rho_pi3k_a':[1000,100000],
        'kd_rho_pi3k_a':[1000,100000],
        'o_rho_pi3k_a':[0,1],
        'k_rho_stat3_a':[1000,10000],
        'kd_rho_stat3_a':[1000,100000],
        'o_rho_stat3_a':[0,1], 

        ## global factors
        '_pSTAT3': [1,100], 
        'k_pstat3_d': [0,1],
        'k_pstat3_b': [0,1],
        'k_pstat3_ub': [0,1], 
        '_pPI3K': [1,100],
        'k_ppi3k_d':[0,1],
        'K_ppi3k_a':[0,1],
        'K_ppi3k_da':[0,1]
    },
    

    'combined':{
        ## regulated production of IL8
        'k_nfkb_il8_p':[1,1000],
        'kd_nfkb_il8_p':[1,100000],
        'o_nfkb_il8_p': [0,1],
        'k_ap1_il8_p': [1,1000],
        'kd_ap1_il8_p':[1,100000],
        'o_ap1_il8_p':[0,1],

        'kd_mg_ikb_d':[0.001,1],
        
        'kd_h3s10_ikb_p' : [0.1,10000],
        'k_h3s10_ikb_p' : [0.1,1000],


        'kd_h3s10_il8_p': [1,10000],
        'k_h3s10_il8_p': [0.1,1000],

        'o_h3s10_ikb': [0,1],
        'o_h3s10_il8': [0,1],

    } 
}
fixed_params = {}
for key,item in fixed_params_p.items():
    fixed_params = {**fixed_params,**item}
