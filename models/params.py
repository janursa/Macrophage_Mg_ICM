fixed_params_p = {
    # 'P1':{
    #     'k_m_i': 0.253, 
    #     'k_m_e': 0.255, 
    #     'k_atp_p0': 9.814,
    #     'k_matp_ub': 0.946,
    #     'k_matp_b': 0.5857,
    #     'k_atp_c': 0.62,
    #     'k_atp_pm' :2.33,
    #     'k_mim_b' :0.95,
    #     'k_mim_ub': 0.04,
    #     'n_matp_b': 0.02,
    #     'n_mim_b': 8.2,
    #     'IM': 88.61,
    #     'ATP': 16
    # },
    # 'P2':{
    #     'k_ntrpm_p': 0.11, 
    #     'n_ntrpm_p': 12.22, 
    #     'k_ntrpm_t': 0.961, 
    #     'n_ntrpm_t': 0.516, 
    #     'k_trpm_s': 0.955, 
    #     'k_m7ck_t1': 0.83, 
    #     'k_m7ck_t2': 0.96, 
    #     'n_m7ck_t': 1.90, 
    #     'k_m7ck_deg': 0.89
    # },
    # 'P3':{
    #     'k_h3s10_p': 0.0008, 
    #     'k_h3s10_up': 0.0007, 
    #     'n_h3s10_p': .8
    # },
    # 'P4':{}
}
free_params_p = {
    'P1':{
        ###---- Mg/Mg_ATP stabilization----###
        # not included
        'ATP':[1,100],
        'Mg_ATP':[5,15],
        'k_atp_p0':[0,10],
        'k_atp_pm':[0,10],
        'k_atp_c':[0,1],
        'IM':[0,100],
        'Mg_IM':[5,15],
        #------
        'k_m_i':[0.001,1],
        'k_m_e':[0.001,1],
        'k_matp_b':[0.001,1],
        'k_matp_ub':[0.001,1],
        'k_mim_b':[0.001,1],
        'k_mim_ub':[0.001,1],
        
        'n_matp_b':[1,10],
        'n_mim_b':[1,10],
        ###----------------------------------###
    },
    'P2':{
        'k_ntrpm_p':[0,1],
        'n_ntrpm_p':[0,30],
        'k_ntrpm_t':[0,1],
        'n_ntrpm_t':[0,10],
        'k_trpm_s':[0,1],
        
        'k_m7ck_t1':[0,1],
        'k_m7ck_t2':[0,1],
        'n_m7ck_t':[0,10],
        'k_m7ck_deg':[0,1],

    },
    'P3':{
        'k_h3s10_p':[0,1],
        'k_h3s10_up':[0,1],
        'n_h3s10_p':[0,10],

    },
    'P4':{
        'kd_ikb_d':[0.001,100],
        'n_ikb_d':[0.1,100],
        'kd_ikb_p': [0.001,100],
        'n_ikb_p': [0.1,30]
    }    
}
fixed_params = {}
for key,item in fixed_params_p.items():
    fixed_params = {**fixed_params,**item}
