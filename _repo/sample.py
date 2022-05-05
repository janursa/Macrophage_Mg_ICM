##-----sample from the original model --------###
if False:
    Zhao_2021 = te.loadSBMLModel(dir_Zhao_model)
    inputs = {'pIKK':[50]}
    targets = {'NFKB_n':[i for i in range(0,500,10)]}
    target_keys = list(targets.keys())
    simulation_duration = 1500
    if len(list(inputs.keys()))>1:
        raise ValueError('the rest of the code is designed for 1 key input')
    results = []
    for key, values in inputs.items():
    #     rr_dict_i={'inputs'}
        for value in values:
            results_i = run_Zhao_model(target_keys=target_keys,params = {key:value},duration=simulation_duration,step=simulation_duration)
            rr_dict = {}
            for target_key,target_time in targets.items():
                results_key = list([results_i[target_key][i] for i in target_time])
                rr_dict[target_key] = {'time':target_time,'values':results_key}
            results.append({'inputs': {key:value},'results':rr_dict})
    with open(dir_samples_Zhao_model, 'w', encoding='utf-8') as f:
        json.dump(results, f, ensure_ascii=False, indent=4)

if False:
    ###------calibrate the model according on the samples taken from the Zhao model ----##
    ## load the samples from the original model
    with open(dir_samples_Zhao_model) as json_file:
        samples = json.load(json_file)
    ## load the model
    model = te.loadSBMLModel(dir_model)
    ## define the free parameters
    free_params = dict( 
        IL8 = [0,100000], 
        K248 = [0,1],  # IL8 regulates TRAF6
        Kd248 = [0,1000], # IL8 regulates TRAF6
    )
    calib_Obj = Calibrate(model = model,free_params=free_params,target=samples,max_iteration=50)
    inferred_params = calib_Obj.optimize()
    print(inferred_params)