import tellurium as te
import numpy as np
from pathlib import Path
import os
import sys
import random
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
from tools.tools import InvalidParams
from tools.dirs import dir_model


class Macrophage:
    def __init__(self,observations):
      self.model = te.loadSBMLModel(dir_model)
      self.obs = observations
    def simulate_studies(self,params):
        """
        Simulate all studies
        """
        studies_results = {}
        for study in self.obs['studies']:
            try:
              results = self.simulate_study(study,params)
            except InvalidParams:
              raise InvalidParams()
            studies_results[study]= results
        return studies_results

    def simulate_study(self,study,params):
        """
        Simulte one study
        """
        measurement_scheme = self.obs[study]['measurement_scheme']
        experiment_period = self.obs[study]['experiment_period']
        IDs = self.obs[study]['IDs']
        results = {}
        for ID in IDs:
            inputs = self.obs[study][ID]['inputs']
            try:
              ID_results = self.simulate(study=study,params=params,inputs=inputs,measurement_scheme=measurement_scheme,experiment_period=experiment_period)
            except InvalidParams:
              raise InvalidParams()
            results[ID] = ID_results

        return results
    # def cost_study(self,study,study_results):
    #     """
    #     calculates cost values for each ID of the given study
    #     """
    #     measurement_scheme = self.observations[study]['measurement_scheme']
    #     errors = {}
    #     for ID, ID_results in study_results.items():
    #         ID_observations = self.observations[study][ID]['expectations']
    #         target_errors = {}
    #         for target in measurement_scheme.keys():
    #             abs_diff =abs(np.array(ID_results[target])-np.array(ID_observations[target]['mean']))
    #             means = [ID_observations[target]['mean'],ID_results[target]]
    #             mean = np.mean(means)
    #             target_error = abs_diff/mean
    #             target_errors[target] = target_error
    #         errors[ID] = target_errors
    #     return errors
    # def cost_studies(self,studies_results):
    #     """
    #     Calculates cost values for all studies
    #     """
    #     errors = {}
    #     for study,study_results in studies_results.items():
    #         errors[study] = self.cost_study(study,study_results)
    #     return errors
    def sim_vs_exp(self,studies_results):
        """
        sort simulation results along side the experiments
        """
        sorted_results = {}
        for study,study_results in studies_results.items():
          sorted_results[study] = {}
          for ID, ID_results in study_results.items():
              for target,target_results in ID_results.items():
                if target not in sorted_results[study]:
                  sorted_results[study][target] = {'sim' : [],'exp':[]}
                sorted_results[study][target]['sim'] += list(target_results)
                sorted_results[study][target]['exp'] += self.obs[study][ID]['expectations'][target]['mean']
                
        return sorted_results
    @staticmethod
    def quantitative_cost_func(exp,sim): # calculates error for quantitative measurements
      diff = abs(exp - sim)
      mean = 0.5*(exp+sim)
      errors = diff/mean
      error = np.mean(errors)
      if error <0:
        error = 1000
        # raise ValueError()
      return error
    @staticmethod
    def limit_cost_func(exp,sim): # calculates error for given limits 
      error_upper = 0
      error_lower = 0
      if max(sim)>exp[1]:
        diff = abs(max(sim)-exp[1])
        mean = exp[1]
        error_upper = diff/mean
      if min(sim)<exp[0]:
        diff = abs(min(sim)-exp[0])
        mean = exp[0]
        error_lower = diff/mean
      error = error_upper + error_lower
      return error

    def calculate_cost(self,results):
      costs = {}
      for study,study_result in results.items():
        
        costs[study] = {}
        for target,target_results in study_result.items():          
          exp,sim = np.array(target_results['exp']),np.array(target_results['sim'])
          ##---- for certain studies, the simulation results should be normalized ----##
          if study == 'R05_mg_f_n' and target == 'Mg':
            error = self.limit_cost_func(exp=exp,sim=sim)
          elif study == 'Q21_nTRPM' or study == 'Q21_TRPM' or study == 'Q21_nM7CK' or study == 'Q21_H3S10':
            if max(sim)>10:
              error = max(sim)
            else:
              sim_8_error = abs(sim[1]/sim[2]-(1/exp[2]))/(1/exp[2])
              mg_dot8_error =  abs(sim[1]/sim[2]-(1/exp[0]))/(1/exp[0])
              error = (sim_8_error + mg_dot8_error)/2
              # error = sim_8_error
              # print(sim_8_error)

          else:
            try:
              error = self.quantitative_cost_func(exp=exp,sim=sim)
            except ValueError:
              print(study,target)
              aa

          costs[study][target] = error
      mean_errors = []
      for study,study_costs in costs.items():
        errors = list(study_costs.values())
        mean_errors.append(np.mean(errors))
      return costs,np.mean(mean_errors)
    def simulate(self,study,params,inputs,measurement_scheme,experiment_period):
      if experiment_period == None:
        raise ValueError('Value of experiment_period is none')
      target_keys = list(measurement_scheme.keys())
      self.model.reset()
      for key,value in params.items():
          self.model[key] = value
      for key,value in inputs.items():
          self.model[key] = value
      try:
        results_raw = self.model.simulate(0,experiment_period+1,experiment_period+1,selections = ['TIME']+target_keys)
      except RuntimeError:
        raise InvalidParams()
      #   raise InvalidParams()
      results ={}
      for key in target_keys:
        measured_time = measurement_scheme[key]
        results[key] = results_raw[key][measured_time]
      return results
    def run(self,params):
      try:
        results = self.simulate_studies(params)
        sim_vs_exp_results = self.sim_vs_exp(results)
        target_costs,mean_cost = self.calculate_cost(sim_vs_exp_results)
      except InvalidParams:
        mean_cost = 1+random.random()
        # mean_cost = 10
      return mean_cost