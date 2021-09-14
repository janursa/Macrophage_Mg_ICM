
import numpy as np
import json
import tellurium as te

# Mg_M = te.loadSBMLModel("Mg_M.xml") 
class MG_MODEL:
    def __init__(self):
        self.create_model()
    def simulate(self,start,end,selections):
        return self.model.simulate(start,end,selections)
    def create_model(self):
        self.originalModel = te.loadSBMLModel("Zhao_2021.xml")
        self.originalModelStr = self.originalModel.getAntimony ()
        self.replacements = {
            # '$ikk_prod => IKK; k191*ikk_prod':'NEMO_IKK + Mg => Mg_NEMO + IKK; k302*Mg*NEMO_IKK- k303*IKK - k304*Mg_NEMO'
            '$ikk_prod => IKK; k191*ikk_prod':'a=>b'

        }
        self.extra_additions = """
            #// Mg initial condition and diffusion
            Mg_e = 0.8;
            Mg = 0.8 ; # internal/physiological Mg
            $Mg_e -> Mg; k301_1*Mg_e - k301_2*Mg ; # Mg diffuses from extra- to intracellular space through345 TRPM7  
            Mg->deg; k302*Mg;
            k301_1 = 0.1; k301_2 = 0.1; k302=0.01;

            #// IKK production
            $ikk_prod => IKK; k303*Mg;
            k303=100; ikk_prod = 1;
            #// params related to the NEMO_IKK + Mg interaction
            #Mg_NEMO = 1000;
            

            #// NEMO_IKK production and degrdation
            #$NEMO_IKK_prod => NEMO_IKK; k308-k310*NEMO_IKK; #production
            #Mg_NEMO -> deg; k309*Mg_NEMO; # degradation
            #NEMO_IKK_prod = 1; NEMO_IKK = 1000;k310=100;k309=.1;k308=100;

            #// TRMP production
            #$TRPM_prod -> TRPM ; k311*Mg^n300 - k312*TRPM;
            #TRPM -> deg; k313*TRPM
            #TRPM_prod = 1; TRPM = 1000; k311 = .1; k312 = .1; k313 = .1;n300=1

            #// activation of M7CKs
            #TRPM7 -> M7CKs; k314*Mg*TRPM7 - k315*M7CKs;
            #M7CKs = 1000; k314 = 100; k315 = 100;

            #// nuclear translocation of M7CKs and vice versa
            #M7CKs -> M7CKs_n; k316*M7CKs*Mg  - k317*M7CKs_n
           # M7CKs_n -> M7CKs ; k318*M7CKs_n - k319*M7CKs
            #M7CKs_n= 100 ; k316=100 ; k317=100; k318=100; k319=100;

            #// activation and phosphorylation of H3S10
            #$H3S10_prod -> H3S10; k320*M7CKs_n - k321*H3S10 ;
            #H3S10 -> deg; k322*H3S10;
            #H3S10 -> pH3S10; k323*M7CKs_n*H3S10-k324*pH3S10
            #H3S10_prod = 1; H3S10 = 100; pH3S10 = 1000; k320=100; k321 = 100; k322 = 100; k323 = 100; k324 = 100; 

            #// IL8 production, extracellular and intracellular transportation
            #$IL8_prod -> IL8; k325*pH3S10 - k326*IL8;
            #IL8 -> deg; k327*IL8;
            #IL8 -> IL8_e; k328*IL8 - k329*IL8_e;
            #IL8_e -> IL8; k330*IL8_e - k331*IL8;
            #IL8 = 100; k325 = 100; k326 = 100; k327 = 100; k328 = 100; k329 = 100; k330 = 100; k331 = 100;

        """
        modified_model = self.modify(self.originalModelStr,self.replacements)
        combined = self.merge(modified_model,self.extra_additions)
        self.model = te.loada(combined)
        self.model.exportToSBML('Mg_M.xml')
    def set(self,key,value):
        self.model[key] = value
    def get(self,key):
        return self.model[key] 
    def reset(self,params = None): # resets the given model and also sets those that cannot be reset by default
        self.model['Mg_e']=0.8
        self.model.reset()
        if params == None:
            pass
        else:
            for key,value in params.items():
                self.model[key] = value
        return self.model
    @staticmethod
    def modify(modelStr,replacements):
        for key,value in replacements.items():
            modelStr = modelStr.replace(key,value)
        return modelStr
    @staticmethod
    def merge(modelStr,sub_model):
        rr = modelStr.split('\nend\n\npad_mac')
        combined = rr[0]+sub_model+'\nend\n\npad_mac is "pad mac"\n'
        return combined
    




