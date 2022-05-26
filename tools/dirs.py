from pathlib import Path
import os
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')

"""
Define relative directories

"""
dir_outputs = os.path.join(main_dir,"outputs")
dir_Zhao_model = os.path.join(main_dir,"models/Zhao_sbml.xml")
dir_Zhao_model_original = os.path.join(main_dir,"models/Zhao_original.xml")
dir_samples_Zhao_model = os.path.join(dir_outputs,'samples.json')
dir_M1_model = os.path.join(main_dir,'models/M1_sbml.xml')
dir_IL8_model = os.path.join(main_dir,'models/IL8_sbml.xml')
dir_IL6_model = os.path.join(main_dir,'models/IL6_sbml.xml')
dir_model = os.path.join(main_dir,'models/combined.xml')
dir_M1_matlab_model = os.path.join(main_dir,'models/M1_matlab_sbml.xml')
dir_activation_stimuli = os.path.join(main_dir,'models/activation_stimuli.json')
dir_LPS_model = os.path.join(main_dir,'models/LPS_sbml.xml')