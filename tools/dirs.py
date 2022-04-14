from pathlib import Path
import os
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')

"""
Define relative directories

"""
dir_outputs = os.path.join(main_dir,"outputs")
dir_Zhao_model = os.path.join(main_dir,"models/Zhao_2021.xml")
dir_samples_Zhao_model = os.path.join(dir_outputs,'samples.json')
dir_model = os.path.join(main_dir,'models/edited_matlab_SBML.xml')
dir_matlab_model = os.path.join(main_dir,'models/matlab_SBML.xml')
dir_calib_output = os.path.join(dir_outputs,'inferred_parameters.json')
